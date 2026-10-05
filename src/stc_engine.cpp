// stc_engine.cpp -- the permutation engine (see stc_engine.h). No R API in this file.
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "stc_engine.h"  // includes stc_fp.h: no FP contraction below
#include "stc_rng.h"

namespace stc {

namespace {

inline double nan_value() { return std::numeric_limits<double>::quiet_NaN(); }

std::string num(double x) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.15g", x);
  return std::string(buf);
}

template <typename T>
void grow(std::vector<T>& v, std::size_t n) {
  if (n < 1) n = 1;  // never empty, so that data() is never a null pointer
  if (v.size() < n) v.resize(n);
}

}  // namespace

const char* engine_status_text(int status) {
  switch (status) {
    case ENG_OK:
      return "ok";
    case ENG_SOURCE_NONFINITE:
      return "the permuted values contain NA, NaN or Inf";
    case ENG_TARGET_NONFINITE:
      return "the values correlated with the surrogates contain NA, NaN or Inf";
    case ENG_SOURCE_CONSTANT:
      return "the permuted values are constant (zero variance)";
    case ENG_TARGET_CONSTANT:
      return "the values correlated with the surrogates are constant (zero variance)";
    case ENG_PLAN_FAILED:
      return "the variogram cannot be computed";
    case ENG_SMOOTHER_FAILED:
      return "the smoother cannot be built for a delta of the grid";
    case ENG_OLS_DEGENERATE:
      return "the variogram of a smoothed permutation is degenerate (lm() gives an NA slope)";
    case ENG_SURROGATE_NONFINITE:
      return "a rescaled permutation has non-finite values";
    case ENG_RSS_NAN:
      return "every residual sum of squares of a permutation is NaN";
    case ENG_COR_NA:
      return "a null correlation is NA";
  }
  return "unknown status";
}

// ---------------------------------------------------------------------------------------------
// Threads
// ---------------------------------------------------------------------------------------------

bool parallel_for(int n_items, int n_threads, const std::function<void(int, int)>& fn,
                  const std::function<bool()>& interrupted) {
  if (n_items <= 0) return true;
  const int nw = std::max(1, std::min(n_threads, n_items));
  std::atomic<int> next(0);
  std::atomic<bool> stop(false);
  std::mutex m;
  std::condition_variable cv;
  int running = nw;
  std::exception_ptr err;

  auto loop = [&](int w) {
    try {
      while (!stop.load(std::memory_order_relaxed)) {
        const int i = next.fetch_add(1, std::memory_order_relaxed);
        if (i >= n_items) break;
        fn(i, w);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lk(m);
      if (!err) err = std::current_exception();
      stop.store(true);
    }
  };
  auto worker = [&](int w) {
    loop(w);
    {
      std::lock_guard<std::mutex> lk(m);
      running--;
    }
    cv.notify_all();
  };

  std::vector<std::thread> threads;
  threads.reserve(nw);
  int started = 0;
  try {
    for (int w = 0; w < nw; w++) {
      threads.emplace_back(worker, w);
      started++;
    }
  } catch (...) {
    // Fewer threads than asked for (std::system_error): the ones that started take all items.
    std::lock_guard<std::mutex> lk(m);
    running -= nw - started;
  }
  if (started == 0) {
    // No thread at all: run on this thread, without interrupt polling.
    loop(0);
    if (err) std::rethrow_exception(err);
    return true;
  }

  bool was_interrupted = false;
  {
    std::unique_lock<std::mutex> lk(m);
    while (running > 0) {
      cv.wait_for(lk, std::chrono::milliseconds(100));
      if (running > 0 && !stop.load()) {
        lk.unlock();
        const bool intr = interrupted ? interrupted() : false;
        lk.lock();
        if (intr) {
          was_interrupted = true;
          stop.store(true);
        }
      }
    }
  }
  for (auto& t : threads) t.join();
  if (err) std::rethrow_exception(err);
  return !was_interrupted;
}

// ---------------------------------------------------------------------------------------------
// Session
// ---------------------------------------------------------------------------------------------

bool Engine::build(const double* x1_, const double* x2_, int n_, const int* ids_, int n_ids,
                   const VariogTable& table, bool plan_ok_, const std::string& plan_reason_,
                   const double* delta_values, int n_deltas, const EngineOptions& options,
                   const std::function<bool()>& interrupted) {
  n = n_;
  x1.assign(x1_, x1_ + n);
  x2.assign(x2_, x2_ + n);
  ids.assign(ids_, ids_ + n_ids);
  vt = table;
  plan_ok = plan_ok_;
  plan_reason = plan_reason_;
  opt = options;
  opt.n_threads = std::max(1, opt.n_threads);
  opt.chunk = std::max(1, opt.chunk);
  deltas.assign(n_deltas, EngineDelta());
  for (int d = 0; d < n_deltas; d++) deltas[d].delta = delta_values[d];
  // The trees and operators depend only on the coordinates and delta: one item per delta.
  return parallel_for(
      n_deltas, opt.n_threads,
      [this](int d, int) {
        EngineDelta& D = deltas[d];
        SmoothTree T;
        D.status = smooth_tree_build(x1.data(), x2.data(), n, D.delta, 300, T);  // maxk = 300
        D.k = T.nnk;
        D.nv = T.nv;
        D.nvm = T.nvm;
        D.depth = T.depth;
        if (D.status == SMOOTH_OK) D.status = smooth_operator_build(T, D.op);
      },
      interrupted);
}

void Engine::set_threads(int n_threads, int chunk) {
  opt.n_threads = std::max(1, n_threads);
  opt.chunk = std::max(1, chunk);
}

// ---------------------------------------------------------------------------------------------
// Tasks and units
// ---------------------------------------------------------------------------------------------

void Engine::define(const double* pool_, int V, const std::vector<int>& task_source,
                    const std::vector<std::vector<int>>& task_grid,
                    const std::vector<std::vector<int>>& task_targets,
                    const std::vector<std::vector<double>>& task_rabs,
                    const std::vector<int>& task_unit, const std::vector<int>& unit_fail_on_dir) {
  if (tasks_defined) throw std::logic_error("the tasks of this session are already defined");
  const int nt = (int)task_source.size();
  const int nu = (int)unit_fail_on_dir.size();
  if ((int)task_grid.size() != nt || (int)task_targets.size() != nt ||
      (int)task_rabs.size() != nt || (int)task_unit.size() != nt) {
    throw std::invalid_argument("task arguments of different lengths");
  }
  for (int t = 0; t < nt; t++) {
    if (task_source[t] < 0 || task_source[t] >= V) throw std::invalid_argument("task source out of range");
    if (task_unit[t] < 0 || task_unit[t] >= nu) throw std::invalid_argument("task unit out of range");
    if (task_grid[t].empty()) throw std::invalid_argument("empty delta grid");
    for (int d : task_grid[t]) {
      if (d < 0 || d >= (int)deltas.size()) throw std::invalid_argument("delta index out of range");
    }
    if (task_rabs[t].size() != task_targets[t].size()) throw std::invalid_argument("one |r| per target");
    for (int v : task_targets[t]) {
      if (v < 0 || v >= V) throw std::invalid_argument("task target out of range");
    }
  }

  // The data vectors: pre-checks that do not depend on the role of a vector.
  pool_cols = V;
  pool.assign(pool_, pool_ + (std::size_t)n * V);
  pool_cor.assign(V, CorTarget());
  pool_status.assign(V, ENG_OK);
  for (int v = 0; v < V; v++) {
    const double* y = &pool[(std::size_t)v * n];
    bool finite = true;
    for (int i = 0; i < n; i++) {
      if (!std::isfinite(y[i])) {
        finite = false;
        break;
      }
    }
    if (!finite) {
      pool_status[v] = ENG_SOURCE_NONFINITE;
    } else {
      bool constant = true;
      for (int i = 1; i < n && constant; i++) constant = (y[i] == y[0]);
      if (constant) pool_status[v] = ENG_SOURCE_CONSTANT;
    }
    cor_target_prepare(y, n, opt.cor_mode, pool_cor[v]);
  }

  tasks.assign(nt, EngineTask());
  fail_b_.reset(new std::atomic<int>[nt]);
  kmax = 0;
  const VariogView view = vt.view();
  const double npairs = (double)(vt.ptr.empty() ? 0 : vt.ptr.back());
  std::vector<double> z(std::max<std::size_t>(1, ids.size()));
  for (int t = 0; t < nt; t++) {
    fail_b_[t].store(INT_MAX);
    EngineTask& T = tasks[t];
    T.source = task_source[t];
    T.grid = task_grid[t];
    T.targets = task_targets[t];
    T.rabs = task_rabs[t];
    T.unit = task_unit[t];
    kmax = std::max(kmax, (int)T.grid.size());
    T.dir_fail_b.assign(T.targets.size(), INT_MAX);
    T.target_status.assign(T.targets.size(), ENG_OK);
    for (std::size_t j = 0; j < T.targets.size(); j++) {
      const int ps = pool_status[T.targets[j]];
      if (ps == ENG_SOURCE_NONFINITE) T.target_status[j] = ENG_TARGET_NONFINITE;
      if (ps == ENG_SOURCE_CONSTANT) T.target_status[j] = ENG_TARGET_CONSTANT;
    }
    T.cost = 0.0;
    for (int d : T.grid) T.cost += (double)deltas[d].op.m * n + 2.0 * npairs + 4.0 * ids.size();

    if (!plan_ok) {
      T.status = ENG_PLAN_FAILED;
      T.message = "the variogram cannot be computed: " + plan_reason;
    } else if (vt.nbins < 2) {
      T.status = ENG_PLAN_FAILED;
      T.message = "the variogram has " + std::to_string(vt.nbins) +
                  " bins with at least 2 pairs; lm() of the variograms needs 2 (try a larger "
                  "maxDistPrctile)";
    } else if (pool_status[T.source] != ENG_OK) {
      T.status = pool_status[T.source];
      T.message = engine_status_text(T.status);
    } else {
      for (int d : T.grid) {
        const EngineDelta& D = deltas[d];
        if (D.status != SMOOTH_OK) {
          T.status = ENG_SMOOTHER_FAILED;
          T.message = "delta = " + num(D.delta) + ": " + smooth_status_message(D.status);
          break;
        }
      }
    }
    if (T.status == ENG_OK) {
      // the target variogram, geoR::variog(X[ids]) of the unpermuted vector
      const double* src = &pool[(std::size_t)T.source * n];
      for (std::size_t r = 0; r < ids.size(); r++) z[r] = src[ids[r]];
      T.tvar.assign(vt.nbins, 0.0);
      variog_eval_tile(view, z.data(), 1, 1, T.tvar.data(), 1);
    }
  }

  units.assign(nu, EngineUnit());
  for (int u = 0; u < nu; u++) units[u].fail_on_direction = unit_fail_on_dir[u] != 0;
  for (int t = 0; t < nt; t++) {
    EngineUnit& U = units[tasks[t].unit];
    U.tasks.push_back(t);
    for (std::size_t j = 0; j < tasks[t].targets.size(); j++) {
      U.dir_task.push_back(t);
      U.dir_target.push_back((int)j);
    }
  }
  for (int u = 0; u < nu; u++) {
    EngineUnit& U = units[u];
    if (U.tasks.empty()) throw std::invalid_argument("a unit has no task");
    U.counts.assign(U.dir_task.size(), 0);
    for (int t : U.tasks) {
      if (tasks[t].status != ENG_OK) {
        fail_unit(U, 0, t, -1, tasks[t].status, tasks[t].message);
        break;
      }
    }
    if (U.state != UNIT_FAILED && U.fail_on_direction) {
      for (std::size_t d = 0; d < U.dir_task.size(); d++) {
        const int st = tasks[U.dir_task[d]].target_status[U.dir_target[d]];
        if (st != ENG_OK) {
          fail_unit(U, 0, U.dir_task[d], (int)d, st, engine_status_text(st));
          break;
        }
      }
    }
  }
  tasks_defined = true;
}

void Engine::fail_unit(EngineUnit& U, int b, int task, int dir, int status, const std::string& msg) {
  U.state = UNIT_FAILED;
  U.status = status;
  U.message = msg;
  U.fail_b = b;
  U.fail_task = task;
  U.fail_dir = dir;
}

// ---------------------------------------------------------------------------------------------
// Batches
// ---------------------------------------------------------------------------------------------

void Engine::ensure_capacity(EngineTask& T, int cap) {
  if (T.cap >= cap) return;
  T.dstar.resize((std::size_t)cap, -1);
  T.nulls.resize((std::size_t)cap * T.targets.size(), nan_value());
  if (opt.keep_surrogates) T.surr.resize((std::size_t)cap * n, nan_value());
  T.cap = cap;
}

void Engine::size_workspaces(const std::vector<int>& task_list, int n_workers) {
  const std::size_t C = (std::size_t)opt.chunk, R = ids.size(), NB = (std::size_t)vt.nbins;
  std::size_t K = 1, summ = 1, nT = 1;
  for (int t : task_list) {
    const EngineTask& T = tasks[t];
    K = std::max(K, T.grid.size());
    std::size_t s = 0;
    for (int d : T.grid) s += (std::size_t)deltas[d].op.m;
    summ = std::max(summ, s);
    nT = std::max(nT, T.targets.size());
  }
  if ((int)ws_.size() < n_workers) ws_.resize(n_workers);
  for (int w = 0; w < n_workers; w++) {
    Workspace& W = ws_[w];
    grow(W.Y, (std::size_t)n * C);
    grow(W.centre, C);
    grow(W.Z, summ * C);
    grow(W.zoff, K);
    grow(W.Xd, R * C);
    grow(W.H, R * C);
    grow(W.G, NB * C);
    grow(W.G2, NB * C);
    grow(W.s0, K * C);
    grow(W.s1, K * C);
    grow(W.rss, K * C);
    grow(W.fail_k, C);
    grow(W.fail_st, C);
    grow(W.surr, (std::size_t)n);
    grow(W.r, nT);
    grow(W.rst, nT);
    grow(W.yptr, nT);
    grow(W.tptr, nT);
  }
}

void Engine::record_task_failure(int task_index, int b, int k, int status) {
  std::lock_guard<std::mutex> lk(fail_mutex_);
  EngineTask& T = tasks[task_index];
  const int fb = fail_b_[task_index].load();
  if (b < fb || (b == fb && k < T.fail_k)) {
    T.fail_k = k;
    T.fail_status = status;
    fail_b_[task_index].store(b);
  }
}

// One work item: permutations b_from + p0 .. b_from + p1 - 1 of one task (dev/engine-spec.md 2.7).
// Every column (permutation) goes through the same operations whatever the other columns are: the
// smoother, variogram, least squares and correlation routines keep a fixed per-column order.
void Engine::process_item(int ti, int p0, int p1, int b_from, Workspace& W) {
  EngineTask& T = tasks[ti];
  const int C = p1 - p0;
  const std::size_t ld = (std::size_t)C;
  const int K = (int)T.grid.size();
  const int R = (int)ids.size();
  const int NB = vt.nbins;
  const VariogView view = vt.view();
  const double* src = &pool[(std::size_t)T.source * n];
  const double* tvar = T.tvar.data();
  double* Y = W.Y.data();
  double* cen = W.centre.data();

  // 1. The permuted vectors X[idx[, b]] as columns of Y (n x C, row-major), each centred with locfit's
  //    parametric component of the permuted vector (a two-pass mean in permuted order, as locfit
  //    computes it on its response). The centre is subtracted here once instead of in every
  //    projection: smooth_project() with the centres would compute the same differences, and
  //    smooth_project_blocked() on the centred columns gives its values bit for bit.
  for (int c = 0; c < C; c++) {
    const int* pc = perm_ + (std::size_t)(p0 + c) * n;
    for (int j = 0; j < n; j++) Y[(std::size_t)j * ld + c] = src[pc[j]];
  }
  for (int c = 0; c < C; c++) cen[c] = smooth_centre(Y + c, ld, n);
  for (int j = 0; j < n; j++) {
    double* yj = Y + (std::size_t)j * ld;
    for (int c = 0; c < C; c++) yj[c] = yj[c] - cen[c];
  }
  for (int c = 0; c < C; c++) {
    W.fail_k[c] = -1;
    W.fail_st[c] = ENG_OK;
  }

  // 2. Every delta of the grid, in grid order (matchingVariograms()'s loop over k).
  std::size_t zo = 0;
  for (int k = 0; k < K; k++) {
    const SmoothOperator& op = deltas[T.grid[k]].op;
    double* Z = W.Z.data() + zo;
    W.zoff[k] = zo;
    zo += (std::size_t)op.m * ld;
    // fitted(locfit(X.randomized ~ lp(long, lat, nn = delta))) on the subsample
    smooth_project_blocked(op, Y, ld, C, Z, ld);
    smooth_interp(op, Z, ld, C, cen, ids.data(), R, W.Xd.data(), ld);
    // variog(X.delta[ids]) and lm(target$v ~ 1 + v)
    variog_eval_tile(view, W.Xd.data(), ld, C, W.G.data(), ld);
    double* s0 = W.s0.data() + (std::size_t)k * ld;
    double* s1 = W.s1.data() + (std::size_t)k * ld;
    for (int c = 0; c < C; c++) {
      double b0, b1;
      const int st = ols_fit(W.G.data() + c, ld, tvar, 1, NB, &b0, &b1);
      if (st != OLS_OK && W.fail_k[c] < 0) {
        W.fail_k[c] = k;
        W.fail_st[c] = ENG_OLS_DEGENERATE;
      }
      s1[c] = std::sqrt(std::fabs(b1));
      s0[c] = std::sqrt(std::fabs(b0));
    }
    // hat = X.delta * sqrt(abs(slope)) + rnorm(N) * sqrt(abs(intercept)), on the subsample
    for (int r = 0; r < R; r++) {
      const std::size_t i = (std::size_t)ids[r];
      const double* xr = W.Xd.data() + (std::size_t)r * ld;
      double* hr = W.H.data() + (std::size_t)r * ld;
      for (int c = 0; c < C; c++) {
        const double e = noise_[((std::size_t)(p0 + c) * noise_k_ + k) * n + i];
        hr[c] = xr[c] * s1[c] + e * s0[c];
      }
    }
    // geoR::variog() fails on non-finite data
    for (int c = 0; c < C; c++) {
      if (W.fail_k[c] >= 0) continue;
      for (int r = 0; r < R; r++) {
        if (!std::isfinite(W.H[(std::size_t)r * ld + c])) {
          W.fail_k[c] = k;
          W.fail_st[c] = ENG_SURROGATE_NONFINITE;
          break;
        }
      }
    }
    // sum((variog(hat[ids])$v - target$v)^2), summed in bin order
    variog_eval_tile(view, W.H.data(), ld, C, W.G2.data(), ld);
    double* rss = W.rss.data() + (std::size_t)k * ld;
    for (int c = 0; c < C; c++) {
      double s = 0.0;
      for (int q = 0; q < NB; q++) {
        const double d = W.G2[(std::size_t)q * ld + c] - tvar[q];
        s += d * d;
      }
      rss[c] = s;
    }
  }

  // 3. delta* (which.min()), the surrogate at all points, and the null correlations.
  const int nT = (int)T.targets.size();
  for (int c = 0; c < C; c++) {
    const int b = b_from + p0 + c;
    double* nul = &T.nulls[(std::size_t)(b - 1) * nT];
    int ks = -1;
    if (W.fail_k[c] >= 0) {
      record_task_failure(ti, b, W.fail_k[c], W.fail_st[c]);
    } else {
      double best = 0.0;
      for (int k = 0; k < K; k++) {
        const double v = W.rss[(std::size_t)k * ld + c];
        if (std::isnan(v)) continue;  // which.min() discards NaN
        if (ks < 0 || v < best) {
          ks = k;
          best = v;
        }
      }
      if (ks < 0) record_task_failure(ti, b, K, ENG_RSS_NAN);
    }
    if (ks < 0) {
      T.dstar[b - 1] = -1;
      for (int t = 0; t < nT; t++) nul[t] = nan_value();
      continue;
    }
    const SmoothOperator& op = deltas[T.grid[ks]].op;
    double* out = opt.keep_surrogates ? &T.surr[(std::size_t)(b - 1) * n] : W.surr.data();
    smooth_interp(op, W.Z.data() + W.zoff[ks] + c, ld, 1, cen + c, nullptr, n, out, 1);
    const double a1 = W.s1[(std::size_t)ks * ld + c];
    const double a0 = W.s0[(std::size_t)ks * ld + c];
    const double* e = noise_ + ((std::size_t)(p0 + c) * noise_k_ + ks) * n;
    for (int i = 0; i < n; i++) out[i] = out[i] * a1 + e[i] * a0;

    int nv = 0;
    for (int t = 0; t < nT; t++) {
      if (T.target_status[t] != ENG_OK) continue;
      W.yptr[nv] = &pool[(std::size_t)T.targets[t] * n];
      W.tptr[nv] = &pool_cor[T.targets[t]];
      nv++;
    }
    cor_with_targets(out, 1, nv, W.yptr.data(), W.tptr.data(), W.r.data(), W.rst.data());
    int q = 0;
    for (int t = 0; t < nT; t++) {
      if (T.target_status[t] != ENG_OK) {
        nul[t] = nan_value();
        continue;
      }
      double r = W.r[q];
      const int st = W.rst[q];
      q++;
      if (st != COR_OK || std::isnan(r)) {
        r = nan_value();
        std::lock_guard<std::mutex> lk(fail_mutex_);
        if (b < T.dir_fail_b[t]) T.dir_fail_b[t] = b;
      }
      nul[t] = r;
    }
    T.dstar[b - 1] = ks;
  }
}

void Engine::scan_unit(EngineUnit& U, int b_from, int b_to, double h, int n_max) {
  const std::size_t nd = U.dir_task.size();
  for (int b = b_from; b <= b_to; b++) {
    for (int t : U.tasks) {
      if (fail_b_[t].load() == b) {
        const EngineTask& T = tasks[t];
        std::string msg = "permutation " + std::to_string(b);
        if (T.fail_k >= 0 && T.fail_k < (int)T.grid.size()) {
          msg += ", delta = " + num(deltas[T.grid[T.fail_k]].delta);
        }
        msg += std::string(": ") + engine_status_text(T.fail_status);
        fail_unit(U, b, t, -1, T.fail_status, msg);
        return;
      }
    }
    if (U.fail_on_direction) {
      for (std::size_t d = 0; d < nd; d++) {
        const EngineTask& T = tasks[U.dir_task[d]];
        if (T.dir_fail_b[U.dir_target[d]] == b) {
          fail_unit(U, b, U.dir_task[d], (int)d, ENG_COR_NA,
                    "permutation " + std::to_string(b) + ": " + engine_status_text(ENG_COR_NA));
          return;
        }
      }
    }
    bool stop = false;
    for (std::size_t d = 0; d < nd; d++) {
      const EngineTask& T = tasks[U.dir_task[d]];
      const int j = U.dir_target[d];
      const double v = T.nulls[(std::size_t)(b - 1) * T.targets.size() + j];
      if (std::fabs(v) >= T.rabs[j]) U.counts[d]++;  // NaN never counts
      if (U.counts[d] >= h) stop = true;
    }
    U.len = b;
    if (stop) {
      U.state = UNIT_STOPPED_H;
      return;
    }
  }
  if (U.len >= n_max) U.state = UNIT_STOPPED_CAP;
}

bool Engine::run_batch(const std::vector<int>& unit_list, int b_from, int b_to, const int* perm,
                       const double* noise, int kn, double h, int n_max,
                       const std::function<bool()>& interrupted) {
  if (!tasks_defined) throw std::logic_error("no tasks defined");
  if (b_from < 1 || b_to < b_from || b_to > n_max) throw std::invalid_argument("invalid batch range");
  const int nb = b_to - b_from + 1;

  // The units to run: re-derive each listed unit's state for this (h, n_max).
  std::vector<int> run_units;
  for (int u : unit_list) {
    if (u < 0 || u >= (int)units.size()) throw std::invalid_argument("unit out of range");
    EngineUnit& U = units[u];
    if (U.state == UNIT_FAILED) continue;
    bool reached = false;
    for (int c : U.counts) reached = reached || (c >= h);
    if (reached) {
      U.state = UNIT_STOPPED_H;
      continue;
    }
    if (U.len >= n_max) {
      U.state = UNIT_STOPPED_CAP;
      continue;
    }
    if (U.len != b_from - 1) {
      throw std::invalid_argument("a unit has " + std::to_string(U.len) +
                                  " permutations; the batch must start at permutation " +
                                  std::to_string(U.len + 1));
    }
    U.state = UNIT_ACTIVE;
    run_units.push_back(u);
  }
  if (run_units.empty()) return true;

  std::vector<int> task_list;
  int kneed = 1;
  for (int u : run_units) {
    for (int t : units[u].tasks) {
      task_list.push_back(t);
      kneed = std::max(kneed, (int)tasks[t].grid.size());
    }
  }
  for (int t : task_list) ensure_capacity(tasks[t], b_to);

  // The legacy noise of the batch, shared by all tasks: block k of permutation b is the k-th
  // rnorm(N) after set.seed(seed + b) under L'Ecuyer-CMRG (block k does not depend on the number of
  // blocks drawn, so tasks with shorter grids read a prefix).
  if (opt.noise_supplied) {
    if (noise == nullptr || kn < kneed) throw std::invalid_argument("the batch needs supplied noise");
    noise_ = noise;
    noise_k_ = kn;
  } else {
    const std::size_t need = (std::size_t)n * kneed * nb;
    if (noise_buf_.size() < need) noise_buf_.resize(need);
    const double seed = opt.seed;
    const int nn = n;
    double* buf = noise_buf_.data();
    const bool ok = parallel_for(
        nb, opt.n_threads,
        [seed, nn, kneed, b_from, buf](int p, int) {
          // as.integer(seed + b): truncation towards zero (the caller checked the range)
          const int s = (int)(seed + (double)(b_from + p));
          legacy_noise(s, nn, kneed, buf + (std::size_t)p * kneed * nn);
        },
        interrupted);
    if (!ok) return false;
    noise_ = noise_buf_.data();
    noise_k_ = kneed;
  }
  perm_ = perm;

  // Work items: (task, sub-chunk of permutations), most expensive tasks first (the order changes no
  // result, only how well the threads are balanced at the end of the batch).
  struct Item {
    int task, p0, p1;
    double cost;
  };
  std::vector<Item> items;
  for (int t : task_list) {
    for (int p0 = 0; p0 < nb; p0 += opt.chunk) {
      const int p1 = std::min(nb, p0 + opt.chunk);
      items.push_back(Item{t, p0, p1, tasks[t].cost * (p1 - p0)});
    }
  }
  std::stable_sort(items.begin(), items.end(),
                   [](const Item& a, const Item& b) { return a.cost > b.cost; });
  const int nw = std::max(1, std::min(opt.n_threads, (int)items.size()));
  size_workspaces(task_list, nw);
  const bool ok = parallel_for(
      (int)items.size(), nw,
      [this, &items, b_from](int i, int w) {
        const Item& it = items[i];
        // A task that already failed at an earlier permutation needs no later ones: skipping them
        // cannot change the smallest failing permutation, nor any result a unit keeps.
        if (b_from + it.p0 > fail_b_[it.task].load(std::memory_order_relaxed)) return;
        process_item(it.task, it.p0, it.p1, b_from, ws_[w]);
      },
      interrupted);
  perm_ = nullptr;
  noise_ = nullptr;
  if (!ok) return false;

  for (int u : run_units) scan_unit(units[u], b_from, b_to, h, n_max);
  return true;
}

}  // namespace stc
