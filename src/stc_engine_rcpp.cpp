// stc_engine_rcpp.cpp -- Rcpp entry points of the permutation engine (src/stc_engine.h), used by
// R/engine.R. Every R object is validated and copied here, on R's main thread, before any worker
// thread starts; the engine itself never touches R. The only R call made while workers run is the
// interrupt check below, on the main thread.
//
// Every export has a dot-prefixed name (not exported by NAMESPACE's exportPattern) and rng = false
// (no GetRNGstate()/PutRNGstate(), which would create a .Random.seed; see src/stc_rcpp.cpp).
#include <Rcpp.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "stc_engine.h"

namespace {

void check_interrupt_fn(void*) { R_CheckUserInterrupt(); }

// TRUE when the user has asked for an interrupt. R_CheckUserInterrupt() jumps out on an interrupt;
// R_ToplevelExec() catches the jump and returns FALSE (the interrupt is then consumed and re-raised
// below with Rcpp's InterruptedException, after the workers have joined). Main thread only.
bool r_interrupted() { return R_ToplevelExec(check_interrupt_fn, nullptr) == FALSE; }

SEXP engine_tag() { return Rf_install("stc_engine_session"); }

stc::Engine* get_engine(SEXP s) {
  if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != engine_tag()) {
    Rcpp::stop("not an engine session (see .stc_engine_new())");
  }
  stc::Engine* e = static_cast<stc::Engine*>(R_ExternalPtrAddr(s));
  if (e == nullptr) {
    Rcpp::stop("the engine session is no longer valid (it does not survive saving and reloading)");
  }
  if (e->magic != stc::Engine::kMagic) Rcpp::stop("not an engine session");
  return e;
}

void finalize_engine(SEXP s) {
  stc::Engine* e = static_cast<stc::Engine*>(R_ExternalPtrAddr(s));
  if (e != nullptr) {
    e->magic = 0;
    delete e;
    R_ClearExternalPtr(s);
  }
}

std::vector<int> int_vector(SEXP x, const char* what) {
  if (TYPEOF(x) != INTSXP) Rcpp::stop("%s must be an integer vector", what);
  Rcpp::IntegerVector v(x);
  for (R_xlen_t i = 0; i < v.size(); i++) {
    if (v[i] == NA_INTEGER) Rcpp::stop("%s must not contain NA", what);
  }
  return std::vector<int>(v.begin(), v.end());
}

int check_threads(int n_threads, int chunk) {
  if (n_threads == NA_INTEGER || n_threads < 1) Rcpp::stop("n_threads must be at least 1");
  if (chunk == NA_INTEGER || chunk < 1) Rcpp::stop("chunk must be at least 1");
  return 0;
}

}  // namespace

// A new engine session (dev/engine-spec.md 3.1). x1 = long = pos[, 2] and x2 = lat = pos[, 1] (locfit's
// order), ids the 0-based variogram subsample, plan the result of .stc_variog_plan() on ids (ok = FALSE
// gives tasks that fail their pre-check with plan$reason), deltas the distinct delta values of all
// tasks: their smoothers are built here, in parallel. cor_mode from .stc_cor_mode(); seed the legacy
// seed (noise of permutation b: set.seed(seed + b) under L'Ecuyer-CMRG); noise_supplied: the noise of
// every batch comes from R instead.
// [[Rcpp::export(name = ".stc_engine_new", rng = false)]]
SEXP stc_engine_new(Rcpp::NumericVector x1, Rcpp::NumericVector x2, Rcpp::IntegerVector ids,
                    Rcpp::List plan, Rcpp::NumericVector deltas, int cor_mode, double seed,
                    bool noise_supplied, int n_threads, int chunk, bool keep_surrogates) {
  const R_xlen_t n = x1.size();
  if (x2.size() != n) Rcpp::stop("x1 and x2 must have the same length");
  if (n < 2 || n > INT_MAX / 2) Rcpp::stop("the engine needs at least 2 points");
  for (R_xlen_t i = 0; i < n; i++) {
    if (!std::isfinite(x1[i]) || !std::isfinite(x2[i])) Rcpp::stop("coordinates must be finite");
  }
  if (ids.size() < 1 || ids.size() > n) Rcpp::stop("invalid subsample");
  for (R_xlen_t r = 0; r < ids.size(); r++) {
    if (ids[r] == NA_INTEGER || ids[r] < 0 || ids[r] >= n) Rcpp::stop("subsample index out of range");
  }
  if (cor_mode == NA_INTEGER || cor_mode < 0 || cor_mode > 3) Rcpp::stop("cor_mode must be 0 to 3");
  if (!std::isfinite(seed)) Rcpp::stop("seed must be finite");
  check_threads(n_threads, chunk);
  if (deltas.size() < 1) Rcpp::stop("no delta values");

  bool plan_ok = plan.containsElementNamed("ok") && Rcpp::as<bool>(plan["ok"]);
  std::string reason = plan.containsElementNamed("reason") ? Rcpp::as<std::string>(plan["reason"]) : "";
  stc::VariogTable vt;
  if (plan_ok) {
    if (!plan.containsElementNamed("ptr") || !plan.containsElementNamed("i") ||
        !plan.containsElementNamed("j") || !plan.containsElementNamed("n")) {
      Rcpp::stop("the plan has no pair table");
    }
    Rcpp::IntegerVector ptr = plan["ptr"], pi = plan["i"], pj = plan["j"];
    Rcpp::NumericVector pn = plan["n"];
    const int nb = (int)pn.size();
    if (ptr.size() != nb + 1 || ptr[0] != 0) Rcpp::stop("invalid pair table: ptr");
    for (int k = 0; k < nb; k++) {
      if (ptr[k + 1] < ptr[k]) Rcpp::stop("invalid pair table: ptr must be non-decreasing");
      if (!(pn[k] > 0)) Rcpp::stop("invalid pair table: n must be positive");
    }
    if (ptr[nb] != pi.size() || pi.size() != pj.size()) Rcpp::stop("invalid pair table: pair count");
    const int nid = (int)ids.size();
    for (R_xlen_t p = 0; p < pi.size(); p++) {
      if (pi[p] < 0 || pi[p] >= nid || pj[p] < 0 || pj[p] >= nid) {
        Rcpp::stop("invalid pair table: a pair index is out of range for the subsample");
      }
    }
    vt.npoints = nid;
    vt.nbins = nb;
    vt.ptr.assign(ptr.begin(), ptr.end());
    vt.pi.assign(pi.begin(), pi.end());
    vt.pj.assign(pj.begin(), pj.end());
    vt.n.assign(pn.begin(), pn.end());
    vt.bin.assign(nb, 0);
  }
  stc::EngineOptions opt;
  opt.n_threads = n_threads;
  opt.chunk = chunk;
  opt.keep_surrogates = keep_surrogates;
  opt.cor_mode = cor_mode;
  opt.seed = seed;
  opt.noise_supplied = noise_supplied;

  std::unique_ptr<stc::Engine> e(new stc::Engine());
  const bool ok = e->build(x1.begin(), x2.begin(), (int)n, ids.begin(), (int)ids.size(), vt, plan_ok,
                           reason, deltas.begin(), (int)deltas.size(), opt, r_interrupted);
  if (!ok) throw Rcpp::internal::InterruptedException();
  SEXP p = PROTECT(R_MakeExternalPtr(e.release(), engine_tag(), R_NilValue));
  R_RegisterCFinalizerEx(p, finalize_engine, TRUE);
  UNPROTECT(1);
  return p;
}

// The smoother of every delta of a session.
// [[Rcpp::export(name = ".stc_engine_deltas", rng = false)]]
Rcpp::List stc_engine_deltas(SEXP session) {
  stc::Engine* e = get_engine(session);
  const int nd = (int)e->deltas.size();
  Rcpp::NumericVector delta(nd);
  Rcpp::IntegerVector status(nd), k(nd), m(nd), nv(nd), nvm(nd), depth(nd);
  Rcpp::CharacterVector message(nd);
  for (int d = 0; d < nd; d++) {
    const stc::EngineDelta& D = e->deltas[d];
    delta[d] = D.delta;
    status[d] = D.status;
    message[d] = stc::smooth_status_message(D.status);
    k[d] = D.k;
    m[d] = D.op.m;
    nv[d] = D.nv;
    nvm[d] = D.nvm;
    depth[d] = D.depth;
  }
  return Rcpp::List::create(Rcpp::_["delta"] = delta, Rcpp::_["status"] = status,
                            Rcpp::_["message"] = message, Rcpp::_["k"] = k, Rcpp::_["m"] = m,
                            Rcpp::_["nv"] = nv, Rcpp::_["nvm"] = nvm, Rcpp::_["depth"] = depth);
}

// Tasks and units of a session (once). pool: N x V data vectors. Per task (0-based indices):
// task_source (pool column), task_grid (list of delta indices), task_targets (list of pool columns),
// task_rabs (list of |r| per target), task_unit. unit_fail_on_dir: one logical per unit.
// [[Rcpp::export(name = ".stc_engine_define", rng = false)]]
Rcpp::List stc_engine_define(SEXP session, Rcpp::NumericMatrix pool, Rcpp::IntegerVector task_source,
                             Rcpp::List task_grid, Rcpp::List task_targets, Rcpp::List task_rabs,
                             Rcpp::IntegerVector task_unit, Rcpp::LogicalVector unit_fail_on_dir) {
  stc::Engine* e = get_engine(session);
  if (e->tasks_defined) Rcpp::stop("the tasks of this session are already defined");
  if (pool.nrow() != e->n) Rcpp::stop("nrow(pool) must equal the number of points");
  const int nt = (int)task_source.size();
  if (task_grid.size() != nt || task_targets.size() != nt || task_rabs.size() != nt ||
      task_unit.size() != nt) {
    Rcpp::stop("task arguments of different lengths");
  }
  std::vector<int> src = int_vector(task_source, "task_source"), unit = int_vector(task_unit, "task_unit");
  std::vector<std::vector<int>> grid(nt), targets(nt);
  std::vector<std::vector<double>> rabs(nt);
  for (int t = 0; t < nt; t++) {
    grid[t] = int_vector(task_grid[t], "task_grid");
    targets[t] = int_vector(task_targets[t], "task_targets");
    Rcpp::NumericVector ra = task_rabs[t];
    rabs[t].assign(ra.begin(), ra.end());
  }
  std::vector<int> fod(unit_fail_on_dir.size());
  for (R_xlen_t u = 0; u < unit_fail_on_dir.size(); u++) {
    if (unit_fail_on_dir[u] == NA_LOGICAL) Rcpp::stop("unit_fail_on_dir must not contain NA");
    fod[u] = unit_fail_on_dir[u] ? 1 : 0;
  }
  // R_xlen_t-safe: pool has n * V values
  e->define(pool.begin(), pool.ncol(), src, grid, targets, rabs, unit, fod);

  Rcpp::IntegerVector status(nt);
  Rcpp::CharacterVector message(nt);
  Rcpp::List target_status(nt);
  for (int t = 0; t < nt; t++) {
    status[t] = e->tasks[t].status;
    message[t] = e->tasks[t].message;
    target_status[t] = Rcpp::IntegerVector(e->tasks[t].target_status.begin(),
                                           e->tasks[t].target_status.end());
  }
  return Rcpp::List::create(Rcpp::_["status"] = status, Rcpp::_["message"] = message,
                            Rcpp::_["target_status"] = target_status, Rcpp::_["kmax"] = e->kmax);
}

// One batch: permutations b_from..b_to for the listed units (0-based), with perm the N x nb matrix of
// 1-based permutation indices (column j is permutation b_from + j - 1), noise NULL or N x K x nb normals
// (sessions with noise_supplied), and the stopping rule h (Inf: none) and cap n_max >= b_to.
// [[Rcpp::export(name = ".stc_engine_run", rng = false)]]
void stc_engine_run(SEXP session, Rcpp::IntegerVector units, int b_from, int b_to,
                    Rcpp::IntegerMatrix perm, SEXP noise, double h, int n_max) {
  stc::Engine* e = get_engine(session);
  if (!e->tasks_defined) Rcpp::stop("the session has no tasks (see .stc_engine_define())");
  if (b_from == NA_INTEGER || b_to == NA_INTEGER || n_max == NA_INTEGER || b_from < 1 ||
      b_to < b_from || n_max < b_to) {
    Rcpp::stop("invalid batch: need 1 <= b_from <= b_to <= n_max");
  }
  if (std::isnan(h) || h < 1) Rcpp::stop("h must be at least 1 (Inf for no stopping)");
  const int n = e->n;
  const int nb = b_to - b_from + 1;
  if (perm.nrow() != n || perm.ncol() != nb) Rcpp::stop("perm must be an N x (b_to - b_from + 1) matrix");
  std::vector<int> ul = int_vector(units, "units");
  for (int u : ul) {
    if (u < 0 || u >= (int)e->units.size()) Rcpp::stop("unit index out of range");
  }
  std::vector<int> P((std::size_t)n * nb);
  const int* pp = perm.begin();
  for (std::size_t q = 0; q < P.size(); q++) {
    const int v = pp[q];
    if (v == NA_INTEGER || v < 1 || v > n) Rcpp::stop("permutation index out of range");
    P[q] = v - 1;
  }
  std::vector<double> noise_copy;
  int kn = 0;
  if (e->opt.noise_supplied) {
    if (TYPEOF(noise) != REALSXP) Rcpp::stop("this session needs the noise of every batch");
    const R_xlen_t len = Rf_xlength(noise);
    if (len == 0 || len % ((R_xlen_t)n * nb) != 0) Rcpp::stop("noise must have N * K * nb values");
    kn = (int)(len / ((R_xlen_t)n * nb));
    const double* nz = REAL(noise);
    noise_copy.assign(nz, nz + len);
  } else if (!Rf_isNull(noise)) {
    Rcpp::stop("this session generates its own noise (noise must be NULL)");
  }
  const bool ok = e->run_batch(ul, b_from, b_to, P.data(), noise_copy.empty() ? nullptr : noise_copy.data(),
                               kn, h, n_max, r_interrupted);
  if (!ok) throw Rcpp::internal::InterruptedException();
}

// The state of every unit: permutations kept (len), state (0 active, 1 stopped by h exceedances,
// 2 stopped at n_max, 3 failed), failure status, message, permutation and task, and the exceedance
// counts of every direction.
// [[Rcpp::export(name = ".stc_engine_units", rng = false)]]
Rcpp::List stc_engine_units(SEXP session) {
  stc::Engine* e = get_engine(session);
  const int nu = (int)e->units.size();
  Rcpp::IntegerVector len(nu), state(nu), status(nu), fail_b(nu), fail_task(nu), fail_dir(nu);
  Rcpp::CharacterVector message(nu);
  Rcpp::List counts(nu);
  for (int u = 0; u < nu; u++) {
    const stc::EngineUnit& U = e->units[u];
    len[u] = U.len;
    state[u] = U.state;
    status[u] = U.status;
    message[u] = U.message;
    fail_b[u] = U.state == stc::UNIT_FAILED ? U.fail_b : NA_INTEGER;
    fail_task[u] = U.fail_task >= 0 ? U.fail_task : NA_INTEGER;
    fail_dir[u] = U.fail_dir >= 0 ? U.fail_dir : NA_INTEGER;
    counts[u] = Rcpp::IntegerVector(U.counts.begin(), U.counts.end());
  }
  return Rcpp::List::create(Rcpp::_["len"] = len, Rcpp::_["state"] = state, Rcpp::_["status"] = status,
                            Rcpp::_["message"] = message, Rcpp::_["fail_b"] = fail_b,
                            Rcpp::_["fail_task"] = fail_task, Rcpp::_["fail_dir"] = fail_dir,
                            Rcpp::_["counts"] = counts);
}

// Results of tasks (0-based) for the permutations their unit keeps (1..len): dstar (1-based grid
// positions, NA where a permutation failed), nulls (len x T), surrogates (N x len, or NULL), the
// pre-check status and message, and per target its pre-check status and first failing permutation.
// [[Rcpp::export(name = ".stc_engine_task_results", rng = false)]]
Rcpp::List stc_engine_task_results(SEXP session, Rcpp::IntegerVector tasks, bool surrogates) {
  stc::Engine* e = get_engine(session);
  const int n = e->n;
  Rcpp::List out(tasks.size());
  for (R_xlen_t q = 0; q < tasks.size(); q++) {
    const int t = tasks[q];
    if (t == NA_INTEGER || t < 0 || t >= (int)e->tasks.size()) Rcpp::stop("task index out of range");
    const stc::EngineTask& T = e->tasks[t];
    const int L = std::min(e->units[T.unit].len, T.cap);
    const int nT = (int)T.targets.size();
    Rcpp::IntegerVector dstar(L);
    for (int b = 0; b < L; b++) dstar[b] = T.dstar[b] >= 0 ? T.dstar[b] + 1 : NA_INTEGER;
    Rcpp::NumericMatrix nulls(L, nT);
    for (int b = 0; b < L; b++) {
      for (int j = 0; j < nT; j++) {
        const double v = T.nulls[(std::size_t)b * nT + j];
        nulls[(std::size_t)j * L + b] = std::isnan(v) ? NA_REAL : v;
      }
    }
    Rcpp::RObject surr;  // R_NilValue unless kept and asked for
    if (surrogates && e->opt.keep_surrogates) {
      Rcpp::NumericMatrix S(n, L);
      if (L > 0) std::copy(T.surr.begin(), T.surr.begin() + (std::size_t)n * L, S.begin());
      surr = S;
    }
    Rcpp::IntegerVector dfb(nT);
    for (int j = 0; j < nT; j++) dfb[j] = T.dir_fail_b[j] == INT_MAX ? NA_INTEGER : T.dir_fail_b[j];
    out[q] = Rcpp::List::create(
        Rcpp::_["dstar"] = dstar, Rcpp::_["nulls"] = nulls, Rcpp::_["surrogates"] = surr,
        Rcpp::_["status"] = T.status, Rcpp::_["message"] = T.message,
        Rcpp::_["target_status"] = Rcpp::IntegerVector(T.target_status.begin(), T.target_status.end()),
        Rcpp::_["dir_fail_b"] = dfb);
  }
  return out;
}

// Change the number of worker threads and the permutations per work item of a session.
// [[Rcpp::export(name = ".stc_engine_set_threads", rng = false)]]
void stc_engine_set_threads(SEXP session, int n_threads, int chunk) {
  stc::Engine* e = get_engine(session);
  check_threads(n_threads, chunk);
  e->set_threads(n_threads, chunk);
}

// Test hook for the thread pool (stc::parallel_for()): n_items items on n_threads workers; item i
// sleeps sleep_ms milliseconds, then records the worker that ran it; item fail_item (0-based; -1 for
// none) throws instead. Returns the worker of every item (NA for items abandoned after a failure), or
// an R error carrying the worker's message; an interrupt while it waits is re-raised in R.
// [[Rcpp::export(name = ".stc_parallel_selftest", rng = false)]]
Rcpp::IntegerVector stc_parallel_selftest(int n_items, int n_threads, int fail_item = -1,
                                          int sleep_ms = 0) {
  if (n_items == NA_INTEGER || n_items < 0) Rcpp::stop("n_items must be non-negative");
  check_threads(n_threads, 1);
  std::vector<int> who((std::size_t)n_items, NA_INTEGER);
  const bool ok = stc::parallel_for(
      n_items, n_threads,
      [&who, fail_item, sleep_ms](int i, int w) {
        if (sleep_ms > 0) std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
        if (i == fail_item) throw std::runtime_error("item " + std::to_string(i) + " failed on a worker");
        who[(std::size_t)i] = w;
      },
      r_interrupted);
  if (!ok) throw Rcpp::internal::InterruptedException();
  return Rcpp::IntegerVector(who.begin(), who.end());
}
