// stc_rcpp.cpp -- thin Rcpp entry points to the engine's building blocks.
//
// Every use of the R API is in this file: the functions validate and copy R objects on the main
// thread, then call the plain C++ routines declared in stc_*.h (which never touch R and can later
// run on worker threads). The R names start with a dot so that NAMESPACE's
// exportPattern("^[[:alpha:]]+") does not export them. They are used by R/engine.R and the tests.
//
// Every export has rng = false: none of them uses R's RNG (the L'Ecuyer-CMRG stream here is a C++
// copy of it), and without it Rcpp wraps each call in GetRNGstate()/PutRNGstate(), which creates a
// time-seeded .Random.seed in the global environment when there is none (dev/engine-spec.md 1.5:
// the global RNG state must be left as it was).
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "stc_rng.h"
#include "stc_smoother.h"
#include "stc_stats.h"
#include "stc_variog.h"

namespace {

inline double na_if_nan(double x) { return std::isnan(x) ? NA_REAL : x; }

// Columns of R matrices are addressed as begin() + column * nrow rather than &X(0, column): Rcpp's
// element access warns ("subscript out of bounds") on a matrix with no rows.

Rcpp::IntegerVector int_vec(const std::vector<int>& v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

Rcpp::NumericVector num_vec(const std::vector<double>& v) {
  return Rcpp::NumericVector(v.begin(), v.end());
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// Variogram
// ---------------------------------------------------------------------------------------------

// Pair table of geoR's binit() loop. x1, x2: coordinates in geoR's column order; lims: geoR's
// bins.lim with the first value set to -1; keep_nugget: geoR's nt.ind. Pair indices are 0-based.
// [[Rcpp::export(name = ".stc_variog_pairs", rng = false)]]
Rcpp::List stc_variog_pairs(Rcpp::NumericVector x1, Rcpp::NumericVector x2, double max_dist,
                            Rcpp::NumericVector lims, bool keep_nugget, int pairs_min = 2) {
  if (x1.size() != x2.size()) Rcpp::stop("x1 and x2 must have the same length");
  stc::VariogTable T;
  const int st = stc::variog_build_table(x1.begin(), x2.begin(), (int)x1.size(), max_dist,
                                         lims.begin(), (int)lims.size(), keep_nugget, pairs_min, T);
  if (st != stc::VARIOG_OK) {
    Rcpp::stop("invalid input for the variogram pair table (2 to 46340 points, ascending limits "
               "starting below 0, pairs_min >= 1)");
  }
  return Rcpp::List::create(
      Rcpp::_["ptr"] = int_vec(T.ptr), Rcpp::_["i"] = int_vec(T.pi), Rcpp::_["j"] = int_vec(T.pj),
      Rcpp::_["bin"] = int_vec(T.bin), Rcpp::_["n"] = num_vec(T.n),
      Rcpp::_["count_all"] = int_vec(T.count_all),
      Rcpp::_["npairs_le_maxdist"] = T.npairs_le_maxdist);
}

// Variograms of the columns of Z (npoints x B) with a pair table (a list with ptr, i, j, n, as
// returned by .stc_variog_pairs() or .stc_variog_plan()). Columns are evaluated in tiles of `tile`
// columns. Returns an nbins x B matrix.
// [[Rcpp::export(name = ".stc_variog_eval", rng = false)]]
Rcpp::NumericMatrix stc_variog_eval(Rcpp::List plan, Rcpp::NumericMatrix Z, int tile = 16) {
  if (!plan.containsElementNamed("ptr") || !plan.containsElementNamed("i") ||
      !plan.containsElementNamed("j") || !plan.containsElementNamed("n")) {
    Rcpp::stop("the plan has no pair table (see plan$ok and plan$reason)");
  }
  Rcpp::IntegerVector ptr = plan["ptr"], pi = plan["i"], pj = plan["j"];
  Rcpp::NumericVector n = plan["n"];
  const int nb = (int)n.size();
  const int M = Z.nrow(), B = Z.ncol();
  if (tile < 1) Rcpp::stop("tile must be at least 1");
  if (ptr.size() != nb + 1 || ptr[0] != 0) Rcpp::stop("invalid pair table: ptr");
  for (int k = 0; k < nb; k++) {
    if (ptr[k + 1] < ptr[k]) Rcpp::stop("invalid pair table: ptr must be non-decreasing");
    if (!(n[k] > 0)) Rcpp::stop("invalid pair table: n must be positive");
  }
  if (ptr[nb] != pi.size() || pi.size() != pj.size()) Rcpp::stop("invalid pair table: pair count");
  for (R_xlen_t p = 0; p < pi.size(); p++) {
    if (pi[p] < 0 || pi[p] >= M || pj[p] < 0 || pj[p] >= M) {
      Rcpp::stop("invalid pair table: a pair index is out of range for nrow(Z)");
    }
  }
  for (R_xlen_t q = 0; q < Z.size(); q++) {
    if (!std::isfinite(Z[q])) Rcpp::stop("NA/NaN/Inf in the data (geoR::variog fails on them too)");
  }
  stc::VariogView view;
  view.nbins = nb;
  view.ptr = ptr.begin();
  view.pi = pi.begin();
  view.pj = pj.begin();
  view.n = n.begin();
  Rcpp::NumericMatrix V(nb, B);
  if (nb == 0 || B == 0) return V;  // nothing to evaluate (and no empty buffers below)
  const int w_max = std::max(1, std::min(tile, B));
  // at least one element each, so that data() is never a null pointer (M can be 0 for a plan without
  // pairs, whose bins are then all empty)
  std::vector<double> zt(std::max<std::size_t>(1, (std::size_t)M * w_max)),
      vt(std::max<std::size_t>(1, (std::size_t)nb * w_max));
  for (int c0 = 0; c0 < B; c0 += w_max) {
    const int w = std::min(w_max, B - c0);
    for (int c = 0; c < w; c++) {
      const double* zc = Z.begin() + (std::size_t)(c0 + c) * M;
      for (int m = 0; m < M; m++) zt[(std::size_t)m * w + c] = zc[m];
    }
    stc::variog_eval_tile(view, zt.data(), (std::size_t)w, w, vt.data(), (std::size_t)w);
    for (int k = 0; k < nb; k++) {
      for (int c = 0; c < w; c++) V(k, c0 + c) = vt[(std::size_t)k * w + c];
    }
  }
  return V;
}

// ---------------------------------------------------------------------------------------------
// Smoother
// ---------------------------------------------------------------------------------------------

// Build locfit's tree for coordinates (x1, x2) = (long, lat) and nn = delta, and optionally the
// factored operator. Never fails for bad delta or capacity: those come back as a status
// (0 ok, 1 bad input, 2 bad delta, 3 N * delta < 2, 4 out of vertex space, 5 internal).
// Wn is returned as an m x N matrix and M as a 0-based CSR list (ptr, col, val). depth is the
// deepest cell of the tree, the recursion depth locfit reaches on the same input.
// [[Rcpp::export(name = ".stc_smoother", rng = false)]]
Rcpp::List stc_smoother(Rcpp::NumericVector x1, Rcpp::NumericVector x2, double delta,
                        int maxk = 300, bool build_operator = true) {
  if (x1.size() != x2.size()) Rcpp::stop("x1 and x2 must have the same length");
  stc::SmoothTree T;
  int st = stc::smooth_tree_build(x1.begin(), x2.begin(), (int)x1.size(), delta, maxk, T);
  // RObject, not SEXP: the objects built in the block below must stay protected after the Rcpp
  // objects that made them go out of scope, until the result list holds them (the allocations that
  // follow, including the result list itself, can run the garbage collector).
  Rcpp::RObject tree, Wn, M;  // R_NilValue unless set
  int m = NA_INTEGER;
  if (st == stc::SMOOTH_OK) {
    Rcpp::NumericMatrix xev(T.nv, 2);
    Rcpp::NumericVector h(T.nv);
    Rcpp::IntegerVector s(T.nv), lo(T.nv), hi(T.nv);
    for (int v = 0; v < T.nv; v++) {
      xev(v, 0) = T.xev[2 * v];
      xev(v, 1) = T.xev[2 * v + 1];
      h[v] = T.h[v];
      s[v] = T.s[v];
      lo[v] = T.lo[v];
      hi[v] = T.hi[v];
    }
    tree = Rcpp::List::create(Rcpp::_["xev"] = xev, Rcpp::_["h"] = h, Rcpp::_["s"] = s,
                              Rcpp::_["lo"] = lo, Rcpp::_["hi"] = hi,
                              Rcpp::_["bbox"] = Rcpp::NumericVector(T.fl, T.fl + 4));
    if (build_operator) {
      stc::SmoothOperator op;
      st = stc::smooth_operator_build(T, op);
      if (st == stc::SMOOTH_OK) {
        Rcpp::NumericMatrix W(op.m, op.n);
        for (int k = 0; k < op.m; k++) {
          for (int j = 0; j < op.n; j++) W(k, j) = op.Wn[(std::size_t)k * op.n + j];
        }
        Wn = W;
        M = Rcpp::List::create(Rcpp::_["ptr"] = int_vec(op.mptr), Rcpp::_["col"] = int_vec(op.mcol),
                               Rcpp::_["val"] = num_vec(op.mval));
        m = op.m;
      }
    }
  }
  return Rcpp::List::create(
      Rcpp::_["status"] = st, Rcpp::_["message"] = std::string(stc::smooth_status_message(st)),
      Rcpp::_["n"] = (int)x1.size(), Rcpp::_["delta"] = delta, Rcpp::_["k"] = T.nnk,
      Rcpp::_["nvm"] = T.nvm, Rcpp::_["nv"] = T.nv, Rcpp::_["depth"] = T.depth, Rcpp::_["m"] = m,
      Rcpp::_["tree"] = tree, Rcpp::_["Wn"] = Wn, Rcpp::_["M"] = M);
}

// Exact replica of fitted(locfit(y ~ lp(x1, x2, nn = delta, deg = 0), kern = "gauss", maxk = maxk)).
// [[Rcpp::export(name = ".stc_smoother_fitted_exact", rng = false)]]
Rcpp::List stc_smoother_fitted_exact(Rcpp::NumericVector x1, Rcpp::NumericVector x2,
                                     Rcpp::NumericVector y, double delta, int maxk = 300) {
  if (x1.size() != x2.size() || y.size() != x1.size()) {
    Rcpp::stop("x1, x2 and y must have the same length");
  }
  stc::SmoothTree T;
  int st = stc::smooth_tree_build(x1.begin(), x2.begin(), (int)x1.size(), delta, maxk, T);
  Rcpp::RObject fitted;  // R_NilValue unless set; RObject keeps it protected (see stc_smoother())
  if (st == stc::SMOOTH_OK) {
    Rcpp::NumericVector f(x1.size());
    st = stc::smooth_fitted_exact(T, y.begin(), f.begin());
    if (st == stc::SMOOTH_OK) fitted = f;
  }
  return Rcpp::List::create(Rcpp::_["status"] = st,
                            Rcpp::_["message"] = std::string(stc::smooth_status_message(st)),
                            Rcpp::_["fitted"] = fitted);
}

// Apply a factored operator (as returned by .stc_smoother()) to the columns of Y (N x B):
// M[rows, ] %*% (Wn %*% (Y - c)) + c, with c the centring constant of each column (locfit's
// parametric component; see smooth_centre() in stc_smoother.h). rows: 1-based data points to
// return, or NULL for all.
// [[Rcpp::export(name = ".stc_smoother_apply", rng = false)]]
Rcpp::NumericMatrix stc_smoother_apply(Rcpp::List op, Rcpp::NumericMatrix Y,
                                       Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue) {
  if (!op.containsElementNamed("Wn") || !op.containsElementNamed("M") || Rf_isNull(op["Wn"]) ||
      Rf_isNull(op["M"])) {
    Rcpp::stop("the smoother has no operator (see its status and message)");
  }
  Rcpp::NumericMatrix W = op["Wn"];
  Rcpp::List Ml = op["M"];
  Rcpp::IntegerVector mptr = Ml["ptr"], mcol = Ml["col"];
  Rcpp::NumericVector mval = Ml["val"];
  stc::SmoothOperator o;
  o.m = W.nrow();
  o.n = W.ncol();
  if (Y.nrow() != o.n) Rcpp::stop("nrow(Y) must equal ncol(Wn)");
  if (mptr.size() != o.n + 1 || mptr[0] != 0 || mptr[o.n] != mcol.size() ||
      mcol.size() != mval.size()) {
    Rcpp::stop("invalid interpolation matrix M");
  }
  for (int i = 0; i < o.n; i++) {
    if (mptr[i + 1] < mptr[i]) Rcpp::stop("invalid interpolation matrix M: ptr");
  }
  for (R_xlen_t q = 0; q < mcol.size(); q++) {
    if (mcol[q] < 0 || mcol[q] >= o.m) Rcpp::stop("invalid interpolation matrix M: column index");
  }
  o.Wn.resize((std::size_t)o.m * o.n);
  for (int k = 0; k < o.m; k++) {
    for (int j = 0; j < o.n; j++) o.Wn[(std::size_t)k * o.n + j] = W(k, j);
  }
  o.mptr.assign(mptr.begin(), mptr.end());
  o.mcol.assign(mcol.begin(), mcol.end());
  o.mval.assign(mval.begin(), mval.end());

  std::vector<int> r0;
  if (rows.isNotNull()) {
    Rcpp::IntegerVector rr(rows.get());
    r0.resize(rr.size());
    for (R_xlen_t q = 0; q < rr.size(); q++) {
      if (rr[q] == NA_INTEGER || rr[q] < 1 || rr[q] > o.n) Rcpp::stop("rows out of range");
      r0[q] = rr[q] - 1;
    }
  }
  const int nrows = rows.isNotNull() ? (int)r0.size() : o.n;
  const int B = Y.ncol();
  std::vector<double> Yt((std::size_t)o.n * B), Z((std::size_t)o.m * B), out((std::size_t)nrows * B);
  std::vector<double> centre(B);
  for (int b = 0; b < B; b++) {
    const double* yb = Y.begin() + (std::size_t)b * o.n;
    centre[b] = stc::smooth_centre(yb, 1, o.n);
    for (int j = 0; j < o.n; j++) Yt[(std::size_t)j * B + b] = yb[j];
  }
  const double* cen = B > 0 ? centre.data() : nullptr;
  stc::smooth_project(o, Yt.data(), (std::size_t)B, B, cen, Z.data(), (std::size_t)B);
  stc::smooth_interp(o, Z.data(), (std::size_t)B, B, cen, rows.isNotNull() ? r0.data() : nullptr,
                     nrows, out.data(), (std::size_t)B);
  Rcpp::NumericMatrix res(nrows, B);
  for (int r = 0; r < nrows; r++) {
    for (int b = 0; b < B; b++) res(r, b) = out[(std::size_t)r * B + b];
  }
  return res;
}

// ---------------------------------------------------------------------------------------------
// RNG
// ---------------------------------------------------------------------------------------------

namespace {
void check_seed(int seed) {
  if (seed == NA_INTEGER) Rcpp::stop("supplied seed is not a valid integer");
}
}  // namespace

// State after set.seed(seed) under L'Ecuyer-CMRG, as unsigned values (.Random.seed[2:7] read as
// unsigned 32-bit integers).
// [[Rcpp::export(name = ".stc_lecuyer_seed", rng = false)]]
Rcpp::NumericVector stc_lecuyer_seed(int seed) {
  check_seed(seed);
  stc::LecuyerCMRG g(seed);
  std::uint32_t s[6];
  g.get_state(s);
  Rcpp::NumericVector out(6);
  for (int j = 0; j < 6; j++) out[j] = (double)s[j];
  return out;
}

// runif(n) after RNGkind("L'Ecuyer-CMRG"); set.seed(seed).
// [[Rcpp::export(name = ".stc_lecuyer_runif", rng = false)]]
Rcpp::NumericVector stc_lecuyer_runif(int seed, int n) {
  check_seed(seed);
  if (n < 0) Rcpp::stop("n must be non-negative");
  stc::LecuyerCMRG g(seed);
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = g.unif_rand();
  return out;
}

// rnorm(n) after RNGkind("L'Ecuyer-CMRG", "Inversion"); set.seed(seed).
// [[Rcpp::export(name = ".stc_lecuyer_rnorm", rng = false)]]
Rcpp::NumericVector stc_lecuyer_rnorm(int seed, int n) {
  check_seed(seed);
  if (n < 0) Rcpp::stop("n must be non-negative");
  stc::LecuyerCMRG g(seed);
  Rcpp::NumericVector out(n);
  g.rnorm(out.begin(), (std::size_t)n);
  return out;
}

// The legacy noise of one permutation: an N x K matrix whose column k is the k-th rnorm(N) after
// set.seed(seed) under L'Ecuyer-CMRG.
// [[Rcpp::export(name = ".stc_legacy_noise", rng = false)]]
Rcpp::NumericMatrix stc_legacy_noise(int seed, int N, int K) {
  check_seed(seed);
  if (N < 0 || K < 0) Rcpp::stop("N and K must be non-negative");
  Rcpp::NumericMatrix out(N, K);
  stc::legacy_noise(seed, N, K, out.begin());
  return out;
}

// ---------------------------------------------------------------------------------------------
// Correlation and least squares
// ---------------------------------------------------------------------------------------------

// cor(X[, b], y) for every column b, as R's cor(X, y) computes it (NA where R gives NA).
// mode: 0 plain (long double), 1 every accumulation fused, 2 CRAN macOS arm64 pattern, 3 double
// accumulators (see stc_stats.h).
// [[Rcpp::export(name = ".stc_cor_cols", rng = false)]]
Rcpp::NumericVector stc_cor_cols(Rcpp::NumericMatrix X, Rcpp::NumericVector y, int mode = 0) {
  const int n = (int)y.size();
  if (X.nrow() != n) Rcpp::stop("nrow(X) must equal length(y)");
  if (mode < 0 || mode > 3) Rcpp::stop("mode must be 0, 1, 2 or 3");
  stc::CorTarget t;
  stc::cor_target_prepare(y.begin(), n, mode, t);
  Rcpp::NumericVector out(X.ncol());
  for (int b = 0; b < X.ncol(); b++) {
    double r;
    const int st = stc::cor_with_target(X.begin() + (std::size_t)b * n, 1, y.begin(), t, &r);
    out[b] = (st == stc::COR_OK) ? r : NA_REAL;
  }
  return out;
}

// lm(y ~ 1 + x)$coefficients by the closed form, with lm()'s rank-deficiency rule.
// status: 0 ok, 1 rank deficient (NA slope), 2 fewer than 2 observations, 3 non-finite input.
// [[Rcpp::export(name = ".stc_ols", rng = false)]]
Rcpp::List stc_ols(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  if (x.size() != y.size()) Rcpp::stop("x and y must have the same length");
  double b0, b1;
  const int st = stc::ols_fit(x.begin(), 1, y.begin(), 1, (int)x.size(), &b0, &b1);
  return Rcpp::List::create(
      Rcpp::_["coefficients"] = Rcpp::NumericVector::create(na_if_nan(b0), na_if_nan(b1)),
      Rcpp::_["status"] = st);
}
