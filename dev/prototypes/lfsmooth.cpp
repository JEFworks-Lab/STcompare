// lfsmooth.cpp -- self-contained re-implementation of the smoother used by
// STcompare::matchingVariograms():
//
//   fit <- locfit::locfit(y ~ locfit::lp(long, lat, nn = delta, deg = 0),
//                         kern = "gauss", maxk = 300)
//   fitted(fit)
//
// for 2-d coordinates, locfit 1.5-9.12 semantics.  Everything here mirrors the
// locfit C sources (file:line refer to locfit_1.5-9.12/src):
//
//  * kernel   W(u) = exp(-(GFACT*u)^2/2), GFACT = 2.5          weight.c:31, local.h:137
//             never truncated (gauss is non-compact)           weight.c:42-46
//  * distance spherical (Euclidean), scale = 1 (lp(scale=FALSE)) lf_nbhd.c:14-44
//  * bandwidth h = max(fixh, k-th smallest distance),
//             k = (int)(n*nn + 1e-12)                          locfit.c:355, lf_nbhd.c:94-111
//  * local constant fit: weighted mean + one Newton step       locfit.c:220-239, m_max.c
//  * parametric component: global mean subtracted at vertices,  pcomp.c:55-123,125-149,185-202
//             added back after interpolation
//  * evaluation structure "tree" (rbox(), cut = 0.8)           ev_atree.c, ev_main.c:191-225
//  * fitted(): descend the tree, pseudo-vertices = mean of the
//             two parents, bilinear interpolation (deg = 0 => no
//             vertex derivatives, hasd = 0)                    ev_atree.c:162-205, ev_interp.c:64-80,159-166
//
// Build: Rcpp::sourceCpp("lfsmooth.cpp")
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <utility>

using namespace Rcpp;

static const double GFACT = 2.5;    // local.h:137
static const double NOSLN = 0.1278433; // local.h:136

// ---------------------------------------------------------------------------
// basic pieces
// ---------------------------------------------------------------------------

// rho() for KSPH, d = 2, sca = 1 (lf_nbhd.c:14-44); written exactly like the C
// (s = 0; s += r0*r0; s += r1*r1) so that FP contraction behaves the same.
static inline double rho2(double u0, double u1) {
  double s = 0.0;
  s += u0 * u0;
  s += u1 * u1;
  return std::sqrt(s);
}

// gauss kernel weight, weightsph() + W() (weight.c:13-31, 77-90)
static inline double wgauss(double di, double h) {
  if (h == 0) return (di == 0.0) ? 1.0 : 0.0;
  double u = std::fabs(di / h);
  return std::exp(-((GFACT * u) * (GFACT * u)) / 2.0);
}

// compbandwid() (lf_nbhd.c:94-111) -- scratch is overwritten
static double compbandwid(std::vector<double>& di, std::vector<double>& scratch,
                          int n, int d, int nnk, double fxh) {
  if (nnk == 0) return fxh;
  double nnh;
  if (nnk < n) {
    scratch.assign(di.begin(), di.begin() + n);
    std::nth_element(scratch.begin(), scratch.begin() + (nnk - 1), scratch.end());
    nnh = scratch[nnk - 1];               // k-th smallest (1-based), as kordstat()
  } else {
    nnh = 0;
    for (int i = 0; i < n; i++) nnh = std::max(nnh, di[i]);
    nnh = nnh * std::exp(std::log(1.0 * nnk / n) / d);
  }
  return std::max(fxh, nnh);
}

// Local-constant gaussian fit at point (px,py) with bandwidth h, reproducing
// reginit() + one max_nr() Newton step with JAC_EIGD scaling (locfit.c,
// m_max.c, m_jacob.c, m_eigen.c).  Returns the coefficient.
static double lc_fit(const double* X0, const double* X1, const double* y, int n,
                     const std::vector<double>& di, double h) {
  static thread_local std::vector<double> wbuf;
  if ((int)wbuf.size() < n) wbuf.resize(n);
  double* w = wbuf.data();
  // weights (only w > 0 enter the design, lf_nbhd.c:266-273)
  double s0 = 0.0, s1 = 0.0;
  for (int i = 0; i < n; i++) {
    w[i] = wgauss(di[i], h);
    if (w[i] > 0) { s1 += w[i] * (1.0 * y[i]); s0 += w[i] * 1.0; }
  }
  if (s0 == 0) return NA_REAL;  // LF_NOPT -- cannot happen for gauss in practice
  double m0 = (s1 - 0.0) / s0;
  // one Newton step: f1 = sum w (y - m0), Z = sum w   (likereg, locfit.c:120-148)
  double f1 = 0.0, Z = 0.0;
  for (int i = 0; i < n; i++) {
    if (w[i] > 0) {
      double z = y[i] - m0;
      f1 += 1.0 * w[i] * (1.0 * z);
      Z  += (w[i] * 1.0) * 1.0 * 1.0;
    }
  }
  // jacob_solve with JAC_EIGD for p = 1 (m_jacob.c:31-80, m_eigen.c:71-100)
  double dg = (Z <= 0) ? 0.0 : 1 / std::sqrt(Z);
  double Zs = Z * (dg * dg);
  double v = f1 * dg;
  double tol = 1.0e-8 * Zs;
  double wv = 0.0; wv += 1.0 * v;
  if (Zs > tol) wv /= Zs;
  double x = 0.0; x += 1.0 * wv;
  double delta = x * dg;
  return m0 + 1.0 * delta;      // coef updated before NR_BREAK (m_max.c:33-36)
}

// parametric component for deg = 0: unweighted global mean + Newton step
static double par_comp(const double* y, int n) {
  double s0 = 0.0, s1 = 0.0;
  for (int i = 0; i < n; i++) { s1 += 1.0 * (1.0 * y[i]); s0 += 1.0 * 1.0; }
  double m0 = s1 / s0;
  double f1 = 0.0, Z = 0.0;
  for (int i = 0; i < n; i++) { double z = y[i] - m0; f1 += 1.0 * 1.0 * (1.0 * z); Z += (1.0 * 1.0) * 1.0 * 1.0; }
  double dg = 1 / std::sqrt(Z);
  double Zs = Z * (dg * dg);
  double v = f1 * dg; double wv = 0.0; wv += 1.0 * v;
  if (Zs > 1.0e-8 * Zs) wv /= Zs;
  double x = 0.0; x += 1.0 * wv;
  return m0 + 1.0 * (x * dg);
}

// ---------------------------------------------------------------------------
// 1. exact direct evaluation at the data points (== locfit(..., ev = dat()))
// ---------------------------------------------------------------------------

//' Direct local-constant gaussian smoother at the data points.
// [[Rcpp::export]]
NumericVector lf_direct(NumericMatrix xy, NumericVector y, double nn, double fixh = 0.0) {
  int n = xy.nrow();
  const double* X0 = &xy(0, 0); const double* X1 = &xy(0, 1);
  int nnk = (int)(n * nn + 1e-12);
  std::vector<double> di(n), scratch(n);
  NumericVector out(n);
  // ev = dat(): fitted() returns (coef - pc) + pc (procvraw/subparcomp, fitp_int + addparcomp)
  double pc = par_comp(y.begin(), n);
  for (int p = 0; p < n; p++) {
    for (int i = 0; i < n; i++) di[i] = rho2(X0[i] - X0[p], X1[i] - X1[p]);
    double h = compbandwid(di, scratch, n, 2, nnk, fixh);
    double c = lc_fit(X0, X1, y.begin(), n, di, h) - pc * 1.0;
    out[p] = (c + pc * 1.0) + 0.0;
  }
  return out;
}

//' Direct smoother applied to the B columns of Y at once (weights computed
//' once per evaluation point; plain normalised weighted means, no Newton step).
// [[Rcpp::export]]
NumericMatrix lf_direct_multi(NumericMatrix xy, NumericMatrix Y, double nn, double fixh = 0.0) {
  int n = xy.nrow(), B = Y.ncol();
  if (Y.nrow() != n) Rcpp::stop("nrow(Y) != nrow(xy)");
  const double* X0 = &xy(0, 0); const double* X1 = &xy(0, 1);
  int nnk = (int)(n * nn + 1e-12);
  std::vector<double> di(n), scratch(n), w(n), acc(B);
  NumericMatrix out(n, B);
  for (int p = 0; p < n; p++) {
    for (int i = 0; i < n; i++) di[i] = rho2(X0[i] - X0[p], X1[i] - X1[p]);
    double h = compbandwid(di, scratch, n, 2, nnk, fixh);
    double s0 = 0.0;
    for (int i = 0; i < n; i++) { w[i] = wgauss(di[i], h); s0 += w[i]; }
    std::fill(acc.begin(), acc.end(), 0.0);
    for (int b = 0; b < B; b++) { const double* yb = &Y(0, b); double a = 0.0;
      for (int i = 0; i < n; i++) a += w[i] * yb[i];
      acc[b] = a; }
    for (int b = 0; b < B; b++) out(p, b) = acc[b] / s0;
  }
  return out;
}

//' Dense direct-evaluation operator L (n x n): smooth(y) = L %*% y
// [[Rcpp::export]]
NumericMatrix lf_direct_operator(NumericMatrix xy, double nn, double fixh = 0.0) {
  int n = xy.nrow();
  const double* X0 = &xy(0, 0); const double* X1 = &xy(0, 1);
  int nnk = (int)(n * nn + 1e-12);
  std::vector<double> di(n), scratch(n), w(n);
  NumericMatrix L(n, n);
  for (int p = 0; p < n; p++) {
    for (int i = 0; i < n; i++) di[i] = rho2(X0[i] - X0[p], X1[i] - X1[p]);
    double h = compbandwid(di, scratch, n, 2, nnk, fixh);
    double s0 = 0.0;
    for (int i = 0; i < n; i++) { w[i] = wgauss(di[i], h); s0 += w[i]; }
    for (int i = 0; i < n; i++) L(p, i) = w[i] / s0;
  }
  return L;
}

// ---------------------------------------------------------------------------
// 2. the adaptive tree ("tree" evaluation structure, ev_atree.c)
// ---------------------------------------------------------------------------

struct ATree {
  int d = 2, vc = 4, nv = 0, nvm = 0, ncm = 0;
  double cut = 0.8;
  double fl[4];                      // ll0, ll1, ur0, ur1  (set_flim, startlf.c:60-88)
  std::vector<double> xev;           // 2*nvm
  std::vector<double> h;             // bandwidth (pseudo: mean of parents)
  std::vector<int> s, lo, hi;        // pseudo flag, parents
  std::unordered_map<long long, int> mid;
  // data
  int n = 0; const double* X0 = nullptr; const double* X1 = nullptr;
  int nnk = 0; double fixh = 0.0;
  std::vector<double> di, scratch;
  bool overflow = false;

  static long long key(int a, int b) { if (a > b) std::swap(a, b); return ((long long)a << 32) | (unsigned int)b; }

  int findpt(int i0, int i1) const {
    auto it = mid.find(key(i0, i1));
    return (it == mid.end()) ? -1 : it->second;
  }
  // vertex bandwidth (the only thing the tree needs from a "fit")
  void vfun(int v) {
    double px = xev[2 * v], py = xev[2 * v + 1];
    for (int i = 0; i < n; i++) di[i] = rho2(X0[i] - px, X1[i] - py);
    h[v] = compbandwid(di, scratch, n, 2, nnk, fixh);
  }
  int split(const int* ce, double* le, const double* ll, const double* ur) const {
    double hmin = 0.0, score[2];
    for (int i = 0; i < vc; i++) {
      double hh = h[ce[i]];
      if ((hh > 0) && ((hmin == 0) | (hh < hmin))) hmin = hh;
    }
    int is = 0;
    for (int i = 0; i < d; i++) {
      le[i] = (ur[i] - ll[i]) / 1.0;
      if (hmin == 0) score[i] = 2 * (ur[i] - ll[i]) / (fl[i + d] - fl[i]);
      else score[i] = le[i] / hmin;
      if (score[i] > score[is]) is = i;
    }
    if (cut < score[is]) return is;
    return -1;
  }
  int newsplit(int i0, int i1, int pv) {
    int i = findpt(i0, i1);
    if (i >= 0) return i;
    if (i0 > i1) std::swap(i0, i1);
    int v = nv;
    if (v == nvm) { overflow = true; return -1; }   // locfit: ERROR("newsplit: out of vertex space")
    lo[v] = i0; hi[v] = i1;
    for (int k = 0; k < d; k++) xev[2 * v + k] = (xev[2 * i0 + k] + xev[2 * i1 + k]) / 2;
    if (pv) { h[v] = (h[i0] + h[i1]) / 2; s[v] = 1; }
    else { vfun(v); s[v] = 0; }
    mid[key(i0, i1)] = v;
    nv++;
    return v;
  }
  void grow(const int* ce, double* ll, double* ur) {
    if (overflow) return;
    double le[2];
    int ns = split(ce, le, ll, ur);
    if (ns == -1) return;
    int tk = 1 << ns, nce[4];
    for (int i = 0; i < vc; i++) {
      if ((i & tk) == 0) nce[i] = ce[i];
      else {
        int i0 = ce[i], i1 = ce[i - tk];
        int pv = (le[ns] < (cut * std::min(h[i0], h[i1])));
        nce[i] = newsplit(i0, i1, pv);
        if (overflow) return;
      }
    }
    double z = ur[ns]; ur[ns] = (z + ll[ns]) / 2;
    grow(nce, ll, ur);
    if (overflow) return;
    ur[ns] = z;
    for (int i = 0; i < vc; i++) nce[i] = ((i & tk) == 0) ? nce[i + tk] : ce[i];
    z = ll[ns]; ll[ns] = (z + ur[ns]) / 2;
    grow(nce, ll, ur);
    ll[ns] = z;
  }
};

// atree_guessnv() (ev_atree.c:16-48)
static void atree_guessnv(double cut, int d, double alp, int mk, int* nvm, int* ncm) {
  *ncm = 1 << 30; *nvm = 1 << 30;
  int vc = 1 << d;
  if (alp > 0) {
    double a0 = (alp > 1) ? 1 : 1 / alp;
    if (cut < 0.01) cut = 0.01;
    double cu = 1;
    for (int i = 0; i < d; i++) cu *= std::min(1.0, cut);
    int nv = (int)((5 * a0 / cu + 1) * vc);
    int nc = (int)(10 * a0 / cu + 1);
    if (nv < *nvm) *nvm = nv;
    if (nc < *ncm) *ncm = nc;
  }
  if (*nvm == 1 << 30) { *nvm = 102 * vc; *ncm = 201; }
  double ifl = mk / 100.0;
  *nvm = (int)(ifl * *nvm);
  *ncm = (int)(ifl * *ncm);
}

static void build_tree(ATree& T, NumericMatrix xy, double nn, int maxk, double cut, double fixh) {
  int n = xy.nrow();
  T.n = n; T.X0 = &xy(0, 0); T.X1 = &xy(0, 1);
  T.nnk = (int)(n * nn + 1e-12); T.fixh = fixh; T.cut = cut;
  T.di.resize(n); T.scratch.resize(n);
  atree_guessnv(cut, 2, nn, maxk, &T.nvm, &T.ncm);
  T.xev.assign(2 * T.nvm, 0.0); T.h.assign(T.nvm, 0.0);
  T.s.assign(T.nvm, 0); T.lo.assign(T.nvm, 0); T.hi.assign(T.nvm, 0);
  // bounding box (set_flim)
  for (int k = 0; k < 2; k++) {
    const double* X = (k == 0) ? T.X0 : T.X1;
    double mx = X[0], mn = X[0];
    for (int j = 1; j < n; j++) { mx = std::max(mx, X[j]); mn = std::min(mn, X[j]); }
    T.fl[k] = mn; T.fl[k + 2] = mx;
  }
  double ll[2] = {T.fl[0], T.fl[1]}, ur[2] = {T.fl[2], T.fl[3]};
  int ce[4];
  for (int i = 0; i < 4; i++) {
    int j = i;
    for (int k = 0; k < 2; ++k) { T.xev[2 * i + k] = (j % 2) ? ur[k] : ll[k]; j >>= 1; }
    ce[i] = i;
    T.vfun(i);
    T.s[i] = 0;
  }
  T.nv = 4;
  T.grow(ce, ll, ur);
}

// Descend the tree for point x, returning the 4 corner "values" as sparse
// combinations of vertex ids plus the terminal cell (ll, ur).  Values are
// represented generically through a functor so the same code serves the
// numeric replica and the operator extraction.
struct Combo { std::vector<std::pair<int, double>> t; };

static inline Combo cmb_avg(const Combo& a, const Combo& b) {
  Combo r; r.t.reserve(a.t.size() + b.t.size());
  for (auto& e : a.t) r.t.push_back({e.first, e.second / 2});
  for (auto& e : b.t) r.t.push_back({e.first, e.second / 2});
  return r;
}
static inline Combo cmb_lin(double hh, double dd, const Combo& f0, const Combo& f1) {
  if (dd == 0) return f0;
  Combo r; r.t.reserve(f0.t.size() + f1.t.size());
  for (auto& e : f0.t) r.t.push_back({e.first, e.second * (dd - hh) / dd});
  for (auto& e : f1.t) r.t.push_back({e.first, e.second * hh / dd});
  return r;
}

// numeric atree_int (ev_atree.c:162-205) + rectcell_interp nc==1 (ev_interp.c:64-80)
static double atree_int_num(const ATree& T, const std::vector<double>& coef, double x0, double x1) {
  double x[2] = {x0, x1};
  double vv[4]; int ce[4];
  for (int i = 0; i < 4; i++) { vv[i] = coef[i]; ce[i] = i; }
  double le[2];
  int ns = 0;
  while (ns != -1) {
    const double* ll = &T.xev[2 * ce[0]]; const double* ur = &T.xev[2 * ce[3]];
    ns = T.split(ce, le, ll, ur);
    if (ns != -1) {
      int tk = 1 << ns;
      double hh = ur[ns] - ll[ns];
      int lo = (2 * (x[ns] - ll[ns])) < hh;
      for (int i = 0; i < 4; i++) if ((tk & i) == 0) {
        int nv = T.findpt(ce[i], ce[i + tk]);
        if (nv == -1) Rcpp::stop("Descend tree problem");
        if (lo) { ce[i + tk] = nv; vv[i + tk] = T.s[nv] ? (vv[i] + vv[i + tk]) / 2 : coef[nv]; }
        else    { ce[i] = nv;      vv[i]      = T.s[nv] ? (vv[i] + vv[i + tk]) / 2 : coef[nv]; }
      }
    }
  }
  const double* ll = &T.xev[2 * ce[0]]; const double* ur = &T.xev[2 * ce[3]];
  for (int i = 0; i < 4; i++) if (vv[i] == NOSLN) return NOSLN;
  for (int i = 1; i >= 0; i--) {
    int tk = 1 << i;
    for (int j = 0; j < tk; j++) {
      double hh = x[i] - ll[i], dd = ur[i] - ll[i];
      vv[j] = (dd == 0) ? vv[j] : (((dd - hh) * vv[j] + hh * vv[j + tk]) / dd);
    }
  }
  return vv[0];
}

static Combo atree_int_combo(const ATree& T, double x0, double x1) {
  double x[2] = {x0, x1};
  Combo vv[4]; int ce[4];
  for (int i = 0; i < 4; i++) { vv[i].t = {{i, 1.0}}; ce[i] = i; }
  double le[2];
  int ns = 0;
  while (ns != -1) {
    const double* ll = &T.xev[2 * ce[0]]; const double* ur = &T.xev[2 * ce[3]];
    ns = T.split(ce, le, ll, ur);
    if (ns != -1) {
      int tk = 1 << ns;
      double hh = ur[ns] - ll[ns];
      int lo = (2 * (x[ns] - ll[ns])) < hh;
      for (int i = 0; i < 4; i++) if ((tk & i) == 0) {
        int nv = T.findpt(ce[i], ce[i + tk]);
        if (nv == -1) Rcpp::stop("Descend tree problem");
        Combo nvv; if (T.s[nv]) nvv = cmb_avg(vv[i], vv[i + tk]); else nvv.t = {{nv, 1.0}};
        if (lo) { ce[i + tk] = nv; vv[i + tk] = nvv; } else { ce[i] = nv; vv[i] = nvv; }
      }
    }
  }
  const double* ll = &T.xev[2 * ce[0]]; const double* ur = &T.xev[2 * ce[3]];
  for (int i = 1; i >= 0; i--) {
    int tk = 1 << i;
    for (int j = 0; j < tk; j++) vv[j] = cmb_lin(x[i] - ll[i], ur[i] - ll[i], vv[j], vv[j + tk]);
  }
  return vv[0];
}

static List tree_to_list(const ATree& T) {
  int nv = T.nv;
  NumericMatrix xev(nv, 2);
  NumericVector h(nv); IntegerVector s(nv), lo(nv), hi(nv);
  for (int v = 0; v < nv; v++) {
    xev(v, 0) = T.xev[2 * v]; xev(v, 1) = T.xev[2 * v + 1];
    h[v] = T.h[v]; s[v] = T.s[v]; lo[v] = T.lo[v]; hi[v] = T.hi[v];
  }
  return List::create(_["xev"] = xev, _["h"] = h, _["s"] = s, _["lo"] = lo, _["hi"] = hi,
                      _["nv"] = nv, _["nvm"] = T.nvm, _["ncm"] = T.ncm,
                      _["overflow"] = T.overflow);
}

//' Build the adaptive tree only (depends on coordinates and nn, not on y).
// [[Rcpp::export]]
List lf_tree(NumericMatrix xy, double nn, int maxk = 300, double cut = 0.8, double fixh = 0.0) {
  ATree T; build_tree(T, xy, nn, maxk, cut, fixh);
  return tree_to_list(T);
}

//' Exact replica of fitted(locfit(y ~ lp(x1, x2, nn = nn, deg = 0), kern = "gauss", maxk = maxk)).
//' If the vertex capacity is exceeded, locfit errors ("newsplit: out of vertex
//' space"); here we error too unless allow_overflow = TRUE (then capacity is
//' unlimited).
// [[Rcpp::export]]
NumericVector lf_tree_fitted(NumericMatrix xy, NumericVector y, double nn, int maxk = 300,
                             double cut = 0.8, double fixh = 0.0, bool allow_overflow = false) {
  ATree T;
  build_tree(T, xy, nn, allow_overflow ? 1000000000 : maxk, cut, fixh);
  if (T.overflow) Rcpp::stop("newsplit: out of vertex space");
  int n = T.n;
  double pc = par_comp(y.begin(), n);
  std::vector<double> coef(T.nv, 0.0), di(n);
  for (int v = 0; v < T.nv; v++) {
    if (T.s[v]) continue;
    double px = T.xev[2 * v], py = T.xev[2 * v + 1];
    for (int i = 0; i < n; i++) di[i] = rho2(T.X0[i] - px, T.X1[i] - py);
    double c = lc_fit(T.X0, T.X1, y.begin(), n, di, T.h[v]);
    coef[v] = c - pc * 1.0;   // subparcomp (pcomp.c:139)
  }
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double th = atree_int_num(T, coef, T.X0[i], T.X1[i]);
    th += pc * 1.0;           // addparcomp (pcomp.c:196)
    th += 0.0;                // base
    out[i] = th;
  }
  return out;
}

//' Factorised linear operator of the default (tree) smoother:
//'   fitted = M %*% (W %*% y),  M: n x nreal (interpolation), W: nreal x n
//'   (row-normalised gaussian weights at the real vertices).
// [[Rcpp::export]]
List lf_tree_operator(NumericMatrix xy, double nn, int maxk = 300, double cut = 0.8, double fixh = 0.0,
                      bool allow_overflow = false) {
  ATree T;
  build_tree(T, xy, nn, allow_overflow ? 1000000000 : maxk, cut, fixh);
  if (T.overflow) Rcpp::stop("newsplit: out of vertex space");
  int n = T.n;
  std::vector<int> col(T.nv, -1); int nreal = 0;
  for (int v = 0; v < T.nv; v++) if (!T.s[v]) col[v] = nreal++;
  NumericMatrix W(nreal, n);
  std::vector<double> di(n), w(n);
  for (int v = 0; v < T.nv; v++) {
    if (T.s[v]) continue;
    double px = T.xev[2 * v], py = T.xev[2 * v + 1];
    double s0 = 0.0;
    for (int i = 0; i < n; i++) { di[i] = rho2(T.X0[i] - px, T.X1[i] - py); w[i] = wgauss(di[i], T.h[v]); s0 += w[i]; }
    for (int i = 0; i < n; i++) W(col[v], i) = w[i] / s0;
  }
  NumericMatrix M(n, nreal);
  for (int i = 0; i < n; i++) {
    Combo c = atree_int_combo(T, T.X0[i], T.X1[i]);
    for (auto& e : c.t) {
      if (col[e.first] < 0) Rcpp::stop("pseudo vertex leaked into combo");
      M(i, col[e.first]) += e.second;
    }
  }
  List out = tree_to_list(T);
  out["W"] = W; out["M"] = M; out["nreal"] = nreal;
  return out;
}

//' Apply the factorised operator to many columns at once: M %*% (W %*% Y)
// [[Rcpp::export]]
NumericMatrix lf_apply_factored(NumericMatrix M, NumericMatrix W, NumericMatrix Y) {
  int n = M.nrow(), r = M.ncol(), B = Y.ncol();
  if (W.nrow() != r || W.ncol() != Y.nrow()) Rcpp::stop("dimension mismatch");
  int nY = Y.nrow();
  NumericMatrix V(r, B), out(n, B);
  for (int b = 0; b < B; b++)
    for (int j = 0; j < nY; j++) { double yj = Y(j, b); if (yj == 0) continue;
      for (int k = 0; k < r; k++) V(k, b) += W(k, j) * yj; }
  for (int b = 0; b < B; b++)
    for (int k = 0; k < r; k++) { double vk = V(k, b); if (vk == 0) continue;
      for (int i = 0; i < n; i++) out(i, b) += M(i, k) * vk; }
  return out;
}
