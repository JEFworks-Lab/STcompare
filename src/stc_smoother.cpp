// stc_smoother.cpp -- locfit's adaptive-tree smoother (see stc_smoother.h).
//
// Ported from dev/prototypes/lfsmooth.cpp; the code mirrors the locfit 1.5-9.12 C sources, and the
// comments cite them as file:line (locfit/src/). Statements are written as in locfit so that the
// rounding is the same.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "stc_smoother.h"  // includes stc_fp.h: no FP contraction below

namespace stc {

namespace {

const double GFACT = 2.5;        // local.h:137
const double NOSLN = 0.1278433;  // local.h:136

// rho() for KSPH, d = 2, sca = 1 (lf_nbhd.c:14-44): s = 0; s += u0*u0; s += u1*u1; sqrt(s).
inline double rho2(double u0, double u1) {
  double s = 0.0;
  s += u0 * u0;
  s += u1 * u1;
  return std::sqrt(s);
}

// Gaussian weight: weightsph() + W() (weight.c:13-31, 77-90). Never truncated.
inline double wgauss(double di, double h) {
  if (h == 0) return (di == 0.0) ? 1.0 : 0.0;
  double u = std::fabs(di / h);
  return std::exp(-((GFACT * u) * (GFACT * u)) / 2.0);
}

// compbandwid() (lf_nbhd.c:94-111), d = 2, fixh = 0: the k-th smallest distance (kordstat), or for
// k >= n the largest distance times (k / n)^(1/2).
double compbandwid(const double* di, std::vector<double>& scratch, int n, int nnk) {
  const double fxh = 0.0;
  if (nnk == 0) return fxh;
  double nnh;
  if (nnk < n) {
    scratch.assign(di, di + n);
    std::nth_element(scratch.begin(), scratch.begin() + (nnk - 1), scratch.end());
    nnh = scratch[nnk - 1];
  } else {
    nnh = 0;
    for (int i = 0; i < n; i++) nnh = std::max(nnh, di[i]);
    nnh = nnh * std::exp(std::log(1.0 * nnk / n) / 2);
  }
  return std::max(fxh, nnh);
}

inline std::uint64_t mid_key(int a, int b) {
  if (a > b) std::swap(a, b);
  return ((std::uint64_t)(std::uint32_t)a << 32) | (std::uint32_t)b;
}

// Local-constant Gaussian fit at a vertex with bandwidth h: reginit() + one max_nr() Newton step with
// JAC_EIGD scaling (locfit.c:120-239, m_max.c, m_jacob.c:31-80, m_eigen.c:71-100). di holds the
// distances from the vertex to the data points; w is scratch (n values).
double lc_fit(const double* y, int n, const double* di, double h, double* w) {
  double s0 = 0.0, s1 = 0.0;
  for (int i = 0; i < n; i++) {
    w[i] = wgauss(di[i], h);
    if (w[i] > 0) {  // only w > 0 enter the design (lf_nbhd.c:266-273)
      s1 += w[i] * (1.0 * y[i]);
      s0 += w[i] * 1.0;
    }
  }
  if (s0 == 0) return std::numeric_limits<double>::quiet_NaN();  // cannot happen for k >= 1
  double m0 = (s1 - 0.0) / s0;
  double f1 = 0.0, Z = 0.0;
  for (int i = 0; i < n; i++) {
    if (w[i] > 0) {
      double z = y[i] - m0;
      f1 += 1.0 * w[i] * (1.0 * z);
      Z += (w[i] * 1.0) * 1.0 * 1.0;
    }
  }
  double dg = (Z <= 0) ? 0.0 : 1 / std::sqrt(Z);
  double Zs = Z * (dg * dg);
  double v = f1 * dg;
  double tol = 1.0e-8 * Zs;
  double wv = 0.0;
  wv += 1.0 * v;
  if (Zs > tol) wv /= Zs;
  double x = 0.0;
  x += 1.0 * wv;
  double delta = x * dg;
  return m0 + 1.0 * delta;  // coefficient updated before NR_BREAK (m_max.c:203-206)
}

// Parametric component for deg = 0 (pcomp.c:55-123): the unweighted mean plus a Newton step.
// y has n values with stride inc.
double par_comp(const double* y, std::size_t inc, int n) {
  double s0 = 0.0, s1 = 0.0;
  for (int i = 0; i < n; i++) {
    s1 += 1.0 * (1.0 * y[i * inc]);
    s0 += 1.0 * 1.0;
  }
  double m0 = s1 / s0;
  double f1 = 0.0, Z = 0.0;
  for (int i = 0; i < n; i++) {
    double z = y[i * inc] - m0;
    f1 += 1.0 * 1.0 * (1.0 * z);
    Z += (1.0 * 1.0) * 1.0 * 1.0;
  }
  double dg = 1 / std::sqrt(Z);
  double Zs = Z * (dg * dg);
  double v = f1 * dg;
  double wv = 0.0;
  wv += 1.0 * v;
  if (Zs > 1.0e-8 * Zs) wv /= Zs;
  double x = 0.0;
  x += 1.0 * wv;
  return m0 + 1.0 * (x * dg);
}

// Tree construction state (atree_start, atree_grow, newsplit).
struct TreeBuilder {
  SmoothTree& T;
  std::vector<double> di, scratch;
  bool overflow = false;

  explicit TreeBuilder(SmoothTree& tree) : T(tree), di(tree.n), scratch(tree.n) {}

  // The only part of a vertex fit that the tree needs: the bandwidth (procv -> nbhd).
  void vfun(int v) {
    const double px = T.xev[2 * v], py = T.xev[2 * v + 1];
    for (int i = 0; i < T.n; i++) di[i] = rho2(T.x1[i] - px, T.x2[i] - py);
    T.h[v] = compbandwid(di.data(), scratch, T.n, T.nnk);
  }

  int add_vertex(double px, double py) {
    T.xev.push_back(px);
    T.xev.push_back(py);
    T.h.push_back(0.0);
    T.s.push_back(0);
    T.lo.push_back(0);
    T.hi.push_back(0);
    return T.nv++;
  }

  // newsplit() (ev_main.c:191-225).
  int newsplit(int i0, int i1, int pv) {
    int i = T.findpt(i0, i1);
    if (i >= 0) return i;
    if (i0 > i1) std::swap(i0, i1);
    if (T.nv == T.nvm) {  // ERROR(("newsplit: out of vertex space"))
      overflow = true;
      return -1;
    }
    const int v = add_vertex((T.xev[2 * i0] + T.xev[2 * i1]) / 2,
                             (T.xev[2 * i0 + 1] + T.xev[2 * i1 + 1]) / 2);
    T.lo[v] = i0;
    T.hi[v] = i1;
    if (pv) {  // pseudo-vertex
      T.h[v] = (T.h[i0] + T.h[i1]) / 2;
      T.s[v] = 1;
    } else {
      vfun(v);
      T.s[v] = 0;
    }
    T.mid[mid_key(i0, i1)] = v;
    return v;
  }

  // A cell of the tree: its corner vertices, its bounds, and its depth (the root cell has depth 1).
  struct Cell {
    int ce[4];
    double ll[2], ur[2];
    int depth;
  };

  // atree_grow() (ev_atree.c:83-124) without recursion. locfit's atree_grow() splits a cell, calls
  // itself for the lower half and then for the upper half. This loop visits the cells in the same
  // depth-first order from an explicit stack on the heap (the upper half is pushed first, so the
  // lower half and everything below it is finished before the upper half is split). The vertices are
  // therefore created in the same order and get the same numbers, and an out-of-vertex-space failure
  // happens at the same vertex. The half-cell bounds are computed with locfit's expressions.
  //
  // Why not recursion: the depth is not bounded by the number of points. With exactly duplicated
  // coordinates and k no larger than the number of copies of a point, the cells around that point
  // keep splitting down to floating-point resolution and then without end (each split adds vertices
  // on top of existing ones) until the vertex capacity is used up. That is about nvm / 2 levels
  // (23,044 at N = 1000, delta = 0.002). locfit itself dies there with a C stack overflow (its frames
  // hold about 128 KB); recursion here overflowed R's 8 MB stack from N = 2000 and a 512 KB thread
  // stack from N = 200. This loop needs about 56 bytes of heap per pending level and ends with
  // SMOOTH_OVERFLOW.
  void grow(const int* ce0, const double* ll0, const double* ur0) {
    std::vector<Cell> stack;
    stack.reserve(64);
    Cell root;
    for (int i = 0; i < 4; i++) root.ce[i] = ce0[i];
    for (int k = 0; k < 2; k++) {
      root.ll[k] = ll0[k];
      root.ur[k] = ur0[k];
    }
    root.depth = 1;
    stack.push_back(root);
    while (!stack.empty()) {
      const Cell c = stack.back();
      stack.pop_back();
      if (c.depth > T.depth) T.depth = c.depth;
      double le[2];
      const int ns = T.split(c.ce, le, c.ll, c.ur);
      if (ns == -1) continue;
      const int tk = 1 << ns;
      int nce[4];
      for (int i = 0; i < 4; i++) {
        if ((i & tk) == 0) {
          nce[i] = c.ce[i];
        } else {
          const int i0 = c.ce[i], i1 = c.ce[i - tk];
          const int pv = (le[ns] < (T.cut * std::min(T.h[i0], T.h[i1])));
          nce[i] = newsplit(i0, i1, pv);
          if (overflow) return;
        }
      }
      // upper half: z = ll[ns]; ll[ns] = (z + ur[ns]) / 2; corners: the new vertices, then the upper ones
      Cell up = c;
      up.depth = c.depth + 1;
      up.ll[ns] = (c.ll[ns] + c.ur[ns]) / 2;
      for (int i = 0; i < 4; i++) up.ce[i] = ((i & tk) == 0) ? nce[i + tk] : c.ce[i];
      // lower half: z = ur[ns]; ur[ns] = (z + ll[ns]) / 2; corners: the lower ones, then the new vertices
      Cell lo = c;
      lo.depth = c.depth + 1;
      lo.ur[ns] = (c.ur[ns] + c.ll[ns]) / 2;
      for (int i = 0; i < 4; i++) lo.ce[i] = nce[i];
      stack.push_back(up);
      stack.push_back(lo);
    }
  }
};

// A linear combination of vertex values (vertex id, weight), for the operator.
typedef std::vector<std::pair<int, double>> Combo;

inline Combo combo_avg(const Combo& a, const Combo& b) {
  Combo r;
  r.reserve(a.size() + b.size());
  for (const auto& e : a) r.push_back(std::make_pair(e.first, e.second / 2));
  for (const auto& e : b) r.push_back(std::make_pair(e.first, e.second / 2));
  return r;
}

// linear_interp(h, d, f0, f1) = ((d - h) f0 + h f1) / d (ev_interp.c:8-12) as weights.
inline Combo combo_lin(double hh, double dd, const Combo& f0, const Combo& f1) {
  if (dd == 0) return f0;
  Combo r;
  r.reserve(f0.size() + f1.size());
  for (const auto& e : f0) r.push_back(std::make_pair(e.first, e.second * (dd - hh) / dd));
  for (const auto& e : f1) r.push_back(std::make_pair(e.first, e.second * hh / dd));
  return r;
}

// atree_int() (ev_atree.c:162-205) with symbolic vertex values. Returns false if the descent fails.
bool descend_combo(const SmoothTree& T, double x0, double x1, Combo& out) {
  const double x[2] = {x0, x1};
  Combo vv[4];
  int ce[4];
  for (int i = 0; i < 4; i++) {
    vv[i].assign(1, std::make_pair(i, 1.0));
    ce[i] = i;
  }
  double le[2];
  int ns = 0;
  while (ns != -1) {
    const double* ll = &T.xev[2 * ce[0]];
    const double* ur = &T.xev[2 * ce[3]];
    ns = T.split(ce, le, ll, ur);
    if (ns != -1) {
      const int tk = 1 << ns;
      const double hh = ur[ns] - ll[ns];
      const int lo = (2 * (x[ns] - ll[ns])) < hh;
      for (int i = 0; i < 4; i++) {
        if ((tk & i) != 0) continue;
        const int nv = T.findpt(ce[i], ce[i + tk]);
        if (nv == -1) return false;  // ERROR(("Descend tree problem"))
        Combo nvv;
        if (T.s[nv]) nvv = combo_avg(vv[i], vv[i + tk]);  // exvvalpv(), nc = 1
        else nvv.assign(1, std::make_pair(nv, 1.0));
        if (lo) {
          ce[i + tk] = nv;
          vv[i + tk] = nvv;
        } else {
          ce[i] = nv;
          vv[i] = nvv;
        }
      }
    }
  }
  const double* ll = &T.xev[2 * ce[0]];
  const double* ur = &T.xev[2 * ce[3]];
  // rectcell_interp(), nc == 1 (ev_interp.c:64-80): coordinate 2 first, then coordinate 1.
  for (int i = 1; i >= 0; i--) {
    const int tk = 1 << i;
    for (int j = 0; j < tk; j++) vv[j] = combo_lin(x[i] - ll[i], ur[i] - ll[i], vv[j], vv[j + tk]);
  }
  out.swap(vv[0]);
  return true;
}

// atree_int() (ev_atree.c:162-205) plus rectcell_interp() for nc == 1 (ev_interp.c:64-80), with the
// vertex coefficients coef. Returns false if the descent fails.
bool descend_value(const SmoothTree& T, const std::vector<double>& coef, double x0, double x1,
                   double* value) {
  const double x[2] = {x0, x1};
  double vv[4];
  int ce[4];
  for (int i = 0; i < 4; i++) {
    vv[i] = coef[i];
    ce[i] = i;
  }
  double le[2];
  int ns = 0;
  while (ns != -1) {
    const double* ll = &T.xev[2 * ce[0]];
    const double* ur = &T.xev[2 * ce[3]];
    ns = T.split(ce, le, ll, ur);
    if (ns != -1) {
      const int tk = 1 << ns;
      const double hh = ur[ns] - ll[ns];
      const int lo = (2 * (x[ns] - ll[ns])) < hh;
      for (int i = 0; i < 4; i++) {
        if ((tk & i) != 0) continue;
        const int nv = T.findpt(ce[i], ce[i + tk]);
        if (nv == -1) return false;
        if (lo) {
          ce[i + tk] = nv;
          vv[i + tk] = T.s[nv] ? (vv[i] + vv[i + tk]) / 2 : coef[nv];
        } else {
          ce[i] = nv;
          vv[i] = T.s[nv] ? (vv[i] + vv[i + tk]) / 2 : coef[nv];
        }
      }
    }
  }
  const double* ll = &T.xev[2 * ce[0]];
  const double* ur = &T.xev[2 * ce[3]];
  for (int i = 0; i < 4; i++) {
    if (vv[i] == NOSLN) {  // locfit's "no solution" sentinel (ev_interp.c:70)
      *value = NOSLN;
      return true;
    }
  }
  for (int i = 1; i >= 0; i--) {
    const int tk = 1 << i;
    for (int j = 0; j < tk; j++) {
      const double hh = x[i] - ll[i], dd = ur[i] - ll[i];
      vv[j] = (dd == 0) ? vv[j] : (((dd - hh) * vv[j] + hh * vv[j + tk]) / dd);
    }
  }
  *value = vv[0];
  return true;
}

}  // namespace

const char* smooth_status_message(int status) {
  switch (status) {
    case SMOOTH_OK:
      return "ok";
    case SMOOTH_BAD_INPUT:
      return "invalid input: the smoother needs at least 2 points with finite coordinates";
    case SMOOTH_BAD_DELTA:
      return "invalid delta: it must be a positive finite number with N * delta below 2^31";
    case SMOOTH_K_TOO_SMALL:
      return "N * delta < 2: the smoothing neighbourhood has fewer than 2 points";
    case SMOOTH_OVERFLOW:
      return "newsplit: out of vertex space (the locfit tree needs more vertices than maxk allows; "
             "exactly duplicated coordinates do this when N * delta is at most the number of copies "
             "of a point)";
    case SMOOTH_INTERNAL:
      return "internal error in the smoothing tree";
  }
  return "unknown status";
}

int smooth_nn_k(int n, double delta) {
  if (!(delta > 0) || !std::isfinite(delta)) return -1;
  const double kd = n * delta + 1e-12;  // locfit.c:355: (int)(lfd->n * nn(sp) + 1e-12)
  if (!(kd < (double)INT_MAX)) return -1;
  return (int)kd;
}

int smooth_vertex_capacity(double delta, double cut, int maxk) {
  // atree_guessnv() (ev_atree.c:16-48), d = 2. Computed in double and capped at 2^30 (locfit's own
  // "unlimited" value) so that a tiny delta cannot overflow an int.
  const double cap = (double)(1 << 30);
  const int vc = 4;
  double nvm = cap;
  if (delta > 0) {
    const double a0 = (delta > 1) ? 1 : 1 / delta;
    if (cut < 0.01) cut = 0.01;
    double cu = 1;
    for (int i = 0; i < 2; i++) cu *= std::min(1.0, cut);
    const double nv = std::floor((5 * a0 / cu + 1) * vc);
    if (nv < nvm) nvm = nv;
  }
  if (nvm >= cap) nvm = 102 * vc;
  const double ifl = maxk / 100.0;
  const double out = std::floor(ifl * nvm);
  if (!(out < cap)) return 1 << 30;
  return (int)out;
}

int SmoothTree::findpt(int i0, int i1) const {
  auto it = mid.find(mid_key(i0, i1));
  return (it == mid.end()) ? -1 : it->second;
}

int SmoothTree::split(const int* ce, double* le, const double* ll, const double* ur) const {
  double hmin = 0.0, score[2];
  for (int i = 0; i < 4; i++) {
    const double hh = h[ce[i]];
    if ((hh > 0) && ((hmin == 0) | (hh < hmin))) hmin = hh;
  }
  int is = 0;
  for (int i = 0; i < 2; i++) {
    le[i] = (ur[i] - ll[i]) / 1.0;
    if (hmin == 0) score[i] = 2 * (ur[i] - ll[i]) / (fl[i + 2] - fl[i]);
    else score[i] = le[i] / hmin;
    if (score[i] > score[is]) is = i;
  }
  if (cut < score[is]) return is;
  return -1;
}

int smooth_tree_build(const double* x1, const double* x2, int n, double delta, int maxk,
                      SmoothTree& T) {
  T = SmoothTree();
  T.status = SMOOTH_BAD_INPUT;
  if (n < 2 || maxk < 1) return T.status;
  for (int i = 0; i < n; i++) {
    if (!std::isfinite(x1[i]) || !std::isfinite(x2[i])) return T.status;
  }
  T.n = n;
  T.delta = delta;
  T.maxk = maxk;
  T.cut = 0.8;
  if (!(delta > 0) || !std::isfinite(delta)) return T.status = SMOOTH_BAD_DELTA;
  T.nnk = smooth_nn_k(n, delta);
  if (T.nnk < 0) return T.status = SMOOTH_BAD_DELTA;
  if (T.nnk < 2) return T.status = SMOOTH_K_TOO_SMALL;
  T.nvm = smooth_vertex_capacity(delta, T.cut, maxk);
  T.x1.assign(x1, x1 + n);
  T.x2.assign(x2, x2 + n);

  // Bounding box of the data (set_flim, startlf.c:60-88).
  for (int k = 0; k < 2; k++) {
    const std::vector<double>& X = (k == 0) ? T.x1 : T.x2;
    double mx = X[0], mn = X[0];
    for (int j = 1; j < n; j++) {
      mx = std::max(mx, X[j]);
      mn = std::min(mn, X[j]);
    }
    T.fl[k] = mn;
    T.fl[k + 2] = mx;
  }

  // atree_start() (ev_atree.c:126-160): the 4 corners, then grow.
  if (T.nvm < 4) return T.status = SMOOTH_OVERFLOW;
  const std::size_t reserve = (std::size_t)std::min(T.nvm, 4096);
  T.xev.reserve(2 * reserve);
  T.h.reserve(reserve);
  T.s.reserve(reserve);
  T.lo.reserve(reserve);
  T.hi.reserve(reserve);
  TreeBuilder B(T);
  double ll[2] = {T.fl[0], T.fl[1]}, ur[2] = {T.fl[2], T.fl[3]};
  int ce[4];
  for (int i = 0; i < 4; i++) {
    int j = i;
    double px[2];
    for (int k = 0; k < 2; ++k) {
      px[k] = (j % 2) ? ur[k] : ll[k];
      j >>= 1;
    }
    const int v = B.add_vertex(px[0], px[1]);
    ce[i] = v;
    B.vfun(v);
    T.s[v] = 0;
  }
  B.grow(ce, ll, ur);
  if (B.overflow) return T.status = SMOOTH_OVERFLOW;
  return T.status = SMOOTH_OK;
}

double smooth_centre(const double* y, std::size_t inc, int n) {
  if (n < 1) return 0.0;
  return par_comp(y, inc, n);
}

int smooth_fitted_exact(const SmoothTree& T, const double* y, double* fitted) {
  if (T.status != SMOOTH_OK) return T.status;
  const int n = T.n;
  const double pc = par_comp(y, 1, n);
  std::vector<double> coef(T.nv, 0.0), di(n), w(n);
  for (int v = 0; v < T.nv; v++) {
    if (T.s[v]) continue;
    const double px = T.xev[2 * v], py = T.xev[2 * v + 1];
    for (int i = 0; i < n; i++) di[i] = rho2(T.x1[i] - px, T.x2[i] - py);
    const double c = lc_fit(y, n, di.data(), T.h[v], w.data());
    coef[v] = c - pc * 1.0;  // subparcomp() (pcomp.c:139)
  }
  for (int i = 0; i < n; i++) {
    double th;
    if (!descend_value(T, coef, T.x1[i], T.x2[i], &th)) return SMOOTH_INTERNAL;
    th += pc * 1.0;  // addparcomp() (pcomp.c:196)
    th += 0.0;       // base (fitted.c:82)
    fitted[i] = th;
  }
  return SMOOTH_OK;
}

int smooth_operator_build(const SmoothTree& T, SmoothOperator& op) {
  op = SmoothOperator();
  if (T.status != SMOOTH_OK) return T.status;
  const int n = T.n;
  std::vector<int> col(T.nv, -1);
  int m = 0;
  for (int v = 0; v < T.nv; v++) {
    if (!T.s[v]) col[v] = m++;
  }
  op.n = n;
  op.m = m;

  // Wn: normalised Gaussian weights at the real vertices.
  op.Wn.assign((std::size_t)m * n, 0.0);
  std::vector<double> w(n);
  for (int v = 0; v < T.nv; v++) {
    if (T.s[v]) continue;
    const double px = T.xev[2 * v], py = T.xev[2 * v + 1];
    double s0 = 0.0;
    for (int i = 0; i < n; i++) {
      const double di = rho2(T.x1[i] - px, T.x2[i] - py);
      w[i] = wgauss(di, T.h[v]);
      s0 += w[i];
    }
    double* row = &op.Wn[(std::size_t)col[v] * n];
    for (int i = 0; i < n; i++) row[i] = w[i] / s0;
  }

  // M: interpolation weights of every data point. Duplicate vertices of a row are merged in the
  // order in which the descent produces them (the same additions as accumulating a dense matrix).
  op.mptr.assign(n + 1, 0);
  op.mcol.reserve((std::size_t)n * 4);
  op.mval.reserve((std::size_t)n * 4);
  Combo terms;
  std::vector<int> rc;
  std::vector<double> rv;
  for (int i = 0; i < n; i++) {
    if (!descend_combo(T, T.x1[i], T.x2[i], terms)) return SMOOTH_INTERNAL;
    rc.clear();
    rv.clear();
    for (const auto& e : terms) {
      const int c = col[e.first];
      if (c < 0) return SMOOTH_INTERNAL;  // a pseudo-vertex left in the expansion
      std::size_t q = 0;
      while (q < rc.size() && rc[q] != c) q++;
      if (q == rc.size()) {
        rc.push_back(c);
        rv.push_back(0.0);
      }
      rv[q] += e.second;
    }
    // ascending column order (insertion sort; a row has at most a handful of entries)
    for (std::size_t a = 1; a < rc.size(); a++) {
      const int cc = rc[a];
      const double vvv = rv[a];
      std::size_t b = a;
      while (b > 0 && rc[b - 1] > cc) {
        rc[b] = rc[b - 1];
        rv[b] = rv[b - 1];
        b--;
      }
      rc[b] = cc;
      rv[b] = vvv;
    }
    for (std::size_t q = 0; q < rc.size(); q++) {
      if (rv[q] == 0.0) continue;
      op.mcol.push_back(rc[q]);
      op.mval.push_back(rv[q]);
    }
    op.mptr[i + 1] = (int)op.mcol.size();
  }
  return SMOOTH_OK;
}

void smooth_project(const SmoothOperator& op, const double* Y, std::size_t ldy, int B,
                    const double* centre, double* Z, std::size_t ldz) {
  for (int k = 0; k < op.m; k++) {
    double* zk = Z + (std::size_t)k * ldz;
    for (int b = 0; b < B; b++) zk[b] = 0.0;
    const double* wk = &op.Wn[(std::size_t)k * op.n];
    for (int j = 0; j < op.n; j++) {
      const double wkj = wk[j];
      const double* yj = Y + (std::size_t)j * ldy;
      if (centre != nullptr) {
        for (int b = 0; b < B; b++) zk[b] += wkj * (yj[b] - centre[b]);
      } else {
        for (int b = 0; b < B; b++) zk[b] += wkj * yj[b];
      }
    }
  }
}

void smooth_interp(const SmoothOperator& op, const double* Z, std::size_t ldz, int B,
                   const double* centre, const int* rows, int nrows, double* out, std::size_t ldo) {
  for (int r = 0; r < nrows; r++) {
    const int i = rows ? rows[r] : r;
    double* o = out + (std::size_t)r * ldo;
    for (int b = 0; b < B; b++) o[b] = 0.0;
    for (int q = op.mptr[i]; q < op.mptr[i + 1]; q++) {
      const double a = op.mval[q];
      const double* zc = Z + (std::size_t)op.mcol[q] * ldz;
      for (int b = 0; b < B; b++) o[b] += a * zc[b];
    }
    if (centre != nullptr) {
      for (int b = 0; b < B; b++) o[b] += centre[b];  // as addparcomp() (pcomp.c:196)
    }
  }
}

}  // namespace stc
