// stc_smoother.h -- the smoother of STcompare's matchingVariograms():
//
//   fit <- locfit::locfit(y ~ locfit::lp(x1, x2, nn = delta, deg = 0), kern = "gauss", maxk = 300)
//   fitted(fit)
//
// with locfit 1.5-9.12 semantics: the default adaptive "tree" evaluation structure (rbox(), cut 0.8),
// a local-constant Gaussian fit at the real vertices, pseudo-vertices that take the mean of their two
// parents, and bilinear interpolation to the data points. Ported from dev/prototypes/lfsmooth.cpp
// (identical() to fitted() on 24 coordinate/delta configurations); see
// dev/investigation/01-locfit-smoother.md (file:line citations refer to locfit's src/).
//
// Coordinates are in locfit's order: x1 = long = pos[, 2] and x2 = lat = pos[, 1], because
// STcompare calls lp(long, lat). The order decides the split dimension on ties and the order of the
// bilinear interpolation, so it matters to the last bit.
//
// Plain C++: no R API, safe on worker threads.
#ifndef STC_SMOOTHER_H
#define STC_SMOOTHER_H

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "stc_fp.h"

namespace stc {

enum SmoothStatus {
  SMOOTH_OK = 0,
  SMOOTH_BAD_INPUT = 1,   // fewer than 2 points, a non-finite coordinate, maxk < 1, or cut <= 0
  SMOOTH_BAD_DELTA = 2,   // delta <= 0, not finite, or so large that n * delta overflows an int
  SMOOTH_K_TOO_SMALL = 3, // k = (int)(n * delta + 1e-12) < 2: locfit errors for k = 0 and kills R
                          // (C stack overflow) for k = 1
  SMOOTH_OVERFLOW = 4,    // more vertices than atree_guessnv() allows: locfit's
                          // "newsplit: out of vertex space"
  SMOOTH_INTERNAL = 5     // tree descent failed (cannot happen for a tree built here)
};

// A short English description of a status, for warnings.
const char* smooth_status_message(int status);

// Nearest-neighbour count k = (int)(n * delta + 1e-12) (locfit.c:354-355), or -1 if delta is not a
// positive finite number or n * delta does not fit in an int.
int smooth_nn_k(int n, double delta);

// Vertex capacity of the tree, locfit's atree_guessnv() (ev_atree.c:16-48) for d = 2:
// floor(maxk / 100 * floor((5 * a0 / cut^2 + 1) * 4)) with a0 = 1 / delta (delta <= 1) or 1. It does
// not depend on the number of points. Values beyond 2^30 are capped at 2^30.
int smooth_vertex_capacity(double delta, double cut, int maxk);

// The adaptive tree (ev_atree.c). It depends only on the coordinates and delta, not on y.
struct SmoothTree {
  int status = SMOOTH_BAD_INPUT;
  int n = 0;
  std::vector<double> x1, x2;    // copies of the coordinates
  double delta = 0.0;
  int nnk = 0;                   // k
  double cut = 0.8;
  int maxk = 300;
  int nvm = 0;                   // vertex capacity
  double fl[4] = {0.0, 0.0, 0.0, 0.0};  // bounding box: min x1, min x2, max x1, max x2 (set_flim)
  int nv = 0;                    // vertices
  int depth = 0;                 // deepest cell visited (root = 1): the recursion depth locfit's
                                 // atree_grow() reaches on the same input, about 128 KB of C stack
                                 // per level in locfit; this build uses no recursion
  std::vector<double> xev;       // 2 * nv: vertex coordinates (x1, x2)
  std::vector<double> h;         // vertex bandwidth; pseudo-vertices: mean of their parents
  std::vector<int> s;            // 1 = pseudo-vertex (no local fit), 0 = real vertex
  std::vector<int> lo, hi;       // parents of a split vertex (lo < hi); 0 for the 4 corners
  std::unordered_map<std::uint64_t, int> mid;  // (lo, hi) -> vertex (findpt)

  // findpt() (ev_main.c:178-185): the vertex at the midpoint of i0 and i1, or -1.
  int findpt(int i0, int i1) const;
  // atree_split() (ev_atree.c:55-78): the split dimension of cell (ce, ll, ur), or -1.
  int split(const int* ce, double* le, const double* ll, const double* ur) const;
};

// Build the tree for n points (x1, x2) and nn = delta, with locfit's vertex capacity for maxk
// (STcompare: 300) and cut = 0.8. Returns T.status. Never aborts: every failure is a status. The
// cells are visited in locfit's depth-first order with an explicit heap stack, not by recursion, so
// the C stack use does not depend on the input (exactly duplicated coordinates with a small k make
// the tree refine without end until the capacity is exhausted: SMOOTH_OVERFLOW, where locfit
// overflows the C stack).
int smooth_tree_build(const double* x1, const double* x2, int n, double delta, int maxk,
                      SmoothTree& T);

// Exact replica of fitted(locfit(...)) for one response y (n values): the local fits at the real
// vertices (weighted mean plus one Newton step, minus the parametric component), the tree descent,
// pseudo-vertex averaging and bilinear interpolation, then the parametric component added back,
// in locfit's operation order. For tests; the engine uses the factored operator.
// Returns SMOOTH_OK, T.status if the tree is unusable, or SMOOTH_INTERNAL.
int smooth_fitted_exact(const SmoothTree& T, const double* y, double* fitted);

// The tree smoother as a linear operator: fitted = M (Wn (y - c)) + c, for a constant c near the
// centre of y (smooth_centre(): locfit's parametric component, which locfit itself subtracts at the
// vertices and adds back after interpolation; pcomp.c).
//   Wn  m x n, row-major (row v contiguous): the normalised Gaussian weights w_vi / sum_i w_vi at
//       real vertex v (m real vertices, in vertex order).
//   M   n x m, CSR: the bilinear interpolation weights of each data point, with pseudo-vertices
//       expanded into their real parents; column indices ascending within a row, exact zeros left out.
// The rows of Wn and M sum to 1 only to rounding, so M (Wn y) without the centring carries an error
// of about eps * |mean(y)|, which is large next to the spread of y when |mean(y)| >> spread (7e-8 of
// range(y) at y = 1e8 + N(0, 1)). With the centring the operator agrees with smooth_fitted_exact()
// and fitted(locfit()) to within 1e-14 * range(y) + 2 * eps * max|fitted|; the second term is the
// final rounding of the added c, which locfit's own result also carries. Measured (different
// summation order): below 3e-16 * range(y) for centred data and below 0.7 * eps * max|fitted| for
// data with large offsets (489 random configurations), and below 1.1e-14 * max|fitted| in all of
// 1,653 random configurations.
struct SmoothOperator {
  int n = 0;
  int m = 0;
  std::vector<double> Wn;     // m * n
  std::vector<int> mptr;      // n + 1
  std::vector<int> mcol;      // nnz
  std::vector<double> mval;   // nnz
};

// Returns SMOOTH_OK, T.status if the tree is unusable, or SMOOTH_INTERNAL.
int smooth_operator_build(const SmoothTree& T, SmoothOperator& op);

// The centring constant of a data vector y (n values with stride inc): locfit's parametric
// component for deg = 0 (pcomp.c), a two-pass mean, computed exactly as smooth_fitted_exact() does.
// Any c close to the centre of y gives the same fitted values up to rounding; this one is
// recommended. Returns 0 for n < 1.
double smooth_centre(const double* y, std::size_t inc, int n);

// Z = Wn (Y - centre) for a block of B columns: Z[k, b] = sum_j Wn[k, j] * (Y[j, b] - centre[b]).
//   Y  n x B, row-major (Y[j * ldy + b]); Z  m x B, row-major (Z[k * ldz + b]).
//   centre  B values, one per column (see smooth_centre()); pass the same array to smooth_interp().
//           nullptr means no subtraction here: either the caller has already subtracted the centres
//           (for example while gathering the columns, which saves the subtraction in this loop; it
//           costs about a third more time at B = 64), or there is no centring at all, which loses
//           accuracy when |mean| >> spread (see above).
// For each (k, b) the sum over j runs in ascending j, so a column's result does not depend on B.
void smooth_project(const SmoothOperator& op, const double* Y, std::size_t ldy, int B,
                    const double* centre, double* Z, std::size_t ldz);

// out = M[rows, ] Z + centre for a block of B columns (the centre is added after the sum, as
// locfit adds its parametric component back).
//   Z  m x B, row-major; out  nrows x B, row-major (out[r * ldo + b]).
//   centre  the B values subtracted from the columns before the projection (by smooth_project() or
//           by the caller), or nullptr if none were.
//   rows  the data points to produce (0-based), or nullptr for all n points in order (nrows = n).
// For each (row, b) the sum runs over the row's entries in ascending column order.
void smooth_interp(const SmoothOperator& op, const double* Z, std::size_t ldz, int B,
                   const double* centre, const int* rows, int nrows, double* out,
                   std::size_t ldo);

}  // namespace stc

#endif  // STC_SMOOTHER_H
