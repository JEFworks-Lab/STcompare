// stc_variog.h -- the binned empirical variogram of geoR 1.9-6,
//
//   geoR::variog(coords = cbind(long, lat), data = z, max.dist = max_dist, option = "bin",
//                messages = FALSE)            (every other argument at its default)
//
// split into a pair table that depends only on the coordinates and max_dist (built once) and an
// evaluator for many data vectors. Plain C++: no R API, safe on worker threads.
//
// The R side (R/engine.R, .stc_variog_plan()) computes what geoR computes with R's dist(): umax, the
// nugget flag and the bin limits (geoR's variogram.R:93-128). Building the table here then repeats
// geoR's binit() pair loop (src/geoR.c:370-422), and evaluation repeats its sums. See
// dev/engine-spec.md, sections 2.3 and 2.4, and dev/investigation/02-geoR-variogram.md.
#ifndef STC_VARIOG_H
#define STC_VARIOG_H

#include <cstddef>
#include <vector>

#include "stc_fp.h"

namespace stc {

enum VariogStatus {
  VARIOG_OK = 0,
  VARIOG_BAD_INPUT = 1  // fewer than 2 or more than 46340 points, bad limits, pairs_min < 1
};

// A view of a pair table (pointers only), so that tables stored in R vectors and in std::vectors can
// both be evaluated.
struct VariogView {
  int nbins = 0;               // number of output bins
  const int* ptr = nullptr;    // nbins + 1 offsets into pi/pj (CSR by output bin)
  const int* pi = nullptr;     // first point of each pair (0-based; pi > pj)
  const int* pj = nullptr;     // second point of each pair (0-based)
  const double* n = nullptr;   // pairs per output bin, as a double (geoR's $n)
};

struct VariogTable {
  int npoints = 0;
  int nbins = 0;                 // output bins
  std::vector<int> bin;          // geoR bin of each output bin: 0 = nugget bin [-1, 1e-12), b >= 1 regular bin b
  std::vector<int> ptr;          // CSR offsets (nbins + 1)
  std::vector<int> pi, pj;       // pairs of the output bins; within a bin in geoR's loop order
  std::vector<double> n;         // pairs per output bin
  std::vector<int> count_all;    // pairs per geoR bin, before the pairs_min and nugget filters
  double npairs_le_maxdist = 0;  // pairs with hypot(dx, dy) <= max_dist

  VariogView view() const {
    VariogView v;
    v.nbins = nbins;
    v.ptr = ptr.data();
    v.pi = pi.data();
    v.pj = pj.data();
    v.n = n.data();
    return v;
  }
};

// Build the pair table: geoR's binit() loop.
//   x1, x2       coordinates in geoR's column order (STcompare: x1 = long, x2 = lat), npoints each.
//   max_dist     the max.dist passed to geoR (STcompare: prctile); a pair is considered only if
//                hypot(dx, dy) <= max_dist.
//   lims         geoR's bins.lim after "bins.lim[1] <- -1": nlims ascending values, lims[0] < 0.
//                Pair (i, j) goes to bin b (0-based) when lims[b] <= d < lims[b + 1]; pairs at or
//                beyond the last limit (umax) are dropped.
//   keep_nugget  geoR's nt.ind: keep bin 0 (the nugget bin) in the output.
//   pairs_min    geoR's pairs.min (2): bins with fewer pairs are left out of the output.
// Pairs are visited with j outer ascending and i > j inner ascending; d = std::hypot(x1[i] - x1[j],
// x2[i] - x2[j]), the same libm hypot() that geoR calls.
int variog_build_table(const double* x1, const double* x2, int npoints, double max_dist,
                       const double* lims, int nlims, bool keep_nugget, int pairs_min,
                       VariogTable& out);

// Variograms of ncol data vectors at once.
//   z   npoints x ncol tile, row-major: value of column c at point m is z[m * ldz + c].
//   v   nbins x ncol, row-major: v[k * ldv + c].
// For every column and bin: acc = 0; for each pair in table order: t = z_i - z_j; t = t * t;
// t = t / 2; acc += t. Then v = acc / n. The pair loop is outside the loop over columns, so the
// inner loop vectorises, and each column's sum is accumulated in the same order as geoR's, so the
// result does not depend on ncol or on how columns are grouped into tiles. Does nothing when
// T.nbins or ncol is 0 (z and v are then not used and may be null).
void variog_eval_tile(const VariogView& T, const double* z, std::size_t ldz, int ncol, double* v,
                      std::size_t ldv);

}  // namespace stc

#endif  // STC_VARIOG_H
