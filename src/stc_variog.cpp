// stc_variog.cpp -- geoR::variog pair table and evaluator (see stc_variog.h).
// Ported from dev/prototypes/variog_fast.cpp, which is bit-identical to geoR 1.9-6.
#include <cmath>
#include <cstddef>
#include <vector>

#include "stc_variog.h"  // includes stc_fp.h: no FP contraction below

namespace stc {

int variog_build_table(const double* x1, const double* x2, int npoints, double max_dist,
                       const double* lims, int nlims, bool keep_nugget, int pairs_min,
                       VariogTable& out) {
  out = VariogTable();
  // 46340 * 46339 / 2 < 2^31: pair offsets fit in an int.
  if (npoints < 2 || npoints > 46340 || nlims < 2 || pairs_min < 1) return VARIOG_BAD_INPUT;
  // lims[0] < 0 guarantees that every distance passes the first limit, so the bin index below is
  // never -1 (geoR sets bins.lim[1] to -1 for this reason; variogram.R:127-128).
  if (!(lims[0] < 0.0)) return VARIOG_BAD_INPUT;
  for (int b = 1; b < nlims; b++) {
    if (!(lims[b] >= lims[b - 1])) return VARIOG_BAD_INPUT;  // also rejects NaN
  }
  const int nb = nlims - 1;
  out.count_all.assign(nb, 0);

  // geoR.c binit(): pair loop, j outer, i > j inner; distance with hypot(); [lo, hi) bins.
  std::vector<int> pi, pj, pb;
  double nle = 0.0;
  for (int j = 0; j < npoints; j++) {
    for (int i = j + 1; i < npoints; i++) {
      const double dx = x1[i] - x1[j];
      const double dy = x2[i] - x2[j];
      const double d = std::hypot(dx, dy);
      if (d <= max_dist) {
        nle += 1.0;
        int ind = 0;
        while (ind < nb && d >= lims[ind]) ind++;
        if (d < lims[ind]) {
          pi.push_back(i);
          pj.push_back(j);
          pb.push_back(ind - 1);
          out.count_all[ind - 1]++;
        }
      }
    }
  }

  // Output bins (variogram.R:141-157): n >= pairs.min, and the nugget bin only if nt.ind.
  std::vector<int> map(nb, -1);
  for (int b = 0; b < nb; b++) {
    if (b == 0 && !keep_nugget) continue;
    if (out.count_all[b] >= pairs_min) {
      map[b] = out.nbins++;
      out.bin.push_back(b);
    }
  }

  // CSR by output bin; within a bin the pairs keep the loop order (= geoR's summation order).
  out.ptr.assign(out.nbins + 1, 0);
  for (int k = 0; k < out.nbins; k++) out.ptr[k + 1] = out.ptr[k] + out.count_all[out.bin[k]];
  out.pi.resize(out.ptr[out.nbins]);
  out.pj.resize(out.ptr[out.nbins]);
  std::vector<int> fill(out.ptr.begin(), out.ptr.end() - 1);
  for (std::size_t p = 0; p < pb.size(); p++) {
    const int k = map[pb[p]];
    if (k < 0) continue;
    out.pi[fill[k]] = pi[p];
    out.pj[fill[k]] = pj[p];
    fill[k]++;
  }
  out.n.resize(out.nbins);
  for (int k = 0; k < out.nbins; k++) out.n[k] = (double)out.count_all[out.bin[k]];
  out.npoints = npoints;
  out.npairs_le_maxdist = nle;
  return VARIOG_OK;
}

namespace {

// TW columns: the pair loop runs once over the table, the inner loop over the TW columns
// vectorises. Each column keeps its own accumulator, summed in table order.
template <int TW>
void eval_fixed(const VariogView& T, const double* z, std::size_t ldz, double* v, std::size_t ldv) {
  for (int k = 0; k < T.nbins; k++) {
    double acc[TW];
    for (int c = 0; c < TW; c++) acc[c] = 0.0;
    for (int p = T.ptr[k]; p < T.ptr[k + 1]; p++) {
      const double* zi = z + (std::size_t)T.pi[p] * ldz;
      const double* zj = z + (std::size_t)T.pj[p] * ldz;
      for (int c = 0; c < TW; c++) {
        double t = zi[c] - zj[c];
        t = t * t;
        t = t / 2.0;
        acc[c] += t;
      }
    }
    const double nk = T.n[k];
    double* vk = v + (std::size_t)k * ldv;
    for (int c = 0; c < TW; c++) vk[c] = acc[c] / nk;
  }
}

// Fewer than 4 columns left.
void eval_small(const VariogView& T, const double* z, std::size_t ldz, int ncol, double* v,
                std::size_t ldv) {
  for (int k = 0; k < T.nbins; k++) {
    double acc[4] = {0.0, 0.0, 0.0, 0.0};
    for (int p = T.ptr[k]; p < T.ptr[k + 1]; p++) {
      const double* zi = z + (std::size_t)T.pi[p] * ldz;
      const double* zj = z + (std::size_t)T.pj[p] * ldz;
      for (int c = 0; c < ncol; c++) {
        double t = zi[c] - zj[c];
        t = t * t;
        t = t / 2.0;
        acc[c] += t;
      }
    }
    const double nk = T.n[k];
    double* vk = v + (std::size_t)k * ldv;
    for (int c = 0; c < ncol; c++) vk[c] = acc[c] / nk;
  }
}

}  // namespace

void variog_eval_tile(const VariogView& T, const double* z, std::size_t ldz, int ncol, double* v,
                      std::size_t ldv) {
  // Nothing to compute. Returning here also keeps the column offsets below (v + c) away from an
  // output buffer that is empty and may be a null pointer, which would be undefined behaviour.
  if (T.nbins <= 0 || ncol <= 0) return;
  int c = 0;
  for (; c + 16 <= ncol; c += 16) eval_fixed<16>(T, z + c, ldz, v + c, ldv);
  for (; c + 4 <= ncol; c += 4) eval_fixed<4>(T, z + c, ldz, v + c, ldv);
  if (c < ncol) eval_small(T, z + c, ldz, ncol - c, v + c, ldv);
}

}  // namespace stc
