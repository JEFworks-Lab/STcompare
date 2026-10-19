// Rcpp re-implementation of geoR 1.9-6
//   variog(coords = coords, data = z, max.dist = max.dist, option = "bin", messages = FALSE)
// (all other arguments at defaults) that precomputes the pair -> bin index once for fixed
// coordinates / max.dist and then evaluates many variograms (columns of an N x B matrix).
//
// Exactness notes (see report):
//  * geoR mixes THREE distance formulas:
//      - prctile  : R dist(cbind(lat, long))  -> sqrt(fma(dlong, dlong, dlat*dlat)) on CRAN arm64 R
//      - umax     : R dist(cbind(long, lat))  -> sqrt(fma(dlat, dlat, dlong*dlong))  (variogram.R:103,132)
//      - binning  : C hypot(dlong, dlat)                                             (geoR.c:386)
//    R's dist() is compiled with FP contraction on arm64 (verified), so we emulate it with std::fma
//    (rdist_fma = 1) or plain sqrt(a*a + b*b) (rdist_fma = 0, e.g. x86_64 builds without FMA).
//  * per-bin sums are accumulated in binit's loop order (j outer, i > j inner) in double.
//  * FP contraction is disabled in this file so that no fma sneaks into the accumulation or the
//    quantile interpolation; explicit std::fma is still honoured.

// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

using namespace Rcpp;

// Euclidean distance exactly as R's stats:::R_euclidean for a 2-column matrix whose columns are
// (c1, c2): dist = 0; dist += d1*d1; dist += d2*d2; sqrt(dist)   [the += may be fused]
static inline double rdist2(double d1, double d2, int use_fma) {
  if (use_fma) {
    double a = d1 * d1;              // fma(d1, d1, 0.0) == round(d1*d1)
    return std::sqrt(std::fma(d2, d2, a));
  } else {
    double a = d1 * d1;
    double b = d2 * d2;
    double s = a + b;
    return std::sqrt(s);
  }
}

// prctile <- quantile(dist(pos), probs = prob)  (type 7), with pos columns in the order used by
// STcompare (cbind(lat_s, long_s)).  Needs O(P) memory, P = N(N-1)/2.
// [[Rcpp::export]]
double rdist_quantile_cpp(NumericMatrix pos, double prob, int rdist_fma = 1) {
  const int n = pos.nrow();
  const double* c1 = &pos(0, 0);
  const double* c2 = &pos(0, 1);
  const R_xlen_t P = (R_xlen_t)n * (n - 1) / 2;
  std::vector<double> d;
  d.reserve(P);
  for (int j = 0; j < n; j++)
    for (int i = j + 1; i < n; i++)
      d.push_back(rdist2(c1[i] - c1[j], c2[i] - c2[j], rdist_fma));
  double index = 1.0 + std::max((double)P - 1.0, 0.0) * prob;   // quantile.default line 37
  double lo = std::floor(index), hi = std::ceil(index);
  R_xlen_t ilo = (R_xlen_t)lo - 1, ihi = (R_xlen_t)hi - 1;
  std::nth_element(d.begin(), d.begin() + ilo, d.end());
  double qlo = d[ilo];
  double qhi = qlo;
  if (ihi != ilo) qhi = *std::min_element(d.begin() + ilo + 1, d.end());
  double qs = qlo;
  if (index > lo && qhi != qlo) {                                // line 44-46
    double h = index - lo;
    double a = (1.0 - h) * qlo;
    double b = h * qhi;
    qs = a + b;
  }
  return qs;
}

// Precompute everything that depends only on coordinates and max.dist.
// coords: N x 2 in the column order passed to geoR::variog (STcompare: cbind(long, lat)).
// umax:   pass NA to compute it here with R-dist emulation, or pass R's value.
// [[Rcpp::export]]
List variog_prep_cpp(NumericMatrix coords, double max_dist, double umax = NA_REAL,
                     int rdist_fma = 1, int pairs_min = 2, int uvec = 13) {
  const int n = coords.nrow();
  const double* xc = &coords(0, 0);
  const double* yc = &coords(0, 1);
  const double nugget_tol = 1e-12;

  // (1) R side: min(u) and umax from R's dist(coords) (variogram.R:103-113,132)
  double min_u = R_PosInf, umax_c = R_NegInf;
  for (int j = 0; j < n; j++)
    for (int i = j + 1; i < n; i++) {
      double d = rdist2(xc[i] - xc[j], yc[i] - yc[j], rdist_fma);
      if (d < min_u) min_u = d;
      if (d < max_dist && d > umax_c) umax_c = d;
    }
  if (ISNAN(umax)) umax = umax_c;
  const bool nt_ind = (min_u < nugget_tol);

  // (2) .define.bins (variogram.R:1023-1028) + bins.lim[1] <- -1 (variogram.R:137)
  std::vector<double> bl(uvec + 1);
  const double step = (umax - 0.0) / (double)uvec;                 // del/n1 in seq.default
  bl[0] = 0.0;
  for (int k = 1; k < uvec; k++) bl[k] = 0.0 + (double)k * step;   // from + seq_len(.) * (del/n1)
  bl[uvec] = umax;                                                 // 'to' appended exactly
  std::vector<double> bins_lim;
  bins_lim.push_back(0.0);
  bins_lim.push_back(nugget_tol);
  for (int k = 0; k <= uvec; k++) if (bl[k] > nugget_tol) bins_lim.push_back(bl[k]);
  const int nbins = (int)bins_lim.size() - 1;
  std::vector<double> mid(nbins);
  for (int m = 0; m < nbins; m++) mid[m] = 0.5 * (bins_lim[m + 1] + bins_lim[m]);
  std::vector<double> lims(bins_lim);
  if (lims[0] < 1e-16) lims[0] = -1.0;

  // (3) C side: binit() pair loop with hypot (geoR.c:380-405); record (i, j, bin)
  std::vector<int> pi, pj, pb;
  pi.reserve((size_t)n * (n - 1) / 6); pj.reserve(pi.capacity()); pb.reserve(pi.capacity());
  R_xlen_t n_le_maxdist = 0;
  for (int j = 0; j < n; j++)
    for (int i = j + 1; i < n; i++) {
      double dx = xc[i] - xc[j];
      double dy = yc[i] - yc[j];
      double dist = std::hypot(dx, dy);
      if (dist <= max_dist) {
        n_le_maxdist++;
        int ind = 0;
        while (ind < nbins && dist >= lims[ind]) ind++;
        if (dist < lims[ind]) { pi.push_back(i); pj.push_back(j); pb.push_back(ind - 1); }
      }
    }
  std::vector<int> cnt(nbins, 0);
  for (int b : pb) cnt[b]++;

  // (4) which bins are returned: n >= pairs.min, nugget bin only if co-located data (variogram.R:151-168)
  std::vector<int> out_bin;  // indices into 0..nbins-1
  std::vector<int> map(nbins, -1);
  for (int b = 0; b < nbins; b++) {
    if (b == 0 && !nt_ind) continue;
    if (cnt[b] >= pairs_min) { map[b] = (int)out_bin.size(); out_bin.push_back(b); }
  }
  const int nout = (int)out_bin.size();

  // (5) CSR: pairs grouped by output bin, original loop order preserved inside each bin
  IntegerVector ptr(nout + 1);
  for (int k = 0; k < nout; k++) ptr[k + 1] = ptr[k] + cnt[out_bin[k]];
  IntegerVector ci(ptr[nout]), cj(ptr[nout]);
  std::vector<int> fill(ptr.begin(), ptr.end() - 1);
  for (size_t p = 0; p < pb.size(); p++) {
    int k = map[pb[p]];
    if (k < 0) continue;
    ci[fill[k]] = pi[p]; cj[fill[k]] = pj[p]; fill[k]++;
  }
  NumericVector u(nout), nn(nout);
  for (int k = 0; k < nout; k++) { u[k] = mid[out_bin[k]]; nn[k] = (double)cnt[out_bin[k]]; }
  NumericVector bl_out(bins_lim.begin() + (nt_ind ? 0 : 1), bins_lim.end());
  return List::create(_["i"] = ci, _["j"] = cj, _["ptr"] = ptr, _["u"] = u, _["n"] = nn,
                      _["bins.lim"] = bl_out, _["umax"] = umax, _["umax_cpp"] = umax_c,
                      _["max.dist"] = max_dist, _["nt.ind"] = nt_ind, _["n.points"] = n,
                      _["npairs_maxdist"] = (double)n_le_maxdist, _["npairs_used"] = (double)pb.size());
}

// One pass per column: v[k, b] = sum_{p in bin k} ((z_i - z_j)^2 / 2) / n_k
// [[Rcpp::export]]
NumericMatrix variog_compute_cpp(List prep, NumericMatrix Z) {
  IntegerVector ci = prep["i"], cj = prep["j"], ptr = prep["ptr"];
  NumericVector nn = prep["n"];
  const int nb = ptr.size() - 1, B = Z.ncol();
  const int* pi = ci.begin(); const int* pj = cj.begin(); const int* pp = ptr.begin();
  NumericMatrix V(nb, B);
  for (int b = 0; b < B; b++) {
    const double* z = &Z(0, b);
    for (int k = 0; k < nb; k++) {
      double acc = 0.0;
      for (int p = pp[k]; p < pp[k + 1]; p++) {
        double d = z[pi[p]] - z[pj[p]];
        double t = d * d;
        t = t / 2.0;
        acc += t;
      }
      V(k, b) = acc / nn[k];
    }
  }
  return V;
}

// Same result, but walks the pair list once for all columns: Zt is t(Z) (B x N, so each point's
// B values are contiguous) and the inner loop over columns vectorises.  Per-column summation
// order is unchanged, so results are identical to variog_compute_cpp.
// [[Rcpp::export]]
NumericMatrix variog_compute_t_cpp(List prep, NumericMatrix Zt) {
  IntegerVector ci = prep["i"], cj = prep["j"], ptr = prep["ptr"];
  NumericVector nn = prep["n"];
  const int nb = ptr.size() - 1, B = Zt.nrow();
  const int* pi = ci.begin(); const int* pj = cj.begin(); const int* pp = ptr.begin();
  const double* zt = Zt.begin();
  NumericMatrix V(nb, B);
  std::vector<double> acc(B);
  for (int k = 0; k < nb; k++) {
    std::fill(acc.begin(), acc.end(), 0.0);
    double* a = acc.data();
    for (int p = pp[k]; p < pp[k + 1]; p++) {
      const double* zi = zt + (R_xlen_t)pi[p] * B;
      const double* zj = zt + (R_xlen_t)pj[p] * B;
      for (int b = 0; b < B; b++) {
        double d = zi[b] - zj[b];
        double t = d * d;
        t = t / 2.0;
        a[b] += t;
      }
    }
    for (int b = 0; b < B; b++) V(k, b) = a[b] / nn[k];
  }
  return V;
}

// Single-column variant that hides FP-add latency: the per-bin accumulation chains are independent,
// so we advance all bins in lock-step (each bin still sums its own pairs in the original order, so the
// result is bit-identical to variog_compute_cpp).
// [[Rcpp::export]]
NumericMatrix variog_compute_il_cpp(List prep, NumericMatrix Z) {
  IntegerVector ci = prep["i"], cj = prep["j"], ptr = prep["ptr"];
  NumericVector nn = prep["n"];
  const int nb = ptr.size() - 1, B = Z.ncol();
  const int* pi = ci.begin(); const int* pj = cj.begin(); const int* pp = ptr.begin();
  std::vector<int> ord(nb), len(nb);
  for (int k = 0; k < nb; k++) { ord[k] = k; len[k] = pp[k + 1] - pp[k]; }
  std::sort(ord.begin(), ord.end(), [&](int a, int b) { return len[a] < len[b]; });
  NumericMatrix V(nb, B);
  std::vector<double> acc(nb);
  for (int b = 0; b < B; b++) {
    const double* z = &Z(0, b);
    std::fill(acc.begin(), acc.end(), 0.0);
    int s = 0;
    for (int q = 0; q < nb; q++) {
      const int end = len[ord[q]];
      for (; s < end; s++) {
        for (int r = q; r < nb; r++) {
          const int k = ord[r];
          const int p = pp[k] + s;
          double d = z[pi[p]] - z[pj[p]];
          double t = d * d;
          t = t / 2.0;
          acc[k] += t;
        }
      }
    }
    for (int k = 0; k < nb; k++) V(k, b) = acc[k] / nn[k];
  }
  return V;
}

// The three distance formulas, for diagnostics (pairs in dist() order).
// cols: (c1, c2) = coords passed to geoR::variog (STcompare: long, lat)
// [[Rcpp::export]]
NumericMatrix pair_distances_cpp(NumericMatrix coords) {
  const int n = coords.nrow();
  const double* x = &coords(0, 0);
  const double* y = &coords(0, 1);
  R_xlen_t P = (R_xlen_t)n * (n - 1) / 2, k = 0;
  NumericMatrix out(P, 4);
  for (int j = 0; j < n; j++)
    for (int i = j + 1; i < n; i++) {
      double dx = x[i] - x[j], dy = y[i] - y[j];
      out(k, 0) = std::hypot(dx, dy);       // binit
      out(k, 1) = rdist2(dx, dy, 1);        // R dist(cbind(c1, c2)), FMA build
      out(k, 2) = rdist2(dy, dx, 1);        // R dist(cbind(c2, c1)), FMA build
      out(k, 3) = rdist2(dx, dy, 0);        // plain sqrt(dx*dx + dy*dy)
      k++;
    }
  return out;
}
