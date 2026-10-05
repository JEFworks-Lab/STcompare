// Exact C++ core of STcompare's viladomatCorrelation() for one direction, given
//  - permuted data Xr (N x B) and noise E (N x K x B) generated in R with the package's RNG stream
//  - factorised smoothers S_k = W_k V_k (W_k: N x m_k interpolation weights, V_k: m_k x N vertex NW weights)
//  - precomputed variogram pair list (1-based i, j, bin) for the subsample ids
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

struct Pairs { std::vector<int> i, j, b; std::vector<double> inv2cnt; std::vector<int> keep; int nb; };

static void variog(const arma::mat& Z, const Pairs& P, arma::mat& out) {
  const int B = Z.n_cols, np = P.i.size();
  arma::mat acc(P.nb, B, arma::fill::zeros);
  for (int c = 0; c < B; c++) {
    const double* z = Z.colptr(c); double* a = acc.colptr(c);
    for (int p = 0; p < np; p++) { double d = z[P.i[p]] - z[P.j[p]]; a[P.b[p]] += d * d * 0.5; }
  }
  out.set_size(P.keep.size(), B);
  for (int c = 0; c < B; c++) for (size_t k = 0; k < P.keep.size(); k++) out(k, c) = acc(P.keep[k], c) * P.inv2cnt[P.keep[k]];
}

// [[Rcpp::export]]
List fastGeneCpp(const arma::mat& Xr, const arma::cube& E, List Ws, List Vs, IntegerVector ids1,
                 IntegerVector pi, IntegerVector pj, IntegerVector pb, IntegerVector cnt, LogicalVector keep,
                 const arma::vec& target) {
  const int N = Xr.n_rows, B = Xr.n_cols, K = Ws.size();
  arma::uvec ids = as<arma::uvec>(ids1) - 1;
  Pairs P; P.nb = cnt.size();
  for (int p = 0; p < pi.size(); p++) { P.i.push_back(pi[p] - 1); P.j.push_back(pj[p] - 1); P.b.push_back(pb[p] - 1); }
  for (int k = 0; k < P.nb; k++) { P.inv2cnt.push_back(cnt[k] > 0 ? 1.0 / cnt[k] : 0.0); if (keep[k]) P.keep.push_back(k); }
  const double tbar = arma::mean(target);
  arma::mat RSS(K, B), A(K, B), C(K, B);
  std::vector<arma::mat> VX(K);
  arma::mat G, GH;
  for (int k = 0; k < K; k++) {
    arma::mat W = as<arma::mat>(Ws[k]); arma::mat V = as<arma::mat>(Vs[k]);
    VX[k] = V * Xr;                                   // m x B
    arma::mat Xds = W.rows(ids) * VX[k];              // n_s x B
    variog(Xds, P, G);
    arma::rowvec gbar = arma::mean(G, 0);
    for (int c = 0; c < B; c++) {
      double sxy = 0, sxx = 0;
      for (arma::uword r = 0; r < G.n_rows; r++) { double gc = G(r, c) - gbar(c); sxy += gc * (target(r) - tbar); sxx += gc * gc; }
      double b1 = sxy / sxx, b0 = tbar - b1 * gbar(c);
      A(k, c) = std::sqrt(std::fabs(b1)); C(k, c) = std::sqrt(std::fabs(b0));
    }
    arma::mat H(Xds.n_rows, B);
    for (int c = 0; c < B; c++) for (arma::uword r = 0; r < ids.n_elem; r++) H(r, c) = Xds(r, c) * A(k, c) + E(ids(r), k, c) * C(k, c);
    variog(H, P, GH);
    for (int c = 0; c < B; c++) { double s = 0; for (arma::uword r = 0; r < GH.n_rows; r++) { double d = GH(r, c) - target(r); s += d * d; } RSS(k, c) = s; }
  }
  IntegerVector dstar(B);
  arma::mat Pm(N, B);
  for (int c = 0; c < B; c++) {
    arma::uword k; RSS.col(c).min(k); dstar[c] = k + 1;
    arma::mat W = as<arma::mat>(Ws[k]);
    arma::vec xd = W * VX[k].col(c);
    for (int r = 0; r < N; r++) Pm(r, c) = xd(r) * A(k, c) + E(r, k, c) * C(k, c);
  }
  return List::create(_["dstar"] = dstar, _["permutations"] = Pm, _["RSS"] = RSS);
}
