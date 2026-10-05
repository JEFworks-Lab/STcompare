// Prototype C++ backend piece for STcompare (std::thread variant, no extra deps).
// Exported to R with dot-prefixed names so that NAMESPACE's
// exportPattern("^[[:alpha:]]+") does not export them.
#include <Rcpp.h>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <cmath>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace {

// Portable parallel-for over [0, n) with dynamic scheduling (atomic counter).
// Worker exceptions are captured and rethrown on the calling (R) thread.
// Workers must not touch the R/Rcpp API.
template <class F>
void parallel_for(int n, int nthreads, F&& f) {
  if (nthreads <= 1 || n <= 1) {
    for (int i = 0; i < n; ++i) f(i);
    return;
  }
  std::atomic<int> next(0);
  std::exception_ptr err = nullptr;
  std::mutex m;
  auto work = [&]() {
    try {
      for (int i; (i = next.fetch_add(1)) < n;) f(i);
    } catch (...) {
      std::lock_guard<std::mutex> g(m);
      if (!err) err = std::current_exception();
      next = n;  // stop handing out work
    }
  };
  std::vector<std::thread> th;
  int nt = std::min(nthreads, n);
  th.reserve(nt);
  for (int t = 0; t < nt; ++t) th.emplace_back(work);
  for (auto& x : th) x.join();
  if (err) std::rethrow_exception(err);
}

constexpr int W = 8;  // columns per tile: 8 doubles = one 64-byte cache line per point

}  // namespace

// Pairs (0-based i > j) within maxdist, binned with geoR::variog's "binit" rule.
// [[Rcpp::export(.variogram_pairs)]]
Rcpp::List variogram_pairs(Rcpp::NumericVector x, Rcpp::NumericVector y,
                           double maxdist, Rcpp::NumericVector lims) {
  int n = x.size(), nbins = lims.size() - 1;
  std::vector<int> pi, pj, pb;
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      double dx = x[i] - x[j], dy = y[i] - y[j];
      double d = std::sqrt(dx * dx + dy * dy);
      if (d <= maxdist) {
        int ind = 0;
        while (ind <= nbins && d >= lims[ind]) ind++;
        if (ind >= 1 && ind <= nbins && d < lims[ind]) {
          pi.push_back(i); pj.push_back(j); pb.push_back(ind - 1);
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("i") = pi, Rcpp::Named("j") = pj,
                            Rcpp::Named("bin") = pb, Rcpp::Named("nbins") = nbins);
}

// S[k, b] = sum over pairs in bin k of (Z[i, b] - Z[j, b])^2, for all columns b.
// [[Rcpp::export(.variogram_sums)]]
Rcpp::NumericMatrix variogram_sums(Rcpp::NumericMatrix Z, Rcpp::IntegerVector pi,
                                   Rcpp::IntegerVector pj, Rcpp::IntegerVector pb,
                                   int nbins, int nthreads = 1) {
  const int n = Z.nrow(), B = Z.ncol();
  const size_t P = pi.size();
  if ((size_t)pj.size() != P || (size_t)pb.size() != P) Rcpp::stop("pair vectors differ in length");
  for (size_t p = 0; p < P; ++p)  // validate on the R thread, before going parallel
    if (pi[p] < 0 || pi[p] >= n || pj[p] < 0 || pj[p] >= n || pb[p] < 0 || pb[p] >= nbins)
      Rcpp::stop("pair index out of range");
  Rcpp::NumericMatrix S(nbins, B);
  // raw pointers only inside workers (no R API)
  const double* zp = Z.begin();
  const int *ip = pi.begin(), *jp = pj.begin(), *bp = pb.begin();
  double* sp = S.begin();
  const int nblk = (B + W - 1) / W;
  parallel_for(nblk, nthreads, [=](int blk) {
    std::vector<double> T((size_t)n * W, 0.0), acc((size_t)nbins * W, 0.0);
    const int b0 = blk * W, wn = std::min(W, B - b0);
    for (int w = 0; w < wn; ++w) {
      const double* z = zp + (size_t)(b0 + w) * n;
      for (int i = 0; i < n; ++i) T[(size_t)i * W + w] = z[i];
    }
    for (size_t p = 0; p < P; ++p) {
      const double* ti = T.data() + (size_t)ip[p] * W;
      const double* tj = T.data() + (size_t)jp[p] * W;
      double* a = acc.data() + (size_t)bp[p] * W;
      for (int w = 0; w < W; ++w) {
        double v = ti[w] - tj[w];
        a[w] += v * v;
      }
    }
    for (int w = 0; w < wn; ++w)
      for (int k = 0; k < nbins; ++k)
        sp[(size_t)(b0 + w) * nbins + k] = acc[(size_t)k * W + w];
  });
  return S;
}

// Safety check: an exception thrown inside a worker must surface as an R error.
// [[Rcpp::export(.test_worker_exception)]]
int test_worker_exception(int nthreads) {
  std::atomic<int> count(0);
  parallel_for(100, nthreads, [&](int i) {
    if (i == 37) throw std::runtime_error("boom in worker");
    count++;
  });
  return count;
}
