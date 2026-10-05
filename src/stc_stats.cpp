// stc_stats.cpp -- R-compatible correlation and closed-form least squares (see stc_stats.h).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "stc_stats.h"  // includes stc_fp.h: no FP contraction below

namespace stc {

namespace {

typedef long double LD;  // R's LDOUBLE (src/include/Defn.h)

inline double nan_value() { return std::numeric_limits<double>::quiet_NaN(); }

// Two-pass mean with accumulator type A (macro MEAN_ in cov.c with A = LDOUBLE).
template <typename A>
double mean_t(const double* x, std::size_t inc, int n) {
  A sum = 0.0;
  for (int k = 0; k < n; k++) sum += x[k * inc];
  A tmp = sum / n;
  if (std::isfinite((double)tmp)) {
    sum = 0.0;
    for (int k = 0; k < n; k++) sum += (x[k * inc] - tmp);
    tmp = tmp + sum / n;
  }
  return (double)tmp;
}

// sum += a * b as cov.c writes it, or fused (only used where A has the size of a double).
template <typename A>
inline A acc_prod(A sum, A a, A b, bool fused) {
  if (fused) return (A)std::fma((double)a, (double)b, (double)sum);
  return sum + a * b;
}

// The standard deviation of y and its mean (MEAN_ and COV_SDEV in cov_na_2()).
template <typename A>
void target_t(const double* y, int n, bool fused, CorTarget& t) {
  t.mean = mean_t<A>(y, 1, n);
  const A yym = t.mean;
  A sum = 0.0;
  for (int k = 0; k < n; k++) {
    const A d = y[k] - yym;
    sum = acc_prod<A>(sum, d, d, fused);
  }
  sum /= (A)(n - 1);
  t.sd = (double)std::sqrt(sum);
}

// cov_na_2() for one column: mean, centred cross products (iterations k >= first_fused fused), the
// standard deviation of x (fused if sd_fused), division and CLAMP.
template <typename A>
int cor_t(const double* x, std::size_t incx, const double* y, const CorTarget& t, int first_fused,
          bool sd_fused, double* r) {
  const int n = t.n;
  const A xxm = mean_t<A>(x, incx, n);
  const A yym = t.mean;
  A sum = 0.0;
  for (int k = 0; k < n; k++) {
    sum = acc_prod<A>(sum, (A)(x[k * incx] - xxm), (A)(y[k] - yym), k >= first_fused);
  }
  double ans = (double)(sum / (A)(n - 1));
  A ss = 0.0;
  for (int k = 0; k < n; k++) {
    const A d = x[k * incx] - xxm;
    ss = acc_prod<A>(ss, d, d, sd_fused);
  }
  ss /= (A)(n - 1);
  const double xsd = (double)std::sqrt(ss);
  if (xsd == 0. || t.sd == 0.) return COR_SD_ZERO;
  ans /= (xsd * t.sd);
  ans = (ans >= 1. ? 1. : (ans <= -1. ? -1. : ans));  // CLAMP
  *r = ans;
  return COR_OK;
}

// Fused accumulation exists only where long double is double.
inline bool fused_mode(int mode) {
  return sizeof(LD) == sizeof(double) && (mode == COR_FMA || mode == COR_FMA_TAIL8);
}

// Two-pass mean in double, for the least squares.
double mean2(const double* x, std::size_t inc, int K) { return mean_t<double>(x, inc, K); }

// ||x - c||, scaled by the largest |x - c| so that it neither overflows nor underflows (as
// LINPACK's dnrm2 does for the norms in dqrdc2).
double scaled_norm(const double* x, std::size_t inc, int K, double c) {
  double amax = 0.0;
  for (int k = 0; k < K; k++) amax = std::max(amax, std::fabs(x[k * inc] - c));
  if (amax == 0.0) return 0.0;
  double s = 0.0;
  for (int k = 0; k < K; k++) {
    const double t = (x[k * inc] - c) / amax;
    s += t * t;
  }
  return amax * std::sqrt(s);
}

}  // namespace

double r_mean(const double* x, std::size_t inc, int n) { return mean_t<LD>(x, inc, n); }

void cor_target_prepare(const double* y, int n, int mode, CorTarget& t) {
  t = CorTarget();
  t.n = n;
  t.mode = (mode >= COR_PLAIN && mode <= COR_DOUBLE) ? mode : COR_PLAIN;
  for (int k = 0; k < n; k++) {
    if (std::isnan(y[k])) t.has_na = true;
  }
  if (t.has_na || n < 2) return;
  if (t.mode == COR_DOUBLE) target_t<double>(y, n, false, t);
  else target_t<LD>(y, n, fused_mode(t.mode), t);
}

int cor_with_target(const double* x, std::size_t incx, const double* y, const CorTarget& t,
                    double* r) {
  *r = nan_value();
  const int n = t.n;
  if (n < 2) return COR_TOO_FEW;  // COV_n_le_1
  bool na = t.has_na;             // find_na_2()
  for (int k = 0; k < n && !na; k++) {
    if (std::isnan(x[k * incx])) na = true;
  }
  if (na) return COR_NA;
  if (t.mode == COR_DOUBLE) return cor_t<double>(x, incx, y, t, n, false, r);
  const bool fused = fused_mode(t.mode);
  int first_fused = n;
  if (fused) first_fused = (t.mode == COR_FMA_TAIL8) ? (n & ~7) : 0;
  return cor_t<LD>(x, incx, y, t, first_fused, fused, r);
}

int ols_fit(const double* x, std::size_t incx, const double* y, std::size_t incy, int K,
            double* b0, double* b1) {
  *b0 = nan_value();
  *b1 = nan_value();
  if (K < 1) return OLS_TOO_FEW;
  for (int k = 0; k < K; k++) {
    if (!std::isfinite(x[k * incx]) || !std::isfinite(y[k * incy])) return OLS_NONFINITE;
  }
  const double ybar = mean2(y, incy, K);
  if (K < 2) {
    *b0 = ybar;
    return OLS_TOO_FEW;
  }
  const double xbar = mean2(x, incx, K);
  // dqrdc2's rank decision for the slope column (tol = 1e-7, work(j, 2) = 1 when ||x|| = 0)
  const double xnorm = scaled_norm(x, incx, K, 0.0);
  const double cnorm = scaled_norm(x, incx, K, xbar);
  if (!(cnorm >= 1e-7 * (xnorm == 0.0 ? 1.0 : xnorm))) {
    *b0 = ybar;
    return OLS_RANK_DEFICIENT;
  }
  double sxx = 0.0, sxy = 0.0;
  for (int k = 0; k < K; k++) {
    const double dx = x[k * incx] - xbar;
    const double dy = y[k * incy] - ybar;
    sxx += dx * dx;
    sxy += dx * dy;
  }
  *b1 = sxy / sxx;
  *b0 = ybar - *b1 * xbar;
  return OLS_OK;
}

}  // namespace stc
