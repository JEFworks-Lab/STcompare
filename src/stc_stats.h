// stc_stats.h -- R-compatible Pearson correlation and the closed-form least squares that replaces
// lm() in matchingVariograms(). Plain C++: no R API, safe on worker threads.
#ifndef STC_STATS_H
#define STC_STATS_H

#include <cstddef>

#include "stc_fp.h"

namespace stc {

// ---------------------------------------------------------------------------------------------
// Correlation: R's cor(x, y) for a column x and a vector y with use = "everything" and
// method = "pearson", i.e. cov_na_2() in R's src/library/stats/src/cov.c (R 4.5), applied to one
// column:
//   - means: two passes in long double (macro MEAN_), rounded to double;
//   - cross products and sums of squares of the centred values, accumulated in long double;
//     sd = sqrt(sum / (n - 1)) in long double, rounded to double;
//   - r = (cross / (n - 1)) / (sd_x * sd_y), clamped to [-1, 1].
// A NaN or NA anywhere in x or y gives NA (status 1); a zero standard deviation gives NA (status 2;
// R also warns); n < 2 gives NA (status 3). Infinite values propagate as in R (NaN).
//
// mode: how the running R build evaluates "sum += a * b" in those loops (R/engine.R, .stc_cor_mode(),
// finds the one that reproduces the running R exactly):
//   COR_PLAIN      long double accumulators, no contraction: R's code as written. Where long double
//                  is wider than double (x86-64: 80-bit x87; arm64 Linux: binary128) there is no
//                  fused operation, and this is R's arithmetic.
//   COR_FMA        long double = double, every accumulation fused: sum = fma(a, b, sum).
//   COR_FMA_TAIL8  CRAN's R 4.5 for macOS arm64 (Apple clang): the cross-product loop is vectorised
//                  in blocks of 8 with separate multiplies, followed by fused scalar iterations for
//                  the last n % 8 values; both sums of squares are fused (read from the stats.so
//                  disassembly, corcov(), and confirmed against cor()).
//   COR_DOUBLE     double accumulators, no contraction, no x87 code. This is R's own arithmetic when
//                  R is built with --disable-long-double (as in CRAN's "noLD" checks), and the
//                  fallback for platforms where the long double code cannot be trusted: Rosetta 2's
//                  x86-64 emulation (Docker on Apple silicon) mis-executes GCC's -O2 x87 code for
//                  COR_PLAIN (the same binary returns r = 1 there and the correct 0.9196 under QEMU;
//                  SSE2 double arithmetic is emulated correctly). Against an R with long double
//                  accumulators it is not exact: its rounding errors grow with n and with badly
//                  centred or badly scaled data, to 2.1e-15 on 2,280 correlations with n <= 5,000
//                  (offsets up to 1e8, scales 1e-150 to 1e150; Linux arm64 and Rosetta).
// Where long double is wider than double, COR_FMA and COR_FMA_TAIL8 behave as COR_PLAIN. Modes
// other than the exact one agree with cor() to about 3e-16 on the same inputs where long double is
// double (macOS arm64).
// ---------------------------------------------------------------------------------------------

enum CorStatus { COR_OK = 0, COR_NA = 1, COR_SD_ZERO = 2, COR_TOO_FEW = 3 };
enum CorMode { COR_PLAIN = 0, COR_FMA = 1, COR_FMA_TAIL8 = 2, COR_DOUBLE = 3 };

// R's two-pass mean (cov.c's MEAN_; also mean.default's C code): long double accumulators.
double r_mean(const double* x, std::size_t inc, int n);

// The parts of the correlation that depend only on y (computed once for many columns).
struct CorTarget {
  int n = 0;
  int mode = COR_PLAIN;
  bool has_na = false;  // ISNAN anywhere in y
  double mean = 0.0;    // ym (two-pass mean)
  double sd = 0.0;      // sqrt(sum((y - ym)^2) / (n - 1)), rounded to double
};

void cor_target_prepare(const double* y, int n, int mode, CorTarget& t);

// r = cor(x, y) in t.mode; x has t.n values with stride incx. Returns a CorStatus; *r is NaN unless
// COR_OK.
int cor_with_target(const double* x, std::size_t incx, const double* y, const CorTarget& t,
                    double* r);

// r[j] = cor(x, ys[j]) for nt targets, each bit-identical to cor_with_target(x, incx, ys[j], *ts[j],
// &r[j]) with the same status in status[j]: the same operations in the same order, except that the
// mean and the standard deviation of x are computed once for all targets (the within-sample engine
// correlates one surrogate with every other gene). All targets must have the same n and mode.
void cor_with_targets(const double* x, std::size_t incx, int nt, const double* const* ys,
                      const CorTarget* const* ts, double* r, int* status);

// ---------------------------------------------------------------------------------------------
// Least squares with intercept, the fit lm(y ~ 1 + x) makes in matchingVariograms() (y = target
// variogram, x = candidate variogram; dev/engine-spec.md, section 2.7):
//   xbar, ybar two-pass means; Sxx = sum((x - xbar)^2); Sxy = sum((x - xbar) (y - ybar));
//   b1 = Sxy / Sxx; b0 = ybar - b1 * xbar.
// Rank deficiency follows lm(): LINPACK dqrdc2 with tol = 1e-7 drops the slope column when its
// norm after removing the intercept, ||x - xbar||, is below 1e-7 times its original norm ||x||
// (1e-7 itself when ||x|| = 0). lm() then reports an NA slope and the intercept-only fit mean(y);
// so does this function (status OLS_RANK_DEFICIENT, *b1 = NaN). With one observation lm() reports
// y[1] and an NA slope (OLS_TOO_FEW); with none it fails (OLS_TOO_FEW, both NaN). A non-finite input
// gives OLS_NONFINITE and NaN for both.
// Agreement with lm(): on variogram pairs as matchingVariograms() fits them, slope within 2.4e-15 and
// intercept within 8.2e-14 relative. Both suffer cancellation (Sxy, and ybar - b1 * xbar), so the
// agreement follows the conditioning of the pair: on random ill-conditioned pairs (x offset by up to
// 100 times its spread, weak correlation) the slope differed by up to 2.0e-11 and the intercept by up
// to 2.6e-9 relative in 20,000 fits. The rank decisions (NA slope) matched lm() in all of them. lm()
// itself is not bit-stable across BLAS builds and memory alignment.
// ---------------------------------------------------------------------------------------------

enum OlsStatus { OLS_OK = 0, OLS_RANK_DEFICIENT = 1, OLS_TOO_FEW = 2, OLS_NONFINITE = 3 };

int ols_fit(const double* x, std::size_t incx, const double* y, std::size_t incy, int K,
            double* b0, double* b1);

}  // namespace stc

#endif  // STC_STATS_H
