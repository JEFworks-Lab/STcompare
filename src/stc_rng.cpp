// stc_rng.cpp -- R-compatible L'Ecuyer-CMRG stream (see stc_rng.h).
#include <cstddef>
#include <cstdint>

#include "stc_rng.h"  // includes stc_fp.h: no FP contraction below

// R's quantile function of the normal distribution (Rmath.h declares qnorm5 as Rf_qnorm5). Declared
// here instead of including Rmath.h, whose macros would leak into this file. See stc_rng.h for why
// calling it from a worker thread is safe.
extern "C" double Rf_qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);

namespace stc {

namespace {
// src/main/RNG.c
const std::int64_t m1 = 4294967087LL;
const std::int64_t m2 = 4294944443LL;
const double normc = 2.328306549295727688e-10;
const std::int64_t a12 = 1403580;
const std::int64_t a13n = 810728;
const std::int64_t a21 = 527612;
const std::int64_t a23n = 1370589;
}  // namespace

LecuyerCMRG::LecuyerCMRG() {
  for (int j = 0; j < 6; j++) s_[j] = 0;
}

LecuyerCMRG::LecuyerCMRG(int seed_value) { seed(seed_value); }

void LecuyerCMRG::seed(int seed_value) {
  // do_setseed(): RNG_Init(RNG_kind, (Int32) seed), Int32 = unsigned int.
  std::uint32_t seed = (std::uint32_t)seed_value;
  // RNG_Init(): initial scrambling ...
  for (int j = 0; j < 50; j++) seed = (69069u * seed + 1u);
  // ... and the LECUYER_CMRG seed loop: every seed below m2.
  for (int j = 0; j < 6; j++) {
    seed = (69069u * seed + 1u);
    while (seed >= m2) seed = (69069u * seed + 1u);
    s_[j] = seed;
  }
}

void LecuyerCMRG::set_state(const std::uint32_t* state) {
  for (int j = 0; j < 6; j++) s_[j] = state[j];
}

void LecuyerCMRG::get_state(std::uint32_t* state) const {
  for (int j = 0; j < 6; j++) state[j] = s_[j];
}

double LecuyerCMRG::unif_rand() {
  // unif_rand(), case LECUYER_CMRG (RNG.c). The seeds are read as unsigned int.
  std::int64_t p1 = a12 * (std::int64_t)s_[1] - a13n * (std::int64_t)s_[0];
  std::int64_t k = p1 / m1;  // R: k = (int)(p1 / m1); |p1 / m1| < 2^21
  p1 -= k * m1;
  if (p1 < 0) p1 += m1;
  s_[0] = s_[1];
  s_[1] = s_[2];
  s_[2] = (std::uint32_t)p1;

  std::int64_t p2 = a21 * (std::int64_t)s_[5] - a23n * (std::int64_t)s_[3];
  k = p2 / m2;
  p2 -= k * m2;
  if (p2 < 0) p2 += m2;
  s_[3] = s_[4];
  s_[4] = s_[5];
  s_[5] = (std::uint32_t)p2;

  return (double)((p1 > p2) ? (p1 - p2) : (p1 - p2 + m1)) * normc;
}

double LecuyerCMRG::norm_rand() {
  // norm_rand(), case INVERSION (snorm.c): two uniforms, because one is not precise enough.
  const double BIG = 134217728; /* 2^27 */
  double u1 = unif_rand();
  u1 = (int)(BIG * u1) + unif_rand();
  return Rf_qnorm5(u1 / BIG, 0.0, 1.0, 1, 0);
}

void LecuyerCMRG::rnorm(double* out, std::size_t n) {
  // rnorm(n, 0, 1) returns mu + sigma * norm_rand() = norm_rand() exactly (nmath/rnorm.c).
  for (std::size_t i = 0; i < n; i++) out[i] = norm_rand();
}

void legacy_noise(int seed, int N, int K, double* out) {
  LecuyerCMRG g(seed);
  for (int k = 0; k < K; k++) g.rnorm(out + (std::size_t)k * N, (std::size_t)N);
}

}  // namespace stc
