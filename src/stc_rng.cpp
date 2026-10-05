// stc_rng.cpp -- R-compatible L'Ecuyer-CMRG stream and the independent streams (see stc_rng.h).
#include <cmath>
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

// ---------------------------------------------------------------------------------------------
// Independent streams
// ---------------------------------------------------------------------------------------------

namespace {
const std::uint64_t kGolden = 0x9e3779b97f4a7c15ULL;  // SplitMix64's increment
const std::uint64_t kTag = 0x5354636f6d706172ULL;     // "STcompar": separates these keys from other uses
const double kTwoPow52 = 4503599627370496.0;          // 2^52

inline std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
}  // namespace

std::uint64_t mix64(std::uint64_t z) {
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

std::uint64_t fnv1a64(const char* s, std::size_t n) {
  std::uint64_t h = 0xcbf29ce484222325ULL;
  for (std::size_t i = 0; i < n; i++) {
    h ^= (std::uint64_t)(unsigned char)s[i];
    h *= 0x100000001b3ULL;
  }
  return h;
}

std::uint64_t stream_key(std::uint32_t seed, const char* name, std::size_t name_len, int direction) {
  std::uint64_t k = mix64(kTag ^ (std::uint64_t)seed);
  k = mix64(k ^ fnv1a64(name, name_len));
  return mix64(k ^ (std::uint64_t)(std::uint32_t)direction);
}

Xoshiro256ss::Xoshiro256ss(std::uint64_t seed) {
  // SplitMix64 (Steele, Lea and Flood 2014): x += golden; output mix64(x)
  for (int j = 0; j < 4; j++) {
    seed += kGolden;
    s_[j] = mix64(seed);
  }
}

std::uint64_t Xoshiro256ss::next() {
  // xoshiro256** 1.0 (https://prng.di.unimi.it/xoshiro256starstar.c)
  const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
  const std::uint64_t t = s_[1] << 17;
  s_[2] ^= s_[0];
  s_[3] ^= s_[1];
  s_[1] ^= s_[2];
  s_[0] ^= s_[3];
  s_[2] ^= t;
  s_[3] = rotl(s_[3], 45);
  return result;
}

std::uint32_t Xoshiro256ss::bounded(std::uint32_t range) {
  // Lemire (2019), "Fast random integer generation in an interval", ACM TOMACS 29(1): the high 32 bits of
  // x * range for a 32-bit x, rejecting the low part below 2^32 mod range, so every value is equally likely.
  std::uint64_t m = (next() >> 32) * (std::uint64_t)range;
  std::uint32_t low = (std::uint32_t)m;
  if (low < range) {
    const std::uint32_t threshold = (std::uint32_t)(0u - range) % range;
    while (low < threshold) {
      m = (next() >> 32) * (std::uint64_t)range;
      low = (std::uint32_t)m;
    }
  }
  return (std::uint32_t)(m >> 32);
}

double Xoshiro256ss::norm() {
  if (has_spare_) {
    has_spare_ = false;
    return spare_;
  }
  double u, v, s;
  do {
    // uniforms in [-1, 1): k * 2^-52 - 1 for k = 0..2^53-1, exact in double
    u = (double)(next() >> 11) / kTwoPow52 - 1.0;
    v = (double)(next() >> 11) / kTwoPow52 - 1.0;
    s = u * u + v * v;
  } while (s >= 1.0 || s == 0.0);
  const double f = std::sqrt(-2.0 * std::log(s) / s);
  spare_ = v * f;
  has_spare_ = true;
  return u * f;
}

std::uint64_t stream_seed(std::uint64_t key, std::uint64_t b, std::uint64_t s) {
  return mix64(mix64(key ^ b) ^ s);
}

void stream_permutation(std::uint64_t key, int b, int N, int* perm) {
  Xoshiro256ss g(stream_seed(key, (std::uint64_t)b, 0));
  for (int i = 0; i < N; i++) perm[i] = i;
  for (int i = N - 1; i > 0; i--) {
    const int j = (int)g.bounded((std::uint32_t)i + 1u);
    const int t = perm[i];
    perm[i] = perm[j];
    perm[j] = t;
  }
}

void stream_noise(std::uint64_t key, int b, int k, int N, double* out) {
  Xoshiro256ss g(stream_seed(key, (std::uint64_t)b, (std::uint64_t)k + 1u));
  for (int i = 0; i < N; i++) out[i] = g.norm();
}

}  // namespace stc
