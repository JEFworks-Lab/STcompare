// stc_rng.h -- the random draws of the engine:
//   - legacy streams (the exported legacy functions): R's L'Ecuyer-CMRG generator and its Inversion
//     normals, reproduced exactly (below);
//   - independent streams (compareSpatial()): one stream per (task key, permutation, sub-stream),
//     generated from xoshiro256** (see "Independent streams" further down).
//
// The legacy code draws its noise inside BiocParallel tasks, where the RNG kind is L'Ecuyer-CMRG
// (dev/engine-spec.md, section 2.6): set.seed(seed + b), then rnorm(N) once per delta. This stream
// gives the same numbers as R:
//
//   RNGkind("L'Ecuyer-CMRG", "Inversion"); set.seed(seed); rnorm(n)
//
// following R's src/main/RNG.c (RNG_Init's initial scrambling and seed loop, unif_rand) and
// src/nmath/snorm.c (norm_rand, INVERSION, BIG = 2^27) in R 4.5.
//
// Thread safety: the generator is a plain value type, and nothing here touches R's RNG state or R
// memory. norm_rand() calls R's own Rf_qnorm5() so that the inversion is R's to the last bit. That is
// safe on worker threads. Its argument is p = u / 2^27 with u = (int)(2^27 * u1) + u2 and both
// uniforms in (0, 1), so p lies in (0, 1]: p > 0 because u >= u2 > 0, but p can be exactly 1, when
// (int)(2^27 * u1) = 2^27 - 1 and u2 > 1 - 2^-27, because (2^27 - 1) + u2 then rounds to 2^27
// (probability about 6e-17 per draw). For p in (0, 1) with mu = 0 and sigma = 1, qnorm5() only reads
// its arguments and calls log() and sqrt(). For p = 1 it returns +Inf from its boundary check
// (R_Q_P01_boundaries) before any computation, exactly as R's own norm_rand() does, so the stream
// still matches R. In both cases it uses no global state, allocates nothing, and none of its warning
// paths (ML_WARNING, only for p outside [0, 1] or invalid mu and sigma) can be reached.
//
// Why not a C++ copy of qnorm5's AS 241 polynomials: it matches R's qnorm() only when compiled with the
// floating-point contraction R itself was built with. On 2e7 arguments drawn like norm_rand()'s, a copy
// compiled without contraction differed from qnorm() in about half of them (by up to 1.1e-15 relative)
// both on CRAN's R for macOS arm64 (clang; -ffp-contract=on or fast matched) and on R for Linux arm64
// (GCC 13; only -ffp-contract=fast matched).
#ifndef STC_RNG_H
#define STC_RNG_H

#include <cstddef>
#include <cstdint>

#include "stc_fp.h"

namespace stc {

class LecuyerCMRG {
 public:
  LecuyerCMRG();
  explicit LecuyerCMRG(int seed);

  // set.seed(seed) under RNGkind("L'Ecuyer-CMRG"): RNG_Init() with Int32 = unsigned int arithmetic,
  // so a negative seed wraps modulo 2^32. The state then equals .Random.seed[2:7] (as unsigned).
  void seed(int seed);
  void set_state(const std::uint32_t* state);  // 6 values: 3 below m1, then 3 below m2
  void get_state(std::uint32_t* state) const;

  double unif_rand();                          // R's unif_rand(), in (0, 1)
  double norm_rand();                          // R's norm_rand() with normal.kind = "Inversion"
  void rnorm(double* out, std::size_t n);      // rnorm(n): n successive norm_rand() values

 private:
  std::uint32_t s_[6];
};

// The legacy noise of one permutation: block k (k = 0..K-1) holds the k-th rnorm(N) after
// set.seed(seed) under L'Ecuyer-CMRG, the noise matchingVariograms() adds for the k-th delta.
// out has N * K values, block k at out + k * N (an N x K column-major matrix).
void legacy_noise(int seed, int N, int K, double* out);

// ---------------------------------------------------------------------------------------------
// Independent streams (compareSpatial(); dev/compare-spatial-spec.md, section 3)
// ---------------------------------------------------------------------------------------------
//
// Every task (a gene and a direction: permute x, or permute y) has a 64-bit key
//
//   key = mix(mix(mix(T ^ seed) ^ fnv1a64(name)) ^ direction)
//
// where mix() is SplitMix64's finaliser (a bijection of 64-bit integers with full avalanche), T a
// fixed tag, seed the 32-bit seed (as unsigned), name the gene name's UTF-8 bytes and direction 1
// (x permuted) or 2 (y permuted). Permutation b (1-based) of the task draws from independent
// sub-streams s = 0, 1, ..., each a xoshiro256** generator (Blackman and Vigna 2018) whose state is
// four SplitMix64 outputs from the seed mix(mix(key ^ b) ^ s):
//   - s = 0: the permutation, Fisher-Yates (Durstenfeld) on 0..N-1, i = N-1 down to 1 swapping i with
//     an index drawn uniformly from 0..i by Lemire's unbiased bounded-integer method (32-bit draws:
//     the upper half of each 64-bit output);
//   - s = 1 + k: the noise added for position k of the delta grid, N standard normals by Marsaglia's
//     polar method (uniforms in [-1, 1) with 53-bit resolution; a pair per accepted point).
// So the draws of (task, b) depend on nothing but (seed, name, direction, b): not on the number of
// threads, the work items, the batches, the gene order or the other genes of the call. A block can be
// regenerated at any time (the engine draws block k* a second time for the full-length surrogate
// instead of keeping every block). The normals use std::log() and std::sqrt() of the platform's C
// library, so they can differ in the last bit between platforms (never between runs); everything else
// is integer arithmetic. Nothing here touches R.
std::uint64_t mix64(std::uint64_t z);
std::uint64_t fnv1a64(const char* s, std::size_t n);
std::uint64_t stream_key(std::uint32_t seed, const char* name, std::size_t name_len, int direction);

class Xoshiro256ss {
 public:
  explicit Xoshiro256ss(std::uint64_t seed);  // state: four SplitMix64 outputs from seed
  std::uint64_t next();
  std::uint32_t bounded(std::uint32_t range);  // uniform in 0..range-1 (Lemire 2019), range >= 1
  double norm();                               // standard normal (Marsaglia's polar method)

 private:
  std::uint64_t s_[4];
  double spare_ = 0.0;
  bool has_spare_ = false;
};

// The seed of sub-stream s of permutation b of the task with key `key`.
std::uint64_t stream_seed(std::uint64_t key, std::uint64_t b, std::uint64_t s);

// Permutation b of a task: perm (N values) receives a permutation of 0..N-1 (sub-stream 0).
void stream_permutation(std::uint64_t key, int b, int N, int* perm);

// The noise of permutation b for grid position k: N standard normals (sub-stream 1 + k).
void stream_noise(std::uint64_t key, int b, int k, int N, double* out);

}  // namespace stc

#endif  // STC_RNG_H
