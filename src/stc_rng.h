// stc_rng.h -- R's L'Ecuyer-CMRG generator and its Inversion normals, reproduced exactly.
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

}  // namespace stc

#endif  // STC_RNG_H
