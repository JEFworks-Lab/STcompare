// Test 4: RNG options for a threaded backend.
//  - dqrng's xoshiro256+ (and its 2-argument seed(seed, stream) = long_jump x stream)
//  - O(1) per-stream seeding by hashing (SplitMix64) the (seed, stream) pair
//  - successive long_jump pre-pass (B states computed once, single-threaded)
//  - boost::random::normal_distribution (ziggurat), as used by dqrng::dqrnorm
//  - per-permutation work item: Fisher-Yates permutation of N + 9*N normals,
//    parallel over permutations with RcppParallel; output must not depend on
//    the number of threads.
// [[Rcpp::depends(dqrng, BH, sitmo, RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <dqrng_distribution.h>
#include <xoshiro.h>
#include <chrono>
#include <cstring>

static double now_sec() {
  using clk = std::chrono::steady_clock;
  return std::chrono::duration<double>(clk::now().time_since_epoch()).count();
}

static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// Lemire's nearly divisionless bounded integer in [0, range)
template <class E>
static inline uint32_t bounded(E& eng, uint32_t range) {
  uint64_t x = eng() >> 32;
  uint64_t m = x * (uint64_t)range;
  uint32_t l = (uint32_t)m;
  if (l < range) {
    uint32_t t = (uint32_t)(-range) % range;
    while (l < t) { x = eng() >> 32; m = x * (uint64_t)range; l = (uint32_t)m; }
  }
  return (uint32_t)(m >> 32);
}

// ---------------------------------------------------------------- micro-benchmarks
// [[Rcpp::export]]
Rcpp::List seeding_costs(double seed_d, int B) {
  uint64_t seed = (uint64_t)seed_d;
  Rcpp::List out;
  // (1) one long_jump
  {
    dqrng::xoshiro256plus e(seed);
    int reps = 2000;
    double t0 = now_sec();
    for (int r = 0; r < reps; ++r) e.long_jump();
    out["long_jump_us"] = 1e6 * (now_sec() - t0) / reps;
  }
  // (2) dqrng 2-arg seeding of streams 0..B-1 (each does `stream` long_jumps)
  {
    double t0 = now_sec(); uint64_t sink = 0;
    for (int s = 0; s < B; ++s) {
      auto g = dqrng::generator<dqrng::xoshiro256plus>(seed, (uint64_t)s);
      sink ^= (*g)();
    }
    out["dqrng_seed_stream_total_ms"] = 1e3 * (now_sec() - t0);
    out["sink1"] = (double)(sink & 0xffff);
  }
  // (3) successive long_jump pre-pass: B states, O(B) jumps total
  {
    double t0 = now_sec();
    dqrng::xoshiro256plus e(seed);
    std::vector<dqrng::xoshiro256plus> states;
    states.reserve(B);
    for (int s = 0; s < B; ++s) { states.push_back(e); e.long_jump(); }
    out["prepass_jump_total_ms"] = 1e3 * (now_sec() - t0);
  }
  // (4) hashed seeding (SplitMix64 inside xoshiro's seed()), O(1) per stream
  {
    double t0 = now_sec(); uint64_t sink = 0;
    for (int s = 0; s < B; ++s) {
      dqrng::xoshiro256plus e(splitmix64(seed ^ splitmix64((uint64_t)s + 1)));
      sink ^= e();
    }
    out["hash_seed_total_ms"] = 1e3 * (now_sec() - t0);
    out["sink2"] = (double)(sink & 0xffff);
  }
  return out;
}

// Single-thread throughput (values/second) of uniform and normal generation.
// [[Rcpp::export]]
Rcpp::NumericVector rng_throughput(double n_d) {
  size_t n = (size_t)n_d;
  Rcpp::NumericVector out(4);
  out.names() = Rcpp::CharacterVector::create(
      "raw_u64_per_s", "uniform01_per_s", "normal_boost_raw_engine_per_s",
      "normal_boost_dqrng_wrapper_per_s");
  dqrng::xoshiro256plus e(42);
  volatile double sink = 0;
  double t0 = now_sec(); uint64_t acc = 0;
  for (size_t i = 0; i < n; ++i) acc ^= e();
  out[0] = n / (now_sec() - t0); sink = (double)(acc & 1);
  t0 = now_sec(); double s = 0;
  for (size_t i = 0; i < n; ++i) s += dqrng::uniform01(e());
  out[1] = n / (now_sec() - t0); sink = s;
  boost::random::normal_distribution<double> nd(0.0, 1.0);
  t0 = now_sec(); s = 0;
  for (size_t i = 0; i < n; ++i) s += nd(e);
  out[2] = n / (now_sec() - t0); sink = s;
  dqrng::random_64bit_wrapper<dqrng::xoshiro256plus> w(42);
  dqrng::normal_distribution nd2(0.0, 1.0);
  t0 = now_sec(); s = 0;
  for (size_t i = 0; i < n; ++i) s += nd2(w);
  out[3] = n / (now_sec() - t0); sink = s;
  (void)sink;
  return out;
}

// ---------------------------------------------------------------- full work items
// Work item b: stream b -> permutation of 0..N-1 (Fisher-Yates) and nd*N normals.
// Output: perms (N x B, int, 1-based) and noise (nd*N x B), or only checksums.
struct PermNoiseWorker : public RcppParallel::Worker {
  uint64_t seed; int N, nd, mode; bool store;
  int* perm; double* noise; double* colsum;
  const std::vector<dqrng::xoshiro256plus>* states;
  PermNoiseWorker(uint64_t seed_, int N_, int nd_, int mode_, bool store_, int* perm_,
                  double* noise_, double* colsum_,
                  const std::vector<dqrng::xoshiro256plus>* states_)
      : seed(seed_), N(N_), nd(nd_), mode(mode_), store(store_), perm(perm_),
        noise(noise_), colsum(colsum_), states(states_) {}
  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> p(N);
    std::vector<double> buf(store ? 0 : (size_t)N * nd);
    boost::random::normal_distribution<double> nd01(0.0, 1.0);
    for (std::size_t b = begin; b < end; ++b) {
      dqrng::xoshiro256plus e = (mode == 0) ? (*states)[b]
          : dqrng::xoshiro256plus(splitmix64(seed ^ splitmix64((uint64_t)b + 1)));
      for (int i = 0; i < N; ++i) p[i] = i;
      for (int i = N - 1; i > 0; --i) std::swap(p[i], p[bounded(e, (uint32_t)(i + 1))]);
      double* z = store ? noise + b * (size_t)N * nd : buf.data();
      double s = 0;
      for (size_t k = 0; k < (size_t)N * nd; ++k) { z[k] = nd01(e); s += z[k]; }
      if (store) for (int i = 0; i < N; ++i) perm[b * (size_t)N + i] = p[i] + 1;
      colsum[b] = s + p[0] + p[N - 1];
    }
  }
};

// mode 0 = pre-pass of successive long_jump states; mode 1 = hashed seeds
// [[Rcpp::export]]
Rcpp::List perm_noise_cpp(double seed_d, int N, int B, int nd, int mode, bool store) {
  uint64_t seed = (uint64_t)seed_d;
  Rcpp::IntegerMatrix perm(store ? N : 0, store ? B : 0);
  Rcpp::NumericMatrix noise(store ? (size_t)N * nd : 0, store ? B : 0);
  Rcpp::NumericVector colsum(B);
  double t0 = now_sec();
  std::vector<dqrng::xoshiro256plus> states;
  if (mode == 0) {
    dqrng::xoshiro256plus e(seed);
    states.reserve(B);
    for (int s = 0; s < B; ++s) { states.push_back(e); e.long_jump(); }
  }
  double t1 = now_sec();
  PermNoiseWorker w(seed, N, nd, mode, store, store ? perm.begin() : nullptr,
                    store ? noise.begin() : nullptr, colsum.begin(), &states);
  RcppParallel::parallelFor(0, B, w, 1);
  double t2 = now_sec();
  return Rcpp::List::create(Rcpp::Named("perm") = perm, Rcpp::Named("noise") = noise,
                            Rcpp::Named("colsum") = colsum,
                            Rcpp::Named("prepass_s") = t1 - t0,
                            Rcpp::Named("parallel_s") = t2 - t1);
}
