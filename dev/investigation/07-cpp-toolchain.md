# STcompare: C++ toolchain and portability investigation

Scope: find out which C++ tooling works on this machine and what will work for users, so we can design a C++ backend for `spatialCorrelation()` and related functions. I did not modify the repository; `git status` is clean. Every script, build and log is under
`/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/toolchain/`
(written as `toolchain/` below).

Machine: Apple M1 Ultra (16 performance + 4 efficiency cores), R 4.5.2, Apple clang 17.0.0 (clang-1700.0.13.5), MacOSX15.5 SDK, Rcpp 1.1.1, RcppArmadillo 15.2.4-1, RcppEigen 0.3.4.0.2, RcppParallel 5.1.11-1, dqrng 0.4.1, BH 1.87.0-1. I also installed RcppThread 2.4.0 into a private library, `toolchain/Rlib`.

Other agents were running throughout. Each timing table lists the `uptime` load averages measured while it ran (20 cores; load ranged from 6 to 20). The scaling numbers are therefore lower bounds on what an idle machine would give.

---

## 0. Key findings and recommendations

| # | Question | Result (evidence below) |
|---|---|---|
| 1 | GEMM speed through R's BLAS | Accelerate dgemm via RcppArmadillo: **480–590 GFLOPS** (default threads); 170–240 GFLOPS with `VECLIB_MAXIMUM_THREADS=1`. Float (`arma::fmat` → `sgemm_`): **1.9–2.2 TFLOPS**, relative error about 2e-6. **R's reference BLAS** (the CRAN macOS default, "recommended for precision") reaches only **5.8–6.0 GFLOPS**, about 80–100× slower. Eigen's built-in kernels are portable and give **31–35 GFLOPS** (double) or **62–72 GFLOPS** (float) per thread. `sgemm_` does not exist in R's reference BLAS, so float GEMM through Armadillo is not portable. |
| 2 | Variogram kernel (N=1000, P=123,987 pairs, B=900 columns) | `geoR::variog` × 900 takes **10.6 s**, and one gene with defaults needs 3,600 calls (**42 s**). With pairs precomputed once: pure R takes 0.72–0.95 s, the naive C++ column loop 133 ms, and a column-blocked tile kernel **31–34 ms** on one thread (about 3.5 G pair-columns/s). With 16 threads it takes **2.8–3.2 ms**: 10.6–12× at load ~10, and every backend gives bitwise-identical output. The C++ results match geoR to 2.2e-16. |
| 3 | OpenMP on this macOS R | `SHLIB_OPENMP_CXXFLAGS` is empty, so a package build from source is serial. Rcpp's `openmp` plugin adds `-Xclang -fopenmp -lomp`, which links **R.framework's bundled libomp**. That libomp lacks `__kmpc_dispatch_deinit`, which Apple clang 17 emits for `schedule(dynamic)`, so the .so **fails to load**; `schedule(static)` works. `-L/usr/local/lib -lomp` is silently ignored because R's link line puts `-L$(R_HOME)/lib` first. **Loading two OpenMP runtimes into one R session** (e.g. CRAN's data.table plus our module) gave **segfaults, deadlocks, or `OMP: Error #15` aborts** in 6 of 8 combinations. |
| 4 | RNG | R `rnorm`: 34 M/s. `dqrng::dqrnorm`: 98–180 M/s. C++ xoshiro256+ with Boost's ziggurat inside threads: 199 M/s per thread, and 1.7–2.1 G normals/s on 16 threads. Building one gene-direction's randomness in R (N=2000, B=1000, 9 deltas) costs **610 ms and 145 MB**. In C++ the same work takes **113 ms on 1 thread or about 9 ms on 16**, with no large allocation, and the output is **bitwise identical for 1, 2, 4, 8 and 16 threads**. Avoid dqrng's 2-argument `seed(seed, stream)` for many streams (it costs O(stream) long-jumps, 408 ms for 1,000 streams); hashed (SplitMix64) seeding takes 0.01 ms. |
| 5 | Package plumbing | I built two prototype packages, one using `std::thread` and one using RcppParallel. Both pass `R CMD check` for everything related to compiled code. The `std::thread` variant needs only `LinkingTo: Rcpp`, `Imports: Rcpp`, and two NAMESPACE lines; it needs no Makevars. RcppParallel adds `SystemRequirements: GNU make`, Makevars and Makevars.win, and `importFrom(RcppParallel, RcppParallelLibs)`; without that import the .so cannot load (`Library not loaded: @rpath/libtbb.dylib`). The repository's NAMESPACE is hand-written (`exportPattern`), so roxygen2 skips it: either edit it by hand or convert it, and converting would stop exporting `assignFill`, `getGenePixelDF` and `threshold`. |
| 6 | Dense smoothing operators | Building all 9 dense N×N operators takes 27 ms, 0.24 s and 0.68 s on one thread for N = 1000/3000/5000, and 6–61 ms on 16 threads. Storing them costs 69 MB, 618 MB and 1.7 GB. One application to N×100 takes 1.0/3.5/12 ms with Accelerate, **0.12/1.1/3.3 s with the reference BLAS**, and 1.8/8.9/39 ms when weights are computed on the fly blockwise with Eigen on 16 threads (no storage, no BLAS dependency). **Bonus:** STcompare's `locfit` smoothing is *exactly* linear in the response (error about 1e-16), so the exact operator can be extracted once per dataset. An exact Gaussian kNN smoother does **not** reproduce locfit's tree-and-interpolation output (max abs difference 0.03–0.2). |

**Recommendations**

1. **Threading: plain `std::thread`** behind a small internal `parallel_for`. It should use an atomic work counter for dynamic scheduling, catch exceptions in workers and re-throw them on the main thread, join with RAII, and poll for interrupts from the main thread. It has no dependencies, no GNU make requirement, no TBB load-order trap, and no persistent pool when R forks. It scaled as well as TBB and OpenMP in interleaved runs (10.9× vs 10.6× vs 11.1× at 16 threads). RcppParallel is an acceptable second choice. **Do not use OpenMP** on macOS. If it is used at all, use only `$(SHLIB_OPENMP_CXXFLAGS)` with `#ifdef _OPENMP` guards, which gives threads on Linux and Windows and serial code on macOS source installs.
2. **Thread count is an explicit argument** passed into C++ (`nThreads`). Default to 1, or at most 2, to follow Bioconductor and CRAN guidance. Do not combine fork-based BiocParallel with threaded or BLAS-heavy C++ work: forked children calling multithreaded Accelerate GEMM crashed in 3 of 3 runs. Make gene-level BiocParallel optional and mutually exclusive with `nThreads > 1`.
3. **Do not make performance depend on R's BLAS**. It is fast only for users who switched to Accelerate, OpenBLAS or MKL. Use Eigen's header-only GEMM inside our own threads by default. Optionally call R's BLAS, from the main thread only, when `La_library()` / `extSoftVersion()["BLAS"]` shows an optimized BLAS. Never call `sgemm_`.
4. **RNG entirely in C++ worker threads.** Seed each work item with a stream keyed by (master seed, gene, direction, permutation). Derive the master seed from the `seed` argument, or from R's RNG once on the main thread. Stop calling `set.seed()` inside package functions: the current code overwrites the user's global RNG state at `R/spatialCorrelation.R:99` and `:242`, and gives every gene the same permutation index sets (§6.4). Prefer a small self-contained xoshiro256++/SplitMix64 plus ziggurat (public-domain reference code) over dqrng, because dqrng is **AGPL-3** and also pulls in BH.
5. **Build the smoothing operators once per dataset pair** and reuse them for all genes and both directions; they depend only on coordinates and delta. Choose between storing them (fastest with a good BLAS, 0.6–1.7 GB at N = 3–5k) and computing them on the fly blockwise (portable, a few tens of ms per gene-direction at 16 threads). If the old results must be reproduced exactly, extract the exact operator from locfit; at N = 3000 that costs roughly 20–70 s per delta on one thread, paid once.

---

## 1. Environment facts (with evidence)

- `R CMD config CXX` → `clang++ -arch arm64 -std=gnu++17`; `CXXFLAGS = -falign-functions=64 -Wall -g -O2`. C++17 is already the default (the package has `Depends: R (>= 4.3.0)`), so `CXX_STD` is not needed.
- `$(R_HOME)/etc/Makeconf` lines 160–162: `SHLIB_OPENMP_CFLAGS =`, `SHLIB_OPENMP_CXXFLAGS =`, `SHLIB_OPENMP_FFLAGS =` (all empty). `R CMD config SHLIB_OPENMP_CXXFLAGS` reports "no information for variable".
- BLAS: `libRblas.dylib -> libRblas.vecLib.dylib`. The symlink is dated Nov 20 2025, while the libraries are dated Oct 31 2025, so the user switched it. `sessionInfo()` reports `BLAS: .../vecLib.framework/.../libBLAS.dylib`. The R for macOS FAQ says: *"Currently the default is to use the R BLAS: this is recommended for precision."* The vecLib shim re-exports all of Accelerate. `dladdr` on Armadillo's `dgemm_`/`sgemm_` resolves both to `/System/Library/Frameworks/Accelerate.framework/.../libBLAS.dylib`. The reference `libRblas.0.dylib` exports `_dgemm_`, `_dgemmtr_` and `_zgemm_`, but **no `_sgemm_`** (`nm -gU`).
- OpenMP runtimes present:
  - `$(R_HOME)/lib/libomp.dylib` (shipped by CRAN R; 1,188 exported text symbols; **no `___kmpc_dispatch_deinit`**)
  - `/usr/local/lib/libomp.dylib` (mac.r-project.org tarball; has it)
  - `/opt/homebrew/opt/libomp` (libomp 23.1.1; has it)
  - `omp.h` only in `/usr/local/include` and Homebrew; R.framework ships no `omp.h`.
- mac.r-project.org/openmp says: *"CRAN R ships with libomp.dylib from here in $R_HOME/lib (corresponding to the Xcode version used on CRAN)"*; *"Xcode 16.3 (Apple clang 1700+) which is incompatible with previous versions"*; *"There is a potential for chaos when more OMP run-times get loaded into one process as they may clash, in fact it is strongly discouraged to use static run-time"*.
- Two of the 240 package `.so` files in the system library link libomp: **RcppArmadillo** and **data.table**. Both are CRAN binaries built 2026 on aarch64 and link `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib`. `data.table::getDTthreads(verbose=TRUE)` reports `_OPENMP 201811` and 10 threads. So CRAN's macOS binaries can be OpenMP-enabled even though source builds on a user's Mac are not.

---

## 2. Test 1: GEMM through RcppArmadillo, R's BLAS, Eigen, and the reference BLAS

Code: `toolchain/01_gemm/gemm_arma.cpp`, `gemm_eigen_ref.cpp`, `run_gemm.R`. Timing is done in C++ around the multiply only. "First" is the first single call; "median" is over 10 further calls.

To measure the reference BLAS I made a copy of `libRblas.0.dylib`, gave it a new install name, re-signed it ad hoc, and loaded it with `dlopen`: `toolchain/refblas/libRrefblas.dylib`.

**(2000×2000) %*% (2000×100)**, 0.8 GFLOP. Load averages 6.4–8.6.

| Method | first call | median | GFLOPS (median) | GFLOPS (first) |
|---|---|---|---|---|
| arma double → dgemm_ (Accelerate, default threads) | 1.72 ms | 1.67 ms | **479** | 466 |
| arma float → sgemm_ (Accelerate) | 0.86 ms | 0.43 ms | **1871** | 928 |
| Eigen double (built-in NEON kernel, 1 thread) | 27.8 ms | 25.9 ms | 31 | 29 |
| Eigen float | 12.9 ms | 13.0 ms | 62 | 62 |
| R `%*%` (Accelerate) | 3.14 ms | 3.27 ms | 245 | 255 |
| **R reference BLAS dgemm_** | 143 ms | 133 ms | **6.0** | 5.6 |

**(1000×3000) %*% (3000×900)**, 5.4 GFLOP.

| Method | first call | median | GFLOPS (median) | GFLOPS (first) |
|---|---|---|---|---|
| arma double → dgemm_ (Accelerate) | 13.1 ms | 9.1 ms | **593** | 413 |
| arma float → sgemm_ | 2.68 ms | 2.42 ms | **2233** | 2016 |
| Eigen double | 157 ms | 156 ms | 34.5 | 34.3 |
| Eigen float | 76 ms | 75 ms | 71.6 | 71.1 |
| R `%*%` | 12.6 ms | 12.0 ms | 451 | 429 |
| **R reference BLAS** | 900 ms | 936 ms | **5.8** | 6.0 |

Effect of Accelerate's thread count (separate processes; `VECLIB_MAXIMUM_THREADS` is read once when Accelerate starts):

| VECLIB_MAXIMUM_THREADS | double, shape A | double, shape B | float, shape A | float, shape B |
|---|---|---|---|---|
| unset (default) | 479 | 593 | 1871 | 2233 |
| 4 | 420 | 536 | 1686 | 1773 |
| 1 | 173 | 241 | 644 | 792 |

Float accuracy (max abs error / max |C|): 1.9e-6 and 2.3e-6.

Notes:
- **Float GEMM is not portable through Armadillo.** Armadillo calls `sgemm_` directly (`armadillo_bits/translate_blas.hpp:63`). The compiled test .so has `_sgemm_` undefined (`nm -u`), and R's reference BLAS does not provide it. With macOS's `-undefined dynamic_lookup`, a user on the default BLAS would get a load-time "symbol not found"; on Linux with R's internal BLAS, the same. It only works here because the vecLib shim re-exports Accelerate. If float is wanted, use Eigen's float kernel.
- Apple's AMX units give R's BLAS on Apple Silicon a large advantage. That advantage disappears for default-BLAS users, which is most macOS and many Linux users.
- Calling R's BLAS from several of our own threads at once causes nested oversubscription with threaded BLAS libraries, and some OpenBLAS builds are not thread-safe. If R's BLAS is used, call it from the main thread with large batched matrices.

## 3. Test 2: variogram-like pairwise kernel

Code: `toolchain/02_vario/vario.cpp`, `run_vario.R` (log `run_vario.log`), pure-R baselines in `pure_r_baselines.R` / `.log`, shared data in `toolchain/common_data.R`.

Data:
- 3,001 unit-grid pixels in a disk, with a random subsample of N = 1000, as `viladomatCorrelation` does (`R/spatialCorrelation.R:253-255`).
- `max.dist` = 25th percentile of pairwise distances, and geoR's bins replicated exactly: `.define.bins` gives 13 regular bins from 0 to umax plus a nugget bin; the `binit` rule assigns bin k when `lims[k] <= d < lims[k+1]`.
- P = **123,987** pairs (24.8% of 499,500).
- Z is a 1000×900 matrix of normals (900 = 100 permutations × 9 deltas).

Correctness:
- Per-bin pair counts are identical to `geoR::variog$n` (645, 2519, …, 14257).
- `v = S/(2n)` matches `geoR::variog$v` to **2.2e-16** on 5 columns.
- All C++ variants give results identical to the column-major reference (relative error 0, or 1.1e-14 for the bin-sorted reordering of sums).

Serial results. Load average 11.3 at start.

| Variant | 900 variograms | G pair-columns/s |
|---|---|---|
| `geoR::variog` × 900 (R loop; 11.73 ms per call, extrapolated from 30 calls) | **10,560 ms** (one gene with defaults = 3,600 calls = **42 s**) | 0.011 |
| (a) naive C: recompute all 499,500 distances and bins per column (what geoR's C code does per call) | 3,835 ms | 0.029 |
| pure R, precomputed pairs, chunked indexing + `rowsum` | 950 ms | – |
| pure R, sparse incidence matrix (Matrix) + `rowsum` | 720 ms | – |
| (b) C++ column-major, precomputed pairs (scatter into `s[bin]`) | 132.6 ms | 0.84 |
| (b2) column-major, pairs sorted by bin, 4 register accumulators | 48.8 ms | 2.29 |
| (c) pair-major on column-major Z (strided inner loop) | 126.3 ms | 0.88 |
| (d) pair-major on transposed Z (contiguous inner loop, transpose included) | 34.2 ms | 3.27 |
| (e) column-blocked tiles, W = 4 / 8 / 16 / 32 / 64 columns | 42.1 / **33.9** / **31.4** / 31.9 / 33.2 ms | 2.65–3.55 |

The block tile layout is N×W row-major; with W = 8, each point's row is one 64-byte cache line.

About two thirds of geoR's per-call cost is R overhead: (a) costs 4.3 ms per column in C, against 11.7 ms per `geoR::variog` call, which runs `dist()` on all 1,000 points every call (line 93 of `deparse(geoR::variog)`: `u <- as.vector(dist(as.matrix(coords)))`, followed by the `.C("binit", …)` pair loop). Precomputing the pairs, which are fixed per gene-direction because `ids` and the coordinates do not change, is the big algorithmic win. C++ adds another 25–30× on one thread.

Thread scaling: the same kernel (W = 8 tiles, 113 column blocks), 7 repetitions each, median. Load average 10.5 throughout (`run_vario.log`).

| threads | RcppParallel, tiles | RcppParallel, column-major kernel | std::thread static | std::thread dynamic |
|---|---|---|---|---|
| 1 | 34.0 ms (1.00×) | 136.4 ms | 34.1 ms | 34.2 ms |
| 2 | 17.1 (1.99×) | 68.3 (2.00×) | 17.2 (1.98×) | 17.1 (2.00×) |
| 4 | 8.71 (3.91×) | 34.0 (4.02×) | 8.60 (3.97×) | 8.75 (3.91×) |
| 8 | 4.68 (7.27×) | 17.8 (7.66×) | 4.61 (7.42×) | 4.60 (7.44×) |
| 12 | 3.56 (9.54×) | 12.4 (10.99×) | 4.39 (7.78×) | 3.20 (10.69×) |
| 16 | **2.96 (11.49×)** | 11.1 (12.32×) | 3.96 (8.62×) | **3.23 (10.59×)** |

Static partitioning loses at 12–16 threads on a loaded machine, where cores are taken by other processes and some threads straggle. Use dynamic scheduling: an atomic counter, or TBB work stealing.

Fixed overheads per call:

| Operation | Time |
|---|---|
| spawn + join 16 empty `std::thread`s | 285 µs |
| RcppParallel `parallelFor` dispatch (16 threads, 113 no-op items) | 96 µs |
| OpenMP empty parallel region, 4 threads | 17–26 µs |
| OpenMP empty parallel region, 16 threads | 95–109 µs |

All are negligible when one dispatch covers ≥ 1 ms of work, so dispatch once per gene-direction or per batch, not per permutation.

## 4. Threading backends compared under identical conditions

Code: `toolchain/09_compare/compare_backends.R`; `toolchain/08_rcppthread/rcppthread.cpp`.

The same kernel ran under five backends, interleaved and in random order within each of 5 rounds, so that changes in load affect all backends equally. The table shows the median of per-round minima over 3 calls. Load average was **~10.7** (results `compare_backends.rds`). Every output was identical to the serial result.

| threads | OpenMP dynamic (/usr/local libomp) | RcppParallel (TBB) | std::thread dynamic | std::thread static | RcppThread 2.4.0 |
|---|---|---|---|---|---|
| 1 | 33.9 ms | 34.0 | 34.0 | 33.8 | 34.3 |
| 4 | 8.72 (3.88×) | 8.84 (3.85×) | 8.73 (3.90×) | 8.79 (3.85×) | 9.27 (3.70×) |
| 8 | 4.55 (7.44×) | 4.60 (7.39×) | 4.59 (7.41×) | 4.62 (7.31×) | 5.18 (6.62×) |
| 12 | 3.20 (10.6×) | 3.49 (9.75×) | 3.23 (10.5×) | 4.03 (8.40×) | 4.12 (8.33×) |
| 16 | **3.05 (11.1×)** | **3.22 (10.6×)** | **3.12 (10.9×)** | 3.89 (8.69×) | 4.68 (7.32×) |

An earlier interleaved run at load ~17–18 (`compare_backends_load17.rds`) gave 6.2×, 6.0×, 5.5×, 4.95× and 4.2× at 16 threads, in the same order. **The backend choice does not matter for performance at this granularity.** It matters for portability, dependencies, fork safety and interrupt handling.

## 5. Test 3: OpenMP

Code: `toolchain/03_openmp/omp_template.cpp`, `omp_build.R`, `run_omp_variant.R`, `clash.R`, `clash_all.sh`; logs `clash_all.log`, `clash_part2.log`. Each variant was built with `sourceCpp(rebuild = TRUE)` using `PKG_CPPFLAGS`/`PKG_LIBS` set in the environment, then loaded in a fresh process.

### 5.1 Build and load matrix (Apple clang 17.0.0)

| Variant | Flags | Linked runtime (`otool -L`) | Result |
|---|---|---|---|
| `// [[Rcpp::plugins(openmp)]]`, `schedule(dynamic,1)` | Rcpp 1.1.1 detects Apple clang and adds `-Xclang -fopenmp` / `-lomp` | R.framework libomp | **load error:** `symbol not found in flat namespace '___kmpc_dispatch_deinit'` |
| same plugin, `schedule(static)` | same | R.framework libomp | works; `_OPENMP = 202011`; 7.1× at 8 threads |
| `-Xclang -fopenmp`, `-L$(R_HOME)/lib -lomp`, dynamic | | R.framework libomp | **load error** (same symbol) |
| same, `schedule(static)` | | R.framework libomp | works; 7.4× at 8 threads |
| `-Xclang -fopenmp`, **`-L/usr/local/lib -lomp`** | | still R.framework libomp: the link line is `... -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -o x.so x.o -L/usr/local/lib -lomp ...` | **load error**; the `-L` is silently shadowed |
| `-Xclang -fopenmp`, full path `/usr/local/lib/libomp.dylib` | | `/usr/local/lib/libomp.dylib` (absolute path; exists only if the user installed the mac.r-project.org tarball) | works; **12.0× at 16 threads** |
| `-Xclang -fopenmp`, full path Homebrew `libomp.dylib` | | `/opt/homebrew/opt/libomp/lib/libomp.dylib` | works; 11.1× at 16 threads |
| `-Xclang -fopenmp`, static Homebrew `libomp.a` | | none (runtime inside our .so) | works alone; 11.6× at 16 threads |

### 5.2 Two OpenMP runtimes in one R process

data.table is a CRAN binary that uses R.framework's libomp. In each run, `fsort` ran on 4 threads and then our kernel ran on 4 threads, in the order listed.

| Load order | Outcome |
|---|---|
| data.table → our module linked to R's libomp | OK |
| data.table → `/usr/local` libomp module | **segfault** (exit 139) |
| `/usr/local` module → data.table | **deadlock**. `sample` shows `fsort → __kmpc_fork_call → __kmp_join_barrier → pthread_cond_wait`. Killed after 10 min. |
| R-libomp module → `/usr/local` module | **deadlock** in `__kmp_hyper_barrier_gather`; killed |
| data.table → Homebrew libomp module | **segfault** |
| data.table → static `libomp.a` module | **segfault** |
| static `libomp.a` module → `/usr/local` module | OK (by luck) |
| `/usr/local` module → Homebrew module | **abort**: `OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized.` (exit 134) |
| `KMP_DUPLICATE_LIB_OK=TRUE`, data.table → `/usr/local` | **segfault** |

### 5.3 What this means per platform

- **macOS source installs (most users who build from source, and Bioconductor devel users):** `$(SHLIB_OPENMP_CXXFLAGS)` is empty, so OpenMP code compiles serially. This is safe, but the "parallel" code is silently single-threaded. Hard-coding `-Xclang -fopenmp -lomp` resolves to R's bundled libomp. That library is matched to CRAN's older Xcode, and with Apple clang 1700+ any `schedule(dynamic)` / `__kmpc_dispatch_deinit` construct fails to load. Linking any other libomp creates a second runtime, which can crash or deadlock as soon as data.table or RcppArmadillo's OpenMP runs in the same session.
- **macOS CRAN/Bioconductor binaries:** these can be OpenMP-enabled if the package's own configure step adds the flags (data.table and RcppArmadillo do), linking R's libomp built with CRAN's matching Xcode. That works for binary users, but the package then behaves differently between binary and source installs. I did not verify how Bioconductor's macOS builders handle it; assume the same constraints as CRAN.
- **Linux:** not tested here. R's Makeconf normally sets `SHLIB_OPENMP_CXXFLAGS = -fopenmp` (GCC with libgomp), so OpenMP works out of the box. Linux has its own caveats: mixing libgomp and libomp, and forking after OpenMP use (libgomp is known to hang in forked children; data.table carries "RestoreAfterFork" logic for that reason).
- **Windows:** not tested here. Rtools (GCC) supports `-fopenmp` with libgomp, and RcppArmadillo enables OpenMP on Windows (`RcppArmadilloConfig.h`). `std::thread` works with Rtools' posix threads. RcppParallel supports TBB on Windows: its NEWS records Windows/TBB support and UCRT-toolchain patches, and the Makevars.win lines are in §7.2.
- **WRE §1.2.1.1:** *"Some compilers (including Apple's for macOS …) have no OpenMP support at all, not even omp.h"*, and code should be guarded with `#ifdef _OPENMP`.

**Recommendation:** use `std::thread` (or RcppParallel) for portable, predictable threading on every platform. If OpenMP is ever added, use only `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)` / `PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)` with every `omp.h` include and `omp_*()` call inside `#ifdef _OPENMP`. My template compiles and runs serially without it. Never hard-code `-lomp` or a libomp path.

## 6. Test 4: RNG

Code: `toolchain/04_rng/rng.cpp`, `run_rng.R`, log `run_rng.log`. Load average 6.3–6.9.

### 6.1 Throughput

| Generator | normals/s |
|---|---|
| R `rnorm` (Mersenne-Twister + inversion) | 33.9 M/s |
| R `runif` | 117 M/s (uniforms) |
| `dqrng::dqrnorm` via R, Xoroshiro128++ / Xoshiro256+ / Xoshiro256++ / pcg64 / Threefry | 134 / **180** / 159 / 178 / 98 M/s |
| C++ xoshiro256+ raw 64-bit / `uniform01` | 333 M/s / 337 M/s (per thread) |
| C++ Boost ziggurat normal (what dqrng uses), raw engine / dqrng wrapper | **199 M/s** / 198 M/s per thread |

### 6.2 One gene-direction's randomness (N = 2000, B = 1000 permutations, 9 deltas)

Option A: generate everything in R and pass it in. Passing a `NumericMatrix` to C++ is zero-copy; the cost is generation plus memory.

| Item | Time | Memory |
|---|---|---|
| R `sample.int(N)` × 1000 permutations | 54.9 ms | 7.6 MB |
| R `rnorm` noise, 18,000 × 1000 | 555 ms | **137 MB** |
| same noise with `dqrng::dqrnorm` | 146 ms | 137 MB |
| **Total in R per gene-direction** | **≈ 610 ms** | **≈ 145 MB** (×2 per gene) |

Option B: generate everything in C++ inside RcppParallel workers. Each work item b is a Fisher–Yates permutation of N (Lemire bounded integers) plus 9N normals, drawn from stream b.

| threads | pre-pass long_jump seeding (store / no store) | hashed seeding (store / no store) | M normals/s (hashed, store) |
|---|---|---|---|
| 1 | 117 / 113 ms | 115 / 112 ms | 157 |
| 2 | 58.9 / 58.3 | 58.3 / 57.1 | 309 |
| 4 | 30.5 / 30.0 | 29.5 / 28.8 | 611 |
| 8 | 15.7 / 15.5 | 14.8 / 14.8 | 1216 |
| 16 | 10.4 / 9.9 | **8.6** / 9.3 | **2100** |

The output is **bitwise identical for 1, 2, 4, 8 and 16 threads** in all modes. I also checked the full matrices for 1 vs 16 threads: every column is a permutation of 1..N, and the noise has mean 0.002 and sd 1.0015.

Option C: R's RNG in a single-threaded pre-pass on the main thread, either for a master seed (negligible cost) or for the permutations (55 ms per gene-direction, as above), with the noise generated in C++.

### 6.3 Seeding costs, B = 1000 streams

| Method | Cost |
|---|---|
| one xoshiro256+ `long_jump()` | 0.79 µs |
| dqrng `generator<xoshiro256plus>(seed, stream)` for streams 0..999 (O(stream) long-jumps each, `dqrng_generator.h:89-91`) | **408 ms total**; impractical if stream indices encode gene × direction × permutation |
| successive long-jump pre-pass, 1,000 states | 0.81 ms |
| hashed seeding, `xoshiro256plus(splitmix64(seed ^ splitmix64(b+1)))` | **0.011 ms** |

### 6.4 Reproducibility implications

- R's RNG and the R API must not be used in worker threads. RcppParallel's documentation: *"The code that you write within parallel workers should not call the R or Rcpp API in any fashion."* Draw any R-RNG values on the main thread before going parallel.
- Keying streams by **work item rather than by thread** (dqrng's vignette example clones per thread) makes results independent of thread count and scheduling. I verified this bitwise.
- A C++ RNG cannot reproduce the current R implementation's random streams, so old p-values will not match bit for bit. Validate distributionally instead, by comparing null-correlation distributions and p-values across many seeds.
- **The current code has RNG side effects worth fixing:**
  - `viladomatCorrelation` calls `set.seed(seed)` (`R/spatialCorrelation.R:242`), and `matchingVariograms` calls `set.seed(seed + i)` (`:99`, called with `seed = seed + i` at `:294`). This overwrites the user's global RNG state.
  - Because `set.seed(seed)` runs at the start of every call, the subsample `ids` (`:255`) and the permutation index order (`sample(X)` at `:287`) are **identical for every gene** of the same N, and the noise streams are identical across genes and directions.
  - A C++ design should include the gene index in the stream key, unless common random numbers across genes are deliberately wanted.
- Cross-platform: xoshiro and SplitMix64 are exact integer algorithms. `std::normal_distribution` output differs between libstdc++ and libc++, so do not use it. Boost's ziggurat is header-only and portable, but its tail path calls `exp`/`log`, so libm differences could change the last bits on rare draws.
- Licences (from the installed DESCRIPTION files):

  | Package | Licence |
  |---|---|
  | dqrng | **AGPL-3** |
  | sitmo | MIT |
  | BH | BSL-1.0 |
  | RcppThread | MIT |
  | RcppParallel | GPL (>= 3) |
  | Rcpp, RcppArmadillo, RcppEigen | GPL (>= 2) |
  | STcompare | GPL-3 |

  GPLv3 §13 permits combining with AGPLv3 code, but that is a governance question for the maintainers. About 50 lines of public-domain xoshiro256++/SplitMix64 plus a ziggurat (or Marsaglia polar) would avoid both the AGPL and the BH dependency.

## 7. Test 5: package plumbing

Code: `toolchain/05_pkg/` contains two complete package copies (DESCRIPTION, NAMESPACE, R/, man/, data/ copied from the repository, plus new files), the private install libraries `lib_stdthread` and `lib_rcppparallel`, `test_pkg.R`, `check_all.sh`, and `{stdthread,rcppparallel}/check.log`.

### 7.1 Minimal changes, `std::thread` variant (recommended)

DESCRIPTION:
```
Imports: BiocParallel, Rcpp, dplyr, ...
LinkingTo: Rcpp
```

NAMESPACE (hand-written, because it uses `exportPattern`):
```
exportPattern("^[[:alpha:]]+")
useDynLib(STcompare, .registration = TRUE)
importFrom(Rcpp, sourceCpp)
```

`R/STcompare-package.R` holds the same tags as roxygen (`#' @useDynLib STcompare, .registration = TRUE`, `#' @importFrom Rcpp sourceCpp`). These reach NAMESPACE only if roxygen manages it; see 7.3.

`src/variogram.cpp` exports functions with **dot-prefixed R names**, for example `// [[Rcpp::export(.variogram_sums)]]`. With `exportPattern("^[[:alpha:]]+")`, an ordinary name would be exported to users. With the prefix, `getNamespaceExports()` lists the same 14 functions as before and none of the wrappers.

No Makevars is needed: C++17 is the default, and no extra libraries are linked.

Then run `Rcpp::compileAttributes()`. It generates `src/RcppExports.cpp`, including `R_init_STcompare`, `R_registerRoutines` and `R_useDynamicSymbols(dll, FALSE)`, and `R/RcppExports.R`.

### 7.2 Additional changes, RcppParallel variant

DESCRIPTION gets `Imports: RcppParallel`, `LinkingTo: Rcpp, RcppParallel` and `SystemRequirements: GNU make`. NAMESPACE gets `importFrom(RcppParallel, RcppParallelLibs)`. The Makevars files follow the RcppParallel documentation:

```
# src/Makevars
PKG_LIBS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "RcppParallel::RcppParallelLibs()")
# src/Makevars.win
PKG_CXXFLAGS += -DRCPP_PARALLEL_USE_TBB=1
PKG_LIBS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript.exe" -e "RcppParallel::RcppParallelLibs()")
```

The resulting .so depends on `@rpath/libtbb.dylib` and `@rpath/libtbbmalloc.dylib` with no rpath. It loads only because RcppParallel's `.onLoad` already loaded TBB and its own DLL with `local = FALSE`. Calling `dyn.load()` on the .so without RcppParallel loaded fails with `Library not loaded: @rpath/libtbb.dylib`, so the `importFrom` line is **mandatory**.

The same trap affects development: a *cached* `sourceCpp` module that uses RcppParallel failed to load in a fresh session until I called `loadNamespace("RcppParallel")` first.

If RcppArmadillo is used, add `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (and `$(SHLIB_OPENMP_CXXFLAGS)` only if OpenMP is used).

### 7.3 roxygen2

With the current hand-written NAMESPACE, `roxygen2::roxygenise()` prints *"Skipping NAMESPACE … It already exists and was not generated by roxygen2"*, so the `@useDynLib` and `@importFrom` tags are ignored.

Converting the NAMESPACE requires an existing roxygen-headed NAMESPACE first: `Rcpp::compileAttributes` errors with *"pkgdir must refer to the directory containing an R package"* when NAMESPACE is missing. After conversion, roxygen writes 11 `export()` lines plus `importFrom(Rcpp,sourceCpp)` and `useDynLib(STcompare, .registration = TRUE)`. **`assignFill`, `getGenePixelDF` and `threshold`**, which `exportPattern` exports today, would no longer be exported.

### 7.4 Build, install and check results

Install time for the `std::thread` variant was 10.3 s wall (4.4 s user); the RcppParallel variant took 23.6 s wall (5.4 s user). Both under load.

`R CMD check --no-manual --ignore-vignettes` on both variants (I did not copy the vignettes):
- "checking foreign function calls", "line endings in C/C++", "compiled code" and, for the RcppParallel variant, "compilation flags in Makevars" are all **OK**.
- The RcppParallel variant adds `INFO: GNU make is a SystemRequirements.`
- All remaining problems already exist in the repository. Fixing them is worthwhile, because Bioconductor submission requires a clean check:
  - ERROR: the `simRanPatternRasts` example calls `assays()` without `SummarizedExperiment::`.
  - WARNING: non-ASCII em dash in `R/iterativePermutations.R:302`.
  - WARNING: `sf` is undeclared, and `library()` calls for ggplot2, gridExtra and patchwork appear in package code.
  - WARNING: Rd `\usage` mismatches in `spatialCorrelationGeneExpWithinSample` and `spatialSimilarity`.
  - NOTE: no `importFrom` for stats/utils (`dist`, `quantile`, `lm`, `rnorm`, `cor.test`, `combn`, …).
  - NOTE: Rd "lost braces".

Runtime check of the installed prototypes (`test_pkg.R`, load ~5):
- `.variogram_sums` matches geoR to 4.4e-16.
- Times were 32.5 ms (1 thread), 8.9 ms (4) and **2.8 ms (16)** for `std::thread`, and 33.4, 9.0 and 3.2 ms for RcppParallel, including the R call overhead.
- Output is identical across thread counts.
- An exception thrown inside a worker comes back as an ordinary R error (`R error: boom in worker`) in both variants. An uncaught exception in a bare `std::thread` would call `std::terminate`, which CRAN policy forbids ("must never terminate the R process").

Compile time for a trivial translation unit with R's flags (best of 3; `toolchain/10_compile/compile_times.sh`):

| Header stack | Compile time |
|---|---|
| Rcpp | 0.90 s |
| + RcppParallel | 0.97 s |
| + RcppThread | 1.02 s |
| RcppArmadillo | 1.64 s |
| RcppEigen | 1.43 s, plus 5 `-Wall` warnings from Eigen's sparse/unsupported headers pulled in by `RcppEigen.h` |
| dqrng + BH | 0.88 s |

### 7.5 Bioconductor and CRAN guidance relevant to C++ and threads

- Bioconductor contributions guide, ch. 17: C/C++ code *"should adhere to the standards and methods described in the System and foreign language interfaces section of the Writing R Extensions manual"*; Rcpp *"allows seamless integration of C++ with R, and is cross-platform"*; *"Use of external libraries whose functionality is redundant with libraries already supported is strongly discouraged"*; bundled third-party code becomes the maintainer's responsibility.
- Bioconductor guide §16.4.6, "Parallel Recommendations":
  - *"Developers should be mindful of supporting the major platforms (Linux, Windows, and Mac) when providing parallel implementations."*
  - *"When using in package development, parallel operations should be set to use either one or two cores by default."*
  - The page now recommends futurize, mirai and furrr and does not mention BiocParallel.
- CRAN policy (Bioconductor's checks follow the same spirit):
  - *"If running a package uses multiple threads/cores it must never use more than two simultaneously"*.
  - compiled code must *"never terminate the R process"* (no `abort`, `exit` or `std::terminate`).
- R's `?parallel::mcfork`: *"It is strongly discouraged to use mcfork and the higher-level functions … in any multi-threaded R process (with additional threads created by a third-party library or package). Such use can lead to deadlocks or crashes."*

What I observed (`toolchain/07_fork/`):
- Our kernels were used in the parent and then in `MulticoreParam(2)` children, which is the pattern STcompare uses today. This ran correctly with `std::thread`, TBB and OpenMP (`/usr/local` libomp).
- With **multithreaded Accelerate dgemm** in the parent followed by dgemm in the forked children, **a child died in 3 of 3 runs** ("1 parallel job did not deliver a result"). With `VECLIB_MAXIMUM_THREADS=1` it succeeded in 3 of 3.
- Today's STcompare (`spatialCorrelation(nThreads = 2)` after a parent dgemm) survives, because its children do no threaded BLAS.
- So a backend that uses BLAS or threads must not run inside fork-based BiocParallel workers.

Interrupts (Ctrl-C), tested with `SIGINT` sent 2 s into a 10 s job (`toolchain/11_interrupt/`):

| Design | Outcome |
|---|---|
| `std::thread` workers, main thread polling `Rcpp::checkUserInterrupt()` every 20 ms, RAII stop-and-join | interrupted after **2.2 s**, R receives an interrupt condition |
| plain join, no polling | not interruptible (ran the full 10 s) |
| RcppThread 2.4.0 `parallelFor` with `RcppThread::checkUserInterrupt()` in workers | raised "C++ call interrupted by the user." only at the end (10 s) |

## 8. Test 6: dense smoothing operators (memory vs compute)

Code: `toolchain/06_operator/op.cpp`, `run_op.R`, log `run_op.log`.

Setup:
- The N grid points closest to the origin, for N = 1000, 3000, 5000.
- Gaussian weights `w_ij = exp(-(2.5·d_ij/h_i)²/2)` (locfit's Gaussian factor), with h_i the distance to the ⌊δN⌋-th neighbour, rows normalised; δ = 0.1, …, 0.9.
- Applied to X of size N×100.
- Methods agree to ≤ 8e-16, and the rows of L sum to 1.
- Load average 12–13.5.

| N | kNN radii for 9 δ (1 / 16 thr) | build 9 operators (1 / 16 thr) | memory, 9 operators (one: 8N² bytes) | apply 1 stored operator, Accelerate dgemm | apply 1 stored operator, **reference BLAS** | on-the-fly, 256-row blocks + Accelerate | on-the-fly, 64-row blocks + Eigen, 1 / 4 / 16 thr |
|---|---|---|---|---|---|---|---|
| 1000 | 32 / 3.8 ms | 27 / 6.1 ms | 69 MB (7.6 MB) | 1.0 ms (200 GFLOPS) | 123 ms (1.6) | 5.8 ms (34.5; weights 53%) | 9.9 / 2.8 / **1.8 ms** (20 / 72 / 112 GFLOPS) |
| 3000 | 335 / 27 ms | 239 / 22 ms | 618 MB (69 MB) | 3.5 ms (508) | 1,143 ms (1.6) | 39.8 ms (45; weights 69%) | 95.6 / 23.1 / **8.9 ms** (19 / 78 / 202) |
| 5000 | 995 / 73 ms | 677 / 61 ms | 1,717 MB (191 MB) | 12.0 ms (417) | 3,325 ms (1.5) | 107.5 ms (46.5; weights 72%) | 277 / 77 / **38.7 ms** (18 / 65 / 129) |

What this shows:
- **Building the operators is cheap.** Under 1 s on one thread even at N = 5000, and ≤ 61 ms on 16 threads. They depend only on coordinates and δ, so they are built **once per dataset pair** and shared by all genes and both directions. Memory, not time, is the cost: 0.6 GB at N = 3000 and 1.7 GB at N = 5000 for 9 operators.
- **Ways to cut memory:**
  - Store only the rows at the variogram subsample `ids`: 1000×N per δ, i.e. 216 MB at N = 3000 or 360 MB at N = 5000 for 9 δ. Then smooth all N rows only for the δ* chosen for each permutation, grouping columns by δ*.
  - Compute weights on the fly in blocks, so nothing is stored. With Eigen on 16 threads this is only about 2.5–3× slower than stored-plus-Accelerate.
  - Use float storage only with Eigen's float kernel, never `sgemm_`.
- **Reference BLAS users** would see 1.1–3.3 s per operator application. That is slower than single-threaded on-the-fly Eigen at 0.1–0.28 s. The `'T'` (transposed) path of the reference dgemm is particularly slow: 1.5 GFLOPS, against 5.9 for `'N'/'N'` in Test 1.

**locfit linearity** (`toolchain/12_locfit_operator/`). On kidney sample A (N = 282), STcompare's smoothing call `fitted(locfit(z ~ lp(x, y, nn = δ, deg = 0), kern = "gauss", maxk = 300))` uses the default `ev = rbox()`: a kd-tree with vertex interpolation. It behaves as follows:
- **Exactly linear:** `sm(2z1 − 3z2) − (2·sm(z1) − 3·sm(z2))` has max error 3e-16 to 9e-16.
- The operator extracted by smoothing the N unit vectors **reproduces locfit to 1.4e-16 to 3.3e-16**. Extraction took 0.3–1.1 s per δ at N = 282.
- An *exact* Gaussian kNN smoother, truncated or not, differs from locfit by up to **0.20 at δ = 0.1** (sd of the fit 0.23), 0.04 at δ = 0.5, and 0.03 at δ = 0.9.
- One locfit call costs 2.6–8.6 ms at N = 1000, 6.4–23 ms at N = 3000, and 11.6–35 ms at N = 5000 (`locfit_timing.R`). Exact extraction would therefore take roughly 20–70 s per δ at N = 3000 on one thread. It is parallelisable across unit vectors and amortised over all genes.
- Design choice for the algorithm team: reproduce the current results exactly by extracting L from locfit, or define the method as the exact kernel smoother (fast to build in C++), which deliberately changes results.

## 9. Projected per-gene cost

These are estimates assembled from the component timings above, not end-to-end measurements.

For N = 3000 pixels, B = 100 permutations, 9 δ, both directions, the current implementation needs about:
- 3,600 × 11.7 ms of `geoR::variog` ≈ **42 s**,
- plus 1,800 locfit calls × about 12 ms ≈ **22 s**,
- plus `lm` and other overhead,

so **more than a minute per gene** on one core.

A C++ design per gene would need:
- smoothing: 2 × (9 applications of the 1000-row subset plus one full application) for the chosen δ,
- 2 × 1,800 variograms,
- 2 × 2.7 M normals.

| Configuration | Smoothing | Variograms | RNG | Total per gene |
|---|---|---|---|---|
| 1 thread, Eigen | ≈ 0.77 s | 0.12 s | 0.03 s | ≈ 0.9 s |
| 16 threads, Eigen | ≈ 70 ms | 12 ms | 3 ms | **≈ 0.09 s** |
| 16 threads, stored operators + Accelerate | ≈ 30 ms | 12 ms | 3 ms | **≈ 0.05 s** |

That is roughly 70× faster on one thread and about 1,000× faster with 16 threads. Gene batching (more columns per GEMM) and the iterative 1,000-permutation round scale linearly in B.

## 10. Files produced

All paths are under `/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/toolchain/`:

| Location | Contents |
|---|---|
| `common_data.R` | shared synthetic data (geoR-exact pairs and bins) |
| `01_gemm/` | `gemm_arma.cpp`, `gemm_eigen_ref.cpp`, `run_gemm.R`, results `gemm_{default,veclib1,veclib4}.rds` |
| `refblas/libRrefblas.dylib` | renamed copy of R's reference BLAS, for `dlopen` benchmarking |
| `02_vario/` | `vario.cpp` (all kernel variants plus RcppParallel / `std::thread`), `run_vario.R`, `run_vario.log`, `vario_results.rds`, `pure_r_baselines.R/.log` |
| `03_openmp/` | `omp_template.cpp`, `omp_build.R`, `run_omp_variant.R`, `clash.R`, `clash_all.sh`, `clash_all.log`, `clash_part2.log`, `omp_*.rds` |
| `04_rng/` | `rng.cpp`, `run_rng.R`, `run_rng.log`, `rng_results.rds` |
| `05_pkg/` | prototype packages `stdthread/STcompare`, `rcppparallel/STcompare` (`src/`, Makevars, DESCRIPTION, NAMESPACE, `R/STcompare-package.R`, generated RcppExports); `install_*.log`; `*/check.log`; `test_pkg.R`; `check_all.sh`; roxygen experiments `roxy_keep/`, `roxy_convert/` |
| `06_operator/` | `op.cpp`, `run_op.R`, `run_op.log`, `op_results.rds` |
| `07_fork/` | `fork_test.R`, `fork_blas.R`, `fork_stcompare.R` and logs |
| `08_rcppthread/` | `rcppthread.cpp`, `run_rcppthread.R`, `rcppthread_results.rds` |
| `09_compare/` | `compare_backends.R`, `compare_backends.rds` (load ~10.7), `compare_backends_load17.rds` |
| `10_compile/` | `compile_times.sh`, `hdr_*.cpp` |
| `11_interrupt/` | `interrupt.cpp`, `interrupt_test.R`, `run_interrupt.sh` |
| `12_locfit_operator/` | `locfit_linearity.R`, `locfit_timing.R` |
| `Rlib/` | private RcppThread 2.4.0 install |
