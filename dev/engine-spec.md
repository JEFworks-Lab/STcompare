# Spec: the C++ engine (legacy-compatible mode)

This is the contract for the C++ engine, which (since 2026-10-04) is the only implementation behind the existing exported functions; the legacy R implementation was removed and survives only as stored outputs (test fixtures and the published results in `inst/extdata`). It is phase 1 of `cpp-backend-plan.md`. The evidence for every semantic rule below is in `investigation/01` (smoother), `02` (variogram), `04` (RNG streams) and `05` (method).

## 1. Goals

1. **Same results as the legacy R implementation** (as pinned by the fixtures and the published results) for every input where the R code succeeded:
   - identical `deltaStarX`/`deltaStarY` for every permutation;
   - identical `pValuePermuteX`/`pValuePermuteY`;
   - `nullCorrelationsX`/`Y` and `permutationsX`/`Y` equal within a relative tolerance of about 1e-9.

   Bit-identical nulls are not a goal. Two things prevent them:
   - R's `lm()` is not bit-stable across BLAS builds.
   - The OLS intercept suffers catastrophic cancellation, which amplifies last-bit differences in the smoothed values or the fit to about 1e-11 in the nulls.
2. **The same NA behaviour as the legacy R code.** Inputs on which the R code returned an NA row get an NA row here too, plus a `warning()` naming the gene and the reason. Where the R code crashes R, this engine returns an NA row and a warning instead.
3. **Speed:**
   - at least 50× faster than the legacy R code per gene on one thread at N = 1000–2000 pixels, default delta grid, B = 100;
   - at least 8× speedup on 16 threads for 64 or more genes;
   - results that never depend on the number of threads.
4. **No R API calls from worker threads:** no `Rcpp` objects, R memory or `R_CheckUserInterrupt`. Everything runs on plain C++ buffers copied before the threads start. The one R entry point that worker threads call is Rmath's `Rf_qnorm5()`, a pure function, for the noise (2.6).
5. **The global RNG state is restored** after every call. This is a side-effect fix only; the legacy streams themselves are reproduced exactly.

Non-goals for this phase: the new main function with independent per-gene streams and other statistical changes (phase 3), and a vectorised `spatialSimilarity()`. The C++ must still keep "where permutations and noise come from" behind one small interface, so that a per-task stream mode can be added later without touching the numerics.

## 2. Legacy semantics to reproduce

Notation for one call of `spatialCorrelation(X, Y, pos, nPermutations = B, deltaX, deltaY, maxDistPrctile, seed)`:

- N = `length(X)`.
- In legacy naming, `lat <- pos[, 1]` and `long <- pos[, 2]`.
- The forward direction permutes X and correlates with Y. The reverse direction permutes Y, with `deltaY`, and correlates with X. Both directions use the same seed.

1. **Subsample and permutation indices** (R, main thread, the session's current `RNGkind()`):
   - `set.seed(seed)`;
   - `ids <- if (N > 1000) sample(N, 1000) else seq_len(N)`;
   - then for b = 1..B, `idx[, b] <- sample.int(N, N)`. This equals `sample(X)`.

   The draws are identical for every gene and both directions. Save `.Random.seed` before and restore it after (or remove it if it did not exist).
2. **Variogram distance** (R):
   - `prctile <- quantile(dist(cbind(lat[ids], long[ids])), probs = maxDistPrctile)` (type 7).
3. **Variogram bins.** The R part:
   - `u <- as.vector(dist(cbind(long[ids], lat[ids])))`;
   - `umax <- max(u[u < prctile])`;
   - `nugget <- min(u) < 1e-12`;
   - `bl <- seq(0, umax, l = 14)`;
   - `bins.lim <- c(0, 1e-12, bl[bl > 1e-12])`;
   - `lims <- bins.lim`, then `lims[1] <- -1`.

   The C++ part:
   - Loop over pairs with `j` outer ascending and `i > j` inner ascending.
   - Compute `d = std::hypot(long_i - long_j, lat_i - lat_j)`.
   - If `d <= prctile`, find the first `ind` with `d < lims[ind]`, scanning `[lo, hi)`. A pair belongs to bin `ind - 1` only if `d < lims[last]`; pairs at or above `umax` are dropped.
   - Bins with fewer than 2 pairs are dropped. The nugget bin is kept only if `nugget`.
   - Within a bin, pairs keep loop order. `dev/prototypes/variog_fast.cpp` already does this bit-identically to `geoR::variog`.
4. **Variogram value** of a vector `z` over `ids`:
   - for each kept bin, accumulate in pair order `t = d*d; t = t/2; acc += t` in double, then divide by `n_bin`;
   - FP contraction must be off in this code.
5. **Smoother** for one delta, on all N points with coordinates `(x1, x2) = (long, lat)` (the `lp(long, lat)` order matters for the tree). This replicates locfit's default `rbox` tree; see `dev/prototypes/lfsmooth.cpp`.
   - **Bandwidth:** `k = (int)(N * delta + 1e-12)`.
   - **Tree:** `cut = 0.8`; vertex capacity from `atree_guessnv(cut, d = 2, alp = delta, maxk = 300)`.
   - **Vertex values:**
     - the kernel is `exp(-((2.5 u)^2) / 2)`, untruncated;
     - the bandwidth at a vertex is the k-th smallest distance (if `k >= N`, use the `nnk >= n` branch of `compbandwid`);
     - pseudo-vertices take the mean of their parents;
     - values reach the data points by bilinear interpolation.
   - **Precompute per delta** a factored operator `fitted = M %*% (Wn %*% y)`:
     - `Wn` is the row-normalised Gaussian weights at the real vertices (m × N, dense);
     - `M` is the interpolation (N × m, sparse, at most a few non-zeros per row, pseudo-vertices folded in).
     - The operator must be applied to centred columns, `fitted = M (Wn (y - c)) + c`, with `c = smooth_centre(y)` (locfit's parametric component), as locfit itself does. Without the centring, the error is about eps·|mean(y)|, which exceeds the 1e-9 null tolerance for inputs whose mean is large relative to their spread.
     - With the centring, agreement with `fitted(locfit(...))` is within 1e-14·range(y) + 2·eps·max|fitted|. The component API is documented in `src/stc_smoother.h`, and the engine passes the column centres to both `smooth_project` (or subtracts them while gathering) and `smooth_interp`.
   - **Errors** become an NA row:
     - `k < 2`: R errors when k = 0 and segfaults when k = 1;
     - vertex overflow: locfit's "newsplit: out of vertex space";
     - `delta <= 0`.
   - **Coordinates:** FP contraction off in the tree and distance code, and the exact coordinate doubles (no rescaling).
6. **Noise** for permutation b. Inside a BiocParallel task, legacy R runs `RNGkind("L'Ecuyer-CMRG")`, then `set.seed(seed + b)`, then for each delta index k in the gene's grid order, `rnorm(N)`. So the noise is the k-th block of N normals of that stream.
   - The normals are always Inversion normals, whatever the session's `normal.kind`. BiocParallel (1.44) runs its tasks under L'Ecuyer-CMRG with Inversion normals and Rejection sampling, whatever `RNGkind()` is in the calling session, so the legacy noise never depended on the session's normal kind, and the engine keeps that definition (tested in `test-engine.R`).
   - The C++ reproduces this stream bit for bit: R's `RNG_Init` scrambling for L'Ecuyer-CMRG, its `unif_rand`, and Inversion `norm_rand` with `BIG = 134217728`, calling R's own `Rf_qnorm5` (`src/stc_rng.h`). For the p in (0, 1] that the generator produces, `Rf_qnorm5` reads only its arguments, allocates nothing and never warns, so worker threads may call it (the exception in 1.4). Calling it, rather than porting its AS241 polynomials to C++, keeps the inversion identical to `rnorm()` on every platform: whether those polynomials are evaluated with fused multiply-adds depends on how R itself was compiled.
   - Noise drawn in R with `rnorm()` (`noise = "R"` in `.stc_engine_correlate()`) exists only as a test hook for the C++ stream; it gives identical results.
   - The block index is the position in the delta grid, not the delta value. `deltaX` and `deltaY` may differ, and so may per-gene grids.
7. **Per permutation b** (both directions):
   - `x_b <- X[idx[, b]]`.
   - For each delta index k:
     1. `xd = S_k x_b`; only the `ids` rows are needed here.
     2. `g = variog(xd[ids])`.
     3. Closed-form OLS of `target = variog(X[ids])` on `g` with intercept: two-pass means, then `b1 = Sxy / Sxx` and `b0 = tbar - b1 * gbar`.
     4. **Degenerate fit:** when `sqrt(sum((g - gbar)^2)) < 1e-7 * sqrt(sum(g^2))` (`lm()`'s rank test, observed threshold between 9e-8 and 1.2e-7), or when there are fewer than 2 bins, `lm()` returns an NA slope. The legacy code then errors in `geoR::variog`, so the whole gene becomes an NA row.
     5. `hat_k = xd * sqrt(abs(b1)) + e_k * sqrt(abs(b0))`.
     6. `rss_k = sum((variog(hat_k[ids]) - target)^2)`, summed sequentially.
   - `k* = ` the first index of the minimum (`which.min`).
   - The surrogate is `hat_{k*}` at all N points.
8. **Null correlation:** `cor(surrogate_b, target_vector)` with R's own algorithm. Means are two-pass in `long double` (R's `MEAN` macro in `cov.c`), followed by centred cross-products and sums of squares, then clamped to [-1, 1].
9. **R side:**
   - `r <- cor.test(X, Y)` as now (for the genes or pairs whose values are all finite, its estimate and p-value are computed all at once with `cor.test()`'s own arithmetic, `.stc_cor_tests()`; identical values);
   - `p <- (sum(abs(null) >= abs(r)) + 1) / (B + 1)`;
   - `deltaStarMedian <- median(deltaStar)`.
10. **NA rows.** When any step fails for a gene (in either direction), every permutation-derived column of that gene is NA, exactly like the legacy `tryCatch`.
    - **Triggers:**
      - a non-finite value in X or Y on the shared pixels;
      - a zero-variance X or Y;
      - `k < 2` or tree overflow for any delta in either grid;
      - a degenerate OLS at any (b, k).
    - `correlationCoef` and `pValueNaive` still come from `cor.test()`, as now.
    - Raise one `warning()` per affected gene. Do not `print()` it.

## 3. Work decomposition and threading

1. **Plan.** One plan per call holds everything that depends only on coordinates, seed and `maxDistPrctile`:
   - `ids`, the pair table, and the smoother operators built lazily per distinct delta value (shared by all genes and both directions);
   - the permutation index matrix (N × B, shared);
   - in legacy mode, noise blocks generated per permutation in C++.
2. **Tasks.** One task is (source vector, delta grid, target vector or vectors):
   - `spatialCorrelationGeneExp`: 2 tasks per gene;
   - within-sample: one task per gene, whose surrogates are correlated with every partner gene. Legacy surrogates do not depend on the partner, so each gene's surrogate set is computed once. This is exact.
3. **Work items.** Permutations are processed in super-chunks, for example 256 at a time, so that the noise memory stays bounded.
   - Within a super-chunk, a work item is (task, sub-chunk of about 16–32 permutations), scheduled dynamically over `nThreads` `std::thread` workers via an atomic counter.
   - A single gene therefore still uses many threads.
   - Each item writes only its own disjoint output slots.
4. **Main thread.** It polls `R_CheckUserInterrupt()` (through a safe wrapper) about every 100 ms and sets a stop flag that workers check. Worker exceptions are caught, stored, and re-thrown as an R error after `join()`.
5. **Determinism.** Every (task, permutation) computation is self-contained, with a fixed operation order. Outputs are identical for any `nThreads`.
6. **Numerics:**
   - no BLAS calls: most Mac users have R's reference BLAS, and float BLAS does not even load there;
   - register-blocked plain C++ loops for `Wn %*% P` (m × N × chunk);
   - variograms evaluated over a tile of columns, with the pair loop outermost so the inner loop over columns vectorises. Per-column summation order is unchanged.
7. **Thread count.** If `BPPARAM` is not NULL, use `BiocParallel::bpnworkers(BPPARAM)`; otherwise use `nThreads`. Never fork in the cpp path.

## 4. R API (revised 2026-10-04: the C++ engine is the only implementation)

The maintainer decided that the legacy R implementation is not kept. There is **no `engine` argument**: the exported functions keep their names, arguments and output structure, and are computed by the C++ engine.
- Functions: `viladomatCorrelation`, `spatialCorrelation`, `spatialCorrelationGeneExp`, `spatialCorrelationGeneExpIterPermutations`, `spatialCorrelationGeneExpWithinSample`.
- `matchingVariograms()` is removed. It was a legacy helper whose interface takes a geoR variogram object.
- The legacy R code paths, and the calls to locfit, geoR and `BiocParallel::bplapply` in them, are deleted. geoR and locfit move to Suggests: the component tests use them as references.
- **Threads:** `nThreads` sets C++ threads. A non-NULL `BPPARAM` is still accepted and sets the thread count through `BiocParallel::bpnworkers()`. Nothing forks.
- **`spatialCorrelationGeneExpIterPermutations`.** Round k+1 extends the promoted genes' nulls with permutations `nPermutations[k] + 1 .. nPermutations[k + 1]` instead of recomputing them. This is exact by the prefix property. Genes with NA p-values are never promoted; the legacy code crashed on them.
- **Internal helpers** live in `R/engine.R`, and their C++ entry points are dot-prefixed so the `exportPattern` does not export them.
- **`verbose`.** One message at the start (genes, directions, permutations, threads) and one at the end (elapsed time).

## 5. Acceptance tests (revised 2026-10-04: lean; nothing re-runs the legacy R code)

- **Main acceptance test: the published analyses** (`bench/validate-published.R`).
  - Rebuild the inputs from the Zenodo and 10x downloads (`data-raw/`, cached outside the repo). Run every published analysis through the exported functions with the authors' parameters, then compare per gene with `inst/extdata`:
    - `correlationCoef` and `pValueNaive`;
    - the number of permutations per gene (the same screening decisions);
    - identical delta*;
    - nulls within 1e-9 relative;
    - identical exceedance counts.
  - Analyses: kidney Visium (iterative and fixed B = 100, 1,046 genes each), MERFISH replicates (affine 483 genes; STalign 483 genes, of which the 122 rows whose stored p-values do not follow from the stored nulls, because the authors patched the table with rows of an earlier run, are compared like the others but counted separately: 116 such rows show it through `pValuePermuteX` and 6 more only through `pValuePermuteY`), brain MERFISH vs Visium (325 genes) and cell types (16).
  - This runs in minutes on 16 threads; the authors needed about 200 CPU-hours.
- **Default test suite:** about a minute, at most 2 threads.
  - Component tests: C++ against geoR, locfit and R's RNG, `cor()` and `lm()`, skipped if geoR or locfit is not installed.
  - Fixture tests through the exported functions against the stored legacy outputs: identical delta*, nulls within 1e-9, p-values computed from the nulls, NA rows. Platform gating stays as it is.
  - Edge cases that must give an NA row and a warning, never a crash:
    - constant gene, NA, `delta * N < 2`, duplicated coordinates;
    - `delta > 1`, a `maxDistPrctile` so small that fewer than 2 bins remain;
    - sparse `dgCMatrix` assays, single-gene inputs.
  - Determinism:
    - `identical()` outputs for 1, 2 and 4 threads and for different batch and sub-chunk sizes;
    - adaptive stopping with h = ∞ equals fixed B;
    - session continuation equals one longer run;
    - `.Random.seed` is unchanged.
- **Opt-in slow tier** (`STCOMPARE_SLOW_TESTS=true`): all realistic-fixture genes, and the calibration jobs.
- **Benchmarks** (`bench/`, excluded from the build): per-gene timings on the fixture genes and thread scaling.

## 6. Adaptive stopping support (added 2026-10-04)

The engine supports adaptive (Besag–Clifford) stopping from the start, so that the phase-3 main function can expose it without reworking the engine. The motivation, the statistical rule and why it is valid are in `cpp-backend-plan.md` §7. The legacy exported functions keep fixed B.

1. **Stopping units.**
   - A unit groups the tasks that stop together. For gene-wise comparisons this is the two directions of one gene, run in lockstep on the same permutation indices.
   - Each unit carries the observed |r|, h (the target number of exceedances) and `n_max` (the cap).
   - Within-sample comparisons, where one gene's surrogates are shared by many pairs, may stay fixed-B in this phase. If they become adaptive later, a gene task stays active while any pair that uses it is active.
2. **Batches.**
   - Permutations are processed in batches of growing size over the active units only. For example the first batch has 64 permutations and each later batch twice as many, capped by `n_max`.
   - Within a batch the work items are those of §3 (task × sub-chunk), dynamically scheduled.
3. **Exact stopping index.**
   - After each batch, the main thread scans each active unit's new nulls in permutation order.
   - In each direction it counts b with |null_b| ≥ |r| until the count reaches h.
   - The unit's stopping index is L = the smallest such permutation index over both directions.
   - If neither direction reaches h by `n_max`, then L = `n_max`.
   - Stopped units are retired. Their stored nulls and delta* are truncated to the first L permutations, and any surrogates computed past L are discarded.
4. **Outputs per unit:**
   - L, the stop reason (`"exceedances"` or `"cap"`) and the exceedance counts of both directions at L;
   - nulls and delta* for permutations 1..L in both directions.

   R computes the p-values:
   - combined: `h / L` if the unit stopped by exceedances, otherwise `(max(bX, bY) + 1) / (n_max + 1)`;
   - per direction: descriptive only.
5. **Invariants, enforced by tests:**
   - adaptive mode with h = ∞ and `n_max = B` is identical to fixed-B mode;
   - results do not depend on the batch schedule or on `nThreads`;
   - the adaptive result equals a brute-force run of `n_max` permutations followed by an offline Besag–Clifford computation on the same sequences.
6. **Random streams in adaptive mode:**
   - Legacy streams: the permutation indices for each batch are generated in R from a saved Mersenne-Twister stream state, so memory stays bounded at any `n_max`. Legacy L'Ecuyer noise is produced per permutation in C++.
   - Phase 3 streams: one counter-based stream per (seed, gene, direction, permutation), generated in C++. The engine reads permutations and noise through the single interface required by §1.
   - Implemented (2026-10-04) as `EngineOptions::rng = RNG_STREAMS` (`src/stc_rng.h`): every task has a 64-bit key from (seed, gene name, direction); the worker that computes permutation b of a task draws its permutation (Fisher-Yates with Lemire's bounded integers) and its noise blocks (polar normals) from xoshiro256** generators seeded from (key, b, sub-stream), through `Engine::item_perm()` and `item_noise()`, which in legacy mode return the shared batch buffers instead. The legacy mode is unchanged: the published validation gives bit-identical results. Two further session options serve `compareSpatial()`: `keep_nulls = false` keeps only the current batch's nulls and delta* (the counts of delta* per grid position are kept in every mode), and `.stc_engine_run()` takes an R progress callback that the main thread calls while the workers run.
