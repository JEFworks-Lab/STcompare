# Plan: a C++ backend for STcompare

Status: proposal, 2026-10-03, based on upstream commit `2983c99`. The evidence for each claim is in `investigation/`, by topic.

## 1. Where the time goes today

`spatialCorrelation()` builds a Viladomat et al. (2014) null for each gene in both directions (permute X, permute Y). It does this B times, by default 100:

1. Shuffle the values.
2. For each delta in the grid (9 by default), smooth with locfit, take a `geoR::variog` variogram, fit OLS against the target variogram, rescale and add noise, then take a second variogram.
3. Keep the delta whose variogram fits best.

At the defaults that is **1,800 locfit fits and 3,600 variogram calls per gene**.

| Shared pixels N | 273 | 1,000 | 2,000 | 4,000 | 8,000 |
|---|---|---|---|---|---|
| s/gene (B = 100, 1 thread) | 8.3 | 39 | 61 | 73 | 94 |

Profile at N = 1,000:

| Component | Share of time |
|---|---|
| `geoR::variog` (includes 26% recomputing all pairwise distances with `dist()` on every call) | 62% |
| locfit + `fitted()` | 22% |
| `lm()` | 8% |
| BiocParallel | 3% |

Garbage collection takes 14–32% of the time. Each gene forks six times, four of them only to pull values out of a list.

Real analyses take hours:

| Analysis | Authors' wall time | Projected single-thread time |
|---|---|---|
| Brain MERFISH vs Visium | 1.8 h on 22 threads | ≈ 34 CPU-h |
| MERFISH replicates | 16.8 h on 20 threads | ≈ 135 CPU-h |

Those projections assume the iterative scheme, in which 45–82% of genes are re-run at 1,000 permutations. This is why the vignettes load precomputed `.RData` files instead of running the analysis.

## 2. The key fact: everything expensive depends only on coordinates

All the expensive work is fixed by the pixel coordinates and the delta grid. None of it depends on the expression values. (All verified; see `investigation/01`, `02` and `05`.)

- **The smoother is a fixed linear operator for each delta.** locfit's default `fitted()` fits on an adaptive `rbox` tree (cut 0.8) and interpolates bilinearly, so S_δ = M_δ · W_δ.
  - W_δ is a dense matrix of m × N vertex weights; m is 9–235 for the default grid and about 1,500–2,000 at δ = 0.01.
  - M_δ is a sparse N × m interpolation matrix with at most 5 non-zeros per row.
  - The C++ replica `prototypes/lfsmooth.cpp` is `identical()` to `fitted()`. Building all 9 operators takes 4–61 ms and at most 24 MB, against 2.7–34 s per gene in locfit.
- **The variogram is a sum over a fixed pair-to-bin table.** The C++ port is bit-identical to geoR and 70–430× faster per call.
- **The OLS step has a closed form**, accurate to 7e-16 relative, and it picks the same delta*.
- **In the legacy code the random draws are also the same for every gene.** One seed is used throughout, so all genes and both directions share the subsample, the permutation indices and the noise.

So one *plan* per dataset pair serves every gene, both directions, every permutation and every iterative round. The plan holds the subsample, the pair table, the operators and, in legacy mode, the random draws.

**Measured.** The prototype engine (`prototypes/fastgene.cpp`, one thread) runs one gene, both directions, at N = 4,984 and B = 100 in **1.1 s, against about 91 s in R**. It gives the same delta* and p-values as the package, with permutations within 5.6e-11. That 1.1 s includes 0.37 s of RNG that a shared plan would pay once for all genes, not once per gene.

**Projected, not yet measured.** Threads over genes scale about 10–11× at 16 threads (`investigation/07`). That should bring the published analyses from hours to minutes.

## 3. Architecture

```
R API (existing exported functions, unchanged signatures, + engine = c("cpp", "R"))
  └─ .stc_plan(pos, seed, maxDistPrctile, deltas)            # once per dataset pair
       ids (legacy: set.seed(seed); sample(N, 1000) if N > 1000)
       prctile and umax computed in R (exact dist()/quantile semantics; see investigation/02)
       variogram pair table on ids (geoR order: j outer, i > j inner; hypot < umax; [lo, hi) bins; n >= 2)
       per-delta operators (M_δ sparse, W_δ dense) from the locfit-tree replica
       legacy draws: permutation index matrix (N x B, Mersenne-Twister) and noise (N x K x B, L'Ecuyer-CMRG,
         set.seed(seed + i)), generated once in R and shared by all genes
  └─ .stc_correlate(plan, X, Y, B, nThreads)                  # all genes in one call
       std::thread pool over (gene, direction) tasks; read-only shared plan; interrupt polling on the main thread
       per task: target variogram -> for each delta: Z = W_δ P (m x B), Xd[ids] = M_δ[ids, ] Z ->
                 variogram -> closed-form OLS -> rescale + noise -> variogram -> RSS ->
                 delta* = first argmin -> full-length surrogate for delta* only -> null r -> p
```

Design rules, with the measurements behind them in `investigation/07`:

- **No dependence on R's BLAS.** Most Mac users run R's reference BLAS, which is 80–100× slower than Accelerate. Float BLAS will not even load against it.
  - The products here are small (m × N × B with small m), so use tiled plain loops or Eigen (header-only) inside our own threads.
- **Threads:** `std::thread` with dynamic scheduling.
  - Never OpenMP on macOS. R's config has no OpenMP flags, `schedule(dynamic)` fails to load, and mixing runtimes crashed in 6 of 8 load orders.
  - Never combine with fork-based BiocParallel. Forked workers running threaded Accelerate crashed in 3 of 3 runs.
- **RNG:**
  - Never call R's RNG from worker threads.
  - Legacy mode generates its draws in R once per plan.
  - A new mode uses a small self-contained counter-based generator, with one stream per (seed, gene, direction, permutation). Results then do not depend on the thread count, the global RNG is left alone, and genes get independent draws. Avoid AGPL dqrng.
- **Exactness:**
  - Compile with `-ffp-contract=off`.
  - Use the exact coordinate doubles: no rescaling, recentring or float32. Rasterized lattices have hundreds of pairwise distances that tie exactly at `max.dist`.
  - Compute `prctile` and `umax` in R so that `dist()`'s FMA behaviour matches whatever R build is running.
- **Memory:** operators take ≤ 24 MB for the default grid at N = 5,000. Small-delta grids at very large N, such as 50 µm MERFISH with about 19k pixels, may need on-the-fly W rows instead of storing W.
- **Iterative rounds:** permutation i depends only on (seed, i), so round k can *extend* round k−1 instead of recomputing it. This is exact and saves about 10%.

## 4. What "correct" means: acceptance tests

The harness is already in the repo: `data-raw/` and `tests/testthat/`. See `data-raw/README.md`.

| Tier | Data | What the C++ engine must match |
|---|---|---|
| 0, kernel | speKidney A–B and A–C (hex lattice, with ties), jittered A–B, `quakes` (irregular), a brain gene with N = 2,170 (subsample path) | variogram u, n, v identical; smoother identical; per-delta RSS within 1e-12; delta* identical; nulls within 1e-10; raw p identical |
| 1, calibration | `simRanPatternRasts` (independent null pairs) plus mixed positives | false-positive rate ≤ nominal; power at ρ = 0.6 (statistical check, for the new RNG mode) |
| 2, realistic | AKI kidney Visium (311 px, 11 deltas) and brain MERFISH vs Visium (2,170 px) subsets, with golden nulls from `inst/extdata` | nulls equal to the golden prefix within 1e-12; delta* identical; engineered negatives satisfy r → −r and nullX → −nullX |
| 3, benchmark | all genes of the AKI, MERFISH and brain analyses (≈ 130 MB of downloads, cached outside the repo) | wall time and end-to-end equality with legacy |

**Platform caveat, found while testing.** The legacy results are not bit-portable across platforms.
- **The cause:** geoR's last-bin edge `umax` comes from R's `dist()`, which uses FMA on macOS arm64 but not on x86-64. Bins are then assigned with libm `hypot`. On lattices, many tied pairs sit exactly at that edge, so one platform may count a pair that another drops.
- **The effect:** one changed bin count moves every null by about 1e-4 and occasionally flips delta*.
- **How the tests cope:** the fixtures were built on macOS arm64 and record a per-coordinate-set "platform signature" (geoR bin counts and edge distances). Exact tests run only where the signature matches; portable tests (invariants, the package against reference kernels on the same machine, loose agreement) run everywhere. The suite passes on macOS, Linux arm64 and Linux x86-64; see `data-raw/README.md`.
- **For the C++ engine:**
  - Compute `prctile` and `umax` in R so that legacy mode matches the running platform's R.
  - Separately, consider a platform-stable binning rule for the new mode. That would be a deliberate, documented method change.

## 5. Phases

| Phase | Content | Exit criterion |
|---|---|---|
| 0 (done) | Investigation, test data (`data-raw/`, `tests/testthat/fixtures/`, 1.2 MB), legacy-reference tests | `devtools::test()` passes against the current R code: 1,618 expectations in about 55 s; slow tiers are opt-in with `STCOMPARE_SLOW_TESTS=true`; mutation testing caught 35 of 36 mutants |
| 1 | Legacy-exact C++ engine behind `engine = "cpp"`, plus `bench/` scripts | Tiers 0 and 2 pass with `engine = "cpp"`; ≥ 50× single-thread speedup per gene; near-linear thread scaling |
| 2 | Correctness fixes (see `investigation/09`), each with a NEWS entry and an updated test | No silent NA rows; no segfaults (validate δ·N ≥ 2); BH applied across genes; `BPPARAM`/`nThreads` honoured; global RNG untouched |
| 3 | Usability: vectorised `spatialSimilarity` (85–90% faster just by hoisting accessors), input validation (pixel coordinates, not just names), tidy results, progress and ETA | Docs agree with code; examples run in < 5 s each |
| 4 | Documentation: "How it works", parameter guide, performance guide, FAQ, pkgdown reference groups | Every new user question in `investigation/06` §"What a new user can't find" has an answer |
| 5 | Hygiene: roxygen-managed NAMESPACE (stop exporting helpers), declare dependencies, move about 32 MB of `inst/extdata` out of the tarball (now 48.6 MB), `R CMD check` clean, CI | 0 errors/warnings; tarball < 10 MB |

## 6. Decisions needed

1. **Default engine semantics.**
   - Recommendation: the existing functions keep legacy-exact results, just faster.
   - Statistical changes go in a new front door (for example `compareSpatial()`) with modern defaults.
   - Alternative: change the old functions' outputs in a major version bump.
2. **The p-value formula.** `extreme/B` gives p = 0 for 573 of 1,046 kidney genes, while the docs promise a minimum of 0.01.
   - Switching to (b+1)/(B+1) changes every p-value.
   - Fix it in the old functions, or only in the new front door?
3. **The BH no-op in `spatialCorrelationGeneExp`.** Recommendation: fix it everywhere. The code contradicts its own documentation and the AKI vignette's text.
4. **Independent random streams per gene** (new mode) versus the legacy shared streams. Shared streams make Monte Carlo errors correlated across genes.
5. **Upstreaming.** Should this be developed as PRs to `JEFworks-Lab/STcompare`? That affects style, scope per PR, and whether behaviour changes need maintainer sign-off.
