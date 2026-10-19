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
| 1 (done 2026-10-04) | The C++ engine is now the only implementation behind the exported functions; the legacy R code and `matchingVariograms()` are removed, and geoR and locfit moved to Suggests. The engine also supports adaptive stopping (§7: growing batches, an active set, exact stopping index) | `bench/validate-published.R` reproduces all six published analyses (3,734,400 stored nulls): identical delta*, nulls within 1.1e-13, BH p bit-identical, in about 2 min on 16 threads against about 26 h of the authors' runs. About 500–650× faster per gene on 1 thread. 2,819 lean tests in about 50 s. ASan/UBSan/TSan clean. `R CMD check --as-cran`: 0 errors |
| 2 (done 2026-10-05) | Correctness fixes (see `investigation/09`), each with a NEWS entry and an updated test: (b+1)/(B+1) p-values everywhere; BH across genes in `spatialCorrelationGeneExp`; `BPPARAM` passed through; NA rows with one warning each instead of crashes (δ·N ≥ 2, duplicated coordinates, fewer than 3 pairs); `spatialSimilarity()`, `savePlots()` and `plotCorrelationGeneExp()` fixes (B12–B36) | No silent NA rows; no segfaults; `BPPARAM`/`nThreads` honoured; global RNG untouched. Met |
| 3 (done 2026-10-05) | `compareSpatial()` (`compare-spatial-spec.md`): independent per-gene random streams, adaptive Besag–Clifford p-values (§7), a combined p = max(pX, pY) with BH, tidy results with methods, input validation (pixel coordinates, not just names), progress and ETA, the vectorised similarity shared with `spatialSimilarity()`. After the acceptance review: rarely detected genes are not tested (`minDetected`, spec §9) | Docs agree with code; examples run in < 5 s each; adaptive p-values are super-uniform on the null calibration tier. Met; the surrogate null itself is too optimistic for sparse genes (§8) |
| 4 (done 2026-10-05) | Documentation: three vignettes (Getting started, How STcompare works, Parameters, performance and reproducibility) and two case-study articles that run the published analyses live in one to two minutes, README, pkgdown menus and reference groups | Every new user question in `investigation/06` §"What a new user can't find" has an answer. Met, except an FAQ |
| 5 (done 2026-10-05) | Hygiene: roxygen-managed NAMESPACE (11 exports), declared dependencies, `inst/extdata` results moved to `bench/published`, `R CMD check` clean | 0 errors/warnings (1 NOTE: CRAN incoming); tarball 2.4 MB. No CI yet |

## 6. Decisions (made 2026-10-04)

1. **Default engine semantics.** The existing functions keep legacy-exact nulls and delta*, just faster. Statistical changes go into a new main function with modern defaults.
2. **The p-value formula.** (b+1)/(B+1), with b counting |null| ≥ |r|, wherever a p-value is computed, including the old functions. Done.
3. **The BH no-op in `spatialCorrelationGeneExp`.** Fixed everywhere. Done.
4. **Random streams.**
   - The old functions keep the legacy shared streams, so published results reproduce.
   - The new main function gives each gene and direction an independent counter-based stream, keyed by gene name.
5. **Upstreaming.** No pull requests yet. The work will go to `JEFworks-Lab/STcompare` as PRs later, so changes are kept separable by PR-sized unit.

## 7. Adaptive permutation testing (added 2026-10-04)

**The idea (Kamil's request).** Stop testing a gene once its p-value is clearly large. Keep testing genes with few exceedances so their small p-values are resolved.
- With b = 90 of B = 100, there is no need to continue.
- With b = 1 of 100, continue to about 1,000 to resolve p ≈ 0.01.
- With b = 0 of 100, the strength is unknown, so push to 1,000 or 10,000.

**The rule.** Besag & Clifford (1991), "Sequential Monte Carlo p-values", *Biometrika* 78:301–304.
- **Stopping:** for each gene, draw permutations until h of them are at least as extreme as the observed statistic, or until a cap n_max.
- **p-value:**
  - if the h-th exceedance occurs at permutation L, then p = h / L;
  - otherwise p = (b + 1) / (n_max + 1), which matches the fixed-B formula.
- **Validity:** this p-value is exactly valid. Under the null, P(p ≤ α) ≤ α at every α, so BH remains applicable. A naive "continue while p looks small" rule that reports (b + 1)/(B + 1) at a data-dependent B does not have this guarantee at every α.
- **Precision:** the relative standard error of p is about 1/√h, which is about 30% at h = 10 and 22% at h = 20.
- **Cost:** a null gene uses on average about h·(1 + ln(n_max / h)) permutations. At h = 10 and n_max = 10⁴ that is about 80, against 100–1,000 today. Significant genes run to n_max.
- **Proposed defaults:** h = 10 and n_max = 10,000, both user-settable. Setting h = ∞ gives fixed B = n_max.

**Both directions in lockstep.** Run the permute-X and permute-Y directions on the same permutation indices b, and stop the gene when either direction reaches h.
- The combined p-value max(pX, pY) is what the package's "both directions significant" rule uses.
- It then equals the maximum of the two standalone Besag–Clifford p-values exactly. When X stops first at L, Y's standalone p is at most h/L.
- It is therefore valid without assuming the directions are independent, and it stops non-significant genes as soon as either direction shows it.
- Per-direction values are reported as descriptive: h/L for the direction that stopped, and an upper bound for the other.

**Exact results with batched computation.**
- The engine computes permutations in batches that grow geometrically (for example 64, 128, 256, …, up to n_max) over the set of still-active genes.
- After each batch it finds each gene's stopping index L exactly from the ordered nulls (the position of the h-th exceedance) and retires stopped genes.
- Permutations computed past L are discarded. Because surrogate b depends only on (seed, gene, direction, b), the results depend neither on the batch schedule nor on the number of threads.

**Multiple testing.**
- Apply BH to the combined p-values.
- Flag genes that reached n_max with b = 0, whose p-value is limited by resolution.
- Warn when 1/(n_max + 1) is too coarse for BH significance at the observed number of genes, and suggest a larger n_max. The minimum p must be at most α·k/G for the k-th ranked gene to be significant.

**Where it lives.**
- It is a statistical change, so it goes into the new main function (phase 3). There it replaces the legacy two-stage scheme of `spatialCorrelationGeneExpIterPermutations()`, which stays unchanged for reproducibility.
- The engine support is built in phase 1, so no rework is needed: grouping tasks into stopping units (both directions of a gene), growing batches, an active set, and exact stopping indices. The engine spec, §6, gives the details.

**Limit of the guarantee.** The p-value is valid with respect to the null distribution that the surrogates
sample. For genes detected in a few pixels, the surrogates' nearly normal values make extreme null
correlations too rare, so p is too small however many permutations are drawn (§8).

**Possible later refinements (opt-in):**
- BH-aware allocation of permutations, spending effort only on genes whose BH decision is still uncertain (MMCTest and QuickMMCTest, Gandy & Hahn).
- Anytime-valid confidence intervals for p.
- Tail approximation for very small p (for example a generalised Pareto fit; Knijnenburg et al. 2009). This is a method change.

## 8. Resolved: the surrogate null of sparse genes (added 2026-10-05, resolved 2026-10-05)

**Resolution.** The first option below, amplitude adjustment, was implemented as `compareSpatial(surrogate =
"remap")` (`dev/surrogate-remap-spec.md`; `EngineTask::remap` in `src/stc_engine.h`) and validated by
`bench/calibrate-surrogates.R`, whose results and recommendation are in `bench/calibration-results.md`: remapped
surrogates are calibrated for independent genes detected in 1 to 100 percent of the pixels on the kidney and
brain grids (gaussian surrogates are anti-conservative below about 10 percent), have the same power on
correlated fields, recover the published genes at least as well, and cost about 7 percent more per permutation.
The maintainer made `"remap"` the default (2026-10-05); `minDetected = NULL` now means 0 pixels (no filter) with
`"remap"` and `sqrt(N)` pixels with `"gaussian"`, which stays available for comparability with the legacy
functions and the published analyses. The other two options were not pursued. The original note follows.

The acceptance review of `compareSpatial()` found that independent sparse genes get p-values that are too
small in the far tail: P(p ≤ 0.001) was 0.008–0.03 for genes with 3–150 nonzero pixels of 311, and 0.003–0.015 for
genes detected in 0.5–10% of 2170 pixels. The cause is the Viladomat surrogate itself (smoothing plus Gaussian noise keeps the
variogram, not the marginal distribution), not the adaptive scheme: the excess kurtosis of the permutation
distribution of r grows roughly as N / (kx·ky) for genes detected in kx and ky pixels, and the surrogates do
not reproduce it. The legacy functions share the problem, hidden by their 100 permutations.

What was done first: `minDetected` (then √N pixels by default) skipped the genes where the effect is largest, and
the documentation said that the p-values of other sparse genes can still be somewhat too small. Options for the
maintainers, each a method change to validate on the calibration tier (the first one was adopted, see above):
- **Amplitude adjustment.** Rank-remap the original values onto each surrogate (as in AAFT surrogates), so that
  every surrogate has exactly the marginal distribution of the data. For a gene without spatial structure this
  becomes the plain permutation test, which is exact.
- **A permutation guard.** Also count exceedances of the plain shuffle (cor(x[π], y), almost free, from the
  permutation already drawn) and report the larger sequential p-value; it is valid wherever either null is.
- **Rank correlation.** Test Spearman's correlation, whose null depends much less on the marginals.

