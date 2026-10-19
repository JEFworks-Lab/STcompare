# STcompare: baseline performance of the pure-R implementation, and where the time goes

Work dir (scripts, raw CSVs, logs, Rprof files):
`/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/profiling/`
(below, `$W` means this directory). Nothing in the repository was modified.

## 0. Summary of findings

* **Cost per gene (spatialCorrelation, 1 thread, default 9 deltas):** about 8 s at N≈275, 17 s at N=487, 39 s at N=1000,
  61 s at N=2000, 73 s at N=4000 and 94 s at N=8000 **for B=100**. Time is exactly linear in B (R² ≥ 0.994),
  so **B=1000 costs 1.3 min (N≈275), 6.4 min (N=1000), 10.2 min (N=2000), 12 min (N=4000) and 15.5 min (N=8000) per gene**.
* **Where it goes (N=1000, B=20, Rprof):** `geoR::variog` 62%, locfit (fit plus `fitted`) 22%, `lm` 8%, BiocParallel 3.5%.
  `sample`, `rnorm`, `cor` and `data.frame` together come to under 0.5%.
  **Only 18–46% of the run time is spent in the three C kernels that do the real numeric work** (geoR `binit`, locfit `slocfit`/`sfitted`).
  The rest is R-level overhead and redundant work:
  * `variog` computes all ~500k pairwise distances with `dist()` on every call, only to take `min` and `max(u[u < max.dist])`.
  * locfit and `lm` spend most of their time in formula and `model.frame` machinery.
  * garbage collection takes 14–32% of the time, and system time 6–18%.
  * every `MulticoreParam()` construction runs a full `gc()` (0.25 s).
* **Everything expensive depends only on the coordinates, not on the values.** This is the key to a C++ backend:
  * The locfit smoother is a fixed linear operator: I checked linearity to 1e-14 for all 9 deltas.
  * The variogram pair-to-bin table is fixed. A precomputed table reproduces `geoR::variog` to 8e-15 and is 5.2× faster even in plain R.
  * Both can be computed once per pixel set and reused for every permutation, delta and gene.
* **Parallelism is over permutations within one gene** and costs 6 `bplapply` forks per gene, 4 of which are trivial list extractions.
  * Speedups at 2/4/8 threads were 2.3/4.1/6.6× at N=1000, B=100, but only 2.0/3.0/3.6× at N=273.
  * A trivial `bplapply` costs 43–54 ms with 4 workers and 59–93 ms with 8.
  * **`spatialCorrelationGeneExp` ignores a user-supplied `BPPARAM`** (confirmed by tracing). In `spatialCorrelationGeneExpIterPermutations`, the
    vignette's `nThreads = 22, BPPARAM = MulticoreParam()` silently runs `detectCores()-2` workers.
* **spatialSimilarity** takes about 20 ms per gene, linear in G up to G=500.
  * At that size the cost is dominated by S4 accessor calls repeated per gene: `assay()` takes 6–7 ms per call.
  * With sparse assays, the per-gene `as.matrix()` of the whole G×N assay adds an O(G²N) term. It is already 41–47% of the time at G=1000–2000, N=4000.
  * Extrapolated to a whole transcriptome (G=20000, N=4000), the run would take about 1.5–2 h instead of seconds.
* **Realistic workloads (N=2000, single thread):**
  * (a) 100 genes at B=100: **1.7 h**.
  * (b) 500 genes with `nPermutations = c(100, 1000)` and 30% promoted: **34 h**. The shipped results show 45–82% promotion, which gives **47–78 h**.
  * (c) the same on 16 threads: **1.8–2.6 h** at 30% promotion and **3.5–4.9 h** at 72% (model range, see §7).

## 1. Environment, method, load conditions

* Apple M1 Ultra: 16 performance cores plus 4 efficiency cores (`sysctl hw.perflevel0/1.physicalcpu` = 16/4). macOS; R 4.5.2 linked to Accelerate vecLib.
* Package versions: locfit 1.5-9.12, geoR 1.9-6, BiocParallel 1.44.0, SEraster 1.2.0, SpatialExperiment 1.20.0,
  SummarizedExperiment 1.40.0, Matrix 1.7.4, bench 1.1.4.
* The package was loaded with `devtools::load_all(<repo>)`.
* Every run used `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1`, so single-thread runs are truly single-threaded.
* Timings use `system.time`. Micro-benchmarks use `bench::mark(filter_gc = FALSE)`, so GC is included.
* Every timing row in the CSVs records the `uptime` load averages taken just before the run.
* **Load from other agents:** the 1-minute load average ranged **4.0–13.7** during my runs (median about 8). It was 6–9 during the parallel-scaling block.
* **Possible bias:**
  * Single-thread runs were most likely scheduled on performance cores, since total load stayed below 16.
  * Repeat-run variability (coefficient of variation) was ≤ 6% for most cells; two B=10 cells reached 12–15%.
  * The same configuration (syn1000, B=100, 1 thread) took a median of 39.3 s in the timing block (load ≈ 11) and 35.2 s in the parallel block (load ≈ 7–9),
    so absolute numbers carry roughly ±10% uncertainty from load.
  * In the 8-worker runs, my 9 processes plus other agents' load approached the 16 performance cores. Some workers may have run on efficiency cores,
    which would understate the 8-thread efficiency. A 16-thread run would certainly contend, so I extrapolated it instead of measuring it.
* `(user + sys) / elapsed` was 0.93–0.99 for serial runs, so little time was lost waiting for a CPU. System time was 6–18% of elapsed (median 10%).

## 2. Test inputs: N at each resolution, and recommended test datasets

`$W/scripts/01_build_inputs.R` builds all inputs and saves them to `$W/inputs.rds`.

**speKidney**, via `SEraster::rasterizeGeneExpression(speKidney, assay_name='counts', resolution=r, fun='mean', square=FALSE)`.
Hexagonal pixels; assay name `pixelval`; dense `matrix`. Note that `names(speKidney)` is `A, C, B`.

| resolution | pixels A / B / C | **shared A∩B (N)** | shared A∩C | cor(A,B) | rasterize time |
|---|---|---|---|---|---|
| 0.20 | 282 / 279 / 287 | **273** | 277 | −0.947 | 0.9 s |
| 0.15 | 459 / 451 / 456 | **418** | 424 | −0.913 | 0.8 s |
| 0.10 | 733 / 743 / 779 | **487** | 527 | −0.873 | 1.2 s |
| 0.07 | 933 / 945 / 983 | **399** | 424 | −0.851 | 2.3 s |
| 0.05 | 1079 / 1090 / 1119 | **259** | 278 | −0.833 | 3.8 s |
| 0.03 | 1168 / 1184 / 1231 | **130** | 121 | −0.844 | 8.2 s |

**speKidney cannot produce large N.** A, B and C are independent point clouds with about 1230–1300 cells each.
Finer pixels raise the per-sample pixel count but shrink the overlap, so N peaks at about 490 (resolution 0.1).

**simRanPatternRasts:** 272–288 pixels per raster (dgCMatrix assay `pixelval`).
Pair 1–7 shares N=275 pixels (r = −0.219); pair 2–3 shares 274 (r = 0.180); pair 4–5 shares 275 (r = 0.151).

**Synthetic data (for N ≥ 1000):** a square grid of exactly N pixels (N = 280, 500, 1000, 2000, 4000, 8000) with two independent smooth fields.
Each field is white noise smoothed with a Gaussian kernel (length scale = side/10, comparable to κ=0.1 on the unit square in R/data.R), plus N(0, 0.3²) noise, plus 10.
On these fields the selected delta (δ*) is mostly 0.1, the bottom of the default grid.

**Recommended test datasets (for regression and C++ equivalence tests):**
* **Fast unit tests** (< 1 s at B=10):
  * kid0.2 A vs B: N=273, strong negative correlation, p_X = p_Y = 0.
  * kid0.2 A vs C: N=277, positive correlation.
  * simRanPatternRasts 1 vs 7: N=275, null case; p_X=0.30 and p_Y=0.28 at B=100 with seed 0.
* **Medium:** kid0.1 (N=487), the largest N speKidney can give.
* **Scaling and benchmarks:** synthetic N = 1000/2000/4000/8000 from `inputs.rds`.
* **Golden outputs:** results are bit-identical across `nThreads` (§5) and deterministic given `seed`, so current R outputs
  (nullCorrelations, deltaStar, p-values) can be stored as golden references.
  Because δ* is a discrete argmin, a numerically different backend can flip choices. Equivalence tests should therefore combine
  exact checks on sub-steps with distributional checks (p-value and null-correlation agreement within tolerance).

## 3. spatialCorrelation wall time vs N and B (1 thread)

Script: `$W/scripts/03_time_spatialCorrelation.R`. Raw data: `$W/timings_spatialCorrelation.csv`; medians: `$W/timings_spatialCorrelation_medians.csv`.
Call: `spatialCorrelation(X, Y, pos, nPermutations = B, nThreads = 1)` with default deltas `seq(0.1, 0.9, 0.1)`, seed 0.
Repetitions: 3 for N ≤ 2000, 2 for N=4000, 1 for N=8000.

| case | N | B=10 | B=50 | B=100 | load1 range |
|---|---|---|---|---|---|
| kid0.2 | 273 | 0.97 s | 3.98 s | 8.27 s | 4.0–4.6 |
| sim1-7 | 275 | 1.07 | 4.12 | 7.65 | 4.7–5.3 |
| kid0.1 | 487 | 1.80 | 8.85 | 17.03 | 4.9–8.6 |
| syn1000 | 1000 | 4.20 | 20.50 | 39.30 | 7.6–11.9 |
| syn2000 | 2000 | 6.34 | 29.98 | 61.08 | 9.9–12.9 |
| syn4000 | 4000 | 7.90 | 36.28 | 72.69 | 8.1–9.8 |
| syn8000 | 8000 | 10.24 | 47.26 | 93.80 | 8.2–10.4 |

Values are median elapsed times. The p-values were sensible: kidney p=0; synthetic and sim cases p_X, p_Y = 0.3–1.0.

### Cost model

Script: `$W/scripts/08_cost_model.R`; output in `$W/logs_08_cost_model.txt`, `$W/cost_model_per_case.csv` and `$W/cost_model_extrapolation.csv`.

**Per case:** elapsed = c0 + c1·B, where c1 is one permutation in both directions (9 locfit fits and 18 variograms per direction).

| N | c1 (ms per permutation, both directions) | R² | **per gene, B=100** | **per gene, B=1000** |
|---|---|---|---|---|
| 273 | 78.9 | 0.995 | 8.0 s | 79 s (1.3 min) |
| 275 | 73.2 | 0.996 | 7.7 s | 74 s |
| 487 | 170.0 | 0.994 | 17.2 s | 170 s (2.8 min) |
| 1000 | 382.4 | 0.998 | 39.1 s | 383 s (6.4 min) |
| 2000 | 610.2 | 0.999 | 61.1 s | 610 s (10.2 min) |
| 4000 | 720.2 | 1.000 | 72.6 s | 721 s (12.0 min) |
| 8000 | 928.4 | 1.000 | 93.8 s | 929 s (15.5 min) |

c0 (fixed cost per call) is 0.06–0.92 s, median 0.41 s. It includes about 0.25 s for the full `gc()` triggered by `MulticoreParam()`; see §4.5.

**Global model** (ms per permutation, both directions):

c1(N) = 51.0 + 0.0530·N + 2.83e-4·min(N,1000)² + 171.6·[N > 1000]

* The N term is locfit, which always runs on all N points.
* The min(N,1000)² term is `geoR::variog`, which runs on a subsample of n_s = min(N, 1000) points (R/spatialCorrelation.R:253-264).
* The step term is real. The 1000-point subsample is in random order, which makes geoR's C pair loop about 60% slower than on grid-ordered points of the same size
  (6.7 ms vs 10.4 ms, `$W/logs_07_algorithmic.txt` (4)), most likely through branch misprediction.

The model fits N ≥ 1000 within 0.5% and N < 1000 within about ±15%.

Predicted cost per gene from the global model:

| N | B=100 | B=1000 |
|---|---|---|
| 3000 | 67 s | 11.1 min |
| 6000 | 83 s | 13.7 min |
| 10000 | 104 s | 17.3 min |
| 20000 | 157 s | 26 min |

**With the AKI vignette's delta grid** (`c(0.01, 0.05, seq(0.1, 0.9, 0.1))`), the small deltas make locfit much slower.
At N=2000, one fit takes 130 ms at δ=0.01, 35 ms at δ=0.05, and 12.7 ms at δ=0.1 (`$W/micro_locfit_small_delta.csv`).
At δ=0.01 locfit builds 1089 vertices, against 81 at δ=0.1. The extra two deltas add about 165 ms of locfit plus about 42 ms of variograms per permutation per direction,
so each permutation costs about 1.7× the default-grid cost at N=2000.

## 4. Where the time goes

### 4.1 Micro-benchmark breakdown of one permutation (one direction, 9 deltas)

Scripts: `$W/scripts/02_micro_components.R` and `02c_decompose.R`. Output: `$W/micro_components.csv`, `$W/micro_decomposition_per_permutation.csv`.
Values are medians in ms. In R/spatialCorrelation.R:107-132, each delta runs 1 `locfit` fit, 1 `fitted`, 2 `geoR::variog` calls, 1 `lm` and 1 `rnorm`.

| N | locfit fit ×9 | fitted ×9 | variog ×18 (of which `dist` alone) | lm ×9 | measured `matchingVariograms` | variog % | locfit % |
|---|---|---|---|---|---|---|---|
| 273 | 13.5 | 4.4 | 13.7 (3.2) | 2.6 | 38.1 | 36% | 47% |
| 487 | 18.4 | 4.2 | 43.2 (8.7) | 2.6 | 78.7 | 55% | 29% |
| 1000 | 31.6 | 6.6 | 116.9 (37.3) | 2.5 | 175.8 | 66% | 22% |
| 2000 | 55.1 | 9.7 | 190.5 (37.1) | 2.8 | 286.3 | 67% | 23% |
| 4000 | 105.7 | 17.5 | 193.2 (40.2) | 2.5 | 334.5 | 58% | 37% |
| 8000 | 201.5 | 30.3 | 197.4 (37.3) | 2.7 | 455.7 | 43% | 51% |

* A single `geoR::variog` call costs 0.76 ms at N=273, 6.5 ms at N=1000 (grid order), and 10.6–11.0 ms on the random 1000-point subsample.
  It allocates 13 MB per call at n=1000.
* A single locfit fit at δ = 0.1 → 0.9 costs 3.4 → 1.0 ms at N=273, 6.9 → 2.1 ms at N=1000, and 48 → 12 ms at N=8000.
* `lm` on about 13 points costs 0.28 ms; `.lm.fit` costs 0.004 ms.
* `rnorm(N)` and `sample(N)` cost under 0.25 ms even at N=8000.

### 4.2 Rprof of the representative run (N=1000 synthetic, B=20, nThreads=1, interval 0.005)

Scripts: `$W/scripts/04_rprof.R` and `04b_attribute_rprof.R`.
Files: `$W/rprof_syn1000_B20.out`, `$W/rprof_syn1000_B20_bytotal.csv` / `_byself.csv` / `_attribution.csv`; full printout in `$W/logs_04_rprof_syn1000_B20.txt`.
Elapsed was 8.26 s, with 7.49 s sampled. Because `nThreads = 1`, BiocParallel 1.44 converts the param to an in-process `SerialParam`
(`BiocParallel:::.bpinit`: `if (bpnworkers(BPPARAM) <= 1L) BPPARAM <- as(BPPARAM, "SerialParam")`), so the profile captures all the work.

`summaryRprof` **by.total**, selected rows (the wrappers system.time, spatialCorrelation, tryCatch, viladomatCorrelation, bplapply and matchingVariograms are all 92–100%):

| function | total s | total % | self % |
|---|---|---|---|
| geoR::variog | 4.625 | 61.75 | 12.08 |
| as.vector (wraps `dist`) | 2.340 | 31.24 | 5.14 |
| dist | 1.965 | 26.23 | 25.63 |
| .C (binit, slocfit and sfitted combined) | 1.635 | 21.83 | 21.83 |
| locfit::locfit | 1.295 | 17.29 | 7.14 |
| unlist / array (wrap `binit` inside variog) | 1.150 | 15.35 | 0.2 |
| model.frame.default | 0.630 | 8.41 | 6.01 |
| lfproc (locfit.raw) | 0.630 | 8.41 | 0.33 |
| lm | 0.590 | 7.88 | 0.00 |
| gc (explicit calls) | 0.465 | 6.21 | 6.21 |
| fitted / fitted.locfit | 0.350 | 4.67 | 0.27 |
| locfit.matrix | 0.265 | 3.54 | 2.47 |
| BiocParallel::MulticoreParam | 0.220 | 2.94 | 0 |
| apply (variog `apply(data, 2, var)`) | 0.140 | 1.87 | 1.67 |
| na.omit.data.frame (inside model.frame) | 0.140 | 1.87 | 1.13 |
| quantile | 0.035 | 0.47 | 0 |
| sample / rnorm | 0.005 each | 0.07 each | – |
| cor, cor.test, output data.frame | not sampled | < 0.07 | – |

**by.self**, top entries: dist 25.63%, .C 21.83%, geoR::variog 12.08%, locfit::locfit 7.14%, gc 6.21%, model.frame.default 6.01%,
as.vector 5.14%, locfit.matrix 2.47%, apply 1.67%, matrix 1.40%, na.omit.data.frame 1.13%.

Of the explicit `gc` time, 3.3 percentage points are a measurement artifact (`system.time(gcFirst = TRUE)`).
The other 2.9 points come from `BiocParallel::MulticoreParam()` → `.snowCoresMax()` → `showConnections()`, which calls `gc()`.

### 4.3 Each sample assigned to exactly one category, at three problem sizes

Scripts: `$W/scripts/04b_attribute_rprof.R` and `04d_combine_attribution.R`; output `$W/rprof_attribution_combined.csv`.
Percentages exclude the `system.time` gcFirst artifact.

| category | N=273, B=50 | **N=1000, B=20** | N=4000, B=10 |
|---|---|---|---|
| variog: `dist()` of all pairs (used only for min/max), incl. the GC its 4 MB allocations trigger | 6.8 | **26.4** | 6.6 |
| variog: other R-level work (`u[u < max.dist]`, `apply(var)`, `array`/`unlist`, `as.matrix`) | 18.5 | **22.2** | 26.3 |
| variog: `.C binit` (the actual binning) | 4.2 | **15.3** | 25.4 |
| locfit fit + fitted: R overhead (formula, model.frame, locfit.matrix) | 33.2 | **15.4** | 11.0 |
| locfit fit: C (slocfit) | 12.2 | **6.5** | 19.4 |
| locfit fitted: C (sfitted) | 1.5 | **0.8** | 1.1 |
| lm() on about 13 points | 15.2 | **8.1** | 4.5 |
| BiocParallel: `MulticoreParam()` → `showConnections()` → full gc | 5.3 | **3.0** | 2.9 |
| BiocParallel: other (SerialParam loop) | 1.3 | **0.6** | 0.8 |
| everything else (sample, rnorm, cor, cor.test, setup dist/quantile, data.frame) | 1.8 | **1.6** | 1.9 |
| **C kernels doing real numerics (binit + slocfit + sfitted)** | **17.9** | **22.6** | **45.9** |

* The "dist" share swings between 7% and 26% because R charges an allocation-triggered GC to whichever function allocated.
* In isolation, `dist()` is 20–32% of a single `variog` call: 2.07 of 6.5 ms at N=1000, and 2.06 of 10.6 ms at N=2000.

### 4.4 Garbage collection and system time

`$W/scripts/04c_gc_time.R`, `$W/logs_04c_gc_time.txt` (B=20, 1 thread):

| N | elapsed | GC time | GC share |
|---|---|---|---|
| 273 | 1.75 s | 0.56 s | 32% |
| 1000 | 7.80 s | 1.69 s | 22% |
| 4000 | 15.02 s | 2.13 s | 14% |

* System time is 6–18% of elapsed. Each `variog` call allocates 13 MB of short-lived vectors.
* The heap after loading the package's imports holds 7.4M cons cells (~400 MB) with plain `library()`, and 7.8M with `load_all`. A full `gc()` costs about 0.25 s in both, so these numbers represent a normal user session.
* The same 20 `matchingVariograms` calls ran 7–11% faster in a forked child than in the parent
  (parent sys time 0.57–0.81 s vs child 0.29–0.33 s; `$W/logs_05b_parent_vs_child.txt`). This explains part of the superlinear 2-thread speedup in §5.

### 4.5 What each hot spot is doing, with source references

* **geoR::variog** (called at R/spatialCorrelation.R:114-117 and :125-128; geoR 1.9-6 source in `$W/ext_src/geoR/`):
  * The R body runs `u <- as.vector(dist(as.matrix(coords)))`, computing all n(n−1)/2 distances (499,500 at n=1000).
    It uses them only for `min(u)` and `umax <- max(u[u < max.dist])`, building 4 MB doubles plus a logical vector of the same length.
  * It then calls `.C("binit")`, which recomputes every distance with `hypot` and bins the pairs (geoR/src/geoR.c:370-417).
  * Also run on every call: `apply(data, 2, var)`, `trend.spatial`, and `array(unlist(lapply(as.data.frame(data), bin.f)))`.
  * **Coordinates, `ids` and `max.dist` are the same for all 3,602 variog calls per gene, and for every gene in a dataset.** Only the data values change.
  * Binning semantics a port must reproduce: geoR 1.9-6 uses left-closed bins `[lims[k], lims[k+1])` (geoR.c:396-397).
    The last limit is `umax`, which is itself an exact pair distance, so pairs at exactly `umax`, and pairs with `umax < d ≤ max.dist`, are always dropped.
    On the grid I tested, 1,225 pairs sit exactly at umax and 1,272 lie in (umax, max.dist].
    The default `uvec = 13` gives 13 bins. Bins with fewer than 2 pairs are dropped: 10–12 bins remained on the small hex and square grids I used, 13 on the larger ones.
  * A pure-R version with a precomputed pair-to-bin table (`findInterval(d, lims)`, then `rowsum((z[I]-z[J])^2, bin) / (2 n_b)`)
    reproduces geoR to 8.4e-15 with identical per-bin counts. It runs in 2.1 ms against 11.1 ms for geoR (5.2×) and allocates 3.8 MB against 13 MB
    (`$W/scripts/07_algorithmic_checks.R` (2)). In C++ the pair loop is about 122k pairs, roughly 0.1 ms.
* **locfit** (R/spatialCorrelation.R:109-112):
  * Each call goes through formula parsing, `model.frame` and `lp()`, `locfit.raw` setup, and `fitted.locfit` → `locfit.matrix` → `model.frame` again.
    At N ≤ 1000 this R overhead is larger than the C fit itself.
  * **The smoother is linear in the response for fixed coordinates and δ:**
    max |S(a·y1 + b·y2) − (a·S y1 + b·S y2)| = 1.4e-14 (N=273) and 7.1e-15 (N=1000, 2000) across all 9 deltas (`$W/logs_07_algorithmic.txt` (1)).
  * So X.delta = S_δ · X.randomized, where S_δ depends only on (pos, δ). It can be precomputed once per pixel set and applied to all B permutations of all genes.
    The kd-tree has only 15–81 vertices for δ ≥ 0.1, so S_δ ≈ H (N×V interpolation) × W (V×N kernel weights).
    Exact equivalence with locfit's tree interpolation would require either reimplementing locfit's evaluation structure
    or building S_δ column by column from locfit itself (N fits per δ, once per dataset).
* **lm** (R/spatialCorrelation.R:120): `lm(target_variog$v ~ 1 + variog.X.delta[[k]]$v)` on about 13 points takes 278 µs.
  The closed-form two-parameter least-squares solution takes 6 µs (46× faster) and gives identical coefficients.
  In the run itself, lm's share (8% at N=1000) is higher than 9 × 0.28 ms would suggest, because GC is charged to its many small allocations.
* **BiocParallel** (R/spatialCorrelation.R:238-240, 290-313 and 507-509):
  * Every `MulticoreParam()` construction calls `.snowCoresMax` → `showConnections()` → `gc()`, about 0.25 s.
  * `spatialCorrelation` constructs one per call; `spatialCorrelationGeneExp` constructs one at line 778 that is never used, then one per gene (§5).
* **Negligible:** `sample` (permutations at :286-288), `rnorm`, `cor` (:320-323), `cor.test` (:534), output `data.frame` construction (:553-577)
  and `dist`/`quantile` for max.dist (:268-269) each take ≤ 0.5%.
  The `dist`/`quantile` step is computed twice per gene with identical inputs, which is redundant but cheap.

## 5. Parallel overhead

Scripts: `$W/scripts/05_parallel.R` and `05c_bpparam_check.R`. Data: `$W/parallel_bplapply_overhead.csv`, `$W/parallel_spatialCorrelation_scaling.csv`;
logs `$W/logs_05_parallel.txt`, `$W/logs_05c_bpparam.txt`.

**Constructor cost:** `MulticoreParam(workers = 1)` 248 ms, `MulticoreParam(workers = 8)` 250 ms, `gc()` 250 ms, `showConnections(all = TRUE)` 245 ms.
Adding 800 MB of large numeric vectors to the heap left the constructor at 249 ms, because the cost is driven by the number of cons cells.

**Trivial `bplapply`** over 100 elements (median ms). The extraction loops mirror R/spatialCorrelation.R:298-313.

| N | workers | `function(i) i` | extract hat.X.delta.star (= :298-304) | extract δ* (= :306-313) | plain `lapply` |
|---|---|---|---|---|---|
| 1000 | 1 (SerialParam fallback) | 9.3 | 7.1 | 7.0 | 0.4 |
| 1000 | 4 | 47.0 | 49.7 | 49.3 | 1.1 |
| 1000 | 8 | 58.8 | 69.7 | 72.8 | 0.4 |
| 4000 | 1 | 5.6 | 7.2 | 5.8 | 1.3 |
| 4000 | 4 | 42.9 | 54.4 | 46.7 | 4.3 |
| 4000 | 8 | 81.0 | 93.5 | 87.1 | 1.5 |

Each gene makes 6 `bplapply` calls (confirmed by trace: `MulticoreParam(workers=n) x6`), 4 of them trivial extraction loops.
With 8 workers that is about 6 × 75 ms of fork overhead plus 0.25 s of constructor gc, roughly **0.7 s per gene of pure overhead**,
and about 0.3 s of it is the extraction loops, which a plain `lapply` would do in under 1.5 ms.

**spatialCorrelation scaling** (median elapsed). Results were **bit-identical to the serial run** for every nThreads (`identical_to_serial = TRUE`).

| case | 1 thread | 2 | 4 | 8 | speedup at 2/4/8 | efficiency at 2/4/8 |
|---|---|---|---|---|---|---|
| kid0.2, N=273, B=100 | 7.82 s | 3.94 | 2.62 | 2.17 | 1.99 / 2.98 / 3.60 | 0.99 / 0.74 / 0.45 |
| syn1000, B=100 | 35.19 s | 15.13 | 8.52 | 5.32 | 2.33 / 4.13 / 6.62 | 1.16 / 1.03 / 0.83 |

**Model of the forked path** (nThreads ≥ 2): T(n) = a + c·imb(n,B)/n + d·n, where imb = ceil(B/n)/(B/n) is the chunk imbalance with the default `bptasks = 0`.

| case | a (serial part per gene) | c (parallelizable work) | d (per worker) | predicted T(16) | speedup at 16 |
|---|---|---|---|---|---|
| syn1000 | 1.97 s | 26.4 s | ≈ 0 | 3.66 s | 9.6× (60% efficiency) |
| kid0.2 | 0.96 s | 5.73 s | 0.059 s | 2.30 s | 3.4× |

For kid0.2, more threads stop helping. The superlinear speedup at 2 and 4 threads comes partly from the 7–11% child-vs-parent effect (§4.4) and partly from load variation.

**BPPARAM handling, confirmed by tracing `BiocParallel::bplapply`:**

| call | what `bplapply` actually received |
|---|---|
| `spatialCorrelationGeneExp(nThreads = 2, BPPARAM = SerialParam())` | MulticoreParam(workers=2) ×6, so **BPPARAM is ignored** |
| `spatialCorrelationGeneExp(nThreads = 1, BPPARAM = MulticoreParam(3))` | MulticoreParam(workers=1) ×6, **ignored** |
| `spatialCorrelationGeneExpIterPermutations(nThreads = 2, BPPARAM = SerialParam())` | SerialParam ×6 (honored) |
| `...IterPermutations(nThreads = 22, BPPARAM = MulticoreParam())` (the vignette call, AKI Rmd:413-422 and inst/scripts) | MulticoreParam(workers=**18**) ×6: nThreads=22 is overridden by `multicoreWorkers()` = detectCores()−2 |
| `spatialCorrelation(nThreads = 2, BPPARAM = SerialParam())` | SerialParam (honored) |

The cause is R/spatialCorrelation.R:830, which passes `BPPARAM = NULL` (and lines 777-779 build a `BPPARAM` that is never used).
Separately, the AKI vignette comment at line 419 ("nThreads = 22, # parallelize genes across threads") is wrong:
genes run sequentially (`lapply` at spatialCorrelation.R:808-845 and iterativePermutations.R:25-55), and only the permutations within one gene and direction run in parallel.

## 6. spatialSimilarity scaling with G

Scripts: `$W/scripts/06_similarity.R` and `06b_similarity_accessors_bigG.R`.
Data: `$W/timings_spatialSimilarity.csv`, `$W/timings_spatialSimilarity_bigG.csv`; logs `$W/logs_06_similarity.txt`, `$W/logs_06b_similarity.txt`.

**Inputs:** G genes made from perturbed copies of the synthetic fields (log-normal multiplicative noise, about 30% zeros) on a shared N-pixel grid.
Both dense `matrix` and sparse `dgCMatrix` assays were tested. SEraster returns a dense matrix for dense input and a dgCMatrix for sparse input;
simRanPatternRasts is stored as dgCMatrix, as real count data usually is.

| G | N=1000 dense | N=2000 dense | N=1000 sparse | N=2000 sparse |
|---|---|---|---|---|
| 10 | 0.24 s | 0.22 s | 0.20 s | 0.18 s |
| 100 | 2.19 s | 2.29 s | 2.15 s | 1.86 s |
| 500 | 9.43 s | 11.79 s | 10.79 s | 14.66 s |

Values are median elapsed times (one repetition for G=500). That works out to about 18–29 ms per gene, roughly linear in G at this size.

**Larger sparse problems (N=4000):**

| G | elapsed | per gene | in `as.matrix.Matrix` | in `rbind` | in `getGenePixelDF` |
|---|---|---|---|---|---|
| 1000 | 40.6 s | 40.6 ms | 41% | 3.2% | 91% |
| 2000 | 101.5 s | 50.8 ms | 47% | 2.9% | 92% |

The cost per gene rises with G×N: that is the quadratic term.

**Hot spots** (Rprof at G=500, N=2000; `$W/rprof_similarity_G500_N2000_{dense,sparse}*`):
* `getGenePixelDF` takes 86–89% of the time (R/packageFunction.R:17-32).
  * It calls `SummarizedExperiment::assay()` twice per gene (lines 25-26).
  * It calls `colnames()` twice through `intersect(colnames(x), colnames(y))` (line 20).
  * It runs `as.matrix()` on the **whole** assay twice per gene.
* Measured accessor costs (`$W/logs_06b_similarity.txt`), dense / sparse:

  | operation | dense | sparse |
  |---|---|---|
  | `assay(x, 1)` on a SpatialExperiment | 5.9 ms | 7.1 ms |
  | `colnames(x)` | 0.3 ms | 0.5 ms |
  | `as.matrix()` of a pre-extracted assay | 0.004 ms | 1.5 ms |
  | `a[gene, pixels]` on a pre-extracted matrix | 0.06 ms | 1.3 ms |
  | `getGenePixelDF()` | 20.0 ms | 22.0 ms |

  The profile shows the accessor time inside S4 dispatch and validity checks (`updateObject` → `rowRanges` → `stopifnot`; `stopifnot` is 55–64% self time).
  **Hoisting the accessors and `as.matrix` out of the loop removes about 85–90% of the time (5–10× faster in pure R) and the quadratic term.**
  After that, the remaining per-gene cost is about 2–3 ms (`threshold`, `quantile`, `rbind`, `data.frame`). A vectorized all-genes implementation would remove most of that too.
* `rbind` growth of the output tables (lines 241-255, 278-297) takes 4.4–4.9% at G=500 and 2.9–3.2% at G=1000–2000.
  It is O(G²) in principle but has a small constant: 0.5 s at G=500 and 2.9 s at G=2000.
* `threshold()` (na.omit, setdiff, data.frame) takes 2–4%, and `quantile` about 2%.
* **Extrapolation to a whole transcriptome** (G=20000, N=4000, sparse). Per-gene densification costs 2·G·N elements at about 1.5–2.1 ns each,
  which is 0.24–0.34 s per gene and 1.3–1.9 h in total. Add about 0.1 h of accessor overhead, plus a 0.6 GB dense temporary twice per gene.
  With the matrices extracted once, the same computation is O(G·N), seconds in R and milliseconds in C++.

## 7. Wall time for realistic workloads

Computed by `$W/scripts/08_cost_model.R`; output in `$W/logs_08_cost_model.txt`.
Basis: N=2000, default deltas, single-thread cost per gene of 61.1 s at B=100 and 610.3 s at B=1000.

**Promotion fractions in the shipped results:** genes rerun at B=1000, counted from `deltaStarX` lengths with `$W/scripts/00_inspect_precomputed.R`.
* kidneyCorrelation: 757/1046 = **72%**
* merfishCorrelation: 397/483 = **82%** (affine version: 386/483 = 80%)
* brainCorrelation: 147/325 = **45%**
* ctCorrelation: 10/16 = 62%

The high rates follow from the screening rule `t <- (alpha / nPermutes) * 100` (iterativePermutations.R:64). At B=100 this is p < 0.05 in both directions,
not the documented `alpha / nPermutations[k]` (lines 90-94). Every strongly patterned gene therefore gets promoted, and round 2 dominates the total.

| workload | single thread | 16 threads (range: conservative to fitted model) |
|---|---|---|
| (a) 100 genes, N=2000, B=100 | **6,108 s = 1.70 h** | 625 s to 485 s (**8–10 min**) |
| (b) 500 genes, N=2000, `nPermutations = c(100, 1000)`, **30% promoted** (150 genes) | **122,081 s = 33.9 h** | **2.6 h to 1.8 h** (c) |
| (b) with 45% promoted (brain) | 46.6 h | 3.4 h to 2.4 h |
| (b) with 72% promoted (kidney) | 69.5 h | 4.9 h to 3.5 h |
| (b) with 82% promoted (MERFISH) | 78.0 h | 5.5 h to 3.9 h |

How the 16-thread range was computed:
* **Conservative:** T16 = a + T1·imb(16,B)/16, using a = 1.97 s per gene measured for the forked path and no superlinear credit.
  This gives 6.2 s per gene at B=100 (9.8×, 61% efficiency) and 40.4 s at B=1000 (15.1×, 94%).
* **Fitted model:** the §5 fork-path model scaled to N=2000.
* Both ignore contention. Sixteen workers plus a parent on a machine with 16 performance cores, while other jobs run, would land some workers on efficiency cores.
  Each `bplapply` waits for its slowest chunk, so real times would probably sit at or above the conservative end.
* **Plausibility check:** the brain vignette reports 1.78 h for 325 genes (45% promoted) with nThreads=22, on the authors' machine and with unknown N.
  The conservative model at 22 threads and N=2000 gives 1.68 h. This is order-of-magnitude agreement only.
* **Memory:** with `returnPermutations = TRUE` and B=1000, the results hold N×B doubles per direction, about 32 MB per gene at N=2000, or about 16 GB for 500 genes.

## 8. Implications for a C++ backend

These are estimates for the next phase, not measurements, except where marked "measured".

1. **Precompute once per pixel set, shared by all genes, permutations and deltas:**
   * the variogram subsample, max.dist, bin limits and pair list. That is about 122k pairs within max.dist at n_s=1000 (24.5% of pairs, measured).
   * each locfit smoothing operator S_δ, which is linear (measured to 1e-14).

   Per gene, compute the target variogram (O(P)).
   Per permutation and δ, compute S_δ·x (O(N·V) with V ≈ 21–81 vertices, or a GEMM over all B permutations at once), two O(P) variograms with the noise term,
   and a closed-form two-parameter least squares (6 µs in R, measured).
   That is roughly 2–3 M flops per permutation per direction against 286 ms today at N=2000, which suggests **100–300× per thread**.
   Threads (RcppParallel) would remove the 0.7–2 s per-gene fork overhead and the 0.25 s `gc()`, so parallel efficiency should be near-linear even at small N.
2. **Parallelize across genes or (gene, permutation) blocks, not just across the permutations of one gene.** At N≈275, 8 forked workers give only 3.6×.
3. **Keep reproducibility:** results are currently deterministic given `seed` and independent of nThreads. A C++ RNG (dqrng is installed) should keep both properties.
   Because δ* is a discrete argmin, exact agreement with the R implementation needs geoR-exact binning (`hypot`, `[lo, hi)`, umax dropped)
   and locfit-exact smoothing (tree interpolation). Otherwise, accept distributional equivalence.
4. **Quick pure-R wins before any C++ work:**
   * hoist accessors and `as.matrix` out of the spatialSimilarity loop (removes 85–90% of its time and the quadratic term);
   * replace the two extraction `bplapply` calls with `lapply` (saves about 0.3 s per gene with 8 workers);
   * construct `MulticoreParam` once (saves 0.25 s per gene);
   * use `.lm.fit` or the closed form instead of `lm` (about 5–8%);
   * replace `geoR::variog` with the precomputed-pair variogram (variog is 50–67% of the time, and the pure-R replacement is 5.2× faster, measured).

## 9. Side findings (correctness and usability, from reading the code)

* **The BH correction in `spatialCorrelationGeneExp` does nothing:** `p.adjust` runs inside the per-gene `lapply` on a one-row data frame
  (R/spatialCorrelation.R:836-842 inside 808-845), and adjusting n=1 p-value returns it unchanged. The IterPermutations version correctly adjusts the full vector (iterativePermutations.R:343-344).
* **`BPPARAM` is ignored** by `spatialCorrelationGeneExp` (line 830). The vignette's `nThreads=22` is overridden by `BPPARAM = MulticoreParam()`. The vignette comment says genes are parallelized; they are not.
* **Global RNG side effects:** `set.seed()` runs inside library functions (R/spatialCorrelation.R:99 and 242).
  Both directions use the same seed, so permutation i of Y uses the same index permutation and the same noise vector as permutation i of X. The two p-values are therefore not independent.
* The **doc and code disagree on the screening threshold**: documented as `alpha/nPermutations[k]`, coded as `alpha/nPermutations[k]*100`.
* `returnPermutations = TRUE` with B=1000 can produce outputs of tens of GB on large gene sets.
* **speKidney cannot provide N > ~490 shared pixels.** Large-N testing needs synthetic data (as built here) or real datasets.

## 10. Files produced

All under `$W`.

**Scripts** (`scripts/`):
* `00_inspect_precomputed.R`: promotion fractions in the shipped results.
* `01_build_inputs.R`: rasterized inputs and synthetic fields → `inputs.rds`.
* `02_micro_components.R`, `02b_locfit_small_delta.R`, `02c_decompose.R`: component micro-benchmarks.
* `03_time_spatialCorrelation.R`: N×B timings.
* `04_rprof.R`, `04b_attribute_rprof.R`, `04c_gc_time.R`, `04d_combine_attribution.R`: profiles.
* `05_parallel.R`, `05b_parent_vs_child.R`, `05c_bpparam_check.R`: parallel overhead, scaling and the BPPARAM trace.
* `06_similarity.R`, `06b_similarity_accessors_bigG.R`: spatialSimilarity.
* `07_algorithmic_checks.R`: locfit linearity, precomputed-pair variogram, lm vs closed form, point-order effect.
* `08_cost_model.R`: cost models and extrapolations.

**Raw and derived data:**
* `timings_spatialCorrelation.csv`, `timings_spatialCorrelation_medians.csv`
* `cost_model_per_case.csv`, `cost_model_extrapolation.csv`
* `micro_components.csv`, `micro_decomposition_per_permutation.csv`, `micro_locfit_small_delta.csv`
* `rprof_*.out` with their `*_bytotal.csv`, `*_byself.csv` and `*_attribution.csv`; `rprof_attribution_combined.csv`
* `parallel_bplapply_overhead.csv`, `parallel_spatialCorrelation_scaling.csv`
* `timings_spatialSimilarity.csv`, `timings_spatialSimilarity_bigG.csv`

**Other:** `logs_*.txt` (console output, including load averages); `ext_src/geoR/` (the geoR 1.9-6 CRAN source, for the `binit` code).
