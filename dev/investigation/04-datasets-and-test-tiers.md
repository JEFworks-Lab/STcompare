# STcompare: datasets used or shipped, and a tiered test-data plan for a C++ backend

Scope: every dataset STcompare ships or references, the precomputed results in `inst/extdata`, how far the published results can be reproduced, and a proposed set of test-data tiers for (i) exact-equivalence tests of a future C++ implementation, (ii) statistical calibration and (iii) realistic regression testing and benchmarking.

Repo: `/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare`, HEAD `2983c99`. Nothing in the repo was modified (`git status` is clean).
Work dir (all scripts, logs, downloads and fixtures): `/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/datasets/` (called `$WD` below).
Environment: macOS, M1 Ultra, R 4.5.2 with vecLib BLAS, locfit 1.5-9.12, geoR 1.9-6, SEraster 1.2.0, BiocParallel 1.44.0. MERINGUE is **not** installed; its `normalizeCounts()` was re-implemented verbatim from `MERINGUE/R/process.R`, and Moran's I was not needed because gene lists come from the stored results. Unless noted otherwise, timings were measured single-threaded on this machine.

---

## 0. Key findings

1. **The published results can be regenerated with the current R code.** For all three real datasets (AKI kidney Visium, MERFISH replicates, MERFISH-vs-Visium brain) the inputs can be rebuilt from public downloads (130 MB in total). The rebuilt inputs reproduce the stored `correlationCoef` to ≤ 2.6e-14. Re-running `spatialCorrelation()` with seed 0 on 1–4 threads on this Mac reproduces the stored null correlations to ≤ 1.0e-15, with identical `deltaStar`. The originals were run on Linux with 20–22 workers. The stored `.RData` files are therefore valid golden references for the current code, with the BH caveats in §2.
2. **The results do not depend on the thread count, and the two RNG streams can be pinned down exactly.**
   - Permutation order: Mersenne-Twister, from `set.seed(seed)` in the calling session (`R/spatialCorrelation.R:242,286-288`). When N > 1000, the 1000-point variogram subsample `sample(N, 1000)` is drawn first (`:253-255`).
   - Gaussian noise: drawn inside BiocParallel tasks, where BiocParallel has switched `RNGkind` to **L'Ecuyer-CMRG**. The code then calls `set.seed(seed + i)` (`:99, :290-295`). An MT-based hand replication differs (max |Δ| = 9.68); an L'Ecuyer-based one matches bit for bit.
   - `MulticoreParam(1)`, `MulticoreParam(2)` and `SerialParam` give identical permutations.
3. **Exact equivalence will be hard in two places.**
   - **Smoother:** `locfit` by default fits on an adaptive kd-tree and interpolates. For N = 2170 it uses only 62 vertices. This differs from an exact Gaussian-kernel smoother evaluated at every point (`ev = dat()`): correlation 0.95–0.99, max |Δ| 0.7–1.2 SD. A new C++ smoother therefore gives *statistically* rather than *numerically* equivalent nulls, unless locfit's C code is wrapped or ported.
   - **Variogram on lattices:** rasterized data sit on a lattice, so pairwise distances tie massively. geoR keeps only pairs with `u < max.dist` (strict). In `speKidney`, 75 pixel pairs lie exactly at `max.dist` (the 25% distance quantile). Translating the coordinates by (1000, −50) moves 50 of those pairs out of the last bin (818 → 768), which changes `deltaStar` and the nulls. Rotating by 90° changes nothing (to 5e-16), and jittered coordinates are translation-invariant. Bit-exact C++ on lattice data needs R's exact distance arithmetic, with no FMA contraction on arm64.
4. **BH correction is a no-op in `spatialCorrelationGeneExp`.** `p.adjust` runs inside the per-gene closure on a length-1 vector (`R/spatialCorrelation.R:836-842`). Evidence: in `kidneyCorrelationNoIter.RData` the stored p equals the raw empirical p for 1046/1046 genes. `spatialCorrelationGeneExpIterPermutations` does apply BH across genes (`R/iterativePermutations.R:343-344`); its stored p equals BH(raw) for every row in the kidney, brain, ct and MERFISH-affine results.
5. **`merfishCorrelation.RData` combines two runs.** Rows that were NA in the 2026-06-23 run were overwritten with rows from a 2026-06-22 run (`inst/scripts/biological-replicates-example.R:210-217`). As a result, 116/483 stored p-values match neither the raw nor the BH-adjusted values; their nulls are still internally consistent. **`brainCorrelation.RData` and `brainCorrelation_1.RData` are byte-identical** (md5 `5661fd72aa53da9ebf788764326808f0`).
6. **Crash: `locfit` segfaults ("C stack overflow") and kills R when 1 ≤ nn·N < 2.** Examples: N = 12–15 with delta 0.1; N = 100–150 with delta 0.01. With the extended delta grid (0.01, …) that the AKI and MERFISH analyses use, any comparison with fewer than 200 shared pixels crashes R; with the default grid, fewer than 20. When nn·N < 1, R raises an error that becomes an NA row instead.
7. **Clearly negative genes exist only in AKI.** AKI has 24 significantly negative genes but only 311 shared pixels. The MERFISH datasets have 1371 and 2170 pixels but no negatives (min r = −0.04). I therefore propose a realistic tier made of two real subsets plus engineered negatives Y′ = max(Y) − Y. These satisfy an exact relation: r → −r and nullX → −nullX (verified to 1e-16).
8. **Current R cost.**

   | Dataset | Shared pixels | Deltas | Per gene, B = 100, both directions, 1 thread |
   |---|---:|---:|---:|
   | simRan pair | 273 | 9 | 6.9 s |
   | AKI | 311 | 11 | 19–20 s |
   | MERFISH replicates | 1371 | 11 | 109 s |
   | Brain | 2170 | 9 | 69 s |

   The two small deltas (0.01, 0.05) cost 35–57% of runtime. Scaling in N is sub-linear: MERFISH at 1371 / 5249 / 19 088 pixels costs 94 / 177 / 457 s per gene.

   Estimated single-thread totals for the published analyses: AKI ≈ 48 CPU-h, MERFISH replicates ≈ 135 CPU-h, brain ≈ 34 CPU-h. The authors reported wall times of 1.78 h on 22 threads for brain and 16.8 h on 20 threads for MERFISH.

---

## 1. Inventory

### 1.1 Built-in data (`data/`)

Commands and output: `$WD/01_inventory_builtin.R` and `.log`.

| | `speKidney` (74.8 KB `.rda`; 947 KB in memory) | `simRanPatternRasts` (598 KB `.rda`; **50.7 MB** in memory, mostly sf geometry and cell-ID lists) |
|---|---|---|
| Structure | Named list `A, C, B` (in that order). Each element is a SpatialExperiment with 1 gene `"Gene"`, a dense `counts` assay (non-integer simulated values) and `colData$sample_id` | Unnamed list of 100 rasterized SpatialExperiments. Each has 1 gene `"1"`, a dense `pixelval` assay, and colData `num_cell, cellID_list, type (hexagon), resolution (0.2), geometry, sample_id` |
| Size | A: 1229 cells, C: 1297, B: 1242 | 272–288 pixels per field (mean 281.7) from 1201–1381 cells each; 305 distinct pixel IDs on one shared hex grid; 101 pixels common to all 100 fields; 264–278 shared pixels per pair |
| Coordinates | x ∈ [0.25, 2.61], y ∈ [0.21, 4.79] | x ∈ [−1.30, 1.20], y ∈ [−2.31, 2.37]; nearest-neighbour spacing 0.2 |
| Values | A: 0–39.8 (median 20.3); B: 0.13–39.4 (median 20.5); C: 25.5–105.6 (median 61.4) | 5.89–14.27 (median 10.0), i.e. G = W + Z + 10 |
| Generation | Not documented ("Simulated for demonstration purposes", `R/data.R:27`) and no code in the repo or its history | Recipe only, in `R/data.R:41-72`: exponential (Matérn ν = 0.5, κ = 0.1) GRF via `MASS::mvrnorm` on 5000 uniform cells, plus N(0, 0.3) noise, plus 10, then a kidney mask. Generation code is not in the repo. "N = 5000 cells" is misleading: about 1240 survive the mask |
| After `rasterizeGeneExpression(res = 0.2, mean, hex)` | A 282, C 287, B 279 pixels; **shared A–B 273, A–C 277**; 1–11 cells per pixel (mean 4.4); rasterizing takes 1.6 s | Already rasterized |
| Controls | r(A,B) = **−0.947** (naive p = 5.7e-136); r(A,C) = **+0.943**. `spatialSimilarity`: S(A,B) = 0.535, S(A,C) = 0 because C is about 3× A | 4950 unordered pairs, all independent nulls. Naive p < 0.05 for 49.8% of pairs; BH(naive) < 0.05 for 42.8% (the vignette says "~43%"). Median |r| 0.118, max 0.674 |

### 1.2 External datasets referenced by vignettes and `inst/scripts`

All URLs returned HTTP 200 (`$WD/03_url_check.log`). Sizes are `Content-Length`. Every file was downloaded to `$WD/raw/` (130 MB total). Zenodo licences and metadata come from the Zenodo API (`raw/zenodo_*.json`).

| File | Bytes | Record (licence) | Used by |
|---|---:|---|---|
| `IL3_filtered_feature_bc_matrix.h5` (AKI, 24 h post-injury) | 15 103 452 | Zenodo 19074288, concept DOI 10.5281/zenodo.17676991 (CC-BY-4.0) | AKI vignette, `visiumKidneySpatialCorrelation.R`, `KidneyNoIter.R` |
| `NL3_filtered_feature_bc_matrix.h5` (sham control) | 13 822 500 | same | same |
| `IL3_tissue_positions.csv` / `NL3_tissue_positions.csv` | 80 786 / 87 431 | same | same |
| `aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz` | 90 127 | Zenodo 19486091 (CC-BY-4.0; also contains the STalign notebook and region one-hots) | same |
| `STalign_S2R2.csv.gz` (MERFISH slice 2, replicate 2) | 12 517 332 | Zenodo 10724029, STalign paper (CC-BY-4.0) | `biological-replicates-example.R` |
| `STalign_S2R3_to_S2R2.csv.gz` (STalign + affine coordinates) | 17 110 794 | same | same |
| `STalign_S2R3_to_Visium.csv.gz` (includes `Pmatch`) | 16 177 162 | same | brain vignette and script |
| `Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz` | 45 925 294 | 10x Genomics "adult-mouse-brain-ffpe-1-standard-1-3-0"; licence not stated in the JS-rendered page (10x public data are usually CC BY 4.0; confirm before redistributing) | brain |
| `Visium_FFPE_Mouse_Brain_spatial.tar.gz` | 10 632 472 | same | brain |
| `STalign_cell_type_transcriptional_correlations.csv.gz` | 5 041 | Zenodo 19582556 (CC-BY-4.0) | brain cell-type section |
| `STalign_S2R3_cell_type_annotations.csv.gz` | 1 244 272 | same | same |
| `STalign_Visium_cell_type_annotations.csv.gz` | 264 260 | same | same |

The roxygen examples also use base R `datasets::quakes`.

### 1.3 Preprocessing, parameters and resulting sizes, per analysis

Each input was rebuilt with `$WD/05_build_aki.R`, `10_build_merfish_replicates.R`, `11_build_brain_merfish_visium.R` and `24_ct_input.R`.

| Analysis (result file) | Inputs → preprocessing | Rasterization | Genes kept | Shared pixels | Call and parameters | Reported time |
|---|---|---|---|---|---|---|
| **AKI** (`kidneyCorrelation.RData`) | IL3 32 285 genes × 2970 barcodes (2965 in tissue); NL3 32 285 × 3229 (3215 in tissue); 40 duplicated symbols. Positions are Visium *array indices*, and the CSV labels are swapped: `array_row` spans 0–127. AKI coordinates come from the STalign affine-only region alignment; both sections are rotated 90° | `rasterizeGeneExpression(res = 5, fun = "sum", square = FALSE)`. Units are index units and **anisotropic**: x′ ≈ 86.6 µm and y′ ≈ 50 µm per unit. About 10–11 spots per pixel. Then CPM | 1046 genes that are SVGs in both conditions (MERINGUE Moran's I, `p.adj == 0`, `filterDist = 10`; 1673 control SVGs, 2318 AKI SVGs) | ctrl 323 px, AKI 314 px, **311 shared** | `…IterPermutations(nPermutations = c(100, 1000), deltas = c(0.01, 0.05, 0.1…0.9), assayName = "CPM", nThreads = 22, BPPARAM = MulticoreParam(), seed = 0)`, BH | not recorded (≈ 48 CPU-h by my estimate) |
| **AKI, no iteration** (`kidneyCorrelationNoIter.RData`) | same | same | same | same | `spatialCorrelationGeneExp(nPermutations = 100, same deltas, nThreads = 22, seed = 0)`; BH requested but has no effect | — |
| **MERFISH replicates** (`merfishCorrelation.RData`) | S2R2: 84 172 cells × 649 genes (166 `Blank-*`). S2R3: 85 958 cells, coordinates `STalign_x/y`. Cell IDs are parsed as doubles (`1.00442548580637e+38`); no collisions were observed | `rasterizeGeneExpression(res = 200)` with defaults: first assay (counts), `fun = "mean"`, **square** pixels. Target 1402 px, source 1393 px; ~60 cells per pixel | 483 (Blank-* controls dropped) | **1371** | `…IterPermutations(deltas = 11 values, nThreads = 20, BPPARAM = MulticoreParam())`, default seed 0, `c(100, 1000)`, BH. **Note:** `MulticoreParam()` means `detectCores() − 2` workers, so `nThreads` is ignored | 16.8 h (`:184`); NA rows patched from a 06-22 run |
| **MERFISH, affine** (`merfishCorrelation_affine.RData`) | S2R3 `affine_x/y` | same | 483 | **1299** | same | 7.65 h (`:563`) |
| **Brain MERFISH vs Visium** (`brainCorrelation.RData`) | MERFISH S2R3→Visium: 85 958 cells, keep `Pmatch > 0.95` → 40 697 cells. Visium FFPE: 19 465 × 2264 spots (11 duplicated symbols), coordinates × `tissue_hires_scalef` (0.0729), rotated. 466 shared genes; libnorm = CPM | `res = 20` (hires px ≈ 80 µm), `fun = "mean"`, hex; then `lognorm = log10(x + 1)`. MERFISH 2627 px (15 cells/px); Visium 2204 px (~1 spot/px) | 325 genes detected in > 1% of pixels in both (SVG filter used only for reporting) | **2170** | `…IterPermutations(assayName = "lognorm", nThreads = 22, seed = 0)`, defaults: 9 deltas, `c(100, 1000)`, BH | 1.78 h (`vignettes/brain-MERFISH-10x-visium.Rmd:316`) |
| **Brain cell types** (`ctCorrelation.RData`) | Visium deconvolved proportions (2264 × 16); MERFISH one-hot labels (39 986 annotated cells with `Pmatch > 0.95`) | `res = 20`, mean, hex. **X = Visium**, Y = MERFISH | 16 cell types | **2174** | defaults, `nThreads = 22` | 6.95 min (`:646`) |
| **simRan nulls** (`simRanPatternResults.RData`) | `simRanPatternRasts`, all **9900 ordered** pairs | already rasterized | 1 | 264–278 | `…IterPermutations(c(100, 1000), nThreads = 22, BPPARAM = MulticoreParam())`, default deltas and seed. **Only `pValuePermuteX` is kept**; 0 is replaced by 0.01 | not recorded |

Reproduction of `correlationCoef` from the rebuilt inputs:

| Dataset | max \|Δr\| | Bit-identical |
|---|---:|---:|
| AKI | 5e-16 | 759/1046 |
| MERFISH (STalign) | 2.0e-14 | 11/483 |
| MERFISH (affine) | 2.6e-14 | 10/483 |
| Brain | 1.05e-15 | 249/325 |
| Cell types | 1.1e-16 | — |

---

## 2. Precomputed results in `inst/extdata` (≈ 36 MB in total)

Sources: `$WD/02_inventory_extdata.R`, `04_crosscheck_extdata.R` and their logs. Every correlation file is a data.frame with 10 columns: `correlationCoef, pValueNaive, pValuePermuteX, pValuePermuteY, deltaStarMedianX/Y, deltaStarX/Y, nullCorrelationsX/Y`. **No file stores permutations.** The nPermutations used per gene was inferred from `length(nullCorrelationsX)`; X and Y lengths always agree with each other and with `deltaStar`.

| File (bytes) | Dims | Genes at 100 / 1000 perms | Stored p vs raw empirical p | Both p < 0.05 (pos / neg) | Last change (git) |
|---|---|---|---|---|---|
| `kidneyCorrelation` (13.0 MB) | 1046 × 10 | 289 / 757 | equals BH(raw), 1046/1046 | 731 (707 / 24), matches vignette lines 453/461 | 2026-05-11/14 |
| `kidneyCorrelationNoIter` (1.77 MB) | 1046 × 10 | 1046 / 0 | **equals raw for 1046/1046 (BH no-op)** | 757 (729 / 28) | 2026-04-28 |
| `merfishCorrelation` (6.67 MB) | 483 × 10 | 86 / 397 | BH-consistent for only 367/483 (composite run) | 397 (397 / 0) | 2026-06-25 |
| `merfishCorrelation_affine` (6.53 MB) | 483 × 10 | 97 / 386 | equals BH(raw), 483/483 | 383 (383 / 0) | 2026-06-29 |
| `brain-…/brainCorrelation` (2.76 MB) | 325 × 10 | 178 / 147 | equals BH(raw), 325/325 | 128 (128 / 0) | 2026-05-26 |
| `brain-…/brainCorrelation_1` (2.76 MB) | **byte-identical duplicate** | | | | 2026-05-26 |
| `brain-…/ctCorrelation` (171 KB) | 16 × 10 | 6 / 10 | equals BH(raw), 16/16 | 10 (10 / 0) | 2026-05-26 |
| `simRanPatternResults` (163 KB) | `cors_df`, 9900 × 5 (`Var1, Var2, cors, corspv, corspv_corrected`) | 13 values come from 1000-perm reruns (p = 0.001–0.009, 0.011) | pX only | empirical p < 0.05: **3.88%**; naive 49.8%; BH(naive) 42.8% | 2026-04-24 |

Notes:

- **Delta grids.** Kidney and MERFISH use 11 deltas (0.01, 0.05, 0.1…0.9). Brain and ct use the default 9.
- **Code provenance.** The last semantic changes to `R/spatialCorrelation.R` and `R/iterativePermutations.R` are `6b27d1d` (seeds, 2026-03-22), `002b151` (BH added to GeneExp, 04-22) and `0a27a4c` (screening threshold becomes `(alpha/nPermutes)*100`, 04-28). So kidney, kidneyNoIter, brain, ct, merfish and merfish_affine were all produced by the code that is current today.
- **`simRanPatternResults` predates the threshold fix.** It was made on 04-24 with `spatialCorrelationGeneExp_test` (commits `895a232`/`400400a`), i.e. before 04-28. Its 1000-perm reruns therefore followed the old rule (`t = alpha/nPermutes = 0.0005`). It is a fine *statistical* reference but cannot be regenerated exactly.
- **Determinism across runs.** `kidneyCorrelation` (iterative, 05-11) and `kidneyCorrelationNoIter` (04-28) agree exactly:
  - identical `correlationCoef` and `pValueNaive`;
  - identical nulls for all 289 genes run at 100 perms (max |Δ| = 0);
  - for all 757 genes run at 1000 perms, the **first 100 nulls and `deltaStar` equal the 100-perm run** (max |Δ| = 0);
  - the set of genes screened into the 1000 round equals "raw pX < 0.05 and raw pY < 0.05" in NoIter (757 = 757).
- **Prefix property.** Permutation i depends only on (seed, i). The 1000-permutation round therefore recomputes the first 100 permutations it already has, about 10% wasted work.

**Can they serve as golden references for the current code?** Yes for `correlationCoef`, `pValueNaive`, nulls, `deltaStar` and raw p, in every file except `simRanPatternResults`. I verified this by recomputing genes on 1–4 threads (`06_aki_golden_check.R`, `12_golden_check_generic.R`):

| Gene (dataset) | Max \|Δ null\| |
|---|---:|
| Lactb2 (AKI) | 4.0e-16 |
| Ech1, Depp1 (AKI) | ≤ 2.5e-16 |
| Cxcr2, Oxgr1 (MERFISH) | ≤ 1.0e-15 |
| Npbwr1, Oprk1, Efemp1 (brain) | ≤ 2.2e-16 |

For BH-adjusted p, use kidney, brain, ct and merfish_affine. Do not use merfishCorrelation (composite) or kidneyNoIter (BH no-op).

---

## 3. Behaviour that the fixtures must capture

Sources: `07_rng_semantics`, `08`, `15`, `17`, `18`, `19–21` logs.

**RNG streams.**
- Permutations: MT with Rejection sampling, `sample(X, N)` ≡ `X[sample.int(N, N)]` after `set.seed(seed)`, preceded by `sample(N, 1000)` when N > 1000.
- Noise: L'Ecuyer-CMRG + Inversion, `set.seed(seed + i)`, then `rnorm(N)` once per delta, in delta order.
- Forward and reverse directions use the same seed and therefore the same index permutations (`:537-548`). Swapping X and Y swaps pX and pY exactly.
- Trap: the exported `matchingVariograms()` gives different results when called directly (MT) than inside the package.

**Metamorphic relations** (B = 20, kidney A–B):

| Transformation | Result |
|---|---|
| Swap X and Y | exact swap |
| X → 3X + 7 | invariant (Δnull ≤ 1e-15, same `deltaStar`) |
| Y → max(Y) − Y | r → −r and nullX → −nullX (2e-16), so pX identical; nullY differs by up to 0.12 (pY only statistically equal) |
| Rotate 90° | invariant (7e-16) |
| Translate (1000, −50) or scale × 37 | **not invariant** (Δnull 0.02–0.12, `deltaStar` changes), caused by lattice ties at `max.dist` (see below); jittered data are invariant |
| Reorder pixels | changes the permutation stream; results only statistically equal |
| B = 40 vs B = 20 | prefix property holds exactly |

**geoR variogram details.**
- 13 nominal bins on [0, umax], with `umax = max(u[u < max.dist])` (strict).
- Bins with fewer than 2 pairs are dropped. Lattices leave 9 bins (AKI) or 10 (kidney); the brain 1000-point subsample keeps all 13.
- In kidney, 75 pairs sit exactly at the quantile (`18_translation_cause.log`).

**locfit.**
- Default evaluation is an adaptive tree (`maxk = 300`) with interpolation; `fitted()` interpolates. Compared with `ev = dat()`: correlation 0.950–0.992, max |Δ| 0.67–1.22 SD.
- For N = 2170 the tree has 62 vertices. The exact evaluation is 10–70× slower in locfit (0.28 s vs 0.004–0.034 s).
- With delta = 0.01 and N = 311, every fit warns "Estimated rdf < 1.0".

**Edge cases** (one process per case):

| Input | Current behaviour |
|---|---|
| Constant or all-zero X | caught error ("NA/NaN/Inf in foreign function call (arg 4)") → all-NA row |
| One NA in X | r is computed but pX = pY = NA |
| Single non-zero pixel | runs |
| Duplicated coordinates | runs |
| N = 30 | runs, with warnings |
| nn·N < 1 | error → NA |
| **1 ≤ nn·N < 2** | **segfault, R session dies** |

Any error inside `spatialCorrelation()` is silently turned into an NA row (`:585-619`). That is how `merfishCorrelation` ended up with NA rows that had to be patched by hand.

**Data quirks.**
- AKI coordinates are anisotropic index units (see §1.3).
- MERFISH cell IDs are stored as doubles.
- On the brain Visium side the median gene is zero in 81% of pixels; some genes exceed 98%.

---

## 4. STexampleData and other public data

`13_stexampledata.log`, `14_st_mouseob.log`. STexampleData 1.18 offers 12 SpatialExperiments: Visium_humanDLPFC, Visium_mouseCoronal, seqFISH_mouseEmbryo, ST_mouseOB, SlideSeqV2_mouseHPC, Janesick breast cancer (Chromium, Visium, Xenium rep1, Xenium rep2), CosMx_lungCancer, MERSCOPE_ovarianCancer and STARmapPLUS_mouseBrain.

All of them are single samples. The only natural pairs are the Janesick Xenium rep1/rep2 and the Visium-vs-Xenium block, and those would first need alignment (STalign) and rasterization. That makes them less convenient than STcompare's own data, which already ship STalign coordinates. `ST_mouseOB` is tiny (15 928 genes × 262 spots, 1.9 MB download, coordinates 7.9–28.0 × 9.0–24.0). It could serve a within-sample test (`spatialCorrelationGeneExpWithinSample`) or a "noisy copy of itself" positive control. **Recommendation:** use the package's own datasets; their inputs are CC-BY-4.0 on Zenodo, and the 10x dataset is public.

---

## 5. Proposed test-data tiers

All four tiers are prototyped in `$WD/build_fixtures.R`; §6 has the build results. Fixtures are plain base-R lists and matrices. They avoid SpatialExperiment, whose sf geometry makes `simRanPatternRasts` 50 MB in memory. That keeps C++ tests light; a tiny helper can rebuild SpatialExperiment objects for end-to-end tests.

### Tier 0 — tiny deterministic "kernel" fixture (exact equivalence)

- **What it contains.** Five cases, each with inputs (X, Y, pos, pixel IDs, delta grid, `maxDistPrctile`, seed, B), a plain-English RNG description, and the following intermediates from the current algorithm:
  - `ids`, `prctile`, and the target variogram (`u, v, n, bins.lim`);
  - the N × B permutation-index matrix;
  - for permutation 1, per delta: locfit tree fit, **exact `ev = dat()` fit**, variogram of the fit, `lm` coefficients, the noise vector, rescaled field, its variogram, RSS, and `delta*`;
  - expected permutations, nulls, `deltaStar` and p for both directions, plus the `spatialCorrelation()` row.

  The fixture also stores `spatialSimilarity` outputs for A–B and A–C, and the edge-case behaviour table from §3.

- **Cases.**
  1. `kidney_AB`: hex lattice, N = 273, negative control r = −0.947.
  2. `kidney_AC`: N = 277, positive control r = +0.943.
  3. `kidney_AB_jitter`: U(±0.02) jitter, no ties, so exact agreement is well-posed.
  4. `quakes_irregular`: base R, 300 irregular points; r = −0.23, pX = 0, pY = 0.4 at B = 10.
  5. `brain_Oprk1_subsample`: N = 2170 > 1000, exercising the subsample path; zero-inflated Y; B = 5; summary intermediates only.

- **Recipe.** `capture_viladomat()` / `make_case()` in `build_fixtures.R` replay `viladomatCorrelation()` step by step. A built-in self-check asserts bit-identity with the package (all TRUE). Essentials:
  ```r
  set.seed(seed); ids <- if (N > 1000) sample(N, 1000) else seq_len(N)
  prctile <- quantile(dist(pos[ids, ]), 0.25)
  target  <- geoR::variog(data = X[ids], coords = pos[ids, 2:1], max.dist = prctile, option = "bin", messages = FALSE)
  perm_index <- vapply(1:B, function(i) sample.int(N, N), integer(N))
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed + i)   # then, per delta: locfit → variog → lm → rnorm(N) → rescale → variog → RSS
  ```
- **Size.** 772 KB (xz). Can be trimmed to about 300–400 KB by keeping full intermediates only for `kidney_AB_jitter` and `quakes`.
- **Runtime (current R).** 26 s to build. The template tests (`$WD/example_test_fixtures.R`, 56 expectations, all passing) take **15 s**.
- **Location.** `tests/testthat/fixtures/kernel_fixture.rds`, with the builder in `data-raw/` (`.Rbuildignore`d).
- **What it validates.** RNG streams, binning (including ties), smoother, rescaling, argmin over delta, null correlations, the p-value formula (`>`, no +1), swap symmetry and `spatialSimilarity`.

**Acceptance criteria for a C++ port.**
- Bit-compatible mode (wrap or port locfit's tree, inject or reproduce R's RNG): equality to about 1e-12 on all cases. Require identical `n` per bin on the lattice cases; that needs R-identical distance arithmetic and `-ffp-contract=off`.
- New-smoother mode: compare the smoother against `fitted_exact_evdat` and validate everything else via Tiers 1–2.

### Tier 1 — null calibration and power

- **Data.** `simRanPatternRasts` is already shipped. `calibration_fixture.rds` (317 KB) adds:
  - a 100 × 305 field matrix (NA where a field lacks a pixel), coordinates and `num_cell`;
  - a table of all 4950 pairs (shared-pixel count, naive r, naive p);
  - the shipped reference (`cors_df`; empirical rate 3.88%);
  - `mix(f_i, f_j, rho) = rho·(f_i − 10) + sqrt(1 − rho²)·(f_j − 10) + 10`. This produces positive controls with known population correlation and the *same* covariance model.
- **Recipe.** `build_fixtures.R tier1`; run with `$WD/22_calibration_subset.R`.
- **Runtime (current R).** 6.9 s per pair (B = 100, 1 thread). All 4950 unordered pairs ≈ **9.5 CPU-h**. The published script runs all 9900 ordered pairs, which is 2× redundant because (j, i) is an exact swap of (i, j), and keeps only pX. 100 null + 50 power pairs took 17 CPU-min (4.4 min on 4 cores).
- **Reference numbers, current code** (`22_calibration_subset.log`):

  | Setting | Naive p < 0.05 | pX < 0.05 | pY < 0.05 | max(pX, pY) < 0.05 |
  |---|---:|---:|---:|---:|
  | Null, 100 pairs | 46% | 2% | 3% | 0% |
  | ρ = 0.2, 25 pairs (mean r 0.25) | — | 44% | — | 28% |
  | ρ = 0.4, 25 pairs (mean r 0.44) | — | 72% | — | 68% |

  pX deciles under the null: 0.14, 0.20, …, 0.89, i.e. close to uniform and slightly conservative. The shipped full run has 3.88%.
- **Acceptance.**
  - Null rejection rate at α = 0.05 at most the binomial upper bound (0.056 for n = 4950).
  - The R and C++ p-value distributions agree (two-sample KS on the same pairs, or a binned χ²).
  - Power at ρ ∈ {0.2, 0.4} within Monte Carlo error of R.
  - `speKidney` A–B and A–C significant.
- **Location.** Compact fixture in `tests/testthat/fixtures/`. The full 4950-pair run is an opt-in slow test (for example `skip_if_not(Sys.getenv("STCOMPARE_SLOW") == "true")`) or a nightly CI job. It becomes cheap once C++ is 50–100× faster (about 6–12 CPU-min).

### Tier 2 — realistic regression and benchmark subset (`realistic_fixture.rds`, 1.16 MB; 554 KB trimmed)

The fixture has two real pairs, each with gene classes, parameters and golden values copied from `inst/extdata` (r, naive p, published final p, all stored nulls, `deltaStar`).

| Pair | Genes × shared pixels | Classes (selection rule) | Deltas | Golden |
|---|---|---|---|---|
| **AKI IL3 vs NL3** (res 5, CPM; X = control, Y = AKI) | **35 × 311** | 10 positive (r 0.40–0.80, the 5 strongest plus a spread); **10 negative** (r −0.35 to −0.47, the most negative of the 24 significant); 10 null (raw p100 > 0.87, \|r\| ≤ 0.012, some zero-inflated, e.g. Upk2 with 94% zeros); 5 borderline (screened into the 1000 round, final BH p ≈ 0.05) | 11 | `kidneyCorrelation` / `NoIter` |
| **Brain MERFISH vs Visium** (res 20, lognorm) | **30 × 2170** (N > 1000, subsample path) | 10 positive (r 0.095–0.57); 10 null (\|r\| ≤ 0.003); 5 sparse-Y (Visium ≥ 98.6% zeros); 5 borderline; plus **5 engineered negatives** Y′ = max(Y) − Y of the top positives | 9 | `brainCorrelation` |

- **Recipe.** Run `05_build_aki.R` and `11_build_brain_merfish_visium.R` (download, align/rotate, rasterize, normalize, exactly as in the scripts), then `build_fixtures.R tier2`. Selection is deterministic from the stored results; the code is in the report's builder.
- **Verified.** Re-running (1 thread, B = 100) matches golden nulls to ≤ 2.5e-16 (Ech1, Depp1, Efemp1). For the engineered negative Slc17a7, r goes 0.5686 → −0.5686 and nullX(flip) = −nullX to 1.1e-16.
- **Runtime (current R, 1 thread).** About 20 s per gene for AKI (35 genes ≈ 12 min at B = 100) and about 69 s per gene for brain (30 genes ≈ 35 min). The full iterative protocol adds roughly 10× per screened gene (AKI 25 genes: +83 min; brain about 20: +3.8 h).
- **Location.** `tests/testthat/fixtures/` (trimmed) for `skip_on_cran()` regression tests. The same file can drive a fast end-to-end vignette and `bench/` scripts.
- **Acceptance.**
  - R path: nulls equal golden to ≤ 1e-12.
  - C++ path: r and naive p exact; p within Monte Carlo error; every positive and negative gene significant and every null gene not (p > 0.2); pX(flip) = pX.

### Tier 3 — full-scale benchmarks (downloaded on demand and cached)

- **Data.** Full AKI (1046 × 311), MERFISH replicates (483 × 1371 STalign, 483 × 1299 affine), brain (325 × 2170), cell types (16 × 2174). All have golden references (§2, with the BH caveats). There is also a scaling ladder from MERFISH at 200 / 100 / 50 µm (1371 / 5249 / 19 088 shared pixels; `25_merfish_scaling.log`).
- **Recipe.** The builder scripts above. Cache inputs with `tools::R_user_dir("STcompare", "cache")` or BiocFileCache. Downloads total 130 MB; the rasterized RDS files are 38 MB (AKI, all genes), 7.5 MB and 9.3 MB.
- **Runtime (current R, single thread).** AKI ≈ 48 CPU-h, MERFISH ≈ 135 CPU-h, brain ≈ 34 CPU-h. Rasterization takes 9–112 s per dataset.
- **Location.** Not in the package. A `bench/` or `data-raw/` script plus a cache directory, run manually or in nightly CI.

---

## 6. Prototype build: what was run and what it produced

```
Rscript build_fixtures.R tier0        # 26 s
Rscript build_fixtures.R tier1 tier2  # 10 s (from the rebuilt RDS inputs)
Rscript build_fixtures.R verify       # 2.4 min, re-runs current code on fixture genes
Rscript -e 'testthat::test_file("example_test_fixtures.R")'   # 15 s, 56 expectations pass
```

| Fixture | Contents | Size on disk |
|---|---|---|
| `fixtures/kernel_fixture.rds` | 5 cases, similarity outputs, edge-case table | 772 KB (per case 158–189 KB) |
| `fixtures/calibration_fixture.rds` | fields 100 × 305, 4950 pairs, reference, `mix()` | 317 KB |
| `fixtures/realistic_fixture.rds` | AKI 35 × 311, brain 30 × 2170 (+ 5 engineered negatives), golden values, selection tables | 1156 KB (data 97 KB + 318 KB; golden 405 KB + 274 KB); 554 KB trimmed to first-100 nulls |

---

## 7. Issues found along the way (relevant to tests and documentation)

1. **BH is a no-op** in `spatialCorrelationGeneExp` (`R/spatialCorrelation.R:836-842`).
2. **User `BPPARAM` is ignored** by `spatialCorrelationGeneExp` (`:830` passes `BPPARAM = NULL`). The iterative function honours it, so scripts that pass `MulticoreParam()` actually run `detectCores() − 2` workers whatever `nThreads` says (18 on this machine).
3. **Screening threshold doc vs code.** The docs say `alpha / nPermutations[k]` (`R/iterativePermutations.R:90-94`); the code uses `(alpha / nPermutes) * 100` (`:64`).
4. **Errors become silent NA rows** (`:585-619`), and the locfit segfault for 1 ≤ nn·N < 2 kills the session. The C++ should validate `nn*N >= 2` and report errors explicitly.
5. **`spatialSimilarity` issues.**
   - The NA branch records `numPixelInThresh = dim(thresh)[1]`, which is always 1 (`R/packageFunction.R:249`).
   - `getGenePixelDF` converts the whole assay to dense for every gene (`:25-26`).
6. **Vignette vs script.** The vignette says `corspv_corrected` is the higher of pX and pY (`vignettes/getting-started-with-STcompare.Rmd:332`), but the script returns pX only (`inst/scripts/simRanPatternSpatialCorrelation.R:29`). The script also computes both (i, j) and (j, i).
7. **`R/data.R` gaps.** `speKidney` has no documented generation and no code. `simRanPatternRasts` has a recipe only, and its "N = 5000" overstates the ~1240 cells that remain.
8. **Packaging.**
   - `inst/extdata` is about 36 MB, including the duplicate `brainCorrelation_1.RData`. Consider hosting results on Zenodo and loading them on demand.
   - The MERFISH script uses hard-coded `~/ST_compare/...` paths (`biological-replicates-example.R:155,186,212,225`).
   - The comments in `inspect_kidney.R` are stale ("360 → 713", "5 → 26").
9. **Inconsistent counts.** The brain vignette says 126/226 SVGs and 89/99 non-SVGs (line 360) while code comments say 230/95. The MERFISH script says 415 SVGs (`:292`) while the `length(svg_int)` comment says 489 (`:162`); these depend on MERINGUE and could not be re-checked here.
10. **AKI coordinates are anisotropic index units.** Physical scaling (86.6 µm and 50 µm per unit) would change the published results, so the fixtures keep the published coordinates.

---

## 8. Files produced

All paths are under `$WD = /private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/datasets/`.

- **Prototype builder:** `build_fixtures.R`, with logs `build_fixtures_tier0.log`, `build_fixtures_tier12.log`, `build_fixtures_verify.log`.
- **Fixtures:** `fixtures/kernel_fixture.rds`, `fixtures/calibration_fixture.rds`, `fixtures/realistic_fixture.rds`.
- **Test template:** `example_test_fixtures.R` (+ `.log`).
- **Inventory and cross-checks:** `01_inventory_builtin.R`, `02_inventory_extdata.R`, `03_url_check.log`, `04_crosscheck_extdata.R`.
- **Input rebuilds:** `05_build_aki.R` → `aki_rast_full.rds`; `10_build_merfish_replicates.R` → `merfish_replicates_rast.rds`; `11_build_brain_merfish_visium.R` → `brain_merfish_visium_rast.rds`; `24_ct_input.R`.
- **Golden checks:** `06_aki_golden_check.R`, `12_golden_check_generic.R`.
- **Behaviour probes:** `07_rng_semantics.R`, `08_aki_warnings.R`, `09_geoR_variog_src.R`, `15_locfit_eval.R`, `16_aki_resolution_scan.R` (AKI at res 4 / 3 / 2.5 / 2 → 469 / 819 / 1165 / 1768 px; negatives weaken to median r −0.27 / −0.22 / −0.16 / −0.10), `17_metamorphic.R`, `18_translation_cause.R`, `19_edge_cases.R`, `20_small_n_crash.R`, `21_locfit_crash_threshold.R`.
- **Calibration and benchmarking:** `22_calibration_subset.R` → `calibration_subset_results.rds`; `23_fixture_sizes.R`; `25_merfish_scaling.R`.
- **Public data surveys:** `13_stexampledata.R`, `14_st_mouseob.R`.
- **Downloads:** `raw/` (130 MB of inputs, Zenodo metadata JSON, `MERINGUE_process.R`, extracted Visium tarballs).
