# STcompare test data

This directory rebuilds the tiered test data in `tests/testthat/fixtures/` and the larger inputs they are cut
from. The fixtures pin down the **current (legacy, pure R) behaviour** of `spatialCorrelation()`,
`viladomatCorrelation()`, `matchingVariograms()`, `spatialCorrelationGeneExp()` and `spatialSimilarity()`. They
are the acceptance harness for a future compiled backend. Every step of the null-generating algorithm is
stored, so a replacement can be checked one kernel at a time.

`data-raw/` is listed in `.Rbuildignore`. Downloads, caches and rasterized inputs are **never** written to the
repository. They go to a cache directory (see below), and the scripts refuse any cache path inside the
repository, including one reached through a symbolic link, before creating anything. The licences and
attribution of the data in the fixtures ship with the tests, in `tests/testthat/fixtures/README.md`.

## Tiers

| Tier | Fixture (xz RDS) | Contents | What it validates | Test file (default runtime) |
|---|---|---|---|---|
| 0 kernel | `kernel_fixture.rds` (305 KB) | Seven small cases (below). Each has its inputs, the variogram subsample `ids`, `max.dist`, the target variogram (`u`, `v`, `n`, `bins.lim`), the N x B permutation-index matrix, the RSS of every delta for every permutation and the arg-min per permutation. For permutation 1 it also stores the per-delta smoother variogram, `lm` coefficients, rescaled-field variogram and RSS. **Complete per-delta vectors** (locfit fit, exact `ev = dat()` fit, noise, rescaled field) are stored for permutation 1, forward direction, of `kidney_AB`, `kidney_AB_jitter` and `quakes_irregular`. In both directions: expected nulls, `deltaStar`, raw p-values, tail counts, a column fingerprint (mean, sd, min, max) of the permuted fields, and the first permuted field when N <= 1000. Also `spatialSimilarity()` outputs for speKidney A-B, A-C and A-C with `foldChange = 2`, and a table of edge-case behaviour | RNG streams; the `max.dist` quantile; geoR binning; the locfit smoother; least squares; rescaling and noise; RSS and arg-min over the delta grid; nulls; the p-value definition; both directions; argument plumbing (seed, `maxDistPrctile`, `deltaX` != `deltaY`, defaults); the N = 1000 subsample boundary; worker-count independence; the `spatialCorrelationGeneExp()` wrapper; `spatialSimilarity()` | `test-reference-kernel.R` (21 s) |
| 1 calibration | `calibration_fixture.rds` (331 KB) | `data(simRanPatternRasts)` as plain matrices: 100 independent simulated fields x 305 hex pixels (NA where absent), coordinates and cells per pixel. All 4950 unordered pairs with shared-pixel count, naive r and p. The shipped reference `inst/extdata/simRanPatternResults.RData`. The mixing recipe (`mix_mu`, `mix_recipe`, and a check vector for the helper `stc_mix()`). 60 test jobs: 40 disjoint null pairs, then the first 20 mixed to rho = 0.6. For these jobs, the current code's results at B = 100: r, naive p, raw pX/pY, tail counts, deltaStar medians, and the 100 x 60 deltaStar matrices of each direction | Type-I error at alpha = 0.05 (one-sided 99.9% binomial bound); power at rho = 0.6; distributional agreement with the stored reference (backend-agnostic); exact reproduction (legacy backend) | `test-calibration.R` (skipped unless `STCOMPARE_SLOW_TESTS=true`; 28 s on 18 workers, about 6.5 CPU-min) |
| 2 realistic | `realistic_fixture.rds` (536 KB) | AKI kidney Visium: 35 genes x 311 pixels, 11 deltas (10 positive, 10 negative, 10 null, 5 borderline). MERFISH-vs-Visium brain: 30 genes x 2170 pixels, 9 deltas (10 positive, 10 null, 5 sparse on the Visium side, 5 borderline), plus 5 engineered negatives defined by a recipe. Golden values copied from `inst/extdata`: r, naive p, the published (BH-adjusted) p for reference, the **first 100 nulls and deltaStar** of each direction, and the selection table rows | Real data end to end: zero inflation, the N > 1000 variogram subsample, small deltas (0.01, 0.05), the prefix property, and the metamorphic relation Y -> max(Y) - Y (r -> -r, nullX -> -nullX) | `test-reference-realistic.R` (28 s; `skip_on_cran()` except AKI Upk2; all 65 genes at B = 100 when `STCOMPARE_SLOW_TESTS=true`, about 160 s on 18 workers) |
| 3 benchmark | not committed | Full rasterized inputs in the cache: AKI (1046 published genes, all 32 285 rasterized), brain (325 genes), brain cell types (16), MERFISH replicates (483 genes; STalign 1371 and affine 1299 shared pixels) | Full-scale regression against every published result, and benchmarking | none (manual or nightly) |

Tier-0 cases (B = 10, deltas 0.1..0.9, `maxDistPrctile = 0.25` and seed 0 unless noted):

| Case | N | What it exercises |
|---|---:|---|
| `kidney_AB` | 273 | `data(speKidney)` A vs B, rasterized with SEraster (resolution 0.2, mean, hexagons). Negative control, r = -0.947. Hexagonal lattice: 75 pixel pairs lie exactly at `max.dist`, and 448 within 4 ulp of geoR's `umax` |
| `kidney_AC` | 277 | A vs C, positive control, r = +0.943 (summary intermediates only) |
| `kidney_AB_jitter` | 273 | `kidney_AB` with coordinates jittered by U(-0.02, 0.02) after `set.seed(7)`: no tied distances. As in every data set, the pair that defines `umax` still sits exactly on the last bin edge (see "Platform dependence") |
| `quakes_irregular` | 300 | `datasets::quakes`, first 300 unique positions (the roxygen example data). The points are irregularly spaced but lie on a 0.01-degree grid, so distances can tie. No download |
| `brain_Oprk1_subsample` | 2170 | Brain gene Oprk1, B = 5. The variogram subsample `sample(N, 1000)` is drawn before the permutations; the Visium side is zero-inflated (summary intermediates only) |
| `kidney_AB_jitter_seed17` | 273 | Inputs of `kidney_AB_jitter` with seed = 17, `maxDistPrctile = 0.3`, `deltaX` = 0.1..0.9 and `deltaY` = c(0.2, 0.5, 0.8); B = 3. Pins the argument plumbing (permutations from `set.seed(17)`, noise from `set.seed(17 + i)`) |
| `brain_Oprk1_N1000` | 1000 | First 1000 pixels of the brain case, deltas c(0.1, 0.3), B = 2. At N = 1000 exactly no subsample is drawn (only when N > 1000) |

## Platform dependence: what "exact" means

The fixtures were built on **macOS arm64** (R 4.5.2, Accelerate BLAS) and encode that machine's arithmetic.

- **Why results depend on the platform.** `geoR::variog()` takes `umax = max(u[u < max.dist])` from R's
  `dist()`. Its C code then counts a pair only if the libm `hypot(dx, dy)` is below `umax`. So in **every**
  data set, not only on lattices, the pair that defines `umax` sits exactly on the last bin edge. Lattice
  coordinates put many mathematically tied pairs there: 448 pairs lie within 4 ulp of `umax` for `kidney_AB`,
  21 for `kidney_AC` and 1210 for the brain subsample (56, 1 and 82 of them equal `umax` exactly in R's
  `dist()`); the jittered and quakes cases have one such pair.
- **How large the effect is.** Whether those pairs are counted depends on the last bit of `dist()` (compiled
  with fused multiply-add on arm64, without it on x86-64) and of the libm `hypot()`. One different bin count
  changes the target variogram. That moves every null correlation by about 1e-4, and some deltaStar values,
  whose nulls then move by up to about 0.2.
- **Smaller differences.** R's `qnorm()` (behind `rnorm()`), long double accumulation in `sum()`, `mean()` and
  `cor()`, and the locfit and `lm()` arithmetic also differ in the last bit between platforms. These move
  results by about 1e-16 relative, which the tolerances absorb.

**How the tests handle it.** Each fixture stores the build machine's signature in `meta$platform_signature`;
the probes are defined in `tests/testthat/helper-platform.R`.

- For every coordinate set the signature holds `max.dist`, `umax`, the number of pairs tied at `umax` and geoR's
  bin counts. The bin counts depend only on the coordinates, so they are those of every variogram of that set.
- The global probes (qnorm and rnorm draws, long double, `dist()` and `hypot()`, `lm()`, locfit) are recorded
  and reported, but they do not gate anything.
- Tests labelled **exact** run only for coordinate sets whose probes match on the running machine. Elsewhere
  they are skipped, and the skip message names the differing probes. `STCOMPARE_EXACT_TESTS=true` forces them,
  and `false` disables them.
- Tests labelled **portable** run everywhere.

Observed with the Bioconductor 3.22 Docker image (R 4.5.2, glibc 2.39, OpenBLAS), with `NOT_CRAN=true`:

| Platform | Kernel sets compared exactly | Realistic | Default-suite result |
|---|---|---|---|
| macOS arm64 (build machine) | 7 of 7 | aki, brain | 1618 expectations, 0 failed (6 slow tests skipped) |
| Linux arm64 | 6 of 7 (not `brain_Oprk1_subsample`) | none: geoR bins aki and brain differently | 1466 expectations, 0 failed (6 slow and 9 exact tests skipped) |
| Linux x86-64 | 3 of 7 (`kidney_AC`, `kidney_AB_jitter`, `kidney_AB_jitter_seed17`) | none | 1109 expectations, 0 failed (6 slow and 17 exact tests skipped) |

In CRAN mode (`NOT_CRAN` unset) the suite passes as well: 787 expectations on Linux x86-64 and 1268 on macOS.
With the slow tests, the calibration's exact block compares 60 of 60 jobs on macOS, 34 on Linux arm64 and 20 on
Linux x86-64; all of them match. Before these changes, the same Linux runs failed 384 (x86-64) and 117 (arm64)
expectations.

Two claims made here earlier were wrong:

- that `kidney_AB_jitter` and `quakes_irregular` "have no ties": `quakes_irregular` bins differently on x86-64;
- that the `inst/extdata` goldens "were produced on Linux": nothing records where they were produced
  (`inst/scripts` used 20-22 workers). They reproduce bit for bit on macOS arm64 and on neither Linux platform,
  so they were most likely computed on macOS arm64 too.

**Exact references for another platform.** Build them with the legacy R code on that platform, into a separate
checkout (`build_test_fixtures.R tier0 tier1`). Tier 2 copies the published values, so its exact layer stays
tied to macOS arm64.

## Rebuilding

Requirements:

- the packages in `DESCRIPTION`, plus `devtools`, `rhdf5` (Bioconductor; for the AKI `.h5` files) and `rjson`
  or `jsonlite`;
- a Unix-alike for parallel parts. On Windows the builder runs serially; `STCOMPARE_BUILD_WORKERS=1` does the
  same elsewhere.

Run every command **from the repository root**:

```sh
# optional: choose the cache root (default: tools::R_user_dir("STcompare", "cache"))
export STCOMPARE_DATA_CACHE=/path/outside/the/repo

Rscript data-raw/download_data.R                   # 13 files, 130 MB; skips files already present with the right md5
Rscript data-raw/build_inputs_aki.R                # ~20 s  -> <cache>/data-raw/inputs/aki_rast.rds (38 MB)
Rscript data-raw/build_inputs_brain.R              # ~40 s  -> brain_merfish_visium_rast.rds (9.3 MB), brain_celltype_rast.rds (0.6 MB)
Rscript data-raw/build_inputs_merfish_replicates.R # ~75 s  -> merfish_replicates_rast.rds (7.5 MB); tier 3 only
Rscript data-raw/build_test_fixtures.R tier0 tier1 tier2 verify
Rscript data-raw/build_test_fixtures.R compare     # rebuild in memory and diff against the committed fixtures
```

**Input builders.** Each `build_inputs_*.R` downloads whatever is missing (through `download_data.R`). It
stops unless the rebuilt input reproduces the published `correlationCoef` of `inst/extdata` to 1e-12.
Observed: AKI 5e-16, brain 1.05e-15, cell types 1.1e-16, MERFISH 2.0e-14 (STalign) and 2.6e-14 (affine).

**Fixture builder.** `build_test_fixtures.R` takes any of `tier0 tier1 tier2 verify` (default: all four),
optionally followed by `compare` or `force`. It writes `tests/testthat/fixtures/{kernel,calibration,realistic}_fixture.rds`.

- **Inputs.** Tier 0 needs the brain input, tier 2 the AKI and brain inputs, tier 1 only the package data.
  Before use, every cached input must carry the expected `meta$script` and reproduce the published `r` and
  naive p of every published gene (1e-12 and 1e-10 relative). A stale or modified input stops the build.
- **Self-checks.** Every tier-0 intermediate is checked for bit-identity against the package (`stopifnot`),
  including every permuted field of `viladomatCorrelation()`.
- **RNG.** The builder sets `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` at start-up and again
  before every `set.seed()`, so a non-default kind from `~/.Rprofile` cannot leak into the fixtures. The kind the
  session started with is recorded in `meta$rng_kind_at_startup`.
- **Writing.** Each new fixture is compared with the committed file, ignoring `meta`. Added, removed and
  changed leaves are listed. The file is rewritten only if its content changed, so a no-op rebuild leaves git
  clean. `force` always rewrites (for example to refresh `meta`). `compare` writes nothing and exits with an
  error if anything differs.
- **Verify.** `verify` re-runs every kernel case at its stored B, every realistic gene at
  B = min(`$STCOMPARE_VERIFY_B`, 100) (default 100) and the 5 engineered negatives at B = min(.., 20). It uses
  the tests' tolerances: r 1e-12, naive p 1e-10 relative, nulls 1e-12, deltaStar identical. It prints the
  effective B, notes coordinate sets whose geoR binning differs from the build machine, and stops on any
  mismatch. The same comparison is available as a test: `STCOMPARE_SLOW_TESTS=true`.
- **Workers.** Parallel parts use `$STCOMPARE_BUILD_WORKERS` processes (default `detectCores() - 2`).

**What `meta` records** in every fixture:

- creation time, script, R version, platform, BLAS and LAPACK, `capabilities("long.double")`, `RNGkind()`;
- package versions (locfit, geoR, SEraster, SpatialExperiment, SummarizedExperiment, BiocParallel, Matrix,
  rhdf5, sf, rearrr, testthat) and `sf::sf_extSoftVersion()` (GEOS, GDAL, PROJ);
- git state:
  - `available`;
  - `in_repository`: TRUE only when the enclosing git work tree *is* this repository;
  - `commit`;
  - `dirty` and `status`: every changed or untracked file under `R`, `NAMESPACE`, `DESCRIPTION`, `data`,
    `data-raw` and `tests`, so untracked files count. Both are NA when git or the repository is unavailable;
- md5 content hashes of `R/*.R`, `NAMESPACE`, `DESCRIPTION`, `data/*.rda`, `data-raw/*.R`, the test helpers and
  the `inst/extdata` result files used. These identify the code even when it is not committed;
- for tiers 0 and 2: the build time, script and download md5s of the cached inputs, and `sources` (record, DOI,
  URL, licence and creators of every download);
- the per-case runtimes, and the platform signature.

Build fixtures only from the legacy R implementation. Once a compiled backend replaces the inner loop, the
fixtures are the reference it is compared with; they are not regenerated from it.

**Cache layout.** The root is `$STCOMPARE_DATA_CACHE` if set, otherwise `tools::R_user_dir("STcompare", "cache")`:

```
<root>/data-raw/            public downloads (md5-checked)
<root>/data-raw/inputs/     rasterized inputs (RDS)
<root>/data-raw/extracted/  unpacked 10x tarballs
```

**Downloads.** `download_data.R` downloads to `<file>.part` and renames the file only after its md5 has been
verified. A cached file with the wrong md5, or a download that does not match, stops with an error naming the
file, the URL and both checksums; a mismatching download is kept as `<file>.md5-mismatch`. To use it from R,
`source("data-raw/download_data.R")`, then call `stc_download(group = "aki")` (groups: `aki`, `brain`,
`celltype`, `merfish`), `stc_download_dir()` or `stc_inputs_dir()`.

## Running the tests

```r
devtools::test()                                     # tiers 0 and 2 plus portable checks, about 50 s
Sys.setenv(STCOMPARE_SLOW_TESTS = "true"); devtools::test()   # also tier 1 and all 65 realistic genes at B = 100 (~4 min on 18 workers)
```

**CRAN mode.** `devtools::test()`, `devtools::check()` and `testthat::test_local()` set `NOT_CRAN=true`. A plain
`R CMD check`, `rcmdcheck::rcmdcheck()` with its default `env` and the Bioconductor build system do not, so
`skip_on_cran()` applies there:

- skipped: tier 2, apart from the AKI Upk2 tests (zero inflation, deltas 0.01 and 0.05), and the
  `nThreads = 2` test;
- still run: all other kernel tests.

A CI job that should cover tier 2 must set `NOT_CRAN=true`.

**Environment variables:**

| Variable | Effect |
|---|---|
| `STCOMPARE_SLOW_TESTS=true` | Runs the calibration tier and the all-genes realistic comparison |
| `STCOMPARE_TEST_WORKERS` | Processes for the slow tests (default `BiocParallel::multicoreWorkers()`). Results do not depend on it: each job resets R's default RNG kinds and uses `BPPARAM = SerialParam()` |
| `STCOMPARE_EXACT_TESTS` | `auto` (default), `true` or `false`; see "Platform dependence" |

**Test labels.** These are defined in `tests/testthat/helper-fixtures.R`. The per-case tests are generated in
loops, so `testthat::test_file(desc = )` cannot select them; run the whole file instead (`filter =`).

| Label | Calls STcompare? | Purpose | Expectations (default suite, macOS) |
|---|---|---|---|
| `exact:` | yes | Stored legacy values at the documented tolerances, on matching platforms | 174 (16 tests) |
| `portable:` | yes | Every platform: STcompare against the reference kernels run on the same machine; prefix property; X/Y symmetry; worker independence; the p-value definition and the `returnPermutations = FALSE` branch; tie semantics (mocked); error paths returning NA; hand-computed `spatialSimilarity()` values; loose agreement with the fixtures | 305 (22 tests) |
| `canary:` | no | Dependency canaries and C++ kernel slots: R's RNG, SEraster, geoR, locfit and `lm()` through the `ref_*()` helpers, compared with the fixture | 816 (9 tests) |
| `fixture integrity:` | no | Consistency of the fixture files | 317 (3 tests) |

Only the `exact` and `portable` expectations are regression coverage of STcompare. With every exported
STcompare function replaced by `stop()`, all of them fail, and the 12 `canary` and integrity tests (1133
expectations) still pass. Of the 36 code mutants of `R/` in the test-strength review's set, the default suite
detects 35. The survivor swaps the lat/long columns, which preserves distances.

**RNG side effects.** Tests set R's default generators only for their own duration (`local_default_rng()`),
and `with_lecuyer()` sets all three kinds. Both restore the caller's kind and seed, so no test depends on the
file order, and running the tests leaves the session's RNG unchanged.

**Warnings.** Only locfit's expected "Estimated rdf < 1.0" warning is muffled (`quiet_locfit()`). Other
warnings, for example from a new backend, are reported by testthat. The error-path test counts and muffles
only the warnings that its inputs are expected to raise.

**Housekeeping:**

- testthat 3.3 creates an empty `tests/testthat/_snaps/` after a filtered run (`devtools::test(filter = ...)`)
  or any run with `CI=true`. It removes the directory only after an unfiltered local run. Empty directories are
  not tracked by git: delete it (`rmdir tests/testthat/_snaps`) or run the full suite. Do not git-ignore
  `_snaps`, because future snapshot files must be committed.
- Loading SparseArray 1.10 (through SEraster, SpatialExperiment and SummarizedExperiment) prints
  "Warning: stack imbalance in ..." once. This is upstream, and it is not an R condition, so testthat does not
  count it. `tests/testthat/setup.R` loads these namespaces before the tests, so a stack imbalance printed
  *during* the tests (for example by a compiled backend with a PROTECT bug) is new.

## Provenance and licences

All 13 downloads are CC BY 4.0. `tests/testthat/fixtures/README.md` is the attribution notice that ships with
the tests: creators, sources, licence and the changes made.

| File | Bytes | md5 | Source | Licence |
|---|---:|---|---|---|
| `IL3_filtered_feature_bc_matrix.h5` | 15 103 452 | `a3ea3a2cb4c3e01403d0116ab580c157` | [Zenodo 19074288](https://zenodo.org/records/19074288), doi:10.5281/zenodo.19074288 (concept doi:10.5281/zenodo.17676991); Clifton, Fan, Rabb | CC-BY-4.0 |
| `NL3_filtered_feature_bc_matrix.h5` | 13 822 500 | `6793a0f7ceeb082bce383748cfc2805a` | same | CC-BY-4.0 |
| `IL3_tissue_positions.csv` | 80 786 | `4e3deabfaa2cd44b3b0231776b3aaf1b` | same | CC-BY-4.0 |
| `NL3_tissue_positions.csv` | 87 431 | `24b23088b7a22a53086f9b62b21fed29` | same | CC-BY-4.0 |
| `aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz` | 90 127 | `adede04898833d385eb52fa2032f2434` | [Zenodo 19486091](https://zenodo.org/records/19486091), doi:10.5281/zenodo.19486091 (concept doi:10.5281/zenodo.17676991); Clifton, Fan, Rabb | CC-BY-4.0 |
| `STalign_S2R3_to_Visium.csv.gz` | 16 177 162 | `c092af84ca82a50539f4fdc0f51e8fdc` | [Zenodo 10724029](https://zenodo.org/records/10724029), doi:10.5281/zenodo.10724029 (STalign; concept doi:10.5281/zenodo.8384018); Clifton, Anant, Aihara, Fan | CC-BY-4.0 |
| `STalign_S2R2.csv.gz` | 12 517 332 | `e104115bb371a92efc06f7a3b908efc2` | same | CC-BY-4.0 |
| `STalign_S2R3_to_S2R2.csv.gz` | 17 110 794 | `73c1969dd2f52c12cbd4870b30f741d1` | same | CC-BY-4.0 |
| `STalign_cell_type_transcriptional_correlations.csv.gz` | 5 041 | `f2a4dc1997ed8ba2344094d534ff64fd` | [Zenodo 19582556](https://zenodo.org/records/19582556), doi:10.5281/zenodo.19582556 (concept doi:10.5281/zenodo.8384018); Clifton, Anant, Aihara, Fan | CC-BY-4.0 |
| `STalign_S2R3_cell_type_annotations.csv.gz` | 1 244 272 | `085e49b9bb8f0693d5854d2013b351dc` | same | CC-BY-4.0 |
| `STalign_Visium_cell_type_annotations.csv.gz` | 264 260 | `897985b732e45a7b4dbee3a64cd357e8` | same | CC-BY-4.0 |
| `Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz` | 45 925 294 | `017fc078b4aa21e1df39027d0bf3ccad` | 10x Genomics, Adult Mouse Brain (FFPE), Space Ranger 1.3.0, published 2021-08-16, <https://www.10xgenomics.com/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0>; file URL `https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/<file>` | CC-BY-4.0 (see below) |
| `Visium_FFPE_Mouse_Brain_spatial.tar.gz` | 10 632 472 | `a78bb08009c4defd791354dd5c292bfa` | same | CC-BY-4.0 |

- **Zenodo licences** were read from the Zenodo API.
- **10x Genomics licence.** The dataset page states "This dataset is licensed under the Creative Commons
  Attribution 4.0 International (CC BY 4.0) license". This was confirmed on 2026-10-03 from the page text as
  indexed by a web search engine; the matching metadata were 2,264 spots under tissue, Space Ranger 1.3.0 and
  2021-08-16. The page itself blocks automated downloads (Vercel checkpoint, HTTP 429), so check it once in a
  browser before a release.
- **md5 values.** For the Zenodo files they equal the checksums published by the Zenodo API.

Download URLs (the ones used by `inst/scripts/` and the vignettes; also in `stc_manifest` in `download_data.R`):

```
https://zenodo.org/records/19074288/files/IL3_filtered_feature_bc_matrix.h5?download=1
https://zenodo.org/records/19074288/files/NL3_filtered_feature_bc_matrix.h5?download=1
https://zenodo.org/records/19074288/files/IL3_tissue_positions.csv?download=1
https://zenodo.org/records/19074288/files/NL3_tissue_positions.csv?download=1
https://zenodo.org/records/19486091/files/aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz?download=1
https://zenodo.org/records/10724029/files/STalign_S2R3_to_Visium.csv.gz?download=1
https://zenodo.org/records/10724029/files/STalign_S2R2.csv.gz?download=1
https://zenodo.org/records/10724029/files/STalign_S2R3_to_S2R2.csv.gz?download=1
https://zenodo.org/records/19582556/files/STalign_cell_type_transcriptional_correlations.csv.gz?download=1
https://zenodo.org/records/19582556/files/STalign_S2R3_cell_type_annotations.csv.gz?download=1
https://zenodo.org/records/19582556/files/STalign_Visium_cell_type_annotations.csv.gz?download=1
https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz
https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/Visium_FFPE_Mouse_Brain_spatial.tar.gz
```

**What is committed.** No downloaded file is committed. The fixtures contain only small derived subsets:

- AKI: CPM of 35 genes on 311 pixels, from the Zenodo AKI records.
- Brain: log10(CPM + 1) of 30 genes, plus Oprk1 in tier 0, on 2170 pixels, from the STalign MERFISH record and
  the 10x Visium dataset.
- Everything else comes from the package itself (`data(speKidney)`, `data(simRanPatternRasts)`,
  `inst/extdata/*.RData`) or from base R (`datasets::quakes`).

**Golden references.** Tier 2 uses `inst/extdata/kidneyCorrelation.RData` and
`inst/extdata/brain-MERFISH-10x-visium/brainCorrelation.RData`. Both were produced by the current code (seed 0).
On macOS arm64, from the rebuilt inputs, the current code reproduces their nulls to 5.4e-15 and their
deltaStar exactly (65 genes x 100 permutations x 2 directions). The machine they were computed on is not
recorded; see "Platform dependence". `simRanPatternResults.RData` predates the 2026-04-28 screening-threshold
change and is used as a statistical reference only.

## RNG semantics (what an exact port must reproduce)

- **Permutation order.** Mersenne-Twister (Inversion, Rejection sampling), seeded by `set.seed(seed)` in the
  calling session. If N > 1000 (not at N = 1000), the variogram subsample `ids <- sample(N, 1000)` is drawn
  first. Then permutation i is `X[sample.int(N, N)]` for i = 1..B in order (`sample(X, length(X))`). Both
  directions use the same seed, and so the same index permutations. Swapping X and Y swaps the two directions
  exactly.
- **Noise.** Permutation i runs inside a BiocParallel task, where `RNGkind()` is L'Ecuyer-CMRG (Inversion).
  `matchingVariograms()` calls `set.seed(seed + i)` and draws `rnorm(N)` once per delta, in delta order.
  Consequence: calling the exported `matchingVariograms()` directly in a Mersenne-Twister session gives
  different noise than inside the package (the tests wrap it in `with_lecuyer()`).
- **Independence of parallelism.** Results do not depend on `nThreads` or on the BiocParallel backend
  (`MulticoreParam(1)`, `MulticoreParam(2)` and `SerialParam()` are bit-identical).
- **Prefix property.** Permutation i depends only on (seed, i). The first k nulls of a B-permutation run equal a
  k-permutation run. This is why tier 2 can test B = 10 against the first 10 of 1000 published nulls.

## Legacy behaviour the fixtures deliberately do not encode

The fixtures store **raw nulls, deltaStar and raw (unadjusted) p-values**, so that fixing the following does
not invalidate them:

- **`p.adjust()` in `spatialCorrelationGeneExp()` is a no-op.** It is applied inside the per-gene loop, to a
  single p-value (`R/spatialCorrelation.R:836-842`). No test checks adjusted p-values of more than one gene. The
  published (BH-adjusted) p of the iterative protocol are stored for reference only.
- **Empirical p-values can be exactly 0** (`extreme / B`, with no +1). The tests compute the expected p-value
  from raw nulls with `stc_empirical_p()` in `tests/testthat/helper-fixtures.R`. If the definition changes,
  update that one function; no fixture rebuild is needed. The tail counts `nExtreme` are stored as well.
- **`spatialSimilarity()` reports `numPixelInThresh = 1`** for a gene skipped by `minPixels` (the row count of
  its one-row summary). The hand-computed test does not check that value.
- **Errors become silent NA rows**, and **locfit segfaults** (the R session dies) when 1 <= delta * N < 2 (for
  example N = 12 with delta 0.1, or N < 200 with delta 0.01). The NA rows for constant X, delta * N < 1 and one
  NA in X are tested; the segfault inputs are documented in `kernel_fixture.rds$edge_cases` and never run.

Further properties of the current code that matter for exactness:

- **Lattice ties and the bin edge.** See "Platform dependence". For the same reason, translating or scaling
  lattice coordinates changes the results, while rotating them by 90 degrees does not.
- **`lm()` is not bit-reproducible** across BLAS libraries or memory alignment. Tests therefore compare
  floating-point values with a tolerance of 1e-12 and RNG draws with 1e-14. They require exact equality only
  for integers, bin counts and deltaStar.
- **locfit's default evaluation** is an adaptive kd-tree with vertex interpolation, not an exact kernel
  smoother. `fitted_exact_evdat` (`ev = dat()`) is stored as the reference for a smoother that is not
  bit-compatible with the tree. The two agree only to a correlation of about 0.95-0.99.

## Acceptance criteria for a compiled backend

- **Kernels, on every machine.** Swap each compiled kernel into the corresponding `ref_*()` helper
  (`helper-fixtures.R`). The portable tests compare STcompare with those helpers on the same machine, so a
  bit-compatible kernel must reproduce R's own result there, whatever the platform.
- **Bit-compatible mode** (locfit's tree ported or wrapped, R's RNG reproduced or injected).
  - All `exact` expectations must pass on the build platform (macOS arm64), together with the portable ones.
  - The fixtures embed the build machine's arithmetic. The binning must match R's `dist()` there (FMA-contracted
    on arm64) and the platform libm `hypot()`.
  - Compiling with `-ffp-contract=off` reproduces x86-64 R, **not** the arm64-built fixtures.
  - On another platform, use exact references built there by the legacy R code (see "Platform dependence");
    otherwise only the portable layer applies.
- **New-smoother mode.**
  - Compare the smoother with `fitted_exact_evdat`.
  - Tier 1 must pass, including "p-values agree in distribution with the stored reference". That test requires:
    - no systematic shift of the 60 jobs' raw p-values against the stored ones (paired Wilcoxon p >= 0.001);
    - mean |dp| <= 0.1;
    - the share of permutations at the smallest and at the largest delta within 0.25 of the reference
      (legacy: 58% and 1% for X, 49% and 1% for Y; about 45% of selected deltas are interior);
    - at least 3 distinct selected deltas.
  - A run of the legacy algorithm with another random stream passes that test, while "always the smallest
    delta", "always the largest delta" and "no smoothing" fail it.
  - On tier 2, every positive and negative gene should stay significant and every null gene should not, and the
    engineered-negative relation must hold.

## Approximate runtimes (macOS arm64, M1 Ultra, R 4.5.2, Accelerate BLAS)

| Step | Time |
|---|---|
| `download_data.R` with everything cached (md5 checks only) | < 1 s |
| `build_inputs_aki.R` / `build_inputs_brain.R` / `build_inputs_merfish_replicates.R` | 20 s / 40 s / 75 s |
| `build_test_fixtures.R tier0` | 30 s |
| `build_test_fixtures.R tier1` (60 jobs at B = 100, 18 workers) | 31 s wall (about 7 CPU-min) |
| `build_test_fixtures.R tier2` | 2 s |
| `build_test_fixtures.R verify` (B = 100, 18 workers) | 170 s wall (about 44 CPU-min) |
| `devtools::test()` (kernel 21 s, realistic 28 s, slow tests skipped) | about 55 s |
| slow tests with `STCOMPARE_SLOW_TESTS=true` (18 workers): calibration / all realistic genes | 28 s / 160 s wall |
| the default suite on Linux x86-64 under emulation / Linux arm64 (Docker, 16 CPUs) | 130 s / 85 s |

For reference, the current R cost per gene at B = 100 (both directions, one thread) is about 7 s for a simulated
pair (N = 273), 19-20 s for AKI (N = 311, 11 deltas) and 65-69 s for brain (N = 2170, 9 deltas).
