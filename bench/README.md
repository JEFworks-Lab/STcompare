# bench/

Scripts that are too slow or need too much data for the test suite. The directory is listed in
`.Rbuildignore`, so nothing here is part of the package.

| File | What it is |
|---|---|
| `published/` | The authors' published results (`.RData`, 30 MB): the reference of `validate-published.R` and of the test fixtures (`data-raw/`). They shipped in `inst/extdata` until version 0.1.0. `published/scripts/` holds the authors' analysis scripts that wrote them (formerly `inst/scripts/`); they need a source checkout, the downloads of `data-raw/` and the GitHub packages MERINGUE and scatterbar, and some paths are the authors' own |
| `validate-published.R` | The main acceptance test of the compiled engine (`dev/engine-spec.md`, section 5): re-runs every published analysis whose results are in `published/` and compares them gene by gene |
| `validation-results.md` | The latest summary written by `validate-published.R`, one section per mode |
| `time-compareSpatial.R` | Run time of `compareSpatial()` with its defaults on the realistic test genes and on the full AKI and brain inputs, with the permutations each gene used |

## validate-published.R

### What it runs

Six analyses, each from inputs rebuilt from the public downloads (`data-raw/`) and with the authors' parameters
(`published/scripts/`, summarised in `dev/investigation/04-datasets-and-test-tiers.md`, section 1.3). The gene lists
are the row names of the stored results. X is the first element of the authors' input list.

| `--analyses` id | Reference (`bench/published/`) | Genes | Shared pixels | Call |
|---|---|---:|---:|---|
| `aki_iter` | `kidneyCorrelation.RData` | 1046 | 311 | `spatialCorrelationGeneExpIterPermutations()`, X = control (NL3), Y = AKI (IL3), `assayName = "CPM"`, deltas `c(0.01, 0.05, 0.1..0.9)`, `nPermutations = c(100, 1000)`, seed 0 |
| `aki_fixed` | `kidneyCorrelationNoIter.RData` | 1046 | 311 | `spatialCorrelationGeneExp()`, same input and deltas, `nPermutations = 100` |
| `merfish_affine` | `merfishCorrelation_affine.RData` | 483 | 1299 | iterative, X = S2R2, Y = S2R3 (affine coordinates), first assay, 11 deltas, `c(100, 1000)` |
| `merfish_stalign` | `merfishCorrelation.RData` | 483 | 1371 | the same with the STalign coordinates |
| `brain` | `brain-MERFISH-10x-visium/brainCorrelation.RData` | 325 | 2170 | iterative, X = MERFISH, Y = Visium, `assayName = "lognorm"`, default deltas and `nPermutations` |
| `celltypes` | `brain-MERFISH-10x-visium/ctCorrelation.RData` | 16 | 2174 | iterative, X = Visium, Y = MERFISH, defaults |

Two modes, selected with `--mode` (or `STCOMPARE_VALIDATE_MODE`, or `validate_mode <- "..."` before `source()`):

- `exported` (default): calls the exported functions exactly as the authors did (`assayName`, delta lists,
  `nPermutations`, `seed`, `adjustMethod = "BH"`), with `nThreads` threads and `BPPARAM = NULL`. The authors
  passed `BPPARAM = MulticoreParam()` in some scripts; with the compiled engine `BPPARAM` only sets the thread
  count, and no result depends on it. The script refuses this mode while the exported functions are still the
  legacy R code (about 200 CPU-hours for these analyses).
- `internal`: calls `.stc_engine_correlate()` directly. Round 1 runs every gene; a gene whose unadjusted
  `(b + 1) / (B + 1)` p-values are both below `(alpha / nPermutations[k]) * 100` (the screening rule of
  `spatialCorrelationGeneExpIterPermutations()`, alpha = 0.05) continues the same engine session up to
  `nPermutations[k + 1]` (`attr(, "state")`; exact by the prefix property). BH across genes at the end.

### What it compares

Per gene, with the published table as the reference:

| Check | Rule |
|---|---|
| `correlationCoef` | \|difference\| ≤ 1e-12 (r is bounded by 1; the per-gene relative difference is reported too) |
| `pValueNaive` | relative difference ≤ 1e-10 |
| permutations | the same number per gene, i.e. the same screening decisions |
| `deltaStarX/Y` | identical for every stored permutation |
| `nullCorrelationsX/Y` | every stored null within 1e-9 × the largest published \|null\| of that gene and direction |
| exceedances | identical counts of \|null\| > \|r\| (what the legacy code counted) and \|null\| ≥ \|r\| (the current definition) |
| final p-values | ours equal `p.adjust((b + 1) / (B + 1), "BH")` recomputed from the **published** nulls across the genes of the analysis (relative tolerance 1e-12) |

The stored p-values are not compared directly: they use the legacy definition `b / B` (strict `>`), and
`kidneyCorrelationNoIter.RData` stores them unadjusted (the legacy `spatialCorrelationGeneExp()` called
`p.adjust()` on one gene at a time). They are only checked for consistency with the published nulls. In
`merfishCorrelation.RData`, which the authors patched with rows of an earlier run
(`published/scripts/biological-replicates-example.R`, lines 210-217), the stored p-values of 122 rows do not follow
from the stored nulls; those rows are compared like the others but reported separately. (They are all but one
of the rows whose stored p-values can reveal this: the other 360 rows have p = 0 under any adjustment, so the
stored values cannot tell which rows came from which run.) Every call is also checked to leave `.Random.seed`
and `RNGkind()` unchanged.

On macOS arm64 the reproduced MERFISH correlations differ from the stored ones by up to 2.6e-14 (10-11 of 483
are bit-identical), because R's `cor()` accumulates in `long double`, which is plain `double` on macOS arm64. The
MERFISH analyses were evidently run on a machine with extended precision: in the Linux arm64 container 467 and
469 of 483 correlations are bit-identical (largest difference 1.1e-16), so the rebuilt inputs are the authors'.
The AKI, brain and cell-type correlations reproduce best on macOS arm64 (to 1e-15 or better).

The checks were verified to catch deliberate errors injected into otherwise exact results: one null changed by
1e-8 relative, one deltaStar changed, one gene not promoted to the second round (its permutations, and the BH
p-values of the other genes, then differ), r + 1e-11 and the naive p-value × (1 + 1e-9) were each reported as
mismatches of the right kind, and nothing else was.

### Running it

```sh
# once: download and rasterize the inputs into the cache (outside the repository; see data-raw/README.md)
Rscript data-raw/download_data.R
Rscript data-raw/build_inputs_aki.R
Rscript data-raw/build_inputs_brain.R
Rscript data-raw/build_inputs_merfish_replicates.R

# then, from the repository root
Rscript bench/validate-published.R                          # exported functions, 16 threads, all analyses
Rscript bench/validate-published.R --mode=internal          # the engine directly
Rscript bench/validate-published.R --threads=8 --analyses=celltypes,brain
```

Options: `--threads=N` (default 16; `STCOMPARE_VALIDATE_THREADS`), `--analyses=id,id`
(`STCOMPARE_VALIDATE_ANALYSES`), `--pkg=DIR` or `--pkg=installed` (default: the package that contains
`bench/`), `--published=DIR` (default: `bench/published`), `--build=install|load_all`, `--out=DIR`,
`--results=FILE`, `--label=TEXT`, and `--reuse` (compare and report the results an earlier run of the same mode
saved in `--out`, without recomputing them). The cache root
is `$STCOMPARE_DATA_CACHE` or `tools::R_user_dir("STcompare", "cache")`.

By default (`--build=install`) the script copies the package sources to a temporary directory and installs them
with `R CMD INSTALL`, so the engine is compiled with R's optimising flags, as users get it. `devtools::load_all()`
compiles with `-g -O0` unless `debug = FALSE`, which makes the engine about 10 times slower (the results are
bit-identical); `--build=load_all` uses `debug = FALSE, recompile = TRUE`. Timings measured on a `load_all()`
build are not representative.

Outputs:

- one per-gene CSV per analysis (every check, the counts, both p-values) and our results as RDS, in `--out`
  (default `<cache>/bench/validate-published/`, never inside the repository);
- `bench/validation-results.md`: the summary tables. The section of the mode that ran is replaced; the file is
  only written when all six analyses ran (or `--results` is given).

The script exits with the summary printed. All six analyses take under 2 minutes on 16 threads of an M1 Ultra
(about 2 minutes 10 seconds including the package installation; see `validation-results.md` for the measured
times). The authors' runs of four of them took about 26 hours of wall time on 20-22 workers (brain 1.78 h,
MERFISH STalign 16.8 h, MERFISH affine 7.65 h, cell types 6.95 min); the kidney runs were not timed.

### Platform

All six analyses reproduce on macOS arm64 (identical deltaStar, nulls to about 1e-13), where the AKI, brain and
cell-type results were computed. On Linux arm64 (Docker, Bioconductor 3.22, GCC 13) the MERFISH and cell-type
analyses reproduce as well, but AKI and brain do not. glibc's `hypot()` rounds a few of
the lattice pairs that sit within a few ulp of the last variogram bin edge differently from macOS, so 4 (AKI)
and 16 (brain) pairs fall out of the last bin. That changes the target variogram, and with it most deltaStar
choices and every null. `geoR::variog()` bins these pairs the same way on Linux, so the legacy R code does not
reproduce these two published analyses there either (see `data-raw/README.md`, "Platform dependence").

## time-compareSpatial.R

Times `compareSpatial()` with its defaults (adaptive p-values with `exceedances = 10` and `nPermutations = 10000`,
the extended delta grid, rank-remapped surrogates with no detection filter) and summarises the permutations per
gene (the `nPermutations` column), how many genes stopped early or reached the limit, and how many have
`padj < 0.05`. Datasets: the realistic test genes
(`aki_fixture`: 35 genes on 311 pixels; `brain_fixture`: 30 genes on 2170 pixels) and, when the cache of
`data-raw/` exists, the full inputs of the published AKI (1046 genes, control vs AKI, assay CPM) and brain
(325 genes, MERFISH vs Visium, assay lognorm) analyses.

```sh
R CMD INSTALL .                                   # the engine compiled with R's optimising flags
Rscript bench/time-compareSpatial.R               # 16 threads, every dataset available
Rscript bench/time-compareSpatial.R --threads=8 --datasets=aki_fixture,brain_fixture --out=/tmp/cs
Rscript bench/time-compareSpatial.R --datasets=aki,brain --nPermutations=1000
```

With the defaults, the run time is dominated by the significant genes, which run to 10000 permutations: on 16
threads of an M1 Ultra the 1046 AKI genes took 118 s (454 genes reached the limit; median 4626 permutations per
gene; 752 significant) and the 325 brain genes 187 s (100 at the limit; median 512; 173 significant). A smaller
`nPermutations` is proportionally faster: with `nPermutations = 1000` they took 18 s and 29 s, with the same
significant genes. (With gaussian surrogates and their sqrt(N) detection filter, the defaults before the
remapped surrogates, the same runs took 106 s and 125 s.)
