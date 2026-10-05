# Validation against the published analyses

Written by `bench/validate-published.R` (see `bench/README.md`). Each analysis whose results are in
`bench/published` is re-run from inputs rebuilt from the public downloads (`data-raw/`) with the authors'
parameters (`bench/published/scripts/`), and compared gene by gene with the stored table:

- `r`: |Δ correlationCoef| ≤ 1e-12; `naive p`: relative difference of pValueNaive ≤ 1e-10;
- `B`: the same number of permutations per gene (the same screening decisions);
- `δ*`: deltaStarX and deltaStarY identical for every stored permutation;
- `nulls`: every stored null correlation within 1e-9 × the largest published |null| of that gene and direction;
- `count >` / `count ≥`: identical numbers of |null| > |r| (the legacy count) and |null| ≥ |r| (the current one);
- `BH p`: our final pValuePermuteX/Y equal p.adjust((b + 1) / (B + 1), "BH") recomputed from the published
  nulls (relative tolerance 1e-12; the stored values used the legacy b / B and are not compared directly).

On macOS arm64 the MERFISH correlationCoef differ from the stored ones by up to 2.6e-14 (10-11 of 483 are
bit-identical) because R's cor() accumulates in long double, which is plain double there: the MERFISH analyses
were evidently computed with extended precision (on Linux arm64, 467 and 469 of 483 are bit-identical, max
difference 1.1e-16), while the AKI, brain and cell-type values reproduce best on macOS arm64. The engine is
not involved in r or the naive p-value.

All six analyses reproduce on macOS arm64, where the AKI, brain and cell-type results were computed. On Linux
arm64 (Docker, Bioconductor 3.22, GCC 13) the MERFISH and cell-type analyses reproduce as well, but AKI and brain
do not: glibc's hypot() puts 4 (AKI) and 16 (brain) lattice pairs that lie within a few ulp of the last
variogram bin edge on the other side of it. geoR::variog() bins them the same way there, so the legacy R
code does not reproduce these two published analyses on Linux either.

<!-- BEGIN validate-published:internal -->
## Mode `internal` — pre-integration: frozen snapshot, engine called directly

Re-reported 2026-10-04 14:01 EDT (`--reuse`) from results computed 2026-10-04 13:52 to 13:53 EDT on DIMKC6JP7VQP1 (Darwin 24.6.0, aarch64-apple-darwin20), R version 4.5.2 (2025-10-31); 16 threads.
Package: STcompare 0.1.0 from /private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/final-validate/mirror/STcompare, installed with R CMD INSTALL (engine compiled with -g -O2); engine/R source fingerprint `37e045c7a01d`.
Call: `Rscript bench/validate-published.R --mode=internal --threads=16 --label="pre-integration: frozen snapshot, engine called directly"`. Per-gene tables: `/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/final-validate/out` (outside the repository).

Mismatching genes by check (0 everywhere = the published analysis is reproduced). Genes = genes compared;
"(+k)": rows of `merfishCorrelation.RData` whose stored p-values do not follow from its stored nulls
(the table was patched with rows of an earlier run; see below). They are compared but counted separately.

| Analysis | Genes | N | B 100/1000 pub | B 100/1000 ours | r | naive p | B | δ* | nulls | count > | count ≥ | BH p | NA rows | max rel null diff | max abs null diff | max abs Δr | max Δr / r |
|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| AKI kidney Visium, iterative (`kidneyCorrelation.RData`) | 1046 | 311 | 289/757 | 289/757 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3e-14 | 1.2e-14 | 5e-16 | 1.2e-15 |
| AKI kidney Visium, fixed B = 100 (`kidneyCorrelationNoIter.RData`) | 1046 | 311 | 1046/- | 1046/- | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.6e-14 | 6.9e-15 | 5e-16 | 1.2e-15 |
| MERFISH replicates, affine (`merfishCorrelation_affine.RData`) | 483 | 1299 | 97/386 | 97/386 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 9.5e-14 | 2.3e-14 | 2.6e-14 | 1.9e-13 |
| MERFISH replicates, STalign (`merfishCorrelation.RData`) | 361 (+122) | 1371 | 86/397 | 86/397 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 8.8e-14 | 8.1e-15 | 2e-14 | 1.6e-12 |
| Brain MERFISH vs Visium (`brainCorrelation.RData`) | 325 | 2170 | 178/147 | 178/147 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 7.8e-14 | 2.8e-14 | 1e-15 | 2.9e-14 |
| Brain cell types (`ctCorrelation.RData`) | 16 | 2174 | 6/10 | 6/10 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6e-15 | 2.4e-15 | 1.1e-16 | 2.6e-16 |

| Analysis | Genes | Wall time | Threads | s / gene | Authors' time | Authors' workers | Speedup | RNG state unchanged |
|---|---:|---|---:|---|---|---|---|---|
| AKI kidney Visium, iterative | 1046 | 13.2 s | 16 | 0.013 | not recorded | 22 |  | yes |
| AKI kidney Visium, fixed B = 100 | 1046 | 2.0 s | 16 | 0.002 | not recorded | 22 |  | yes |
| MERFISH replicates, affine | 483 | 39.6 s | 16 | 0.082 | 7.65 h | 20 (MulticoreParam()) | 697× | yes |
| MERFISH replicates, STalign | 483 | 38.9 s | 16 | 0.080 | 16.8 h | 20 (MulticoreParam()) | 1557× | yes |
| Brain MERFISH vs Visium | 325 | 12.5 s | 16 | 0.039 | 1.78 h | 22 | 511× | yes |
| Brain cell types | 16 | 0.9 s | 16 | 0.057 | 6.95 min | 22 | 455× | yes |

Speedup = the authors' reported wall time / ours (their 20-22 workers on an unrecorded machine, our threads here).

Screening rounds (ours):

- `aki_iter`: round 1 (B = 100): 1046 genes, 757 promoted; round 2 (B = 1000): 757 genes
- `aki_fixed`: round 1 (B = 100): 1046 genes
- `merfish_affine`: round 1 (B = 100): 483 genes, 386 promoted; round 2 (B = 1000): 386 genes
- `merfish_stalign`: round 1 (B = 100): 483 genes, 397 promoted; round 2 (B = 1000): 397 genes
- `brain`: round 1 (B = 100): 325 genes, 147 promoted; round 2 (B = 1000): 147 genes
- `celltypes`: round 1 (B = 100): 16 genes, 10 promoted; round 2 (B = 1000): 10 genes

Stored p-values that do not follow from the stored nulls:

- `merfish_stalign`: the stored p-values of 122 rows are not BH(b / B) of the published nulls: 116 through pValuePermuteX (the number in dev/investigation/04) and 6 more through pValuePermuteY only (their raw pX is 0). They are 122 of the 123 rows whose stored p-values can show this at all: the other 360 rows have no exceedance in either direction, so their stored p is 0 under any adjustment, and Cxcr2 has the largest raw p-values, which BH leaves unchanged. The stored values therefore cannot tell which rows were copied from the earlier run (`bench/published/scripts/biological-replicates-example.R`, lines 210-217); they show that the stored BH adjustment was computed over other raw p-values than the published nulls give. These rows are compared like the others but not counted in the table: 122 match on every check (r, naive p, permutations, deltaStar, nulls, both counts, BH p), 0 differ; max relative null difference 1.1e-13.

No mismatches.

<!-- END validate-published:internal -->

<!-- BEGIN validate-published:exported -->
## Mode `exported` — remap default (compareSpatial surrogate = "remap", minDetected = NULL -> 0 pixels); legacy functions unchanged

Run 2026-10-05 11:43 EDT on DIMKC6JP7VQP1 (Darwin 24.6.0, aarch64-apple-darwin20), R version 4.5.2 (2025-10-31); 16 threads.
Package: STcompare 0.1.0.9000 from /private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/remap-default/mirror/STcompare, installed with R CMD INSTALL (engine compiled with -g -O2); engine/R source fingerprint `19f18547621d`.
Call: `Rscript bench/validate-published.R --mode=exported --threads=16 --label="remap default (compareSpatial surrogate = "remap", minDetected = NULL -> 0 pixels); legacy functions unchanged"`. Per-gene tables: `/Users/ks38/Library/Caches/org.R-project.R/R/STcompare/bench/validate-published` (outside the repository).

Mismatching genes by check (0 everywhere = the published analysis is reproduced). Genes = genes compared;
"(+k)": rows of `merfishCorrelation.RData` whose stored p-values do not follow from its stored nulls
(the table was patched with rows of an earlier run; see below). They are compared but counted separately.

| Analysis | Genes | N | B 100/1000 pub | B 100/1000 ours | r | naive p | B | δ* | nulls | count > | count ≥ | BH p | NA rows | max rel null diff | max abs null diff | max abs Δr | max Δr / r |
|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| AKI kidney Visium, iterative (`kidneyCorrelation.RData`) | 1046 | 311 | 289/757 | 289/757 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3e-14 | 1.2e-14 | 5e-16 | 1.2e-15 |
| AKI kidney Visium, fixed B = 100 (`kidneyCorrelationNoIter.RData`) | 1046 | 311 | 1046/- | 1046/- | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.6e-14 | 6.9e-15 | 5e-16 | 1.2e-15 |
| MERFISH replicates, affine (`merfishCorrelation_affine.RData`) | 483 | 1299 | 97/386 | 97/386 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 9.5e-14 | 2.3e-14 | 2.6e-14 | 1.9e-13 |
| MERFISH replicates, STalign (`merfishCorrelation.RData`) | 361 (+122) | 1371 | 86/397 | 86/397 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 8.8e-14 | 8.1e-15 | 2e-14 | 1.6e-12 |
| Brain MERFISH vs Visium (`brainCorrelation.RData`) | 325 | 2170 | 178/147 | 178/147 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 7.8e-14 | 2.8e-14 | 1e-15 | 2.9e-14 |
| Brain cell types (`ctCorrelation.RData`) | 16 | 2174 | 6/10 | 6/10 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6e-15 | 2.4e-15 | 1.1e-16 | 2.6e-16 |

| Analysis | Genes | Wall time | Threads | s / gene | Authors' time | Authors' workers | Speedup | RNG state unchanged |
|---|---:|---|---:|---|---|---|---|---|
| AKI kidney Visium, iterative | 1046 | 13.7 s | 16 | 0.013 | not recorded | 22 |  | yes |
| AKI kidney Visium, fixed B = 100 | 1046 | 1.9 s | 16 | 0.002 | not recorded | 22 |  | yes |
| MERFISH replicates, affine | 483 | 39.2 s | 16 | 0.081 | 7.65 h | 20 (MulticoreParam()) | 704× | yes |
| MERFISH replicates, STalign | 483 | 38.8 s | 16 | 0.080 | 16.8 h | 20 (MulticoreParam()) | 1558× | yes |
| Brain MERFISH vs Visium | 325 | 12.7 s | 16 | 0.039 | 1.78 h | 22 | 506× | yes |
| Brain cell types | 16 | 0.9 s | 16 | 0.058 | 6.95 min | 22 | 451× | yes |

Speedup = the authors' reported wall time / ours (their 20-22 workers on an unrecorded machine, our threads here).

Stored p-values that do not follow from the stored nulls:

- `merfish_stalign`: the stored p-values of 122 rows are not BH(b / B) of the published nulls: 116 through pValuePermuteX (the number in dev/investigation/04) and 6 more through pValuePermuteY only (their raw pX is 0). They are 122 of the 123 rows whose stored p-values can show this at all: the other 360 rows have no exceedance in either direction, so their stored p is 0 under any adjustment, and Cxcr2 has the largest raw p-values, which BH leaves unchanged. The stored values therefore cannot tell which rows were copied from the earlier run (`bench/published/scripts/biological-replicates-example.R`, lines 210-217); they show that the stored BH adjustment was computed over other raw p-values than the published nulls give. These rows are compared like the others but not counted in the table: 122 match on every check (r, naive p, permutations, deltaStar, nulls, both counts, BH p), 0 differ; max relative null difference 1.1e-13.

No mismatches.

<!-- END validate-published:exported -->
