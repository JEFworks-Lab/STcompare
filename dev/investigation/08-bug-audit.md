# STcompare: correctness, statistics, API, packaging and docs audit

Scope: every file in `R/`, `NAMESPACE`, `DESCRIPTION` and `man/*.Rd`, plus a skim of `vignettes/*.Rmd` and `inst/scripts/*.R`.
Every finding was reproduced in R 4.5.2 (macOS arm64, Bioconductor 3.22) with
`devtools::load_all()`. The repository was not modified. All work is in the work directory:

`WD = /private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/bugs/`

Line numbers are 1-indexed and refer to the repository at commit `2983c99`.

---

## 1. Summary

| Severity | Count | IDs |
|---|---|---|
| critical | 1 | B01 |
| high | 3 | B02, B03, B04 |
| medium | 14 | B05 to B18 |
| low | 19 | B19 to B37 |

The most important problems:

1. **B01 (critical).** `spatialCorrelationGeneExp()` never corrects for multiple testing. `p.adjust()` runs inside the per-gene `lapply`, on a single p-value, so it does nothing. Even so, `adjustMethod = "BH"` is documented, and the kidney vignette says this function "returns two adjusted empirical p-values". In my test, `adjustMethod = "BH"`, `"bonferroni"` and `"none"` returned identical values.
2. **B02 (high).** The empirical p-value is `sum(|r_null| > |r_obs|) / B`. It has no +1 and uses a strict `>`, so it can be exactly 0, and BH leaves a 0 at 0. In a null simulation (B = 19), 12% of independent pairs got pX = 0 and 6% got both pX and pY = 0. All of those stay "significant" after BH at any FDR. The docs say the smallest p-value is 1/B (0.01).
3. **B03 (high).** `spatialCorrelationGeneExpIterPermutations()` crashes with "subscript out of bounds" whenever any gene has NA p-values, for example a gene that is all-zero in one sample, or a single NA value. The NA leaks into the rerun gene list. All first-round work is lost.
4. **B04 (high).** `locfit` dies with "segfault from C stack overflow", killing the whole R session, whenever `delta * N` falls in [1, 2). With the vignette's recommended deltas (`c(0.01, 0.05, seq(0.1, 0.9, 0.1))`), any gene with 100 to 199 shared pixels kills R. With the default grid, 10 to 19 pixels does. With `nThreads > 1` the forked worker dies instead, and the gene silently gets NA.
5. **R CMD check fails** (1 ERROR, 3 WARNINGs, 2 NOTEs). The ERROR is in the `simRanPatternRasts` example. A second broken example (`spatialCorrelationGeneExpWithinSample`) is hidden behind it.

---

## 2. R CMD build / check

Commands:
```
rsync -a --exclude .git --exclude docs "<repo>/" WD/pkgcopy/
R CMD build --no-build-vignettes WD/pkgcopy            # -> STcompare_0.1.0.tar.gz, 48,623,451 bytes
_R_CHECK_FORCE_SUGGESTS_=false BIOCPARALLEL_WORKER_NUMBER=4 BIOCPARALLEL_WORKER_MAX=4 \
  R CMD check --no-manual --ignore-vignettes --no-build-vignettes STcompare_0.1.0.tar.gz
```
(The `BIOCPARALLEL_*` variables only cap examples at 4 workers, to respect the shared machine.)
Full log: `WD/check/check_full.log`, `WD/check/STcompare.Rcheck/00check.log`. Wall time: 52 s.

**Status: 1 ERROR, 3 WARNINGs, 2 NOTEs** (plus INFO: installed size 33.8 Mb; `extdata` 32.3 Mb, `data` 1.2 Mb).

- **ERROR (examples).** The `simRanPatternRasts` example (`R/data.R:97`) calls `assays()` without a namespace: `Error in assays(simRanPatternRasts[[1]]) : could not find function "assays"`. The check stops at the first failing example, so I ran the remaining examples one by one with `WD/examples/run_one.R`:
  - `spatialCorrelationGeneExpWithinSample`: **ERROR**, `unable to find an inherited method for function 'spatialCoords' for signature 'x = "list"'`.
  - `spatialCorrelation`: OK, but took **3.64 min**.
  - `spatialCorrelationGeneExp` (9.9 s), `spatialSimilarity`, `viladomatCorrelation`, `speKidney`: OK.
  - `linearRegression`, `matchingVariograms`, `pixelClass`, `plotCorrelationGeneExp`, `savePlots` ran OK inside the check. After `savePlots`, the check had to detach `patchwork`, `gridExtra` and `ggplot2`.
- **WARNING: non-ASCII characters** in `R/iterativePermutations.R`. This is the em dash at line 302.
- **WARNING: dependencies in R code.** `'::' import not declared from: 'sf'`; `'library' or 'require' calls not declared from: 'ggplot2' 'gridExtra' 'patchwork'`; `'library' or 'require' calls in package code`.
- **WARNING: Rd \usage.** `spatialCorrelationGeneExpWithinSample.Rd` has an undocumented argument `delta`, and documents `deltaX`/`deltaY`, which are not in `\usage`. `spatialSimilarity.Rd` has an undocumented argument `verbose`.
- **NOTE: R code possible problems.** No visible global function definitions for `cor cor.test dist fitted lm median na.omit quantile rnorm p.adjust.methods combn unit plot_layout`. No visible bindings for the NSE variables `x y X Y XGexp YGexp color fill pValuePermuteX`. The check suggests `importFrom("stats", ...)` and `importFrom("utils", "combn")`.
- **NOTE: Rd files.** 89 instances of `Lost braces in \itemize; \value handles \item{}{} directly` across 11 Rd files, plus `savePlots.Rd:47` `"{gene_name}.pdf"` lost braces.
- **Not exercised here but confirmed.** Under `--as-cran` or Bioconductor checks (`_R_CHECK_LIMIT_CORES_=TRUE`), the examples that use `nThreads = 5` would also error: `BiocParallel workers must be <= 2 was (5)` (`WD/12_misc.R`).
- **Vignettes were not built.** They need `MERINGUE` and `scatterbar`, which are not installed and not declared, and they download data from Zenodo or 10x.

---

## 3. Findings table

| ID | Sev | Category | Title | Location |
|---|---|---|---|---|
| B01 | critical | statistical | `p.adjust` applied per gene in `spatialCorrelationGeneExp` is a no-op, so there is no multiple-testing correction | R/spatialCorrelation.R:837 |
| B02 | high | statistical | Empirical p-value `extreme/B` with strict `>` can be 0, and BH cannot fix it; docs say min is 0.01 | R/spatialCorrelation.R:326 |
| B03 | high | correctness | IterPermutations crashes when any gene has NA p-values (NA leaks into rerun list) | R/iterativePermutations.R:65 |
| B04 | high | correctness | `locfit` C-stack-overflow segfault kills R when `delta*N` is in [1,2); no guard | R/spatialCorrelation.R:109 |
| B05 | medium | correctness | Error handler references `corDF`, which does not exist when `cor.test` fails, giving "object 'corDF' not found" | R/spatialCorrelation.R:591 |
| B06 | medium | api | User-supplied `BPPARAM` ignored by `spatialCorrelationGeneExp` | R/spatialCorrelation.R:830 |
| B07 | medium | docs-code-mismatch | Rerun threshold `(alpha/nPermutes)*100` vs documented `alpha/nPermutations[k]` | R/iterativePermutations.R:64 |
| B08 | medium | reproducibility | `set.seed()` inside functions overwrites the user's global RNG stream | R/spatialCorrelation.R:242 |
| B09 | medium | correctness | `deltaX`/`deltaY` given as a numeric vector are silently used one element per gene | R/spatialCorrelation.R:827 |
| B10 | medium | correctness | Pixels matched by name only: separately rasterized inputs silently mis-paired; stale `spatialCoords` rownames crash | R/spatialCorrelation.R:785 |
| B11 | medium | correctness | Gene missing from second object crashes the whole run (no gene intersection) | R/spatialCorrelation.R:822 |
| B12 | medium | statistical | One NA/NaN silently disables the permutation test (NA p-values); `spatialSimilarity` errors on NA | R/spatialCorrelation.R:534 |
| B13 | medium | docs-code-mismatch | `spatialCorrelationGeneExpWithinSample` documents `deltaX`/`deltaY` but takes `delta`; 1-gene input fails | R/spatialCorrelation.R:977 |
| B14 | medium | packaging | R CMD check ERROR: `simRanPatternRasts` example uses unqualified `assays()` | R/data.R:97 |
| B15 | medium | packaging | `spatialCorrelationGeneExpWithinSample` example passes a list and assay "A", so it errors | R/spatialCorrelation.R:969 |
| B16 | medium | correctness | `savePlots` non-geometry branch plots sample 1 in both expression panels | R/visualizationFunctions.R:395 |
| B17 | medium | api | Assay used by `spatialSimilarity` not stored; `savePlots` never forwards `assayName`, so the wrong assay is plotted | R/visualizationFunctions.R:422 |
| B18 | medium | packaging | `savePlots` calls `library()` on undeclared gridExtra/patchwork, attaching them and masking `BiocGenerics::combine` | R/visualizationFunctions.R:351 |
| B19 | low | correctness | `numPixelInThresh = dim(thresh)[1]` (always 1) for skipped genes | R/packageFunction.R:249 |
| B20 | low | api | `getGenePixelDF(assayName = assayName)` gives a recursive-default error when `assayName` is omitted | R/packageFunction.R:17 |
| B21 | low | api | `spatialCorrelation` documents 1 x N matrix input but errors on it | R/spatialCorrelation.R:512 |
| B22 | low | packaging | Undeclared dependencies/imports (sf, stats, utils, vignette packages); `dplyr` needs >= 1.1.0 | DESCRIPTION:26 |
| B23 | low | packaging | `exportPattern("^[[:alpha:]]+")` exports internal helpers; NAMESPACE not roxygen-managed | NAMESPACE:1 |
| B24 | low | packaging | Non-ASCII em dash in R code (check WARNING) | R/iterativePermutations.R:302 |
| B25 | low | packaging | Rd: undocumented `verbose`; 89 "Lost braces" NOTEs | R/packageFunction.R:162 |
| B26 | low | packaging | 48.6 MB tarball / 33.8 MB installed; byte-identical duplicate RData | inst/extdata/brain-MERFISH-10x-visium/brainCorrelation_1.RData:1 |
| B27 | low | packaging | Other example defects (`nThreads = 5` under core limit, 3.6 min example, undefined `negCorrelation`, wrong sample in savePlots example) | R/spatialCorrelation.R:760 |
| B28 | low | docs-code-mismatch | Getting-started vignette: "BH-corrected" and "max of pX/pY" claims are false for the precomputed results | vignettes/getting-started-with-STcompare.Rmd:332 |
| B29 | low | docs-code-mismatch | Kidney vignette: `nThreads` "parallelize genes" (it is permutations); claims `spatialCorrelationGeneExp` returns adjusted p | vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:419 |
| B30 | low | docs-code-mismatch | Assorted doc errors (0.001 vs 1e-4 pseudo-count, "both" vs OR, %, copy-pasted description, N = 5000) | R/packageFunction.R:51 |
| B31 | low | reproducibility | Same seed for every gene and both directions; p depends on source column order; stage 2 recomputes stage 1 | R/spatialCorrelation.R:831 |
| B32 | low | api | `matchingVariograms`: `i` unused, default `seed = 0`, so identical noise in every permutation as in its own example | R/spatialCorrelation.R:97 |
| B33 | low | api | `plotCorrelationGeneExp` drops negative values, mishandles NA in one direction, prints unrounded p | R/spatialCorrelation.R:1145 |
| B34 | low | correctness | `spatialSimilarity` with negative values: proportions sum > 1 and NA pixel IDs; no shared pixels gives NaN | R/packageFunction.R:269 |
| B35 | low | api | Input-validation gaps (3-element list, assay order, `nPermutations` 0/1) and errors `print()`ed to stdout | R/spatialCorrelation.R:586 |
| B36 | low | performance | `getGenePixelDF` densifies the full assay twice per gene (quadratic); `rbind` in loop | R/packageFunction.R:25 |
| B37 | low | performance | Redundant `bplapply` extraction passes, per-call `MulticoreParam`; hard-coded undocumented `N_s = 1000` | R/spatialCorrelation.R:298 |

---

## 4. Detailed findings with evidence

Shared fixtures (made by `WD/01_make_data.R`):
- `WD/rastKidney.rds`: `speKidney` rasterized at resolution 0.2. Note that the list order is A, C, B. Each object has 279 to 287 pixels and a dense `matrix` assay `pixelval`.
- `WD/nullGenes.rds`: two 6-gene SpatialExperiments on 231 shared pixels, built from 12 independent `simRanPatternRasts` fields. Gene `g1` is made strongly correlated (r = 0.998). The other genes are independent nulls.

### B01 (critical, statistical): `spatialCorrelationGeneExp` performs no multiple-testing correction
`R/spatialCorrelation.R:836-842`: `p.adjust()` is called on `output$pValuePermuteX`, which holds a single gene's value, inside the per-gene `lapply` (lines 808-845). For one p-value, BH, Bonferroni and the other methods all return the input unchanged. Yet the `adjustMethod` doc (lines 715-718) says it corrects the final columns. The kidney vignette (`vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:430-431`) says "`STcompare::spatialCorrelationGeneExp` returns two adjusted empirical p-values". Also, `adjustMethod` is not validated up front, so an invalid value errors only after the first gene's permutations. (`spatialCorrelationGeneExpWithinSample` has no correction at all, and does not claim one.)

Evidence (`WD/02_padjust_bpparam.R`, log `WD/02_padjust_bpparam.log`), 6 genes, B = 20:
```
  gene pX_BH pX_none pX_bonferroni pX_BH_expected pX_bonf_expected
1   g1  0.00    0.00          0.00           0.00                0
2   g2  0.55    0.55          0.55           1.00                1
3   g3  0.25    0.25          0.25           0.75                1
...
identical(BH, none): TRUE  identical(bonferroni, none): TRUE
```
Fix: collect all genes, then call `p.adjust()` once per column after `do.call(rbind, ...)`, as `spatialCorrelationGeneExpIterPermutations` does at `R/iterativePermutations.R:343-344`. Validate `adjustMethod` with `match.arg`.

### B02 (high, statistical): empirical p-value can be exactly 0; docs claim min = 1/B
`R/spatialCorrelation.R:326-327`: `extreme <- sum(abs(cor.global) > abs(cor.global.obs)); p.value.global <- extreme / B`. With no +1 and a strict `>`, a gene with no exceedances gets p = 0. Under an exchangeable null this happens with probability 1/(B+1), and `p.adjust` keeps 0 at 0. Those genes stay significant at any FDR, whatever the number of tests. The docs say "Default is 100, such that the smallest p-value is 0.01" (lines 168-169, 357-358, 646-647, 865-866). The standard valid estimator is `(1 + #{|r_b| >= |r_obs|}) / (B + 1)` (Phipson and Smyth 2010).

Evidence:
- Package examples (`WD/check/STcompare.Rcheck/STcompare-Ex.Rout`): kidney A-B and A-C give `pValuePermuteX 0`, `pValuePermuteY 0`.
- Null simulation (`WD/04_pvalue_zero.R`, log `.log`, raw results `WD/04_null_results.rds`). 50 disjoint independent `simRanPatternRasts` pairs, B = 19, default delta grid, 23 s on 4 cores:
```
fraction pX == 0: 0.12   fraction pY == 0: 0.06   fraction both == 0: 0.06
expected P(p==0) under an exact null with strict '>': 0.05
BH-adjusted max(pX,pY) for pairs with p == 0: 0 0 0
With the +1 rule those pairs would get BH-adjusted: 0.8333333 0.8333333 0.8333333
Generic illustration: 10,000 tests, one p = 0 and the rest = 1 -> BH keeps 0
```
- Shipped results `inst/extdata/kidneyCorrelation.RData` (1046 genes; 289 genes finished at B = 100 and 757 at B = 1000): 444 genes have adjusted pX == 0, and 414 have both adjusted p-values == 0.

In the iterative workflow the impact on false discoveries is small, because zeros at B = 100 are re-estimated at B = 1000. With `spatialCorrelationGeneExp` at B = 100, roughly 1% of null genes per direction get p = 0, and no correction can remove them.

### B03 (high, correctness): IterPermutations crashes when any gene has NA p-values
`R/iterativePermutations.R:63-67`: `keep <- pX < t & pY < t` is NA for a gene whose p-values are NA, and `rownames(results_df)[keep]` then returns `NA`. Round 2 calls `.run_spatial_correlation_iteration(genes = c(..., NA))`, which indexes the assay with `NA`, giving "subscript out of bounds". All of round 1 is discarded. NA p-values come from `spatialCorrelation`'s catch-all handler. That happens for genes that are constant (for example, all zeros in one sample's shared pixels, which is common for sparse data), for any NA or NaN value (B12), and for forked workers that crash (B04).

Evidence (`WD/06_iter.R`, `WD/06b_iter_na.R`):
```
.get_genes_to_repermute on a results table where one gene has NA p-values:
[1] "g1" NA
round with genes: g1,g2,g3
round with genes: g1,NA
ERROR: subscript out of bounds
Same data through spatialCorrelationGeneExp (non-iterative) completes with an NA row: ... g3 NA NA NA NA
```
Fix: `keep <- which(pX < t & pY < t)` (or `%in% TRUE`), and handle constant genes explicitly.

### B04 (high, correctness): `locfit` segfault (C stack overflow) kills R for small `delta * N`
`R/spatialCorrelation.R:109-111` calls `locfit::locfit(... lp(long, lat, nn = delta[k], deg = 0), kern = "gauss", maxk = 300)` without checking `delta[k] * N`. When the nearest-neighbour fraction covers between 1 and 2 points, locfit overflows the C stack. R prints `Error: segfault from C stack overflow` and the process exits. `tryCatch` cannot catch it, and neither can the package's handler.

Evidence (`WD/03b_segfault.R`, `WD/03c_segfault_big.R`, log `WD/03c_segfault_big.log`). Every call ran in a fresh process:

| N | nn | result | N | nn | result |
|---|---|---|---|---|---|
| 8 | 0.1 | OK | 8 | 0.2 | **crash** |
| 10, 12, 15 | 0.1 | **crash** | 20, 30 | 0.1 | OK |
| 100 to 199 | 0.01 | **crash** | 200, 201, 250 | 0.01 | OK |
| 30, 35, 39 | 0.05 | **crash** | 40, 41, 60 | 0.05 | OK |
| 1000, 1500 | 0.001 | **crash** | 2000, 2001 | 0.001 | OK |

Through the package API (`WD/03d_edge.R`) with the kidney vignette's grid `c(0.01, 0.05, seq(0.1, 0.9, 0.1))`:
```
N = 150, nThreads = 1 -> Error: segfault from C stack overflow / Execution halted   (R session dies)
N = 150, nThreads = 2 -> <simpleError ... wrong args for environment subassignment>; pValuePermuteX/Y = NA (silent)
N = 210, nThreads = 1 -> OK
```
The original `03_errors.R` run also died this way at N = 10 with the default deltas.
Fix: drop or clamp any delta with `delta * N < 2` (or a safer minimum) and warn. Validate `0 < delta <= 1`.

### B05 (medium, correctness): error handler references `corDF`, which does not exist when `cor.test` fails
`R/spatialCorrelation.R:531-620`: `corDF <- cor.test(...)` is the first statement in the `tryCatch`. If `cor.test` itself errors, the handler (lines 591-592 and 605-606) evaluates `corDF$estimate` and throws a new error, `object 'corDF' not found`. That error propagates out of `spatialCorrelationGeneExp` and aborts the run. Triggers include fewer than 3 shared pixels, an all-NA vector, no shared pixels, and SpatialExperiments without colnames (so `shared_pixels` is NULL).

Evidence (`WD/03_errors.R`, log `WD/03_errors.log`):
```
--- n = 2 observations ---          <simpleError in cor.test.default(...): not enough finite observations>  [1] "ERROR: object 'corDF' not found"
--- all-NA X ---                    same
--- no shared pixels ---            same
--- unnamed pixels (spatialCorrelationGeneExp on SPEs without colnames) --- same
```

### B06 (medium, api): `BPPARAM` ignored in `spatialCorrelationGeneExp`
`R/spatialCorrelation.R:830` passes `BPPARAM = NULL` to `spatialCorrelation()`. That always builds `MulticoreParam(workers = nThreads)`. The `BPPARAM` built at lines 777-779 is never used. A user's `SnowParam` (needed on Windows), `SerialParam` or tuned `MulticoreParam` is silently replaced, despite the docs (lines 699-703). `spatialCorrelationGeneExpWithinSample` (line 1033) and the iterative function both pass it through correctly.

Evidence (`WD/02_padjust_bpparam.R`, via `trace(viladomatCorrelation)`):
```
User passed SnowParam(workers = 3); viladomatCorrelation received: "MulticoreParam workers = 1" "MulticoreParam workers = 1"
Contrast: spatialCorrelationGeneExpWithinSample with SerialParam() passes it through: "SerialParam workers = 1" ...
```

### B07 (medium, docs-code-mismatch): rerun threshold is 100x the documented one
The code (`R/iterativePermutations.R:64`) uses `t <- (alpha / nPermutes) * 100`. The roxygen (lines 90-94) says a gene is carried forward when both p-values are "less than `alpha / nPermutations[k]`". The vignette says "if a gene has an empirical p-value less than 0.05, creates 1000 permutations". The code matches the vignette only when `nPermutations[1] == 100`.

Evidence (`WD/06_iter.R`):
```
nPermutes=   10: code threshold=0.5   -> rerun: p0,p004,p03,p2,p45 | documented 0.005  -> rerun: p0,p004
nPermutes=  100: code threshold=0.05  -> rerun: p0,p004,p03      | documented 0.0005 -> rerun: p0
nPermutes= 1000: code threshold=0.005 -> rerun: p0,p004          | documented 5e-05  -> rerun: p0
```
With `nPermutations = c(10, 20)` the code reruns every gene with p < 0.5.

### B08 (medium, reproducibility): `set.seed()` inside package functions overwrites the global RNG
`matchingVariograms` (line 99) and `viladomatCorrelation` (line 242) call `set.seed()` on the global stream. Any randomness the user draws afterwards is fixed by the package, not by the user's seed.

Evidence (`WD/05_rng.R`):
```
runif(3) after call with user seed 123: 0.5116576 0.08492106 0.9947705
runif(3) after call with user seed 999: 0.5116576 0.08492106 0.9947705
identical -> user's seed has no effect on downstream randomness: TRUE
```
Fix: use `withr::with_seed()`/`local_seed()`, or save and restore `.Random.seed`. Better still, generate per-permutation streams (for example L'Ecuyer or dqrng) without touching global state.
Verified as **not** an issue: the docs say results are reproducible "regardless of parallelization back-end", and that holds. Null correlations were identical for `SerialParam`, `MulticoreParam(1)`, `MulticoreParam(2)` and `MulticoreParam(4)`.

### B09 (medium, correctness): numeric-vector `deltaX`/`deltaY` silently means one delta per gene
`spatialCorrelation()` takes `deltaX` as a numeric vector. `spatialCorrelationGeneExp()` (line 827) and the iterative helper (`R/iterativePermutations.R:44-45`) index it with `deltaX[[i]]`. So a vector such as `c(0.1, 0.5, 0.9)` makes gene 1 use only 0.1 and gene 2 use only 0.5, with no warning. With more genes than elements, the run errors part-way. The list length is never checked.

Evidence (`WD/03_errors.R`):
```
deltaX = c(0.1, 0.5, 0.9) for 2 genes -> unique deltaStarX per gene: [[1]] 0.1  [[2]] 0.5
deltaX = c(0.1, 0.5) for 3 genes      -> "ERROR: subscript out of bounds"
```
Fix: accept a numeric vector as "same grid for all genes" (wrap it in `rep(list(v), G)`), and require a list of length G otherwise.

### B10 (medium, correctness): pixels matched by name only, never by coordinates
`R/spatialCorrelation.R:785-787` (and the same pattern in `iterativePermutations.R:266-268`, `plotCorrelationGeneExp`, `savePlots` and `pixelClass`) pairs pixels by `rownames(spatialCoords())`, then takes `pos` from the first object only. SEraster names pixels by grid index. If the two objects are rasterized in separate calls, their grids differ, so the same name can refer to different locations. The result is silently wrong. A second problem: `rownames(spatialCoords())` are used to index assay columns, but SpatialExperiment does not update them on `colnames<-`, so the two can diverge.

Evidence (`WD/08b_pixels.R`, log `WD/08b_pixels.log`). Here `speKidney$B` was restricted to y > 1.5, then rasterized either separately or jointly with A:
```
separately rasterized: shared names = 166 ; names with different coordinates = 166 ; max offset = 4.291
jointly rasterized: shared names = 205 ; names with different coordinates = 0
correlationCoef: jointly rasterized = -0.9503 (n = 205) ; separately rasterized = -0.2909 (n = 166) -- no warning issued
colnames: px170 px172 px173  rownames(spatialCoords): pixel170 pixel172 pixel173
spatialCorrelationGeneExp ERROR: subscript out of bounds        (spatialSimilarity, which uses colnames, works)
```
(When the bounding boxes happen to be almost identical, as with the full `speKidney` A and B, the separate grids coincide by luck: `WD/08_pixels.log`.)
Fix: match on `colnames`, check that coordinates agree for shared names (error on a mismatch), and document that objects must be rasterized in one `rasterizeGeneExpression(list(...))` call.

### B11 (medium, correctness): a gene missing from the second object aborts the run
`spatialCorrelationGeneExp` loops over `rownames(source)` and indexes `target` by name (`R/spatialCorrelation.R:821-822`; also `R/iterativePermutations.R:38-39`). A gene absent from `target` throws "subscript out of bounds" outside the `tryCatch`, after the earlier genes have already been computed. `spatialSimilarity` uses `intersect(rownames(x), rownames(y))` (`R/packageFunction.R:193`).

Evidence (`WD/03_errors.R`):
```
1: g1
2: g2
[1] "ERROR: subscript out of bounds"
contrast: spatialSimilarity uses intersect(rownames) -> [1] "g1" "g3"
```

### B12 (medium, statistical): a single NA/NaN silently disables the permutation test
`cor.test` (line 534) drops incomplete pairs, so `correlationCoef` and `pValueNaive` are still reported. But `locfit` fails on the NA inside `viladomatCorrelation`, so both empirical p-values come back NA, with the error only `print()`ed. NaN arises naturally from CPM normalization of a zero-count pixel. In `spatialSimilarity`, `quantile()` (`R/packageFunction.R:221,226`) errors before `threshold()`'s `na.omit` (line 57) is reached, unless the user supplies `t1` and `t2`.

Evidence (`WD/03_errors.R`, `WD/07_similarity_plots.R`):
```
--- single NA in X ---  <simpleError in FUN(X[[i]], ...): NA/NaN/Inf in foreign function call (arg 4)>
    correlationCoef pValueNaive pValuePermuteX pValuePermuteY
cor      -0.1606388  0.01473639             NA             NA
spatialSimilarity with one NA: "ERROR: missing values and NaN's not allowed if 'na.rm' is FALSE"
```
Constant (zero-variance) genes also return an all-NA row, after a printed BiocParallel error. In that row `deltaStarX` is logical `NA`, while normal rows hold numeric vectors.

### B13 (medium, docs-code-mismatch): `spatialCorrelationGeneExpWithinSample` argument mismatch
The signature (line 977) has `delta`, while the roxygen (lines 868-889) documents `deltaX` and `deltaY`. R CMD check raises the Rd usage WARNING. Calling it with the documented argument fails: `unused argument (deltaX = list(0.5, 0.5))`. A 1-gene object fails in `combn()` with `n < m`. The `@return` section lists `permutationsX` but not `permutationsY`. (`WD/12_misc.R`)

### B14 (medium, packaging): R CMD check ERROR in the `simRanPatternRasts` example
`R/data.R:97`: `assays(simRanPatternRasts[[1]])$pixelval[1, 1:5]`. SummarizedExperiment is only imported, not attached, so the function is not found: `could not find function "assays"`. This is the check's single ERROR. Fix: `SummarizedExperiment::assay(simRanPatternRasts[[1]], "pixelval")[1, 1:5]`.

### B15 (medium, packaging): `spatialCorrelationGeneExpWithinSample` example cannot run
Lines 969-972 pass `input = rastKidney`, which is a list of 3 objects, and `assayName = "A"`, which is not an assay name. The function expects one SpatialExperiment, and the assay is called `pixelval`. Result (`WD/examples/ex_spatialCorrelationGeneExpWithinSample.log`): `unable to find an inherited method for function 'spatialCoords' for signature 'x = "list"'`. This failure is hidden in R CMD check behind B14. Even with valid input, the object has only 1 gene, so `combn` would fail (B13).

### B16 (medium, correctness): `savePlots` non-geometry branch plots sample 1 twice
`R/visualizationFunctions.R:394-395` builds panel b's `dfb` from the coordinates of `rastGexp[[2]]` but colours it with `assay(rastGexp[[1]], ...)`. Neither panel has a title in this branch.

Evidence (`WD/07_similarity_plots.R`; SEraster objects with the `geometry` column removed):
```
panel b colour values identical to sample A: TRUE ; identical to sample B: FALSE
cor(A,B) for this gene = -0.947 (so panels a and b should look opposite)
```

### B17 (medium, api): plots can show a different assay from the one that was classified
`spatialSimilarity` does not store `assayName` in `$parameters` (`R/packageFunction.R:304-308`). `linearRegression`/`pixelClass` default to assay 1. `savePlots` accepts `assayName` but does not forward it: `pixelClass(spatialSimilarity, gene)` at line 362, `linearRegression(input = spatialSimilarity, gene = gene)` at line 422, and `SEraster::plotRaster(...)` with `assay_name = NULL` at lines 365 and 374. plotRaster's source shows `if (is.null(assay_name)) mat <- SummarizedExperiment::assay(input)`, i.e. the first assay.

Evidence (`WD/07_similarity_plots.R`), with `spatialSimilarity(..., assayName = "lognorm")`:
```
linearRegression(sL, 'Gene') x-values come from assay 'pixelval': TRUE ; from 'lognorm': FALSE
range of plotted x: 4.737711 37.39558  vs lognorm range: 0.7587387 1.584281
savePlots(..., assayName='lognorm') panel 4 x-values from 'pixelval': TRUE
```
The points are coloured by the lognorm classification but placed at pixelval coordinates, so the fold-change guide lines are meaningless.

### B18 (medium, packaging): `savePlots` attaches packages with `library()`
`R/visualizationFunctions.R:351-353` call `library(ggplot2)`, `library(gridExtra)` and `library(patchwork)`. gridExtra and patchwork are not in DESCRIPTION, so `remotes::install_github()` does not install them, and `savePlots` then errors. When they are present, they are attached to the user's search path, and gridExtra masks `BiocGenerics::combine`/`Biobase::combine`. The code relies on that attachment, because `plot_layout` (line 424) and `unit` (lines 370-416) are used without a namespace.

Evidence (`WD/07_similarity_plots.R`):
```
Attaching package: 'gridExtra'
The following object is masked from 'package:BiocGenerics': combine
newly attached: package:patchwork package:gridExtra
```
R CMD check: "'library' or 'require' calls in package code" (WARNING).

### B19 (low, correctness): `numPixelInThresh` is always 1 for skipped genes
`R/packageFunction.R:249` uses `numPixelInThresh = dim(thresh)[1]`, but `thresh` is the 1-row data.frame returned by `threshold()`. It should be `nrow(threshDF)`. The non-skip branch (line 286) is correct.

Evidence (`WD/07_similarity_plots.R`; the gene has 5% non-zero pixels):
```
  gene percentSimilarity numPixelInThresh numPixelOutThresh
2 Sparse                NA                1               247
true number of pixels passing threshold for 'Sparse': 26  (reported: 1 )
```
(1 + 247 is not equal to 273 shared pixels.)

### B20 (low, api): `getGenePixelDF` default argument is self-referential
`R/packageFunction.R:17` has `assayName = assayName`. Calling the exported function without `assayName` gives `promise already under evaluation: recursive default argument reference or earlier problems?` (`WD/07_similarity_plots.R`). The default should be `assayName = 1`.

### B21 (low, api): documented "1 x N matrix" input fails
The roxygen at lines 346-350 says X and Y may be "a 1 x N numeric vector or matrix". `data.frame(X = X, ...)` (lines 512-515) fails before the `tryCatch`: `arguments imply differing number of rows: 0, 231` (`WD/03_errors.R`).

### B22 (low, packaging): dependency declarations
- `sf::` is used in `pixelClass` (line 182) but sf is not declared (check WARNING). gridExtra and patchwork are discussed under B18.
- `stats` functions (`cor`, `cor.test`, `dist`, `fitted`, `lm`, `median`, `na.omit`, `quantile`, `rnorm`, `p.adjust.methods`) and `utils::combn` are used unqualified with no imports (check NOTE).
- `dplyr::case_when(.default = )` is used in `plotCorrelationGeneExp` (line 1131). `.default` arrived in dplyr 1.1.0 (dplyr NEWS.md line 656, section "dplyr 1.1.0"), but DESCRIPTION has `dplyr` with no version.
- The vignettes use `MERINGUE`, `scatterbar`, `rhdf5`, `rjson`, `Matrix`, `BiocGenerics`, `patchwork` and `gridExtra`, none of which are in Suggests. MERINGUE (GitHub-only) and scatterbar are not installed here, so those vignettes cannot be built.
- `class` is in Suggests but is not used anywhere.

### B23 (low, packaging): everything is exported
`NAMESPACE` contains only `exportPattern("^[[:alpha:]]+")` and has no roxygen header. `roxygen2::roxygenise()` therefore leaves it alone, and I confirmed `man/` and `NAMESPACE` are unchanged after regenerating. As a result, `getGenePixelDF`, `threshold` and `assignFill`, which have no `@export` tag, are exported. That includes the very generic name `threshold`.

### B24 (low, packaging): non-ASCII character in R code
`R/iterativePermutations.R:302` contains an em dash (`<e2><80><94>`) in a message string. This causes the check WARNING.

### B25 (low, packaging): Rd problems
`verbose` in `spatialSimilarity` (`R/packageFunction.R:162`) is undocumented (check WARNING). The `@return` sections use `\itemize{\item{a}{b}}`, which produces 89 "Lost braces" NOTE lines across 11 Rd files. `savePlots.Rd:47` has unescaped braces.

### B26 (low, packaging): package size
The source tarball is 48.6 MB and the installed size is 33.8 MB (`extdata` 32.3 MB, `data` 1.2 MB); Bioconductor's limit is 5 MB. `inst/extdata/brain-MERFISH-10x-visium/brainCorrelation_1.RData` is byte-identical to `brainCorrelation.RData` (both MD5 `5661fd72aa53da9ebf788764326808f0`, 2.76 MB each). `images/` (1.3 MB) and the vignette figure directories are not excluded in `.Rbuildignore`.

### B27 (low, packaging): other example defects
- The `spatialCorrelationGeneExp` example (lines 760-761) uses `nThreads = 5`. Under `_R_CHECK_LIMIT_CORES_=TRUE` (as set by `--as-cran`) this errors with `BiocParallel workers must be <= 2 was (5)` (`WD/12_misc.R`).
- The `spatialCorrelation` example (100 permutations, then delta grids of 18 and 25 values) takes 3.64 min. It should be in `\donttest{}` or use a small `nPermutations`.
- The `spatialCorrelationGeneExpIterPermutations` example (inside `\dontrun`) assigns `corr` but prints `negCorrelation` (line 224): `object 'negCorrelation' not found` (`WD/06_iter.R`).
- The `savePlots` example (`R/visualizationFunctions.R:344-345`) computes the similarity on `list(rastKidney$A, rastKidney$B)` but passes the 3-element `rastKidney`, whose order is A, C, B. Panel 2 is therefore kidney C: its title is `C` (`WD/07_similarity_plots.R`). Line 340 also contains a stray `#' #'`.

### B28 (low, docs-code-mismatch): getting-started vignette, simRanPattern section
- `vignettes/getting-started-with-STcompare.Rmd:332` says `corspv_corrected` is "chosen to be the higher of pValuePermuteY and pValuePermuteX". The generating script returns `results$pValuePermuteX` only (`inst/scripts/simRanPatternSpatialCorrelation.R:29`).
- The plot title (line 366) and text call these values "BH-corrected empirical p-values". They are not BH-corrected across the 9,900 pairs, because each pair is a separate 1-gene call.
- The script also replaces zeros with 0.01 (line 68), which the vignette does not mention.
- The vignette text names `spatialCorrelationGeneExpIterPermutation` (missing "s", line 317).

Evidence (`WD/13_vignette_simran.R`):
```
naive p < 0.05 (uncorrected): 0.498 ; after BH (vignette's ~43%): 0.428
empirical p < 0.05 as stored (vignette's ~4%, labelled 'BH-corrected'): 0.039
empirical p < 0.05 after actually applying BH across the 9900 pairs: 0
```
The qualitative conclusion stands, but the comparison mixes BH-corrected naive p-values with uncorrected empirical ones.

### B29 (low, docs-code-mismatch): kidney vignette
Line 419 says `nThreads = 22, # parallelize genes across threads`. Genes are in fact processed sequentially (`lapply` in `R/iterativePermutations.R:25`); parallelism is over the B permutations of one gene (`bplapply(1:B, ...)`, `R/spatialCorrelation.R:290`). Line 430 says `spatialCorrelationGeneExp` returns adjusted p-values, which is false (B01).

### B30 (low, docs-code-mismatch): assorted documentation errors
- `threshold()` doc says zeros become 0.001 (`R/packageFunction.R:51`); the code uses 0.0001 (lines 69-71).
- `numPixelInThresh` is documented as "above the threshold in both experiments" (line 132), but the code keeps `x > t1 | y > t2` (line 60).
- `percentSimilarity` and the dissimilarity columns are described as percentages but are proportions in [0, 1].
- The vignette says "|log2(y/x)| < b", but the code uses `<=` (line 269).
- The `plotCorrelationGeneExp` description (`R/spatialCorrelation.R:1049-1054`) is copied from `WithinSample` ("Function to calculate Pearson's correlation between rows from one SpatialExperiment ...").
- The `viladomatCorrelation` `nullCorGlobal` doc says "between the permutations and X"; it is Y.
- The `simRanPatternRasts` doc says "Each dataset consists of N = 5000 simulated cells" (`R/data.R:37`), but the datasets contain 1,201 to 1,381 cells (sum of `num_cell`). Its colData list omits `type` and `resolution`.
- The variogram subsampling to 1,000 points (lines 251-264) is not documented anywhere.

### B31 (low, reproducibility/statistical): same seed for every gene and both directions
Every gene, and both directions, calls `viladomatCorrelation(seed = seed)` (lines 537-548 and 831). So `set.seed(seed)` then `sample(X)` produces the same permutation index vectors, and `seed + i` the same noise, for every gene and for permute-X versus permute-Y. Two consequences:
- The p-values of different genes share Monte Carlo error, which induces positive dependence. BH tolerates that, but it is undocumented.
- p-values depend on the source's column order.

The iterative stage-2 run also regenerates the stage-1 permutations exactly, which is wasted compute.

Evidence (`WD/05_rng.R`, `WD/08_pixels.R`, `WD/06_iter.R`):
```
permutation #1 index vector identical for gene g2 (permute X) and gene g3 (permute X): TRUE
permutation #1 index vector identical for g2 permute-X and g2 permute-Y: TRUE
shuffling SOURCE columns: correlationCoef identical: TRUE ; pValuePermuteX: 0 0.6 -> 0 0.4   (B = 5)
first 10 null correlations of B=20 run identical to B=10 run: TRUE
```
Shuffling the target's columns gives identical results, so name-based indexing is robust to target column order.

### B32 (low, api): `matchingVariograms` reuses one noise vector
The exported `matchingVariograms(..., i, seed = 0)` ignores `i` and calls `set.seed(seed)`. Called as in its own example (lines 89-92, no `seed`), every permutation receives the same Gaussian noise vector. (`viladomatCorrelation` avoids this by passing `seed + i`.)
Evidence (`WD/14_misc2.R`): `recovered noise vectors for permutation i=1 and i=2 identical: TRUE`.

### B33 (low, api): `plotCorrelationGeneExp` problems
- `xlim(0, max)` and `ylim(0, max)` (lines 1145-1146) silently drop negative values. With z-scored values, only 113 of 231 points were drawn.
- `case_when(pY > pX ~ pY, .default = pX)` returns pX when pY is NA. That is not the stated "greater of the two".
- p_E is printed unrounded (`p_E = 0.333333333333333`).
- A gene missing from the results table gives `r = NA p_E = NA`, with no error.
- `guides(color = guide_legend(override.aes = ...))` and `labs(fill = "Data")` refer to aesthetics that are not mapped, which produces "Ignoring unknown labels" in the check output.

(`WD/09_plotcor.R`)

### B34 (low, correctness): `spatialSimilarity` edge cases
With negative values (scaled data or residuals), `log2(y/x)` is NaN. Subsetting a data.frame with an NA logical index (lines 269-276) returns NA rows, which are counted in all three classes.

Evidence (`WD/16_similarity_negative.R`):
```
percentSimilarity + percentDissimilarityX + percentDissimilarityY = 1.052 (should be 1)
NA pixel IDs inside similarPixelID: 6 of 225
```
With no shared pixels, `percentSimilarity` is `NaN` with no warning (`WD/12_misc.R`).

### B35 (low, api): input validation gaps; errors printed to stdout
- Input with more than 2 SpatialExperiments is accepted; extra elements are ignored.
- With the default numeric `assayName = 1`, objects whose assays are in different orders are compared across different assays with no warning (`pixelval` vs `log`, r = 0.9974).
- `nPermutations = 0` produces a printed BiocParallel error and NA; `nPermutations = 1` gives p = 1.
- `adjustMethod` (in `GeneExp`) and the delta range are not validated.
- `spatialCorrelation` catches every error and `print(cond)`s it (line 586). That output cannot be silenced or detected programmatically, and the resulting NA row has different column types (`deltaStarX` logical instead of numeric vector).

(`WD/12_misc.R`, `WD/14_misc2.R`, `WD/03_errors.R`)

### B36 (low, performance): quadratic densification in `spatialSimilarity`
`getGenePixelDF` runs `as.matrix(assay(x))` on the whole assay, twice per gene (`R/packageFunction.R:25-26`). The output is grown with `rbind` inside the loop (lines 241 and 278).

Evidence (`WD/11b_densify.R`, `WD/11_similarity_perf.R`):
```
G= 3200: as.matrix(full assay)[gene,] 0.077 s/call (98 MB dense) vs assay[gene,] 0.0020 s/call -> 8.2 min vs 0.21 min per run
G=12800: 0.117 s/call (391 MB dense) vs 0.0070 s/call -> 49.9 min vs 2.99 min per run
G=800 x 4000 pixels: spatialSimilarity 27.5 s (34.4 ms/gene, up from 23 ms/gene at G=100)
```

### B37 (low, performance): parallel overhead and hard-coded subsample
- `viladomatCorrelation` uses two extra `bplapply` calls just to pull elements out of a list (lines 298-313).
- `spatialCorrelation` and `viladomatCorrelation` create a fresh `MulticoreParam` per gene and direction (lines 238-240 and 507-509), each forking and tearing down workers.
- Measured overhead (`WD/17_bplapply_overhead.R`): 0.05 s (1 worker) to 0.10 s (4 workers) per call, versus about 0 for `vapply`. That is 18 to 35 min per 10,000 genes.
- The variogram subsample `N_s <- 1000` (line 253) is hard-coded and undocumented, while locfit always fits all N points. I checked this path for consistency (the target and permuted variograms use the same `ids`, and it runs correctly at N = 1,500 with integer counts, `WD/14_misc2.R`), so it is a design and performance note, not a correctness bug.

---

## 5. Candidate issues checked and **not** confirmed as bugs

| Candidate | Result | Evidence |
|---|---|---|
| Results depend on parallel back-end | No: identical nulls for Serial, MC(1), MC(2), MC(4) | `WD/05_rng.log` |
| `locfit maxk = 300` truncating the fit at large N or small nn | No: no warnings; vertex counts and fits identical with `maxk = 1e5` (N up to 10,000, nn down to 0.01) | `WD/10_locfit_maxk.log` |
| Duplicated coordinates | Works; geoR prints "co-locatted data found, adding one bin at the origin" | `WD/03d_edge.R dup` |
| Sparse `dgCMatrix` and `DelayedMatrix` assays | Work in `spatialCorrelationGeneExp` and `spatialSimilarity` | `WD/04_*`, `WD/08_pixels.log`, `WD/15_delayed.R` |
| Target object with columns in a different order | Name-based indexing gives identical results | `WD/08_pixels.log` |
| Variogram subsample when N > 1000 vs locfit on all N | Internally consistent (same `ids` for target and permuted variograms); only a performance/doc point (B37) | code review, `WD/14_misc2.log` |
| `returnPermutations = TRUE` in IterPermutations (list-column row assignment) | Works (231x20, 231x10, 231x10 matrices) | `WD/14_misc2.log` |
| `make.names()` in the kidney script and vignette | Intentional (matches MERINGUE's renaming); 0 names changed in `kidneyCorrelation.RData` | `WD/18_kidney_names.R` |
| `man/` out of date relative to roxygen | No: `roxygenise()` on a copy produced identical `man/` and `NAMESPACE` | `WD/pkgcopy_roxy` |
| "stack imbalance in '::'" warnings during rasterization | Come from SEraster/BiocParallel on this R build, not from STcompare (STcompare has no compiled code) | check Rout |

---

## 6. Notes relevant to a C++ rewrite

- Fix the statistics before porting, so that the port is not tested against wrong reference values:
  - BH across genes (B01).
  - `(1 + #>=)/(B + 1)` p-values (B02).
  - NA-safe rerun selection (B03).
  - Delta guard `delta * N >= 2` (B04).
- A native smoother should reject or clamp neighbourhoods smaller than 2 points, since that region crashes locfit itself.
- Do the RNG per permutation, with independent streams (dqrng or L'Ecuyer substreams) and without touching `.Random.seed` (B08, B31, B32). Consider sharing permutations across genes deliberately (it is cheaper), but document it.
- Validate inputs once at the top:
  - two SpatialExperiments;
  - genes intersected;
  - pixels matched by colnames, with a coordinate check;
  - assay present by name in both objects;
  - deltas in (0, 1];
  - `nPermutations >= 1`;
  - `adjustMethod` checked with `match.arg`.
- Handle constant, NA and too-few-pixel genes explicitly, with a `status` column instead of `print()` (B05, B12, B35).
- Regression fixtures: `WD/nullGenes.rds` (6 genes, 231 pixels: 1 strongly correlated and 5 null) and `WD/rastKidney.rds` are small and fast. `WD/04_null_results.rds` holds null p-values (B = 19) for calibration checks.

## 7. Files produced

Scripts (all in `WD`):
- `00_env.R`
- `01_make_data.R`
- `02_padjust_bpparam.R`
- `03_errors.R`, `03b_segfault.R`, `03c_segfault_big.R`, `03d_edge.R`
- `04_pvalue_zero.R`
- `05_rng.R`
- `06_iter.R`, `06b_iter_na.R`
- `07_similarity_plots.R`
- `08_pixels.R`, `08b_pixels.R`
- `09_plotcor.R`
- `10_locfit_maxk.R`
- `11_similarity_perf.R`, `11b_densify.R`
- `12_misc.R`
- `13_vignette_simran.R`
- `14_misc2.R`
- `15_delayed.R`
- `16_similarity_negative.R`
- `17_bplapply_overhead.R`
- `18_kidney_names.R`
- `examples/run_one.R`

Logs: the matching `*.log` file next to each script; `examples/ex_*.log`.

Data: `rastKidney.rds`, `nullGenes.rds`, `04_null_results.rds`.

Check output: `check/STcompare_0.1.0.tar.gz`, `check/check_full.log`, `check/STcompare.Rcheck/`.

Package copies: `pkgcopy/` (for build and check) and `pkgcopy_roxy/` (roxygen regeneration test).
