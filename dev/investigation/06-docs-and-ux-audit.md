# STcompare: documentation and user-experience audit

**What I looked at:** `README.md`, `DESCRIPTION`, `NAMESPACE`, `_pkgdown.yml`, the four vignettes, the roxygen blocks in `R/*.R`, `man/*.Rd`, and the built site in `docs/`. The repository is at commit `2983c99`, and `docs/pkgdown.yml` says the site was last built on 2026-05-28.
**Ground rules:** This was an investigation only. I did not change anything in the repository. All scripts and outputs are in my work directory (Appendix C).
**Environment:** macOS on an M1 Ultra, R 4.5.2, Bioconductor 3.22, roxygen2 7.3.3, BiocParallel 1.44.0. I used at most 4 cores.

**How I checked things:**
1. I read every documentation source and compared each claim against the code.
2. I ran `R CMD build` and `R CMD check --no-examples --ignore-vignettes` on a copy of the package.
3. I regenerated the man pages with roxygen2 and diffed them against the committed ones.
4. I ran each Rd example in its own R process, with `_R_CHECK_LIMIT_CORES_=TRUE`.
5. I wrote small scripts that test what the documentation says the code does.
6. I inspected the precomputed results shipped in `inst/extdata`.
7. I timed how long one gene takes.

Wherever a claim rests on a run, the output file is cited as `out/NN_*.txt`.

---

## 1. Executive summary

At a high level the science is explained well. The README figure and the contrast with differential gene expression (A and B have the same mean but different patterns; A and C have the same pattern but different means) make the purpose clear quickly, and the three tutorials show realistic end-to-end use. Beyond that, the documentation has four kinds of problems:

- **(a)** several inaccuracies that change statistical conclusions;
- **(b)** broken examples and broken links on the public site;
- **(c)** almost no guidance on choosing parameters, expected runtime or reading the results;
- **(d)** an API with inconsistent conventions and failure modes that make no noise.

The most important findings:

1. **`adjustMethod` does nothing in `spatialCorrelationGeneExp()`.**
   - `p.adjust()` runs inside the per-gene `lapply`, so it only ever sees one p-value (`R/spatialCorrelation.R:836-842`). The docs say it corrects "the final `pValuePermuteX` and `pValuePermuteY` columns separately" (`:715-718`).
   - With 6 genes, the output from `"BH"`, `"bonferroni"` and `"none"` is identical (`out/06`, section C).
   - The shipped `kidneyCorrelationNoIter.RData`, made with this function, has every p-value on the 1/100 grid, which means none were adjusted (`out/01`).
   - `spatialCorrelationGeneExpIterPermutations()` does adjust correctly (`R/iterativePermutations.R:343-344`).
2. **"The smallest p-value is 0.01" is false; the code returns p = 0.**
   - The claim appears at `R/spatialCorrelation.R:169, 358, 647, 866`. The code is `p = #(|r_null| > |r_obs|)/B` (`:326-327`).
   - The Getting Started output shows `pValuePermuteX = pValuePermuteY = 0`.
   - p = 0 for 573 of 1,046 genes in `kidneyCorrelationNoIter`, and for 444 of 1,046 genes in `kidneyCorrelation` even after 1,000 permutations and BH (`out/01`).
3. **The documented screening threshold of the iterative test is 100 times stricter than the code.**
   - The docs say `alpha / nPermutations[k]` (`R/iterativePermutations.R:90-94`), which is 5e-4 for round 1. The code uses `(alpha / nPermutes) * 100` (`:64`), which is 0.05 (`out/06` H).
   - Under the code's rule, 72% of AKI genes (757 of 1,046) were rerun at 1,000 permutations. Counting round 1, the shipped analyses cost 5.5 to 9.2 times a single 100-permutation round (`out/12` prints 5.1–8.4× because it counts only the final round).
4. **`BPPARAM` is silently ignored by `spatialCorrelationGeneExp()`.**
   - The function hard-codes `BPPARAM = NULL` when it calls the engine (`R/spatialCorrelation.R:830`). Passing a character string as `BPPARAM` raises no error (`out/06` D).
   - This contradicts the docs, which say BPPARAM overrides `nThreads` (`:690-697`, and 4 other copies).
   - Separately, `nThreads` starts forked *processes* that work on *permutations within one gene*. Genes are always processed one at a time. The AKI vignette's comment "parallelize genes across threads" is wrong.
5. **The Getting Started figure on null calibration is mislabeled.**
   - It is titled "BH-corrected … empirical p-value (pE)" and the text says the value was "chosen to be the higher of `pValuePermuteY` and `pValuePermuteX`" (`vignettes/getting-started-with-STcompare.Rmd:332, 366`).
   - The script that produced the data uses only `pValuePermuteX`, with no BH correction (`inst/scripts/simRanPatternSpatialCorrelation.R:22-29`). All values sit on the 1/1000 grid.
   - If BH were actually applied, the "~4%" would be 0% (`out/01`).
6. **Several failures produce no error or warning.** All verified in `out/06`, `out/07` and `out/09`:
   - Inputs on mismatched pixel grids are compared anyway. Same pixel ID, coordinates up to 4.16 units apart, `r = 0.238`, no warning.
   - Engine errors are `print()`ed and turned into `NA`.
   - If a gene fails in round 1, `spatialCorrelationGeneExpIterPermutations()` crashes in round 2 with `subscript out of bounds`.
   - The functions reset the global random-number generator (RNG): `runif(1)` after a call is the same whatever seed the user set.
   - A numeric `deltaX` is quietly used one element per gene.
7. **Examples:**
   - 3 of the 12 example blocks fail. Two of those failures are visible on the live site (`spatialCorrelationGeneExpWithinSample`, `simRanPatternRasts`).
   - The `spatialCorrelationGeneExp` example fails whenever cores are limited, as under `--as-cran` and on Bioconductor builders (`nThreads = 5`).
   - The `spatialCorrelation` example takes **219 s**.
   - The iterative-function example is entirely inside `\dontrun` and prints an object that doesn't exist (`out/05`).
8. **`R CMD check` gives 3 WARNINGs and 2 NOTEs even with examples and vignettes skipped** (`out/02`).
   - About 90 "Lost braces" in `\value` sections make help pages render as "‘correlationCoef’Pearson's correlation coefficient."
   - Roxygen markdown is off, so backticks show up literally.
   - `man/spatialCorrelationGeneExpIterPermutations.Rd` is no longer managed by roxygen because of a leading blank line (`out/03_roxygenise.txt`).
9. **Front door:**
   - 8 links on the site are broken because the URLs are quoted inside markdown link syntax (`href="%22https://…%22"`).
   - The README says alignment is done "using SEraster"; it should be STalign.
   - `DESCRIPTION` has no `biocViews`, so `remotes::install_github()` won't add the Bioconductor repositories that some dependencies need.
   - The Title says "Comparision", and `citation("STcompare")` prints that typo instead of the Bioinformatics paper.
10. **No runtime guidance anywhere.** Measured with default settings: about 9 s per gene for 273 pixels on 1 core, about 3 s on 4 workers, and about 44 s, 73 s and 86 s per gene for 1,000, 2,500 and 5,000 pixels (`out/08`). The authors' own runs took between 7 minutes and 16.8 hours with 20–22 workers (vignette and script comments).
11. **The δ grid often runs out at the low end.** In the shipped brain results, 58% of X-direction permutations (and 85% for cell types) picked the smallest δ in the grid (`out/12`). Users are never told to look for this, and nothing warns them.
12. **Interpretation guidance lives only in vignette prose**, never in the function docs or the output. This covers using max(pX, pY), what δ* means, and that `foldChange` is on the log2 scale (1 means 2-fold).

The prioritized action list is in section 9.

---

## 2. What a new user meets today

- **README / home page:** purpose, an overview figure, two one-paragraph test descriptions, an install command, three tutorial links and a citation. There is no runnable quick start, no statement of the exact input requirements, and no runtime note. It names the wrong alignment tool.
- **Navbar:** Install, Tutorials (3) and Functions. There is no GitHub link, no issue tracker, no news page, and no conceptual or how-to guides.
- **Install article:** one `remotes::install_github()` line. It doesn't mention the Bioconductor dependencies.
- **Function reference:** one flat list called "All functions" with 17 entries.
  - Internal helpers (`assignFill`, `getGenePixelDF`, `threshold`, `matchingVariograms`, `viladomatCorrelation`) sit next to the main entry points.
  - Seven core entries have the function name as their title, e.g. "spatialCorrelationGeneExp — spatialCorrelationGeneExp", so the index doesn't say what anything does.
- **Tutorials:**
  - *Getting Started* uses simulated data. It is quick to read but contains statistical mislabels (§3.1).
  - *AKI (Visium)* and *Brain (MERFISH vs Visium)* depend on Zenodo downloads, MERINGUE (GitHub-only), scatterbar, rhdf5 and `devtools::load_all()`. Their expensive steps are `eval=FALSE` and load precomputed `.RData` files.
- **Reference pages on the live site show three things that would confuse anyone:**
  - two example blocks ending in errors (`spatialCorrelationGeneExpWithinSample`, `simRanPatternRasts`);
  - an example that prints the full permutation vector (`matchingVariograms`);
  - plots with "Ignoring unknown labels" and "Removed n rows" noise.

---

## 3. Accuracy: where the docs disagree with the code

Each item gives the documentation location, the code location and the evidence.

### 3.1 Statistical claims (largest impact)

| # | Documentation says | Code does | Evidence | Fix |
|---|---|---|---|---|
| S1 | "Default is `100`, such that the smallest p-value is 0.01" (`R/spatialCorrelation.R:169, 358, 647, 866`) | `extreme <- sum(abs(cor.global) > abs(cor.global.obs)); p <- extreme / B` (`:326-327`), so p = 0 is possible. The strict `>` also drops ties. | Getting Started output shows `pValuePermuteX = 0, pValuePermuteY = 0` for both pairs (`docs/articles/getting-started-with-STcompare.html`; `out/06` B). 573/1,046 (X) and 686/1,046 (Y) zeros in `kidneyCorrelationNoIter`; 444/531 zeros in `kidneyCorrelation` after 1,000 permutations and BH; 90/80 in `brainCorrelation` (`out/01`). The vignette even works around it ("Avoid log(0)", `getting-started…Rmd:340-342`), and the simulation script replaces 0 with 0.01 (`simRanPatternSpatialCorrelation.R:68`). | Use `(extreme + 1) / (B + 1)` with `>=` (Phipson & Smyth 2010), so the minimum is 1/(B+1). Document the floor and its consequences for BH. If the formula is kept, at least correct the sentence. |
| S2 | `adjustMethod`: "multiple-testing correction method passed to `stats::p.adjust()` for the final `pValuePermuteX` and `pValuePermuteY` columns separately" (`R/spatialCorrelation.R:715-718`). The AKI vignette says "`spatialCorrelationGeneExp` returns two adjusted empirical p-values" (`acute-kidney…Rmd:430`). | `p.adjust()` is applied to a 1-row data.frame inside the per-gene `lapply` (`R/spatialCorrelation.R:836-842`), so it never changes anything. `@return` (`:727-730`) doesn't say the values are adjusted, while the iterative function's `@return` does (`R/iterativePermutations.R:180-183`). | 6-gene test: `identical(BH, none) == TRUE` and `identical(bonferroni, none) == TRUE`. Manual BH changes g5 from 0.80 to 0.95 (`out/06` C). `kidneyCorrelationNoIter` p-values are 100% on the 1/100 grid (`out/01`). An invalid method is only caught *after* the first gene has been computed (`out/07` S). | Move the adjustment after the `rbind`, check `adjustMethod` up front, and add a NEWS entry flagged as a behaviour change. |
| S3 | Iterative screening: "a gene is carried forward only when both … are less than `alpha / nPermutations[k]`" (`R/iterativePermutations.R:90-94`) | `t <- (alpha / nPermutes) * 100` (`:64`), i.e. 0.05 at k = 1 and 0.005 at k = 2 of `c(100, 1000, 10000)`. Screening uses *unadjusted* p-values, but the final table is adjusted. The AKI vignette describes the code's behaviour ("less than 0.05", `:383`). | With toy p of 0.03, 0.04 and 0.0004, both genes are carried forward. The docs would carry forward neither the 0.03/0.04 gene nor anything above 5e-4 (`out/06` H). 757/1,046 kidney genes, 147/325 brain genes and 397/483 MERFISH genes were rerun at 1,000 (`out/01`). Including round 1, that is 8.2×, 5.5× and 9.2× the cost of one 100-permutation round (`out/12` prints 7.5×, 5.1× and 8.4× because it counts only the final round). | Decide on the rule, document it as a formula with a worked example, and add an `nPermutations` column. The authors currently infer it from `length(deltaStarX)` (`inst/scripts/inspect_kidney.R:44-48`). |
| S4 | Getting Started: `corspv_correct` "is the permuted p-value for that pair (chosen to be the higher of pValuePermuteY and pValuePermuteX)" (`:332`). Plot title: "BH-corrected Analytical p-value (pA) vs BH-corrected empirical p-value (pE)" (`:366`). Text: "~43% … vs ~4%" (`:396-397`). | The generating script returns only `results$pValuePermuteX` (`inst/scripts/simRanPatternSpatialCorrelation.R:29`). Each pair is run as a 1-gene input, so the empirical values are not BH-adjusted. The vignette applies BH only to the naive values (`:338`). Each unordered pair is counted twice (9,900 ordered pairs = 4,950 × 2). | `corspv_corrected` is 100% on the 1/1000 grid. Naive p < 0.05 for 42.8% after BH and 49.8% raw. Empirical p < 0.05 for 3.88% raw and **0% after BH** (`out/01`). | Relabel the figure ("raw empirical vs BH naive"), or better, compare like with like (raw vs raw and BH vs BH). Use max(pX, pY) if that is what the text says. Deduplicate the pairs. |
| S5 | "The default for minQuantile is 0.05, meaning by default the pixels at the bottom 5% are removed" (`getting-started…Rmd:418`) | A pixel is kept if `x > t1 | y > t2` (`R/packageFunction.R:60`), so it is removed only when it is low in *both* datasets. | 0 of 273 pixels removed for A vs B, and 5 of 277 for A vs C (`out/07` I3). | Say "pixels whose expression is at or below the 5% quantile in **both** datasets are excluded". |
| S6 | `foldChange`: "Fold-change threshold … Default is 1" (`R/packageFunction.R:116-117`). Vignette: "foldChange is the number of fold that are considered similar. The default is 1 fold" (`getting-started…Rmd:421`). Same vignette: "similar if \|log2(y/x)\| < b … b = 1 … less than two-fold" (`:405-406`). | The band is `log2(y/x) >= -foldChange & <= foldChange` (`R/packageFunction.R:269`). So `foldChange = 1` means up to and including 2-fold, and `foldChange = 2` means 4-fold. | `log2(2/1) = 1` counts as similar (`out/06` J). | Rename to `log2FoldChange` or `maxLog2FC` (keep `foldChange` as a deprecated alias), state "1 = 2-fold" everywhere, and use ≤. |
| S7 | Which p-value to report: only the Getting Started prose recommends "the higher of `pValuePermuteY` or `pValuePermuteX`" (`:230`). The AKI and brain vignettes require both < 0.05 (`acute…Rmd:431`). | None of the function docs say this. `plotCorrelationGeneExp()` computes the max internally (`R/spatialCorrelation.R:1127-1131`) but never returns it. | Example outputs on the site: quakes gives pX = 0.03 and pY = 0.13 with the default δ grid, but pY = 0.02 with a different grid (`docs/reference/spatialCorrelation.html`). Users are left to guess. | Return a combined `pValue = pmax(pX, pY)` and a single `padj`. Put the rule and the reasoning in `@return` and in a FAQ. |
| S8 | Iterative `@return`: "each row reflects the last permutation round in which that gene was evaluated" (`R/iterativePermutations.R:172-175`) | p-values computed with different B are mixed in one table, with no column saying which B was used. The first 100 permutations of the 1,000-permutation rerun are, by construction, identical to round 1 (same seed and order), so round 1 is thrown away rather than reused. | `out/01` (nPermutations per gene obtained by counting `deltaStarX`) | Add `nPermutations` and, ideally, a Monte Carlo SE column. Document. |

### 3.2 Parallelism and reproducibility claims

| # | Documentation says | Code does | Evidence |
|---|---|---|---|
| P1 | "If `BPPARAM` argument is not `NULL`, the `BPPARAM` argument would override `nThreads` argument" (`R/spatialCorrelation.R:695-697`, repeated in 4 other places) | `spatialCorrelationGeneExp()` builds a `BPPARAM` (`:777-779`) that is never used, and then calls `spatialCorrelation(..., BPPARAM = NULL, ...)` (`:830`). | `BPPARAM = "this is not a BiocParallelParam"` gives **no error** (`out/06` D). |
| P2 | `nThreads`: "Number of threads … We recommend setting this argument to be the number of cores available" (`:690-697`). AKI vignette: `nThreads = 22, # parallelize genes across threads` together with `BPPARAM = BiocParallel::MulticoreParam()` (`acute…Rmd:419-420`) | `MulticoreParam` forks *processes*. The parallel loop is over *permutations within one gene* (`R/spatialCorrelation.R:290-295`); genes run one at a time (`:808`, `R/iterativePermutations.R:25`). When BPPARAM is given, `nThreads` is ignored (so the vignette's 22 has no effect). On Windows, `MulticoreParam()` warns "not supported on Windows, use SnowParam()" and sets `workers = 1` (checked in BiocParallel 1.44 source). Because `spatialCorrelation()` builds a new param for every gene (`:507-509`), Windows users get that warning every time a gene is processed. Two further `bplapply` calls fork workers just to `cbind` results (`:298-313`). | Code reading; `out/06` D |
| P3 | `seed`: "Seed for the random number generator used to generate noise in the variogram matching step. Ensures reproducibility … regardless of parallelization back-end" (`:710-713`) | `set.seed(seed)` at the start of every `viladomatCorrelation()` call (`:242`) and `set.seed(seed + i)` in every `matchingVariograms()` call (`:99`). The seed also controls the shuffle order and the variogram subsample (`:255, 286-288`). The **same index permutations** are reused for the X and Y directions (`:537-548` pass the same seed) and for **every gene** (`:831`). The user's global RNG stream is reset. | After a call, `runif(1)` gives 0.462357 whether the user set seed 1 or seed 2 (`out/06` E). Identical permutation indices (`out/06` F). |

### 3.3 Argument and return documentation that is wrong or missing

| # | Location (docs ↔ code) | Problem | Evidence |
|---|---|---|---|
| D1 | `spatialCorrelationGeneExpWithinSample`: `@param deltaX`/`deltaY` (`R/spatialCorrelation.R:868-889`) vs the argument `delta` (`:977`) | The arguments in the docs don't exist and the real one isn't documented. `R CMD check` WARNING. The example passes a *list* of SPEs plus `assayName = "A"` (`:969-972`) and errors on the live site. With a 1-gene SPE the function fails with a cryptic `n < m` from `combn` (`:1007`). | `out/02`, `out/05`, `out/06` O, `docs/reference/spatialCorrelationGeneExpWithinSample.html` |
| D2 | `spatialSimilarity`: `verbose` (`R/packageFunction.R:162`) | Not documented. `R CMD check` WARNING. | `out/02` |
| D3 | `percentSimilarity`: "Percentage of similar pixels"; `percentDissimilarityX/Y`: "Percentage…" (`R/packageFunction.R:126-128`) | These are proportions in [0, 1] (`:270-276`). The three add up to 1. | `out/06` I |
| D4 | `numPixelInThresh`: "Number of pixels above the threshold in both experiments" (`:132`) | It is pixels above the threshold in *either* experiment (`:60`). For genes that fail `minPixels` it is always **1**, because `dim(thresh)[1]` measures a 1-row data.frame (`:249`). | NA case reports 1, true count 16 (`out/07` I2) |
| D5 | `minPixels`: "If less than this percentage of pixels…" (`:112-115`) | A proportion of the shared pixels (`:239`). "Percentage" and the default of 0.1 contradict each other. | Code |
| D6 | `threshold()`: "Zero values are approximated to 0.001" (`:51`) | The code uses `0.0001` (`:69-71`). The pseudocount is fixed and not adjustable, although it changes `log2(y/x)` for pixels with zero expression. | `out/06` M |
| D7 | `getGenePixelDF()`: return columns listed as pixel, y, x (`:11-14`); default `assayName = assayName` (`:17`) | The actual order is pixel, x, y. The self-referencing default errors whenever `assayName` isn't given ("promise already under evaluation"). | `out/06` L |
| D8 | `pixelClass`: "classifying pixels into three categories" (`R/visualizationFunctions.R:124-126`) | There are four categories (`:136-142, 297-301`). | Code |
| D9 | `plotCorrelationGeneExp`: `@description` is copied from the within-sample function ("Function to calculate Pearson's correlation between rows from one SpatialExperiment…", `R/spatialCorrelation.R:1049-1054`) | It is a plotting function. `ylim(0, …)` (`:1145-1146`) silently drops negative values (e.g. scaled data), and `fill = "Data"` has no matching aesthetic, which causes "Ignoring unknown labels" (`:1158`; visible on the site). | Site output |
| D10 | `linearRegression`: "Generates linear regression plot" (`R/visualizationFunctions.R:4-8`) | No regression is fit. It is a scatter plot with y = x and ±log2FC lines, and the axes are cut at the 95% quantile, so points are dropped with warnings (`:68-71, 106-107`). | Site output: "Removed 14 rows" |
| D11 | `savePlots`: "Each gene generates a four-panel figure that is saved as a PDF file … 300 DPI" (`:309-311, 332-334`) | Files are written only if `filePath` is given (default `FALSE`, `:349`), and DPI means nothing for PDF. In the branch without geometry, panel 2 plots dataset **1** again (`:394-395`). `assayName` isn't passed to `pixelClass()` or `linearRegression()` (`:362, 422`). The function calls `library(ggplot2/gridExtra/patchwork)` inside package code (`:351-353`). | `out/02` |
| D12 | `viladomatCorrelation` `@return nullCorGlobal`: "correlation coefficients between the permutations and X" (`R/spatialCorrelation.R:196-198`) | It is between permutations of X and **Y** (`:323`). | Code |
| D13 | `matchingVariograms`: `long` = "y-coordinates", `lat` = "x-coordinates" (`:10-12`); `delta` = "percentage of neighbors"; `i` = "the ith permutation" | The names are the reverse of the usual convention (harmless only because distances are isotropic). δ is the nearest-neighbour *fraction* used as locfit's adaptive bandwidth. `i` is never used in the function body (`:96-140`). | Code |
| D14 | `deltaX`/`deltaY`: "`list`: List of single numerics or list of numeric vectors … length … same as the number of rows" (`:649-662`) | Never checked. A plain numeric vector means gene *i* silently gets the single value `deltaX[[i]]`. A list that is too short fails with a bare "subscript out of bounds". | `out/06` G |
| D15 | `input`: "observations at the same coordinate location in both datasets should have the same row names … use `SEraster::rasterizeGeneExpression()`" (`:633-642`) | Correlation matches pixels by `rownames(spatialCoords())` and takes coordinates **only from X** (`:785-787`). Similarity matches by `colnames()` (`R/packageFunction.R:20`). Coordinates are never compared. Correlation assumes every gene of X exists in Y and fails with "subscript out of bounds" if one doesn't (`out/06` P). Similarity uses the intersection of genes (`R/packageFunction.R:193`). The docs never say to rasterize **both objects in one call**. | Mismatched grids: same pixel ID but coordinates up to 4.16 units apart, `r = 0.238` on 137 "shared" pixels, no warning (`out/09`). With `speKidney` the bounding boxes are nearly identical, so separate rasterization happens to line up (`out/09b`). |
| D16 | `spatialSimilarity` `@return parameters`: "A list of input parameters used in the computation" (`R/packageFunction.R:139-141`) | It stores only `foldChange`, `minPixels` and the **entire input objects** (`:304-308`). It does not record `t1`, `t2`, `minQuantile` or `assayName`, so the plotting helpers can't know which assay was used. | Result is 1,072.8 KB for one gene, versus 1,026.3 KB for the two inputs (`out/06` I) |
| D17 | `simRanPatternRasts` example uses `assays(...)` without a namespace prefix (`R/data.R:97`) | Fails with "could not find function 'assays'". Visible on the live site. | `out/05`, site |
| D18 | Iterative example: `corr <- …` and then prints `negCorrelation` (`R/iterativePermutations.R:218-224`) | Prints an undefined object. Wrapped in `\dontrun`, so nothing ever catches it. | `out/05` (0.0 s) |
| D19 | `plotCorrelationGeneExp` `@param geneName`: "specifiying … in both a SpatialExperiments" (`:1064-1065`); `@return` "ggplot grob" (`:1072`) | Typo. It returns a ggplot object, not a grob. | — |

### 3.4 README and vignette prose that disagrees with reality

- **README:**
  - "datasets first must be first aligned using `SEraster`" (`README.md:22`). Alignment is done with STalign (a Python package); SEraster rasterizes. "are coincide" is a grammar slip.
  - The site home page repeats the same text.
- **Install instructions:**
  - `require(remotes); remotes::install_github(...)` appears in `README.md:26-29`, `Install.Rmd:29-31` and the Getting Started install section without `build_vignettes = FALSE` (`:39-41`). It assumes `remotes` is already installed.
  - Because `DESCRIPTION` has no `biocViews:`, `remotes` doesn't add the Bioconductor repositories (`remotes:::is_bioconductor <- function(x) !is.null(x$biocviews)`), so on a fresh system SpatialExperiment, SummarizedExperiment, BiocParallel and SEraster won't resolve.
  - Recommend `BiocManager::install("JEFworks-Lab/STcompare")` and adding `biocViews:`.
- **Getting Started:**
  - "To visualize the results of comparing A and C" introduces a chunk that plots A vs B (`:435-453`).
  - The function is misnamed `spatialCorrelationGeneExpIterPermutation` (`:317`).
  - The text refers to a test called "`SpatialSimilarity`" (`:167`).
  - Numbers are hard-coded (`-0.9472813`, `0.9431531`, `S=0.535`) instead of using inline R.
  - The figure statistics are wrong (S4 above).
- **AKI vignette:**
  - Says Getting Started showed STcompare "works on both simulated and real ST data" (`:37`). Getting Started uses only simulated data.
  - Cites the bioRxiv preprint (`:39`), while the README cites the 2026 Bioinformatics paper.
  - Refers to `spatialCorrelationGeneExp` when the analysis used the iterative function (`:430`).
  - `compute_dissimilarity(x = rast$ctrl[gene, ], y = rast$AKI[gene, ])` (`:556-564`) uses list names that don't exist. `rast` is `list(AKI_ctrl, AKI_aki)`, so `rast$ctrl` is `NULL` and `rast$AKI` is `NULL` (an ambiguous partial match). Both percentages become `NaN`, so "gene2 = smallest difference" is really just the first negatively correlated gene (reproduced in `out/10`). The columns it tries to compute already exist as `ss$similarityTable$percentDissimilarityX/Y`.
  - `deltaList` is defined twice (`:405, 412`).
  - The correlation and similarity tables are joined by position, not by gene name (`:517-519`). The order happens to agree, but this is fragile.
  - `devtools::load_all()` runs inside the vignette (`:32`).
- **Brain vignette:**
  - The prose says "56% (126/226) … 90% (89/99)" (`:360`). The rendered outputs on the site are **0.513 (118/230)** and **0.895 (85/95)**.
  - The comment says "324 genes" (`:286`); the output is 323.
  - The sentence at `:707` stops mid-way ("…we expect to see where the cell type is located in").
  - The comment says "brainCorrelation is saved at" just before loading `ctCorrelation` (`:657`).
  - Joins are by position (`:444-446`).
  - Exactly nine cell-type plots are hard-coded (`pltlist[[1]]…[[9]]`, `:716-729`).
  - `c(plt1, plt2)` is used on ggplot objects (`:713`); it should be `list()`.
  - `names(annot) <- names(df)` assigns column names where row names were meant (`:68-69`).
- **Citations disagree:**
  - The README cites the Bioinformatics 2026 paper, with dos Santos Peixoto.
  - The vignettes cite bioRxiv 2025, without him.
  - `citation("STcompare")` prints an auto-generated package citation with the "Comparision" typo, because there is no `inst/CITATION`.

---

## 4. Clarity for a new user

| Question a new user has | Where it is answered now | Gaps | What to do |
|---|---|---|---|
| **What is it for?** | README overview, figure, DGE contrast; the Getting Started introduction | The "Spatial Fold Change" description is one vague sentence. "Person correlation assumes that each sample is independent…" mixes up the coefficient with its p-value. Nothing says when to use which test, or that they are complementary (pattern vs magnitude). | Add a two-row "Which test answers which question?" table: correlation tests whether the *pattern* agrees and gives a p-value; similarity measures whether *magnitudes* at matched pixels agree and gives a descriptive score with no p-value. |
| **What input does it need?** | `@param input` (`R/spatialCorrelation.R:633-642`); "Input formatting" in Getting Started (`:172-180`) | It never says that both objects must be **rasterized together in one SEraster call**, which pixel identifier is used (`colnames`, equal to `rownames(spatialCoords)`), that the default assay is "the first one" (SEraster calls it `pixelval`), that genes must exist in both, or that coordinates are taken from X. Nothing in the code checks any of this. | Add an "Input checklist" box to the README, the Getting Started article and `@param input`. Back it with a `checkComparable()` validator (§7.5). |
| **How do I get there from raw data?** | Links to an STalign notebook and a SEraster formatting guide (both broken by the quoting). The AKI vignette loads coordinates that were aligned in advance. | Nothing says STalign is Python, or shows how the aligned coordinates come back into R. There is no single raw → aligned → rasterized → normalized walkthrough. Normalization varies between vignettes (CPM after sum-rasterization; log10 of library-normalized mean) and is never discussed as a choice. | Add a short "From raw data to STcompare input" article with code, covering alignment, joint rasterization and normalization. |
| **How do I choose parameters?** | Resolution: a good qualitative paragraph (`acute…Rmd:244`). δ grid: "here we add two smaller values" (`:381`). nPermutations: "more is more accurate" (`:383`). Similarity: code comments only (`getting-started…Rmd:413-421`). | Nothing quantitative. There is no advice on how many pixels to aim for, how B relates to the smallest achievable p-value and to the number of genes, how to diagnose the δ grid, or what `maxDistPrctile` does. It never says that **similarity must be computed on linear-scale values**: the brain vignette computes `log2(y/x)` on `log10(x+1)` values (`brain…Rmd:413`), which is a ratio of logs, not a fold change. | Write a parameter guide (§8.2). Back it with outputs that flag δ* hitting the edge of the grid. 58% (brain X direction) and 85% (cell types) of permutations hit the low end of the default grid (`out/12`), while the authors' own AKI and MERFISH runs extended the grid to 0.01 and 0.05. |
| **How do I read the output?** | The bullet list in Getting Started (`:227-236`); the Rd `@return` sections | It doesn't explain why pX and pY differ, or which one to report (only in prose). It doesn't say what δ* means or what a δ* at the edge implies, what the null-correlation columns are for, that p can be 0, or that `pValuePermute*` is adjusted in one function and not in the other. Similarity S, `percentDissimilarityX` (X higher) and `percentDissimilarityY` (Y higher) aren't defined in plain language. | Return a combined p-value and `padj`, add an "Interpreting results" section and FAQ entries, and rename the columns (§7.3). |
| **How long will it take, and how do I parallelize?** | Only inline comments (`# Time difference of 1.777397 hours`, `brain…Rmd:316`; `6.951119 mins`, `:646`). Scripts mention 16.8 h and 7.65 h (`biological-replicates-example.R:184, 563`). | No cost model. Nothing says genes are processed one at a time while only permutations run in parallel, nothing about Windows, nothing about memory with `returnPermutations`, and no way to resume an interrupted run. | Write a performance guide (§8.3) with the measured numbers below, plus `estimateRuntime()`. |

**Measured runtime** (default settings: B = 100, 9 δ values, both directions; `out/08`):

| Shared pixels N | 1 core | 4 workers |
|---|---|---|
| 273 (`speKidney` A vs B) | 9.1 s/gene | 3.0 s/gene |
| 1,000 (synthetic) | ~44 s/gene | — |
| 2,500 | ~73 s/gene | — |
| 5,000 | ~86 s/gene | — |

- The synthetic rows were extrapolated ×10 from B = 10.
- B = 1,000 costs about 10× more.
- The iterative default multiplied total cost by 5.5–9.2× in the shipped analyses (`out/12`).
- As a worked example, 1,000 genes at 1,000 pixels with `c(100, 1000)` is roughly 44 s × 8 ≈ 6 minutes per gene on 1 core, or about 4 days (about 5 h with 20 workers if scaling were perfect). That is the same order of magnitude as the authors' 16.8 h for 483 genes on 20 workers.

---

## 5. Typos, grammar, broken links and examples

### 5.1 `R CMD check` on a copy (no examples, no vignettes): 3 WARNINGs and 2 NOTEs (`out/02`)

- **WARNING:** non-ASCII character in R code. The em dash is inside a `message()` string (`R/iterativePermutations.R:302`).
- **WARNING:** undeclared dependencies.
  - `sf::` is used but not declared (`R/visualizationFunctions.R:182`).
  - `library(ggplot2)`, `library(gridExtra)` and `library(patchwork)` are called inside package code (`:351-353`).
  - `gridExtra` and `patchwork` aren't declared at all.
- **WARNING:** Rd `\usage`. `spatialCorrelationGeneExpWithinSample` has `delta` undocumented while `deltaX`/`deltaY` are documented; `spatialSimilarity` has `verbose` undocumented.
- **NOTE:** undefined globals.
  - `stats` functions are not imported (`cor`, `cor.test`, `dist`, `fitted`, `lm`, `median`, `na.omit`, `quantile`, `rnorm`, `p.adjust.methods`), nor `utils::combn`.
  - There are unbound NSE variables (`x`, `y`, `color`, `fill`, `XGexp`, …), plus `unit` and `plot_layout`.
  - NAMESPACE is hand-written (`exportPattern("^[[:alpha:]]+")`), so it exports every helper, including one with the generic name `threshold`.
- **NOTE:** about 90 "Lost braces in \itemize; \value handles \item{}{} directly".
  - Effect on the help pages: in `?spatialCorrelationGeneExp` each name runs straight into its description, e.g. "‘correlationCoef’Pearson's correlation coefficient."
  - Also, `Roxygen: list(markdown = TRUE)` is not set, so backticks render literally ("B is \`nPermutations\`").
- **INFO:** installed size 33.8 MB, of which `extdata` is 32.3 MB. `R CMD build` produces a **48.6 MB** tarball. Bioconductor's limit is 5 MB.
- **roxygen2:** `roxygenise()` reports "Skipping spatialCorrelationGeneExpIterPermutations.Rd … not generated by roxygen2". The committed file starts with a blank line, so future roxygen edits to that function will never reach the man page (`out/03_roxygenise.txt`). All other man pages match their sources (`out/03_man_diff.txt` is empty).

### 5.2 Examples

Each was run in its own process with `_R_CHECK_LIMIT_CORES_=TRUE` (`out/05`):

| Rd | Result | Time | Notes |
|---|---|---|---|
| spatialCorrelation | OK | **218.7 s** | Default settings on quakes (N = 998), run twice, the second time with 18 and 25 δ values. Prints p = 0.03 / 0.13 / 0.02. Should be under 5 s; with B = 19 and 3 δ values it takes 4.0 s, and on `speKidney` 1.2 s (`out/11`). |
| spatialCorrelationGeneExp | **ERROR** under a core limit | 4.1 s | `nThreads = 5` exceeds the 2-worker cap. BiocParallel stops: "workers must be <= 2 was (5)". |
| spatialCorrelationGeneExpWithinSample | **ERROR** (also on the live site) | 4.4 s | Passes a list where one SPE is expected, and `assayName = "A"`. |
| simRanPatternRasts | **ERROR** (also on the live site) | 0.1 s | `assays` not found. |
| spatialCorrelationGeneExpIterPermutations | Everything in `\dontrun` | 0.0 s | Prints undefined `negCorrelation`. |
| plotCorrelationGeneExp | OK | 23.9 s | Recomputes two correlations with B = 100. |
| linearRegression / pixelClass / savePlots | OK | about 4 s each | Unnamed input lists give blank axis labels ("Expression of pixel in "; `out/06` K). `savePlots` example has a stray `#'` (`R/visualizationFunctions.R:340`). Several call `MulticoreParam()`, i.e. all cores. |
| spatialSimilarity | OK | 4.0 s | Prints nothing, so the reference page shows no output. |
| matchingVariograms / viladomatCorrelation | OK | about 1 s | `matchingVariograms` prints whole permutation vectors (around 26,000 per value against a data range of 40–680). Permutations are not on the data's scale (also `out/06` Q), and that is never documented. |

### 5.3 Broken links: 8 on the live site, all from quoted URLs

Markdown such as `[text]("https://…")` turns into `href="%22https://…%22"`:

- `getting-started…Rmd:71, 176, 179` (×2)
- `acute…Rmd:219`
- `brain…Rmd:162, 184` (×2)
- `inst/scripts/biological-replicates-example.R:23` (a comment, but the same pattern)

### 5.4 Typos and grammar (file:line)

- **DESCRIPTION**
  - `:3` "Comparision"
- **README**
  - `:18` "Person correlation"
  - `:22` "must be first aligned … are coincide"
- **Getting Started**
  - `:47` "Person correlation"
  - `:123` "comparision"
  - `:378` "signficant"
  - `:443` "in which gene expression each pixels falls"
  - `:483` "neighbors’" (stray curly apostrophe)
- **AKI vignette**
  - `:240` "refered", "share unit system"
  - `:242` "a a"
  - `:244` "larger enough"
  - `:252` "the the"
  - `:278` "is to be able to"
  - `:361` "Filter lowly expression genes"
  - `:363` "zero expression of a gene expression"
  - `:379, 404` "parameters that controlling"
  - `:380` "ofautocorrelation"
  - `:514` "Similiarty"
  - `:624` "eache"
- **Brain vignette**
  - `:41` "This tutorial demonstrate"
  - `:215` "lowly expression"
  - `:437` "showed that the compare the correlation…"
  - `:441` "Similiarty"
  - `:460` "significantly positively SVGs"
- **R sources**
  - `R/spatialCorrelation.R:37`, `:737-738`, `:950-951` "permuation"
  - `:151` "matrix of with"
  - `:160, 381, 673, 892` and `R/iterativePermutations.R:126` "max distance in when"
  - `:651, 870` "should the same as"
  - `:836` and `R/iterativePermutations.R:342` "seperately"
  - `:1064` "specifiying"
- **Brain vignette metadata**
  - `VignetteIndexEntry` is the slug "brain-MERFISH-10x-visium" (`brain…Rmd:5`), not a title.
- **Install.Rmd**
  - Has an empty `author:`.
  - Copy-pasted `fig.path` (`:21`).

### 5.5 Packaging that affects users

- **Vignette dependencies aren't declared in `Suggests`:** MERINGUE (GitHub only), scatterbar, rhdf5, rjson, patchwork, gridExtra, Matrix, BiocGenerics and devtools. MERINGUE and scatterbar aren't even installed here, so the two case-study vignettes can't be rebuilt.
- **`class` is declared in `Suggests` but isn't used in `R/` or the vignettes.**
- **Vignettes depend on the network (Zenodo, 10x) and use `cache = TRUE`.**
- **`devtools::load_all()` runs inside vignettes and scripts.**
- **Scripts hard-code personal paths** (`~/ST_compare/...`, `~/github/STcompare/...`).
- **`biological-replicates-example.R:208-223`** shows the authors patching rows that came back `NA` with rows from another run (made with different settings). That is the result of errors being turned silently into `NA`. It is also statistically awkward, because BH-adjusted values from two different runs end up mixed.

---

## 6. API ergonomics

### 6.1 Inventory and naming

| Concept | Current names | Inconsistency |
|---|---|---|
| Main entry points | `spatialCorrelation(X, Y, pos)`, `spatialCorrelationGeneExp(input)`, `spatialCorrelationGeneExpIterPermutations(input)`, `spatialCorrelationGeneExpWithinSample(input)`, `spatialSimilarity(input)` | Three correlation functions with overlapping jobs. The iterative one is the one the authors actually use (all shipped real-data results). Names are long, and "GeneExp" carries no meaning. |
| Low-level engines | `viladomatCorrelation(data = N×4 matrix)`, `matchingVariograms(X.randomized, long, lat, delta, target_variog, prctile, ids, i, seed)` | Exported, but the arguments use dot and snake case. `prctile` is a *distance*, while `maxDistPrctile` is a *probability*. |
| Helpers exported by accident | `getGenePixelDF`, `threshold`, `assignFill` | `exportPattern` exports everything. `threshold` is a generic name likely to clash with other packages. |
| Plots | `plotCorrelationGeneExp(speList, spatialCorrelation, geneName, assayName)`, `linearRegression(input, gene, assayName)`, `pixelClass(input, gene, assayName)`, `savePlots(geneNames, spatialSimilarity, rastGexp, assayName, filePath = FALSE)` | Four names for the pair of objects (`input`, `speList`, `rastGexp`, plus `input` meaning a *result*). Arguments named after functions (`spatialCorrelation`, `spatialSimilarity`). `gene` vs `geneName` vs `geneNames`. |
| Assay | `assayName` (STcompare) vs `assay_name` (SEraster, used in the same tutorials) | The default is "first assay" and is not recorded in results. |
| Smoothing grid | `deltaX`/`deltaY` (list per gene) vs `delta` (WithinSample: list per gene; engines: numeric) | The same idea takes three different shapes. |
| Parallelism | `nThreads` and `BPPARAM` | Two knobs; one is silently ignored in one function (P1). `nThreads` doesn't mean threads. |
| Thresholds | `t1`, `t2`, `minQuantile`, `minPixels` (a proportion), `foldChange` (log2) | Names hide the units and scales. |
| Output columns | `pValuePermuteX/Y` (adjusted in one function, not in another), `deltaStarMedianX/Y`, list-columns `deltaStarX/Y`, `nullCorrelationsX/Y` (B×1 matrices inside list-columns), `percentSimilarity` (a proportion), `percentDissimilarityX/Y` | The same column name means different things depending on the function. List-columns print as `0.2, 0.2....`. |
| Defaults | `verbose = TRUE` for correlation, `FALSE` for similarity; `seed = 0` | Inconsistent defaults. |

### 6.2 Inputs

- **"A list of exactly two SpatialExperiments" is the only input format.**
  - It doesn't accept two objects (`x`, `y`) or a single SPE with two `sample_id`s, which is the idiomatic Bioconductor container.
  - It doesn't accept the 3-element list SEraster returns (the user has to subset it).
  - The order of `speKidney` is A, C, B (`out/06` A), which catches out anyone who selects by position.
- **No input checking:**
  - list length;
  - classes;
  - that the assay exists;
  - genes present in both;
  - **coordinates agreeing for shared pixel IDs**;
  - constant or all-zero genes;
  - the structure of `deltaX`;
  - the `adjustMethod` value (in `spatialCorrelationGeneExp`);
  - the minimum number of pixels.
- **What failures look like:** see the table in 6.4 below.

### 6.3 Outputs

- **Correlation result:** a `data.frame` with list-columns and `AsIs` classes. The gene name lives only in `rownames` (lost under dplyr/tibble), and no columns record which settings were used.
  - For one gene the row is named "cor" (`spatialCorrelation`).
  - `WithinSample` appends `first`/`second` columns, and its row names are arbitrary.
- **Similarity result:** a list containing a data.frame with six list-columns of pixel IDs, a long-format list-column of log ratios, and `parameters$input`, a copy of both input objects (memory, `out/06` I). The two results (correlation, similarity) can't be joined except by hand on gene name, and the vignettes join them by position.
- **No classes or methods:** no S3/S4 result class, and no `print`, `summary`, `plot`, `as.data.frame`, `head` or `topGenes`.
- **Nothing is written back** to `rowData()` or `metadata()`.

### 6.4 Progress, errors and robustness

- **Progress:** one `message()` per gene (`1: Gene`, or `nPermutations=100 | gene 3/325: X`). There is no progress bar, no ETA and no summary at the end. Similarity is quiet by default.
- **Errors:**
  - Engine errors are caught and `print(cond)`ed, and a row of `NA` is returned (`R/spatialCorrelation.R:585-619`). The user can't catch, suppress or log these as conditions, and there is no status or error column.
  - If `cor.test` itself fails, the handler uses `corDF` before it has been assigned.
  - In the iterative function, an `NA` row in round 1 becomes a gene called `NA` in round 2 and **crashes the whole run** (`out/07` R). A single all-zero gene can lose hours of work.
  - Nothing is saved along the way, and an interrupted run can't be resumed.
- **Error messages users actually see:**

  | Situation | Message |
  |---|---|
  | Gene missing from Y | "subscript out of bounds" |
  | A 1-gene SPE passed to the within-sample function | "n < m" |
  | A list passed to the within-sample function | "unable to find an inherited method … 'spatialCoords' for signature 'list'" |
  | Constant gene | A printed `bplist_error` blob plus "NA/NaN/Inf in foreign function call (arg 4)" |
  | `getGenePixelDF()` without `assayName` | "promise already under evaluation" |

- **RNG:** the global stream is reset (P3).
- **Side effects:**
  - `savePlots()` attaches three packages.
  - The examples use `MulticoreParam()`, i.e. all cores.

### 6.5 Tests, methods, namespace

- **Tests:** there is no `tests/` directory. Many of the problems above would have been caught by minimal testthat checks: BH applied across genes, p-values in (0, 1], the example code, an unchanged RNG, and rejection of mismatched grids.

---

## 7. Proposal: a "comfortable" user-facing design

### 7.1 Principles

- One obvious entry point.
- Validate early and explain clearly.
- Return tidy, self-describing results, with one row per gene, no list-columns in the main table, and the settings recorded.
- Make the defaults reflect what the authors actually do (iterative B and an extended δ grid).
- Make it fast to try: an estimate first, a progress bar, and the option to resume.
- Be reproducible without side effects.
- Keep backwards compatibility through thin wrappers.

### 7.2 Front door (sketch)

```r
res <- compareSpatial(
  x, y,                                   # two SpatialExperiments on ONE shared grid; also accepts list(x, y)
                                          # or a single SPE plus sample_ids = c("ctrl", "aki")
  assay         = "pixelval",             # recorded in the result; accepts assay_name= as an alias
  genes         = NULL,                   # default: shared genes passing minDetection
  minDetection  = 0.05,                   # fraction of shared pixels with expression > 0 in BOTH samples
  tests         = c("correlation", "similarity"),
  nPermutations = c(100, 1000),           # iterative by default; a scalar means one round
  rerunAlpha    = 0.05,                   # documented screening rule (on combined raw p)
  delta         = c(0.01, 0.05, seq(0.1, 0.9, 0.1)),   # what the authors used for AKI and MERFISH
  maxDistQuantile = 0.25,
  log2FC        = 1,                      # similarity band: |log2(y/x)| <= 1, i.e. within 2-fold
  BPPARAM       = BiocParallel::SerialParam(),   # one knob; documented; Windows-safe
  seed          = 1L,                     # local streams; the global RNG is left alone
  progress      = interactive(),
  checkpoint    = NULL                    # optional directory: write per-gene results and resume
)
```

Supporting functions:

- **`checkComparable(x, y)`** is a pre-flight report. It covers shared pixels, coordinate agreement for shared IDs (an error if they differ), gene overlap, constant or low-detection genes, the assays present, the pixel-count range, and an estimated runtime.
- **`estimateRuntime(x, y, ...)`** times 1–2 permutations for a few genes and extrapolates for the chosen B, iterative schedule and BPPARAM. `compareSpatial()` prints this estimate before it starts.
- **`prepareComparison(speList, resolution, fun, normalize = c("cpm", "none"))`** is an optional wrapper. It runs SEraster on the whole list in one call (so the grid is guaranteed to be shared), adds CPM, and runs `checkComparable()`.

### 7.3 Result object

`STcompareResults`: either an S4 class extending `S4Vectors::DataFrame` or a classed tibble. The main table has one row per gene and no list-columns:

| Column | Meaning |
|---|---|
| `gene` | Gene (row) name, stored as a real column |
| `nPixels` | Shared pixels used |
| `correlation` | Pearson r |
| `pValueNaive` | `cor.test` p-value (assumes independent pixels; for reference only) |
| `pValuePermX`, `pValuePermY` | Empirical p from permuting X and from permuting Y, computed as `(b+1)/(B+1)` |
| `pValue` | `pmax(pValuePermX, pValuePermY)`: the conservative p to report |
| `padj` | BH (or the chosen method) across genes, applied to `pValue` |
| `nPermutations` | B actually used for this gene |
| `deltaStarX`, `deltaStarY` | Median δ* |
| `deltaEdgeX`, `deltaEdgeY` | Fraction of permutations whose δ* was at the edge of the grid (warn when > 0.5) |
| `similarity` | Fraction of thresholded pixels within ±log2FC |
| `fractionHigherX`, `fractionHigherY` | Fraction with log2(y/x) < −log2FC (X higher) and > log2FC (Y higher) |
| `nPixelsSimilarity`, `thresholdX`, `thresholdY` | Similarity bookkeeping |
| `status`, `message` | `"ok"`, `"skipped: constant"`, `"error: …"` (never crash, never fail silently) |

Diagnostics such as null correlations, δ* per permutation, variogram-fit residuals and pixel classes go in `metadata(res)$details` and are retrieved with `details(res, gene)`. The settings (`call`, package version, engine, seed, assay, grid, B schedule) go in `metadata(res)$settings`.

Methods:

| Method | What it does |
|---|---|
| `print()` / `show()` | A short header, then counts by direction at `padj < 0.05`, then the top genes |
| `summary()` | |
| `plot(res)` | Correlation vs similarity, coloured by significance; or a volcano-style plot |
| `plotGene(res, gene, type = c("raster", "scatter", "classes", "null"))` | Replaces `plotCorrelationGeneExp`, `linearRegression`, `pixelClass` and `savePlots`; uses the stored assay and names |
| `topGenes()` | |
| `as.data.frame()` | |
| `addToRowData(spe, res, prefix = "STcompare.")` | |

### 7.4 Integration with SpatialExperiment

- Accept a single SPE with a `sample_id` column (the container SEraster/SpatialExperiment users already have after `cbind`).
- Optionally write per-gene results into `rowData()` of either input, with a prefix such as `vs_<other sample>.`.
- Record the assay name, so plots never fall back to "first assay" by mistake.

### 7.5 Validation and error handling

- Error when the coordinates of shared pixel IDs differ: "x and y were not rasterized on the same grid; rasterize them together: `SEraster::rasterizeGeneExpression(list(x, y), ...)`".
- Skip constant genes with a status, instead of `NA` plus a printed blob.
- Turn per-gene errors into `status`/`message` and a single summary `warning()` ("3 genes failed; see `res$status`"). Never `print()`.
- Check `adjustMethod`, `delta`, `nPermutations` and the assay before any computation starts.

### 7.6 Reproducibility without clobbering the RNG

- Wrap all randomness in `withr::with_preserve_seed()`, or use `dqrng`/L'Ecuyer streams.
- Derive one independent stream per (gene, direction) from `seed`, e.g. with `parallel::nextRNGStream` or a dqrng stream id equal to the gene index × 2 + direction. Results then don't depend on the number of workers, the chunking or the gene order, and `.Random.seed` is left untouched.
- Document that permutation p-values are Monte Carlo estimates, and add a `pValueMCSE` column if useful.

### 7.7 Parallelism

- Parallelize **over genes** in the outer loop. This parallelizes better than over permutations, and is what the AKI vignette comment already assumes. Alternatively, flatten (gene × permutation-chunk) tasks.
- Expose a single `BPPARAM`, with `SerialParam()` or `bpparam()` as the default. Keep `nThreads` as a deprecated alias.
- Explain that on Windows you need `SnowParam` (or the C++ threads below).
- With a C++ backend, use RcppParallel threads inside each gene's permutation loop, and say clearly which layer does what.

### 7.8 Progress

- Use `cli::cli_progress_bar` (or `progressr`, so it works across BiocParallel workers) with an ETA.
- Print a one-line summary at the end: genes tested, skipped and failed; elapsed time.
- Offer `verbose = 0/1/2`.

### 7.9 Backwards compatibility

- Keep the existing exported functions as thin wrappers over the new engine, with `lifecycle::deprecate_soft()` notes. They keep their current column names, but the documented bugs get fixed (BH across genes, BPPARAM pass-through, the threshold matching its docs) and are listed under **"Behaviour changes"** in a new `NEWS.md`.
- Argument aliases with deprecation messages: `assayName` → `assay`, `nThreads` → `BPPARAM`, `foldChange` → `log2FC`, `deltaX`/`deltaY` → `delta` (still accepting per-gene lists).
- `as.data.frame(res, legacy = TRUE)` returns the old column layout for downstream scripts.
- Don't export internal helpers any more (`@keywords internal`, a roxygen-managed NAMESPACE). Keep deprecated exports for one release.
- Bump the version (0.2.0), add an `inst/CITATION`, and add snapshot tests comparing old and new wrappers. The tests should check statistical equivalence: identical r and similarity, and permutation p-values within Monte Carlo error once RNG streams change.

### 7.10 Interaction with a C++ backend

- Keep the R-level API and documentation independent of the engine.
- Provide `engine = c("cpp", "R")` for validation.
- Document numerical equivalence: a locfit-equivalent kernel smoother and the same variogram binning.
- State that RNG streams differ from versions before 0.2.0. Offer `legacy_rng = TRUE` to reproduce old numbers exactly through the R engine.
- Adaptive stopping could be added: Besag–Clifford sequential p-values, stopping after h exceedances. This does what the iterative scheme tries to do, but more cheaply and with a formal justification. It is easy to explain in the docs.

---

## 8. Documentation plan

### 8.1 New article: "How STcompare works"

The article should take the reader through the algorithm in order, with one figure per step:

1. **Inputs.** Two samples on a shared pixel grid; X is gene g in sample 1 and Y is gene g in sample 2, at the shared pixels.
2. **Observed statistic.** r = cor(X, Y). Explain why `cor.test` is anti-conservative under spatial autocorrelation, using the shipped null simulation: 49.8% of raw naive p-values are below 0.05 on independent fields (`out/01`).
3. **Null for "permute X":** for b = 1…B:
   - (a) shuffle X across pixels, which destroys autocorrelation;
   - (b) smooth the shuffled values with a Gaussian kernel whose adaptive bandwidth covers the nearest δ·N pixels (locfit, local-constant);
   - (c) compute the empirical variogram, with at most 1,000 sampled pixels and bins up to the 25th percentile of pairwise distances;
   - (d) regress X's variogram on it, giving scale β₁ and nugget β₀;
   - (e) form X̂ = √|β₁|·X_δ + √|β₀|·Z (the location is not preserved, which is why returned permutations look shifted);
   - (f) keep the δ* that minimizes the residual sum of squares between the variograms;
   - (g) compute r_b = cor(X̂_δ*, Y).
4. **Empirical p** = (#{|r_b| ≥ |r|} + 1)/(B + 1). Do the same with Y permuted, report max(pX, pY), and adjust across genes.
5. **Iterative refinement.** Give the exact rule and its cost.
6. **Similarity.** Thresholds (OR rule), pseudocount, log2 ratio, the three classes, S. Stress that the values must be on a linear scale.
7. **Diagnostics.** δ* histogram (edge hits), null distribution against r, variogram fit, scatter of X vs Y with classes.

### 8.2 Parameter guide

| Parameter | Default | Controls | How to choose | Diagnostic |
|---|---|---|---|---|
| SEraster `resolution` | — | Pixel size, hence N and runtime | Larger than the alignment error and smaller than the structures of interest. Roughly 300–3,000 shared pixels is a practical range (runtime table §4). Rasterize **jointly**. | Shared pixel count; cells per pixel |
| `fun` + normalization | mean | Pixel value | `sum` then CPM (AKI), or `mean` of library-normalized values. Correlation tolerates log; similarity needs a linear scale. | Total-expression raster |
| Gene filter | none | Wasted runtime; NA results | Detected in at least 5% (AKI) or 1% (brain) of pixels in **both** samples. Constant genes can't be tested. | `status` column |
| `nPermutations` | 100 / `c(100, 1000)` | Resolution of p (min p = 1/(B+1)) | With G genes and BH at α, genes can only become significant if p is well below α; B = 100 limits you to p ≥ 0.0099. Use the iterative default, or more rounds (`c(100, 1000, 10000)`) for large G. | Fraction of genes at the floor |
| Screening α | 0.05 | Which genes are rerun | State the rule; report how many genes were rerun and the cost | Count of reruns |
| δ grid | `seq(0.1, 0.9, 0.1)` | Candidate smoothing | Include 0.01 and 0.05 for fine or high-resolution patterns (the authors did for AKI and MERFISH) | `deltaEdgeX/Y` (58–85% at the low edge in the shipped brain results) |
| `maxDistPrctile` | 0.25 | Variogram range used for matching | 0.25 is typical; lower for large, heterogeneous sections | Variogram plot |
| `seed` | 0 | Reproducibility | Any value; record it | — |
| `t1`/`t2`/`minQuantile` | NULL / 0.05 | Which pixels count for similarity | Use absolute thresholds on the assay's scale if known. The default removes pixels that are low in *both* samples. | `numPixelOutThresh` |
| `minPixels` | 0.1 | Minimum fraction of pixels kept | — | NA similarity |
| `foldChange` | 1 | Similarity band, **log2** | 1 means within 2-fold; 0.585 means within 1.5-fold | Class fractions |

### 8.3 Performance and parallelization guide

- **Cost model.** Per gene: 2 directions × B × |δ| × (one locfit fit on N pixels plus two variograms on min(N, 1000) pixels). Iterative schedules multiply this by 1 + (rerun fraction × B₂/B₁). That was 5.5–9.2× in the shipped analyses.
- **Measured numbers** from §4, plus the authors' own runs: 1.8 h, 7 min, 16.8 h and 7.65 h.
- **How parallelism works today.** Forked processes over permutations, genes one at a time, the Windows fallback, and why `nThreads` = all cores can oversubscribe a shared node.
- **Memory.** `returnPermutations = TRUE` keeps N × B doubles per gene and direction, e.g. 40 MB for N = 5,000 and B = 1,000.
- **Recipe:**
  1. filter genes;
  2. dry-run 10 genes with B = 100;
  3. estimate;
  4. run with the iterative schedule;
  5. on HPC, split genes with `BatchtoolsParam`, or over array jobs with `genes =` subsets, and combine.

### 8.4 FAQ (suggested entries)

1. *Why do `pValuePermuteX` and `pValuePermuteY` differ, and which do I report?*
   Each direction preserves the autocorrelation of only one sample. Report the max (conservative) and adjust across genes.
2. *Why is my p-value 0?*
   This is the Monte Carlo floor. Report < 1/(B+1), or use the (b+1)/(B+1) estimate.
3. *Why are some results NA?*
   Constant genes, too few pixels, or an engine error. Check the status column.
4. *My two samples were rasterized separately. Is that OK?*
   No. Rasterize them together, then check with `checkComparable()`.
5. *How long will it take?*
   See the table, or run `estimateRuntime()`.
6. *Most δ* values are at the edge of the grid. What does it mean?*
   The best smoothing lies outside the grid. Extend the grid.
7. *Can I compare more than two samples?*
   Yes, as pairwise calls. Mention multiple testing across pairs.
8. *Can I compare cell-type proportions?*
   Yes (see the brain vignette): put the proportions in an assay.
9. *Can I compare genes within one sample?*
   `spatialCorrelationGeneExpWithinSample()`.
10. *Do I need STalign?*
    Only for sections that aren't already aligned. STalign is Python.
11. *Which assay should I use for similarity?*
    A linear-scale normalized one (CPM), not log.
12. *Is the similarity score a test?*
    No. It is descriptive; there is no p-value.
13. *Square or hexagonal pixels?*
    Hexagons fit Visium spots better; either works.
14. *Does it work on Windows?*
    Yes, but parallelism needs `SnowParam`.

### 8.5 `_pkgdown.yml`: reference and article grouping

```yaml
url: https://jef.works/STcompare/
template:
  bootstrap: 5
repo:
  url:
    home: https://github.com/JEFworks-Lab/STcompare/
navbar:
  structure:
    left:  [articles, reference, news]
    right: [search, github]
articles:
- title: Get started
  navbar: ~
  contents: [Install, getting-started-with-STcompare]
- title: Concepts and guides
  navbar: Guides
  contents: [how-it-works, from-raw-data, choosing-parameters, performance, faq]
- title: Case studies
  navbar: Case studies
  contents: [acute-kidney-injury-10x-visium-rasterized, brain-MERFISH-10x-visium]
reference:
- title: Compare two samples
  desc: Main entry points (pattern correlation with empirical p-values; magnitude similarity).
  contents: [spatialCorrelationGeneExpIterPermutations, spatialCorrelationGeneExp, spatialSimilarity]
- title: Visualize results
  contents: [plotCorrelationGeneExp, linearRegression, pixelClass, savePlots]
- title: Other comparisons and low-level engines
  contents: [spatialCorrelationGeneExpWithinSample, spatialCorrelation, viladomatCorrelation, matchingVariograms]
- title: Example data
  contents: [speKidney, simRanPatternRasts]
```

The helpers (`getGenePixelDF`, `threshold`, `assignFill`) should get `@keywords internal` so pkgdown hides them. Give every topic a descriptive `@title`, for example "Spatial correlation test between two rasterized samples". Add `@family` and `@seealso`, enable roxygen markdown, and rewrite the `\value` lists with `\describe{}` or markdown lists to remove the "Lost braces".

### 8.6 README rewrite (outline)

1. One-sentence purpose, plus the figure.
2. **Which test answers which question** (two-row table).
3. **What you need:** two sections, aligned (STalign) and rasterized **together** (SEraster), with matching pixel IDs.
4. **Install:**
   ```r
   install.packages("BiocManager")
   BiocManager::install("JEFworks-Lab/STcompare")
   ```
   This resolves the Bioconductor dependencies.
5. **Quick start.** Runs in about 10 s; 9.1 s measured for the correlation on one core:
   ```r
   library(STcompare)
   data(speKidney)                                   # simulated sections A, B, C (1 gene)
   rast <- SEraster::rasterizeGeneExpression(speKidney, assay_name = "counts",
             resolution = 0.2, fun = "mean", square = FALSE)   # rasterize jointly
   cor_AB <- spatialCorrelationGeneExp(list(A = rast$A, B = rast$B))
   cor_AB[, c("correlationCoef", "pValuePermuteX", "pValuePermuteY")]
   sim_AB <- spatialSimilarity(list(A = rast$A, B = rast$B))
   sim_AB$similarityTable[, c("gene", "percentSimilarity")]
   ```
6. **Learn more:** tutorials, How it works, Choosing parameters, Performance, FAQ.
7. **Citation** (Bioinformatics 2026, from `inst/CITATION`) and **Getting help** (issues link).

### 8.7 Vignette hygiene

- Remove `devtools::load_all()` and `cache = TRUE`.
- Declare `Suggests`, and add `Remotes:` for MERINGUE, or precompute the Moran's I results the way the correlation results already are.
- Replace hard-coded numbers with inline `` `r ...` ``.
- Join tables on gene name.
- Fix the AKI `rast$ctrl` / `rast$AKI` names, and the stale brain numbers and truncated sentence.
- Make the Getting Started null-calibration figure compare like with like.
- Add a short "What you will learn / time to run" box to each case study.
- Move the large `.RData` files (32 MB) to Zenodo or ExperimentHub, and download them on demand with caching (BiocFileCache).

---

## 9. Prioritized action list (impact vs effort)

Effort: S is under 2 h, M is half a day to 2 days, L is more than 2 days. "Code" means a small code change is needed to make the docs true.

| # | Action | Impact | Effort | Type | Where |
|---|---|---|---|---|---|
| 1 | Apply `p.adjust` across genes in `spatialCorrelationGeneExp`, check `adjustMethod` up front, and add a NEWS note on the behaviour change | **Very high** (changes published conclusions) | S | code + docs | `R/spatialCorrelation.R:836-842, 715-718` |
| 2 | Fix the p-value floor: adopt (b+1)/(B+1) with `>=`, or at least correct "smallest p-value is 0.01" in 4 places | Very high | S | code/docs | `:326-327, 169, 358, 647, 866` |
| 3 | Make the iterative screening rule and its docs agree; add an `nPermutations` column | High | S | code/docs | `R/iterativePermutations.R:64, 90-94` |
| 4 | Pass `BPPARAM` through; rewrite the `nThreads` docs (processes, over permutations, Windows) | High | S | code/docs | `R/spatialCorrelation.R:830`; 5 docs blocks |
| 5 | Fix the Getting Started null-calibration figure and text (labels, max vs X only, BH) and the "bottom 5%" sentence | High | S | docs | `getting-started…Rmd:317-397, 418` |
| 6 | Fix the broken examples (WithinSample, simRanPatternRasts, the iterative example), keep each under 5 s (B = 19, 3 δ values), and remove `nThreads = 5` and `MulticoreParam()` from examples | High | S | docs | Rd examples |
| 7 | Fix the 8 quoted links; correct the README's SEraster → STalign; fix "Comparision"; add `biocViews`, `BugReports` and the GitHub URL; install through BiocManager; add `inst/CITATION` | High (first impression) | S | docs/packaging | README, DESCRIPTION, vignettes |
| 8 | Clear the `R CMD check` WARNINGs and NOTEs: WithinSample/verbose docs, declare `sf`, remove `library()` calls, ASCII only, `importFrom(stats, …)`, roxygen-managed NAMESPACE, fix `\value` braces, enable markdown, restore the IterPermutations Rd header | Medium-high | S–M | packaging | R/, man/, DESCRIPTION |
| 9 | Correct the similarity docs (log2 `foldChange`, proportions, "either" not "both", `minPixels`, 1e-4 pseudocount) and fix the `numPixelInThresh = 1` bug and the `savePlots` panel-2 bug | Medium-high | S | code + docs | `R/packageFunction.R`, `R/visualizationFunctions.R:395` |
| 10 | Fix the vignette errors: AKI `rast$ctrl`, stale brain numbers, truncated sentence, joins by position, wrong A/C label | Medium | S | docs | vignettes |
| 11 | Input validation with clear messages (shared-grid coordinates, genes, assay, constant genes, `delta` shape) | **High** | M | code | entry points |
| 12 | Error handling: no `print()`; status/message columns; one summary warning; fix the iterative round-2 crash on `NA` | High | M | code | `R/spatialCorrelation.R:585-619`, `R/iterativePermutations.R:63-67` |
| 13 | Combined `pValue` = max(pX, pY) plus `padj`; document "which p to report" in `@return` | High | S–M | code + docs | outputs |
| 14 | Leave the global RNG alone; use streams per gene and direction; document | Medium-high | M | code + docs | `R/spatialCorrelation.R:99, 242` |
| 15 | Progress bar with ETA, `estimateRuntime()`, a summary at the end | Medium-high | M | code | — |
| 16 | pkgdown reference and article grouping, descriptive titles, internal helpers hidden, GitHub link in the navbar | Medium | S | docs | `_pkgdown.yml`, roxygen |
| 17 | New articles: How it works, From raw data, Choosing parameters, Performance, FAQ; README quick start | **High** | M–L | docs | vignettes/ |
| 18 | Default δ grid including 0.01 and 0.05, plus an edge-of-grid warning and `deltaEdge` columns | Medium | S–M | code + docs | defaults |
| 19 | Front door `compareSpatial()` with a result class, methods, rowData helpers and `checkComparable()`; legacy wrappers with lifecycle notes | High (comfort) | L | code + docs | new |
| 20 | testthat suite: examples as tests, BH across genes, p in (0, 1], RNG unchanged, mismatched grids rejected, similarity arithmetic, calibration on a small subset of `simRanPatternRasts` | High (prevents regressions) | M | tests | `tests/` |
| 21 | Slim the package: tarball 48.6 MB → under 5 MB (move `extdata` out); make the vignettes buildable (Suggests/Remotes, precomputed objects) | Medium | M | packaging | inst/, vignettes |
| 22 | C++ backend (separate workstream). Documentation needs: a numerical-equivalence note, an `engine` option, RNG change notes | High (speed) | L | code + docs | — |

**Suggested order.** Items 1–10 together are about two days of work and remove every known inaccuracy and broken example. Next, 11–18 make the package safe and much easier to understand. Then build 19–20 on that base, alongside the C++ engine (22).

---

## Appendix A: Key evidence excerpts

```
# out/06 C - adjustMethod is a no-op (6 genes, nPermutations = 20)
        r pX_BH pX_none pX_bonf manual_BH manual_bonf
g5 -0.056  0.80    0.80    0.80      0.95           1
identical(BH output, 'none' output): TRUE
identical(bonferroni output, 'none' output): TRUE

# out/06 D - BPPARAM silently ignored
BPPARAM = <character string> -> no error (BPPARAM silently ignored)

# out/06 E - global RNG reset
runif(1) after call, user seed 1: 0.462357   user seed 2: 0.462357   identical: TRUE

# out/06 H - iterative threshold
genes carried to round 2 (alpha=0.05, nPermutations[1]=100): p~0.03, p~0.0004
documented threshold alpha/nPermutations[k] =  5e-04 ; code threshold =  0.05

# out/07 I2 - numPixelInThresh in the minPixels-failure branch
  gene percentSimilarity numPixelInThresh numPixelOutThresh
1 Gene                NA                1               257
actual pixels passing (x > t1 | y > t2): 16

# out/07 R - iterative run with one all-zero gene
nPermutations=20 | gene 2/3: NA
ERROR: subscript out of bounds

# out/09 - grids not shared, nothing detected
separately rasterized: shared pixel IDs = 137 ; max |coordinate difference| for same ID = 4.157
spatialCorrelationGeneExp on mismatched grids -> no error/warning; r = 0.238 computed on 137 'shared' pixels

# out/01 - shipped results
kidneyCorrelationNoIter: pValuePermuteX zero = 573 of 1046; fraction on 1/100 grid: 1 (=> unadjusted)
simRanPatternResults: corspv_corrected on 1/1000 grid: 1; share < 0.05: 0.0388; share BH(empirical) < 0.05: 0

# out/12 - delta* at the edge of the grid
brainCorrelation dir X: permutations at grid min (0.10) = 58.3%; genes with >50% at min = 43.1%
ctCorrelation    dir X: permutations at grid min (0.10) = 84.9%
```

## Appendix B: Runtime numbers

- **`out/08`:** `speKidney` A vs B (273 pixels), defaults: 9.1 s per gene on 1 core, 3.0 s on 4 workers. Synthetic data extrapolated from B = 10 to B = 100 on 1 core: N = 273 → 12 s, 1,000 → 44 s, 2,500 → 73 s, 5,000 → 86 s per gene.
- **`out/11`:** fast example settings take 1.2 s on `speKidney` (B = 19, 3 δ values) and 4.0 s on quakes.
- **`out/05`:** the current `spatialCorrelation` example takes 218.7 s.
- **Authors' notes:** brain, 325 genes, 22 workers: 1.78 h. Cell types, 16 rows: 6.95 min. MERFISH replicates, 483 genes, 20 workers: 16.8 h. The same with affine alignment: 7.65 h.

## Appendix C: Files produced

Work directory: `/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/docs/`

| Script | Output | Purpose |
|---|---|---|
| `scripts/01_inspect_extdata.R` | `out/01_inspect_extdata.txt` | Shipped results: p-value grid, zeros, B per gene, δ* |
| (`R CMD build`/`check`) | `out/02_rcmdcheck_static.txt`, `pkgcopy/STcompare.Rcheck/00check.log`, `pkgcopy/STcompare_0.1.0.tar.gz` | Static check, 48.6 MB tarball |
| (roxygen2 rerun) | `out/03_roxygenise.txt`, `out/03_man_diff.txt` (empty) | Rd drift; unmanaged Rd and NAMESPACE |
| `scripts/04_extract_examples.R` | `examples/*.R` | Rd examples as scripts |
| `scripts/05_run_one_example.R`, `scripts/05_run_all_examples.sh` | `out/05_examples.txt` | Example status and timing under a core limit |
| `scripts/06_verify_behaviors.R` | `out/06_verify_behaviors.txt` | BH no-op, BPPARAM, RNG, δ misuse, threshold rule, similarity docs, labels, helper errors, constant gene, scale of permutations |
| `scripts/07_verify_more.R` | `out/07_verify_more.txt` | `numPixelInThresh` bug, the 5% claim, iterative crash, late `adjustMethod` failure, print output |
| `scripts/08_runtime_estimates.R` | `out/08_runtime_estimates.txt` | Per-gene runtime |
| `scripts/09_separate_rasterization.R`, `scripts/09b_separate_rasterization_noshift.R` | `out/09_*.txt`, `out/09b_*.txt` | Mismatched-grid hazard |
| `scripts/10_aki_vignette_names.R` | `out/10_aki_vignette_names.txt` | AKI `rast$ctrl` bug |
| `scripts/11_fast_example_timing.R` | `out/11_fast_example_timing.txt` | Fast example settings |
| `scripts/12_delta_boundary.R` | `out/12_delta_boundary.txt` | δ* at grid edges; cost of the iterative scheme |
| `scripts/13_similarity_AC.R` | `out/13_similarity_AC.txt` | Checks the Getting Started similarity claims (they hold: S = 0.535 for A–B; A–C has 100% of pixels higher in C) |
