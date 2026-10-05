# Spec: `compareSpatial()`, the new main function

This is phase 3 of `cpp-backend-plan.md`. The legacy functions keep their legacy semantics (`engine-spec.md`). This function is where the statistical improvements go, per the maintainer's decisions in plan §6–7: independent random streams, adaptive Besag–Clifford p-values, one combined p-value with BH, tidy output and a progress display.

## 1. Signature

```r
compareSpatial(x, y = NULL, assay = 1, genes = NULL,
               tests = c("correlation", "similarity"),
               nPermutations = 10000, exceedances = 10,
               delta = c(0.01, 0.05, seq(0.1, 0.9, 0.1)), maxDistPrctile = 0.25,
               foldChange = 1, minQuantile = 0.05, minPixels = 0.1,
               adjustMethod = "BH", seed = 0L,
               nThreads = getOption("STcompare.nThreads", 1L),
               progress = interactive(), verbose = TRUE, keepNulls = FALSE)
```

- **`x`, `y`:** two `SpatialExperiment` objects with matched pixels, for example from one `SEraster::rasterizeGeneExpression(list(x, y), ...)` call. Alternatively `x` is a list of two such objects and `y` is `NULL`; list names, or "x" and "y", become the sample labels.
- **`assay`:** a name or index, used for both objects, or a length-2 vector.
- **`genes`:** defaults to the genes present in both objects, with a message when some are dropped. Explicit genes that are missing from either object are an error.
- **`tests`:** which tests to run. `"similarity"` alone runs no permutations.

## 2. Input validation (errors with actionable messages)

- **Shared pixels.** They are `intersect(colnames(x), colnames(y))`; at least 3 are needed. For every shared pixel, the two `spatialCoords()` must agree to within 1e-6 of the coordinate range. Otherwise stop, with a message to rasterize both samples in one SEraster call, because pixel names alone are not enough. This fixes audit finding B10 for this function.
- **Gene names.** Unique row names are required.
- **Per-gene problems become a `"failed"` row** with a `message` and NA statistics, and one summary warning lists the affected genes. They do not stop the run. Examples: non-finite values on the shared pixels, a constant gene, a tree overflow.
- **Delta grid.**
  - Keep only deltas with `floor(N * delta + 1e-12) >= 2` for the number of shared pixels N.
  - Report dropped deltas in one message.
  - Error if none remain.
  - The default grid is the extended one the authors used for the kidney and MERFISH analyses. The default grid of the legacy functions often left delta* at its smallest value (investigation/05, /06).

## 3. Correlation test

- **Observed statistic:** r = Pearson correlation on the shared pixels. `pNaive` comes from `cor.test()`.
- **Null:** the Viladomat surrogates of `engine-spec.md` §2: the locfit-tree smoother, geoR variogram matching, closed-form OLS and delta search. Both directions run, X permuted and Y permuted.
- **Random streams.** Each (seed, gene name, direction, permutation b) has its own independent stream: permutation indices and Gaussian noise are generated in C++ by a counter-based generator keyed by those four values. Results therefore do not depend on `nThreads`, the batch schedule, the gene order, or which other genes are in the call.
  - The variogram subsample for N > 1000 is part of the per-dataset plan. It is drawn once from `seed`, exactly as the legacy plan draws it, and the user's global RNG state is restored.
- **Adaptive stopping (Besag–Clifford, plan §7).** Both directions of a gene run in lockstep. The gene stops at the first permutation L at which either direction has `exceedances` (h) nulls with |null| ≥ |r|, or at `nPermutations` (n_max).
  - **Combined p-value:** `p = h / L` if the gene stopped by exceedances, otherwise `(max(bX, bY) + 1) / (L + 1)`. This equals the maximum of the two directions' standalone sequential p-values, so it is valid.
  - **Per-direction values:** `pX = (bX + 1)/(L + 1)` and `pY` likewise, at the common L. They are descriptive; the help page says so.
  - **Fixed B:** `exceedances = Inf` gives fixed B = `nPermutations` for every gene, with p = max(pX, pY).
- **Adjustment:** `padj = p.adjust(p, adjustMethod)` across genes. NA rows are excluded from the count.
- **Diagnostics:**
  - `deltaStarMedianX`/`Y`;
  - `deltaGridEdge`: TRUE when, in either direction, more than half of the permutations chose the smallest or the largest delta of the grid. A message summarises how many genes are flagged, with advice to extend the grid.
  - A message when genes reached n_max with b = 0 and `1 / (n_max + 1)` could limit BH significance at the observed number of genes.

## 4. Similarity test

- **Semantics:** the same as `spatialSimilarity()`: per-gene thresholds at the `minQuantile` quantile, pixels kept if above the threshold in either sample, zeros replaced by 1e-4, similarity = the share of kept pixels with |log2(y/x)| ≤ `foldChange`, and NA when fewer than a `minPixels` share of the pixels pass.
- **Implementation:** vectorised over genes (an internal helper shared with the reworked `spatialSimilarity()` in phase 5). It must reproduce `spatialSimilarity()`'s numbers.
- **Columns:** `similarity`, `dissimilarityX`, `dissimilarityY`, `nPixelsSimilarity`, `thresholdX`, `thresholdY`. Pixel ID lists are not included; `spatialSimilarity()` still provides them.

## 5. Output

- **Shape:** a `data.frame` with class `c("STcompareResult", "data.frame")` and one row per gene. Row names are the genes, and there are no list-columns.
- **Columns, in order:**
  - `gene`, `nPixels`
  - correlation: `r`, `pNaive`, `p`, `padj`, `pX`, `pY`, `nPermutations` (L), `stop` (`"exceedances"`, `"limit"` or `"failed"`), `deltaStarMedianX`, `deltaStarMedianY`, `deltaGridEdge`
  - similarity: `similarity`, `dissimilarityX`, `dissimilarityY`, `nPixelsSimilarity`, `thresholdX`, `thresholdY`
  - `status` (`"ok"` or `"failed"`), `message`

  Columns of tests that were not run are absent.
- **Attributes:**
  - `params`: every argument value after defaults, the delta grid actually used, the sample labels and the assay;
  - `call` and `runtime`;
  - with `keepNulls = TRUE`, `details`: a per-gene list of `nullX`, `nullY`, `deltaStarX` and `deltaStarY`.
- **Methods:**
  - `print()`: a header with genes, pixels, tests and permutation settings, the counts with `padj < 0.05` split by the sign of r, and how many genes stopped early, followed by the first rows;
  - `summary()`: a small summary object;
  - `as.data.frame()`: drops the class.

## 6. Progress display and messages

- **Progress line.** With `progress = TRUE` (default `interactive()`), a single line updated at most about twice a second, for example:
  `compareSpatial: 61% | genes done 640/1046 | 1.9M permutations | elapsed 0:12 | ETA 0:08`.
  - It is driven by the engine's batch loop, plus a main-thread callback inside long batches, so that it advances smoothly in fixed and adaptive mode.
  - It ends with a final newline. No R API is called from worker threads.
  - When `progress` is off, there is no progress output at all.
- **Messages.** `verbose` controls the informative messages: dropped genes, dropped deltas, flagged genes and the final timing.

## 7. Threads and reproducibility

- `nThreads` sets the number of C++ threads. The default is `getOption("STcompare.nThreads", 1L)`, so users can set it once.
- Results depend only on the data, the parameters and `seed`, and are deterministic on a given platform. The normal generator may differ in the last bits across platforms.

## 8. Tests (lean)

- **Invariance:**
  - `identical()` results for 1 and 2 threads, for different batch schedules, for reversed gene order and for gene subsets;
  - the global RNG state is unchanged.
- **Adaptive stopping:**
  - it equals an offline Besag–Clifford computation over a fixed-length sequence from the same streams;
  - `exceedances = Inf` gives fixed B.
- **Calibration** (slow tier, or default if fast): on independent `simRanPatternRasts` pairs, the share of p < 0.05 is not above a one-sided binomial bound.
- **Controls:** speKidney A vs B (negative) and A vs C (positive) are significant with modest `nPermutations`.
- **Errors and NA rows:**
  - a coordinate mismatch is an error;
  - the dropped-delta message;
  - a failed gene gives an NA row and a warning.
- **Similarity:** the columns equal `spatialSimilarity()` on the fixtures.
- **Methods:** `print()` and `summary()` work, and `progress = TRUE` does not break the result.

## 9. Amendments after the acceptance review (2026-10-05)

- **Rarely detected genes.** The surrogates' values are close to normally distributed, so the correlation of
  two genes detected in a few pixels has a null distribution with much heavier tails than the surrogates
  give, and its p-value is far too small (on the whole AKI raster, genes detected in one pixel per section were
  significant). New argument `minDetected` (after `maxDistPrctile`): the share of the shared pixels at which a
  gene must be detected in each sample, where detected means above the gene's lowest value. Default `NULL`:
  `ceiling(sqrt(N))` pixels, which bounds the heuristic departure `N / (kx * ky)` by 1. Genes below it, and
  constant genes, get `status = "skipped"` and `stop = "skipped"`, `NA` p-values, their `r`, `pNaive` and
  similarity, and one message (no warning); `padj` counts only the tested genes. The help page and the
  vignettes say that sparse genes above the threshold can still get p-values that are somewhat too small.
- **Negative values.** The similarity columns of a gene with negative values are `NA`, with one warning and a
  note in `message`; the correlation test is unaffected (section 4 kept the old counting of opposite signs).
- **Methods.** `[` keeps the class and attributes for a selection of rows and gives a plain data frame for a
  selection of columns; `print()` and `summary()` treat a result without its essential columns as a data
  frame, and count skipped genes.
- **Progress.** The line leaves out the elapsed time, then is cut, where it is wider than the console; its
  count is the permutations kept before the batch plus those computed in it, and its final count equals the
  `nPermutations` total. Fixed B with fewer permutations than a batch works.
- **Interrupts.** The engine's interrupt check runs under Rcpp's unwind protection, so an R error raised there
  (a time limit) stays an R error.
- **Pixel order (section 7).** The shared pixels stay in the column order of `x`, as in the legacy functions,
  so that the variogram subsample of runs with more than 1000 pixels is the legacy one; reordering the pixels
  of `x`, or swapping `x` and `y`, gives other permutations. This is documented.
