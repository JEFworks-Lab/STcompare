# Changelog

## STcompare 2.0.0

### Changes

- [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
  and
  [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  return the unadjusted (raw) empirical p-values in `pValuePermuteX` and
  `pValuePermuteY`, `(b + 1) / (B + 1)` for each direction, and two new
  columns right after them, `pValuePermuteXAdj` and `pValuePermuteYAdj`,
  hold these p-values adjusted across genes with `adjustMethod` (each
  direction on its own). Earlier versions had no raw columns:
  [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  returned the adjusted values under the names `pValuePermuteX` and
  `pValuePermuteY`, while the released 0.1.0
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
  called [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) on one
  gene at a time, so the values it returned under these names were in
  fact unadjusted (see “Changes to results”). Code that called a gene
  significant when both p-values of the iterative function were below
  0.05 must now read `pValuePermuteXAdj` and `pValuePermuteYAdj`; the
  published tables of the iterative analyses in `bench/published/` keep
  the adjusted values under the old names
  (`kidneyCorrelationNoIter.RData`, from
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md),
  holds unadjusted values). The screening between the rounds of
  [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  has always used the raw p-values. A p-value column carries the suffix
  `Adj` exactly when it is adjusted, in every function:
  [`spatialCorrelation()`](https://jef.works/STcompare/reference/spatialCorrelation.md),
  [`viladomatCorrelation()`](https://jef.works/STcompare/reference/viladomatCorrelation.md)
  and
  [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md)
  never adjust and keep their columns.
  [`plotCorrelationGeneExp()`](https://jef.works/STcompare/reference/plotCorrelationGeneExp.md)
  shows the greater of the two adjusted p-values of a legacy result, now
  labelled `p_adj`; a table saved before 2.0.0 has no `Adj` columns, so
  its `pValuePermuteX` and `pValuePermuteY` are shown with the `p_E`
  label even when the iterative function had adjusted them. The nulls,
  delta stars, `correlationCoef` and `pValueNaive` are unchanged.
- The columns of the new
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
  have the names of the legacy functions where the quantity is the same
  (`correlationCoef`, `pValueNaive`, `pValuePermuteX`, `pValuePermuteY`,
  `percentSimilarity`, …), plus `pValuePermute` and `pValuePermuteAdj`;
  the pre-release development versions on GitHub used `r`, `pNaive`,
  `pX`, `pY`, `p`, `padj`, `similarity`, `dissimilarityX`,
  `dissimilarityY`, `nPixelsSimilarity`, `thresholdX` and `thresholdY`.
- `matchingVariograms()` is removed (see “A compiled engine replaces the
  R implementation”).
- Empirical p-values are `(b + 1) / (B + 1)` instead of `b / B`, so they
  are never 0 (see “Changes to results”).

### New main function: `compareSpatial()`

- `compareSpatial(x, y)` compares every gene of two samples rasterized
  onto one grid, with both tests in one call: the spatial correlation
  test (the Viladomat null of
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md))
  and the spatial similarity of
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).
  It returns one row per gene in a plain table (no list columns) of
  class `"STcompareResult"`, with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods; the settings, the delta grid used, the call and the run time
  are attributes, and `keepNulls = TRUE` keeps the null correlations and
  delta stars of every permutation.
- The columns have the names of the legacy functions where the quantity
  is the same, so that the two kinds of result read alike:
  `correlationCoef`, `pValueNaive`, `pValuePermuteX` and
  `pValuePermuteY` as in
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md),
  and `percentSimilarity`, `percentDissimilarityX`,
  `percentDissimilarityY`, `numPixelInThresh`, `t1` and `t2` as in
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).
  The combined permutation p-value is `pValuePermute`, and
  `pValuePermuteAdj` is its adjustment across genes, the column to use
  to call a gene significant. The help page describes every column in
  order, says what each p-value is and where it differs from the legacy
  functions (`pValuePermuteX` and `pValuePermuteY` are sequential
  p-values at the common stopping point here, and fixed-B p-values in
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md),
  which adjusts each direction on its own, `pValuePermuteXAdj` and
  `pValuePermuteYAdj`, and has no combined p-value), documents the
  attributes, and has a paragraph “Which p-value should I report?”.
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) name the columns
  they count.
- Statistical changes compared with the legacy functions, which are
  unchanged:
  - Adaptive p-values (Besag and Clifford 1991): the two directions of a
    gene (x permuted, y permuted) run in lockstep until either has
    `exceedances = 10` null correlations at least as extreme as the
    observed one, or until `nPermutations = 10000`. A gene that is
    clearly not significant stops after a few dozen permutations; a
    significant gene runs to `nPermutations`, so p-values down to 1 /
    10001 are resolved.
  - One p-value per gene: `exceedances / L` if the gene stopped early
    after `L` permutations, otherwise `(max(bX, bY) + 1) / (L + 1)`.
    This is the larger of the two directions’ sequential p-values and is
    valid without assuming that the directions are independent
    (`pValuePermute`). It is adjusted across genes
    (`adjustMethod = "BH"`, `pValuePermuteAdj`). `exceedances = Inf`
    gives a fixed number of permutations with
    `pValuePermute = max(pValuePermuteX, pValuePermuteY)`.
  - Independent random streams: the permutations and the noise of each
    gene, direction and permutation come from their own stream
    (xoshiro256\*\*, keyed by the seed, the gene name and the
    direction), so the results do not depend on the number of threads,
    the gene order or the other genes in the call. The legacy functions
    share one stream across genes and directions, so that the published
    results reproduce.
  - The default delta grid is `c(0.01, 0.05, seq(0.1, 0.9, 0.1))`, the
    grid of the authors’ kidney and MERFISH analyses; deltas too small
    for the number of shared pixels (`floor(N * delta) < 2`) are dropped
    with a message. `deltaGridEdge` flags genes whose delta star is
    mostly at an end of the grid.
  - Rank-remapped surrogates (`surrogate = "remap"`, the default): the
    values of every surrogate are replaced by the gene’s own values in
    the surrogate’s rank order (the amplitude adjustment of AAFT
    surrogates, Theiler et al. 1992), so that each surrogate has exactly
    the gene’s distribution of values, zeros included, and keeps the
    spatial arrangement of the Viladomat surrogate. The Viladomat
    surrogates alone are smoothed and mixed with Gaussian noise, so
    their values are close to normally distributed; for a gene detected
    in few pixels, whose correlation is decided by the pixels where both
    samples detect it, they give p-values that are far too small: in
    simulations on the kidney grid, independent genes detected in 10 of
    311 pixels had P(p \<= 0.001) = 0.03, Benjamini-Hochberg found 9
    false discoveries among 1800 independent sparse genes, and on the
    whole AKI raster genes detected in a single pixel of each section
    were significant. For a gene without spatial structure the remapped
    surrogates are plain permutations of its values, so the test is the
    exact permutation test whatever the gene’s distribution. On the
    kidney and brain grids, remapped surrogates kept P(p \<= alpha) at
    or below alpha for independent genes detected in 1 to 100 percent of
    the pixels, had the same power on correlated fields, recovered the
    published genes at least as well (728 of the 731 kidney genes
    against 723 with gaussian surrogates; 124 of the 128 brain genes
    against 122), and cost about 7 percent more per permutation
    (`bench/calibrate-surrogates.R`, results in
    `bench/calibration-results.md`). In this mode a null correlation
    that equals the observed one up to rounding (within a relative 1e-9)
    counts as an exceedance, as in the exact test: many rearrangements
    of a sparse gene give exactly the observed correlation, and
    [`cor()`](https://rdrr.io/r/stats/cor.html) rounds each of them
    differently.
  - `surrogate = "gaussian"` keeps the surrogates as the smoothing and
    the noise leave them: the surrogates of the legacy functions and of
    the published analyses, for comparability with them.
  - `minDetected`: a gene must be detected in at least this share of the
    shared pixels of each sample to be tested for correlation (a gene is
    detected at a pixel where its value is above its lowest value; `0`
    tests every gene that is not constant). The default `NULL` means 0,
    no filter, with `surrogate = "remap"`, whose p-values were
    calibrated down to genes detected in a single pixel, and `sqrt(N)`
    of the `N` shared pixels (18 of 311, 47 of 2170) with
    `surrogate = "gaussian"`, which removes the genes whose p-values are
    worst (sparse genes above that filter can still get p-values that
    are somewhat too small). The number of pixels asked for is recorded
    in `attr(result, "params")$minDetectedPixels`, and
    [`print()`](https://rdrr.io/r/base/print.html) and
    [`summary()`](https://rdrr.io/r/base/summary.html) name the
    surrogate mode.
- Input checks with actionable errors: the shared pixels are matched by
  name and their coordinates must agree (pixels of samples rasterized
  separately are no longer silently mis-paired); gene names must be
  unique; genes are those present in both samples unless given. A gene
  with missing values or whose permutations fail gets a `"failed"` row
  with the reason in `message`, and one warning lists such genes. A gene
  that is constant in a sample (or detected in fewer pixels than
  `minDetected` asks for) gets a `"skipped"` row (its similarity is
  still computed), and one message counts such genes, so comparing a
  whole transcriptome does not warn about the genes that are never
  detected.
- The similarity of a gene with negative values is `NA`, with one
  warning that lists such genes: a fold change needs values that are not
  negative
  ([`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  stops on them).
- `progress = TRUE` (the default in interactive sessions) shows one
  progress line with the share of the work done, the genes finished, the
  permutations so far, the elapsed time and an estimate of the remaining
  time. The line fits the width of the console, and its final count of
  permutations is the total of the result. `nThreads` defaults to
  `getOption("STcompare.nThreads", 1L)`.
- Selecting rows of a result with `[` keeps its class and settings;
  selecting columns gives a plain data frame.

### Visualization

- `plotGene(result, gene, x, y)` draws one standardized figure for one
  gene of a
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
  result, ready to print or to save with
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html):
  the gene’s expression in the two samples at native resolution
  (optional, from the objects before rasterization,
  `native = list(...)`), on the shared grid of pixels (with the pixel
  geometry of SEraster, or tiles), the similarity class of every pixel
  (similar, higher in either sample, below the threshold; the legend
  gives the share of each class, so the share of similar pixels is
  `percentSimilarity`), and a scatter plot of the expression at matched
  pixels coloured by class with the line y = x, the fold-change band and
  the thresholds (the figure’s caption is the key), subtitled with the
  gene’s statistics: `correlationCoef`, `pValuePermute` (and the number
  of permutations), `pValuePermuteAdj`, `pValueNaive` and
  `percentSimilarity`, with the reason when a statistic is `NA`. The
  panels are lettered A to D (one letter for each pair of maps, each map
  titled with the sample name), the cells of the native maps are points
  sized for their number; the colour scales are viridis, the class
  colours are distinguishable under colour-vision deficiencies, a pair
  of maps shares one colour scale (`sharedScale = FALSE` gives each map
  its own), `transform = "log1p"` puts the maps and the scatter axes on
  a log1p scale, pixels that one sample lacks are drawn in light grey in
  every map, and `nativeUnit` names the points of the native maps
  (“cells” or “spots”). Every text, legend key and point is scaled from
  `baseSize` (10 points by default; 7 for a journal figure 180 mm wide),
  and the figure is a patchwork object, so it can be changed with
  patchwork and ggplot2. The function warns when the objects are not
  those the result was computed on (another number of shared pixels, or
  `x` and `y` swapped).
- `plotOverview(result)` draws one figure for all the genes of a result:
  `correlationCoef` against `-log10(pValuePermuteAdj)` with the
  significant genes coloured by the sign of the correlation, a histogram
  of `percentSimilarity`, and `correlationCoef` against
  `percentSimilarity`; `label` names genes to label and `baseSize`
  scales the figure. Panels whose columns are absent (a result of one
  test) are left out.
- The tutorials use these figures: Getting started explains the panels
  and how to save and customize a figure, and the case studies show
  [`plotOverview()`](https://jef.works/STcompare/reference/plotOverview.md)
  and [`plotGene()`](https://jef.works/STcompare/reference/plotGene.md)
  on the published analyses.
- patchwork, which arranges the panels, is now imported (it was
  suggested).

### A compiled engine replaces the R implementation

- [`viladomatCorrelation()`](https://jef.works/STcompare/reference/viladomatCorrelation.md),
  [`spatialCorrelation()`](https://jef.works/STcompare/reference/spatialCorrelation.md),
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md),
  [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  and
  [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md)
  are now computed by compiled code (C++ through Rcpp). Their names,
  arguments and output structure are unchanged.
  - The method is unchanged: the code reimplements locfit’s
    local-constant Gaussian kernel smoother (its adaptive k-d tree with
    vertex interpolation) and
    [`geoR::variog()`](https://rdrr.io/pkg/geoR/man/variog.html)’s
    binned variogram, and it draws the same permutations and noise as
    before.
  - It reproduces the published analyses (formerly in `inst/extdata`):
    on the test genes (35 AKI kidney genes on 311 pixels, 30 brain genes
    on 2170 pixels, 100 permutations), delta\* is identical for every
    permutation, and the null correlations agree within 3e-14 relative.
    The empirical p-values are therefore the same.
  - Speed, per gene with 100 permutations in both directions (macOS
    arm64, M1 Ultra), against the R implementation on one core (19-20 s
    for AKI and 65-69 s for brain):
    - one thread: 0.029 s (AKI, 11 deltas) and 0.13 s (brain, 9 deltas),
      about 650 and 500 times faster;
    - 16 threads: 0.0035 s and 0.016 s per gene.
  - `nThreads` is now the number of threads of the compiled code. A
    non-`NULL` `BPPARAM` only sets that number
    (`BiocParallel::bpnworkers(BPPARAM)`); no BiocParallel back-end is
    started and nothing is forked. Results do not depend on the number
    of threads.
  - [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
    computes all genes in one call, and
    [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
    extends the genes it carries forward from `nPermutations[k] + 1`
    instead of recomputing their first permutations. The results are the
    same as a fresh run with the larger number of permutations.
  - [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md)
    computes the permutations of every gene once and correlates them
    with all the other genes, instead of permuting both genes again for
    every pair, and it builds its result table without looping over the
    pairs, whose number grows with the square of the number of genes:
    600 genes (179,700 pairs) with 20 permutations take about 2.6 s on 8
    threads of an M1 Ultra.
  - The global random number generator state is left unchanged. The
    permutations are still drawn with `set.seed(seed)` under the
    session’s RNG kinds, and the noise of permutation b with
    `set.seed(seed + b)` under L’Ecuyer-CMRG.
- A user interrupt stops the compiled code within about 0.1 s, and an
  error that R raises while it runs, such as the time limit of
  [`setTimeLimit()`](https://rdrr.io/r/base/setTimeLimit.html), is an R
  error that [`tryCatch()`](https://rdrr.io/r/base/conditions.html) can
  catch. `nThreads` must be a whole number.
- Failures no longer print errors. Where the R implementation returned
  an NA row (for example a constant gene or a missing value), there is
  still an NA row, plus one
  [`warning()`](https://rdrr.io/r/base/warning.html) per gene that gives
  the reason. Inputs that crashed R or failed with an unrelated error
  also give an NA row and a warning:
  - `1 <= delta * N < 2`, where locfit killed the R session;
  - exactly duplicated coordinates with a small delta, where locfit
    overflowed the C stack;
  - fewer than 3 complete pairs, where
    [`cor.test()`](https://rdrr.io/r/stats/cor.test.html) failed and the
    error handler raised “object ‘corDF’ not found”.
  - [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
    no longer stops when a gene has an NA row: such genes are not
    carried forward to the next round.
- An NA row now has the same column types as other rows (numeric
  p-values, list columns with `NA` elements).
  [`viladomatCorrelation()`](https://jef.works/STcompare/reference/viladomatCorrelation.md)
  returns `NA` elements with a warning instead of an error, and when Y
  is constant or has missing values it returns the permutations with
  `NA` nulls and p-value, as before, plus a warning.
- `verbose = TRUE` prints one message when the computation starts and
  one when it ends, instead of one per gene.
- `matchingVariograms()` is removed. It was the R helper of the old
  implementation, and its interface took a geoR variogram object.
- geoR and locfit are no longer needed at run time; they moved to
  Suggests (the tests use them as references).
- Installing STcompare from source needs a C++17 compiler. A source
  install stops with an explanatory error if the compiler flags include
  `-ffast-math`, `-Ofast` or similar unsafe floating-point options (for
  example from `~/.R/Makevars`): the compiled code reproduces R’s
  arithmetic exactly.

### Changes to results

- Empirical p-values are now `(b + 1) / (B + 1)`.
  - `b` is the number of null correlations whose absolute value is at
    least the absolute value of the observed correlation, and `B` is the
    number of permutations.
  - Previously they were `b / B` with a strict `>`, so they could be
    exactly 0. The smallest possible p-value is now `1 / (B + 1)`.
  - Affected functions:
    [`viladomatCorrelation()`](https://jef.works/STcompare/reference/viladomatCorrelation.md),
    [`spatialCorrelation()`](https://jef.works/STcompare/reference/spatialCorrelation.md),
    [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md),
    [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
    and
    [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md).
  - Null correlations and delta stars are unchanged. With the default
    `alpha = 0.05`,
    [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
    screens the same genes into each round as before. With other values
    of `alpha`, a gene at the screening boundary can be decided
    differently. When `100 * alpha` is a whole number, a gene is now
    carried forward exactly when it has fewer than `100 * alpha`
    exceedances in each direction; the earlier comparison of `b / B`
    with the threshold also kept some genes with exactly `100 * alpha`,
    depending on floating-point rounding (for example `alpha = 0.07`
    with 1000 permutations). Otherwise a gene with exactly
    `floor(100 * alpha)` exceedances may no longer be carried forward:
    with `alpha = 0.053` and 100 permutations, a gene can have at most 4
    exceedances instead of 5.
- [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
  now applies `adjustMethod` across all genes, separately to
  `pValuePermuteX` and `pValuePermuteY`, as documented, and returns the
  result in the new columns `pValuePermuteXAdj` and `pValuePermuteYAdj`
  (see “Changes”). Previously it was applied to one gene at a time,
  which left every p-value unadjusted.

### Bug fixes and minor improvements

- [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
  checks `adjustMethod` before computing anything.

- The documentation of
  [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  now describes the screening rule the code applies: a gene is carried
  forward when both unadjusted p-values are below
  `100 * alpha / nPermutations[k]`.

- The examples of
  [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md)
  and
  [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  run (they passed a list where one object was needed, and printed an
  undefined object).

- [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  computes all genes at once, with the helper that also computes the
  similarity columns of
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md),
  and densifies each assay once instead of once per gene: on the
  published inputs it takes 0.2 s instead of 3.9 s for the 325 brain
  genes, and 3.5 s instead of 25 minutes for the 32,285 genes of the
  sparse AKI raster. Its results are identical, except that:

  - a gene with fewer than `minPixels` kept pixels now reports the
    number and the IDs of its kept pixels in `numPixelInThresh` and
    `pixelIDInThresh` (it reported 1 and `NA`);
  - the list columns are always of class `"AsIs"` (they were not when
    the first gene had no score);
  - `parameters` also records the assay used (`assayName`),
    `minQuantile`, `t1` and `t2`.

- [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  stops with an error that names the genes when values are missing,
  infinite or negative (a fold change needs values that are not
  negative; missing values gave an unclear error from
  [`quantile()`](https://rdrr.io/r/stats/quantile.html), or, with `t1`
  and `t2` given, were counted as below the threshold, and negative
  values counted a pixel as similar and dissimilar at once), and when
  the two objects share no pixel name (the similarity was `NaN` without
  a warning). `t1` and `t2` must be single numbers, `minQuantile` and
  `minPixels` between 0 and 1, and `foldChange` not negative.
  `verbose = TRUE` prints one message instead of one per gene.

- [`linearRegression()`](https://jef.works/STcompare/reference/linearRegression.md),
  [`pixelClass()`](https://jef.works/STcompare/reference/pixelClass.md)
  and
  [`savePlots()`](https://jef.works/STcompare/reference/savePlots.md)
  use the assay that
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  classified when `assayName` is not given (they used the first assay,
  so the scatter plot could show other values than the ones classified).

- [`savePlots()`](https://jef.works/STcompare/reference/savePlots.md):

  - without pixel geometries, panel 2 shows the second sample (it showed
    the first one again);
  - uses `assayName` for every panel (the expression panels of
    rasterized objects always showed the first assay);
  - no longer attaches ggplot2, gridExtra and patchwork to the search
    path (patchwork is imported, and gridExtra is not used);
  - the expression panels have the sample names as titles, and their
    colour bars the same size, in both branches.

- [`plotCorrelationGeneExp()`](https://jef.works/STcompare/reference/plotCorrelationGeneExp.md):

  - draws pixels with negative values (the axes started at 0, so they
    were dropped);
  - shows the greater of the two p-values, or `NA` if either is `NA` (it
    showed `pValuePermuteX` when `pValuePermuteY` was `NA`), rounded to
    3 significant digits;
  - also takes a
    [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
    result, and shows its `pValuePermuteAdj` (labelled `p_adj`);
  - stops with a clear error for a gene that is not in the results, and
    no longer prints “Ignoring unknown labels”.

- [`pixelClass()`](https://jef.works/STcompare/reference/pixelClass.md)
  no longer prints “Coordinate system already present”. It needs the sf
  package (now suggested) for rasterized objects.

- [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md)
  stops with a clear error when the object has no row names.

- Documentation:

  - [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md):
    the proportions are called proportions, `numPixelInThresh` counts
    the pixels above the threshold in either object (not both), the
    `verbose` argument is documented, and the details explain the
    thresholds and the fold-change band (both ends included);
  - [`plotCorrelationGeneExp()`](https://jef.works/STcompare/reference/plotCorrelationGeneExp.md)
    and
    [`linearRegression()`](https://jef.works/STcompare/reference/linearRegression.md)
    describe what they draw (the former had the description of
    [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md));
  - [`spatialCorrelation()`](https://jef.works/STcompare/reference/spatialCorrelation.md):
    `X` and `Y` are vectors (a matrix with one row or one column is
    accepted);
  - `simRanPatternRasts`: each dataset keeps 1201 to 1381 of 5000
    simulated cells, on 272 to 288 pixels, and its `colData` columns are
    listed; `speKidney` is ordered A, C, B;
  - return values are described with `\describe{}` lists (their item
    names were lost in the help pages), and typos are fixed;
  - the examples run in under 5 seconds each with 2 threads, without
    starting BiocParallel workers, and the
    [`savePlots()`](https://jef.works/STcompare/reference/savePlots.md)
    example plots the samples it compared;
  - [`?STcompare`](https://jef.works/STcompare/reference/STcompare-package.md)
    describes the package, and `citation("STcompare")` gives the
    Bioinformatics paper.

### Package size, namespace and dependencies

- The authors’ precomputed results (`inst/extdata/*.RData`, 32 MB) are
  no longer part of the package. They are kept in the source repository,
  in `bench/published/`, as the reference of the validation script and
  the test fixtures; the byte-identical `brainCorrelation_1.RData` was
  removed. The source tarball is now 2.4 MB with the built vignettes
  (2.1 MB without; it was 50 MB) and the installed package 2.6 MB (34
  MB). `system.file("extdata", ...)` no longer finds these files.
- The authors’ analysis scripts (`inst/scripts/`) are no longer
  installed with the package: they needed GitHub-only packages and the
  authors’ file paths. They are kept in the source repository with the
  results they wrote, in `bench/published/scripts/`.
- The NAMESPACE is generated by roxygen2. Only the documented
  user-facing functions are exported: `getGenePixelDF()` and
  `assignFill()` are internal now, and `threshold()` is removed
  ([`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  no longer uses it). The
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  of
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
  results are registered as S3 methods instead of being exported as
  functions.
- DESCRIPTION: R 4.5 is required, as SEraster is in Bioconductor from
  release 3.21, which requires R 4.5; the title is spelled correctly;
  `biocViews` lets
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  and
  [`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html)
  find the Bioconductor dependencies; `URL` and `BugReports` point to
  GitHub; `SystemRequirements: C++17`; ggplot2 (\>= 3.5.0) and patchwork
  are required; sf and Matrix are suggested, class and reshape2 (no
  longer used) are not, and `Config/Needs/website` lists the packages
  that only the case-study articles of the website need.

### Documentation

- The tutorials run the analyses live with
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
  instead of loading precomputed results. “Getting started” uses the
  built-in data, including a null-calibration check on 50 independent
  pairs of `simRanPatternRasts` that replaces
  `simRanPatternResults.RData`. The case studies “Acute kidney injury
  (10x Visium)” and “Comparison of MERFISH and Visium for mouse brain”
  download and cache their inputs and run in one to two minutes on 8
  threads.
- New articles: “How STcompare works” (the tests step by step, with
  figures from the built-in data, and why the surrogates take the gene’s
  own values) and “Parameters, performance and reproducibility”
  (choosing the settings, run times, threads, seeds, and the legacy
  functions compared with
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)).
  Their formulas render offline (MathML), and the figures of all
  tutorials have alternative text.
- The case studies no longer need MERINGUE or scatterbar. The spatially
  variable genes of the published analyses ship in
  `inst/extdata/vignette-aki-svg-genes.txt` and
  `vignette-brain-svg-genes.txt`, and `data-raw/vignette_gene_lists.R`
  in the source repository regenerates them.
- The README has the installation with BiocManager (and remotes, which
  BiocManager needs to install from GitHub), the R and C++17 compiler
  requirements, a quick start, links to every tutorial, and corrected
  descriptions: alignment is done with STalign and rasterization with
  SEraster, and it is “Pearson’s”, not “Person”, correlation.
- The website has Tutorials and Articles menus, a grouped function
  reference with descriptive titles, and a changelog. The case studies
  and the installation page are pkgdown articles, not package vignettes.
