# STcompare 0.1.0.9000

## New main function: `compareSpatial()`

* `compareSpatial(x, y)` compares every gene of two samples rasterized onto one grid, with both tests in one
  call: the spatial correlation test (the Viladomat null of `spatialCorrelationGeneExp()`) and the spatial
  similarity of `spatialSimilarity()`. It returns one row per gene in a plain table (no list columns) of class
  `"STcompareResult"`, with `print()`, `summary()` and `as.data.frame()` methods; the settings, the delta grid
  used, the call and the run time are attributes, and `keepNulls = TRUE` keeps the null correlations and
  delta stars of every permutation.
* Statistical changes compared with the legacy functions, which are unchanged:
  * Adaptive p-values (Besag and Clifford 1991): the two directions of a gene (x permuted, y permuted) run in
    lockstep until either has `exceedances = 10` null correlations at least as extreme as the observed one,
    or until `nPermutations = 10000`. A gene that is clearly not significant stops after a few dozen
    permutations; a significant gene runs to `nPermutations`, so p-values down to 1 / 10001 are resolved.
  * One p-value per gene: `exceedances / L` if the gene stopped early after `L` permutations, otherwise
    `(max(bX, bY) + 1) / (L + 1)`. This is the larger of the two directions' sequential p-values and is valid
    without assuming that the directions are independent. It is adjusted across genes (`adjustMethod = "BH"`).
    `exceedances = Inf` gives a fixed number of permutations with `p = max(pX, pY)`.
  * Independent random streams: the permutations and the noise of each gene, direction and permutation come
    from their own stream (xoshiro256**, keyed by the seed, the gene name and the direction), so the results
    do not depend on the number of threads, the gene order or the other genes in the call. The legacy
    functions share one stream across genes and directions, so that the published results reproduce.
  * The default delta grid is `c(0.01, 0.05, seq(0.1, 0.9, 0.1))`, the grid of the authors' kidney and MERFISH
    analyses; deltas too small for the number of shared pixels (`floor(N * delta) < 2`) are dropped with a
    message. `deltaGridEdge` flags genes whose delta star is mostly at an end of the grid.
  * Rarely detected genes are not tested for correlation (`minDetected`): by default, a gene must be detected
    in at least `sqrt(N)` of the `N` shared pixels of each sample (18 of 311, 47 of 2170). The surrogates'
    values are close to normally distributed, so for genes detected in a few pixels, whose correlation is
    decided by the pixels where both samples detect them, they give p-values that are far too small (on the
    whole AKI raster, genes detected in a single pixel of each section were significant). Sparse genes above
    the threshold can still get p-values that are somewhat too small; the documentation says so.
* Input checks with actionable errors: the shared pixels are matched by name and their coordinates must agree
  (pixels of samples rasterized separately are no longer silently mis-paired); gene names must be unique;
  genes are those present in both samples unless given. A gene with missing values or whose permutations fail
  gets a `"failed"` row with the reason in `message`, and one warning lists such genes. A gene that is constant
  or rarely detected in a sample gets a `"skipped"` row (its similarity is still computed), and one message
  counts such genes, so comparing a whole transcriptome does not warn about the genes that are never detected.
* The similarity of a gene with negative values is `NA`, with one warning that lists such genes: a fold change
  needs values that are not negative (`spatialSimilarity()` stops on them).
* `progress = TRUE` (the default in interactive sessions) shows one progress line with the share of the work
  done, the genes finished, the permutations so far, the elapsed time and an estimate of the remaining time.
  The line fits the width of the console, and its final count of permutations is the total of the result.
  `nThreads` defaults to `getOption("STcompare.nThreads", 1L)`.
* Selecting rows of a result with `[` keeps its class and settings; selecting columns gives a plain data frame.

## A compiled engine replaces the R implementation

* `viladomatCorrelation()`, `spatialCorrelation()`, `spatialCorrelationGeneExp()`,
  `spatialCorrelationGeneExpIterPermutations()` and `spatialCorrelationGeneExpWithinSample()` are now computed
  by compiled code (C++ through Rcpp). Their names, arguments and output structure are unchanged.
  * The method is unchanged: the code reimplements locfit's local-constant Gaussian kernel smoother (its
    adaptive k-d tree with vertex interpolation) and `geoR::variog()`'s binned variogram, and it draws the same
    permutations and noise as before.
  * It reproduces the published analyses (formerly in `inst/extdata`): on the test genes (35 AKI kidney genes on 311
    pixels, 30 brain genes on 2170 pixels, 100 permutations), delta\* is identical for every permutation, and
    the null correlations agree within 3e-14 relative. The empirical p-values are therefore the same.
  * Speed, per gene with 100 permutations in both directions (macOS arm64, M1 Ultra), against the R
    implementation on one core (19-20 s for AKI and 65-69 s for brain):
    * one thread: 0.029 s (AKI, 11 deltas) and 0.13 s (brain, 9 deltas), about 650 and 500 times faster;
    * 16 threads: 0.0035 s and 0.016 s per gene.
  * `nThreads` is now the number of threads of the compiled code. A non-`NULL` `BPPARAM` only sets that number
    (`BiocParallel::bpnworkers(BPPARAM)`); no BiocParallel back-end is started and nothing is forked. Results
    do not depend on the number of threads.
  * `spatialCorrelationGeneExp()` computes all genes in one call, and
    `spatialCorrelationGeneExpIterPermutations()` extends the genes it carries forward from
    `nPermutations[k] + 1` instead of recomputing their first permutations. The results are the same as a fresh
    run with the larger number of permutations.
  * `spatialCorrelationGeneExpWithinSample()` computes the permutations of every gene once and correlates them
    with all the other genes, instead of permuting both genes again for every pair, and it builds its result
    table without looping over the pairs, whose number grows with the square of the number of genes: 600 genes
    (179,700 pairs) with 20 permutations take about 2.6 s on 8 threads of an M1 Ultra.
  * The global random number generator state is left unchanged. The permutations are still drawn with
    `set.seed(seed)` under the session's RNG kinds, and the noise of permutation b with `set.seed(seed + b)`
    under L'Ecuyer-CMRG.
* A user interrupt stops the compiled code within about 0.1 s, and an error that R raises while it runs, such
  as the time limit of `setTimeLimit()`, is an R error that `tryCatch()` can catch. `nThreads` must be a whole
  number.
* Failures no longer print errors. Where the R implementation returned an NA row (for example a constant gene
  or a missing value), there is still an NA row, plus one `warning()` per gene that gives the reason. Inputs that
  crashed R or failed with an unrelated error also give an NA row and a warning:
  * `1 <= delta * N < 2`, where locfit killed the R session;
  * exactly duplicated coordinates with a small delta, where locfit overflowed the C stack;
  * fewer than 3 complete pairs, where `cor.test()` failed and the error handler raised "object 'corDF' not
    found".
  * `spatialCorrelationGeneExpIterPermutations()` no longer stops when a gene has an NA row: such genes are not
    carried forward to the next round.
* An NA row now has the same column types as other rows (numeric p-values, list columns with `NA` elements).
  `viladomatCorrelation()` returns `NA` elements with a warning instead of an error, and when Y is constant or
  has missing values it returns the permutations with `NA` nulls and p-value, as before, plus a warning.
* `verbose = TRUE` prints one message when the computation starts and one when it ends, instead of one per gene.
* `matchingVariograms()` is removed. It was the R helper of the old implementation, and its interface took a
  geoR variogram object.
* geoR and locfit are no longer needed at run time; they moved to Suggests (the tests use them as references).
* Installing STcompare from source needs a C++17 compiler. A source install stops with an explanatory error if
  the compiler flags include `-ffast-math`, `-Ofast` or similar unsafe floating-point options (for example from
  `~/.R/Makevars`): the compiled code reproduces R's arithmetic exactly.

## Changes to results

* Empirical p-values are now `(b + 1) / (B + 1)`.
  * `b` is the number of null correlations whose absolute value is at least the absolute value of the
    observed correlation, and `B` is the number of permutations.
  * Previously they were `b / B` with a strict `>`, so they could be exactly 0. The smallest possible p-value
    is now `1 / (B + 1)`.
  * Affected functions: `viladomatCorrelation()`, `spatialCorrelation()`, `spatialCorrelationGeneExp()`,
    `spatialCorrelationGeneExpIterPermutations()` and `spatialCorrelationGeneExpWithinSample()`.
  * Null correlations and delta stars are unchanged. With the default `alpha = 0.05`,
    `spatialCorrelationGeneExpIterPermutations()` screens the same genes into each round as before.
    With other values of `alpha`, a gene at the screening boundary can be decided differently. When
    `100 * alpha` is a whole number, a gene is now carried forward exactly when it has fewer than
    `100 * alpha` exceedances in each direction; the earlier comparison of `b / B` with the threshold also
    kept some genes with exactly `100 * alpha`, depending on floating-point rounding (for example
    `alpha = 0.07` with 1000 permutations). Otherwise a gene with exactly `floor(100 * alpha)` exceedances
    may no longer be carried forward: with `alpha = 0.053` and 100 permutations, a gene can have at most 4
    exceedances instead of 5.

* `spatialCorrelationGeneExp()` now applies `adjustMethod` across all genes, separately to
  `pValuePermuteX` and `pValuePermuteY`, as documented. Previously it was applied to one gene at a time,
  which left every p-value unadjusted.

## Bug fixes and minor improvements

* `spatialCorrelationGeneExp()` checks `adjustMethod` before computing anything.

* The documentation of `spatialCorrelationGeneExpIterPermutations()` now describes the screening rule the code
  applies: a gene is carried forward when both unadjusted p-values are below `100 * alpha / nPermutations[k]`.

* The examples of `spatialCorrelationGeneExpWithinSample()` and `spatialCorrelationGeneExpIterPermutations()`
  run (they passed a list where one object was needed, and printed an undefined object).

* `spatialSimilarity()` computes all genes at once, with the helper that also computes the similarity columns
  of `compareSpatial()`, and densifies each assay once instead of once per gene: on the published inputs it
  takes 0.2 s instead of 3.9 s for the 325 brain genes, and 3.5 s instead of 25 minutes for the 32,285 genes of
  the sparse AKI raster. Its results are identical, except that:
  * a gene with fewer than `minPixels` kept pixels now reports the number and the IDs of its kept pixels in
    `numPixelInThresh` and `pixelIDInThresh` (it reported 1 and `NA`);
  * the list columns are always of class `"AsIs"` (they were not when the first gene had no score);
  * `parameters` also records the assay used (`assayName`), `minQuantile`, `t1` and `t2`.
* `spatialSimilarity()` stops with an error that names the genes when values are missing, infinite or negative
  (a fold change needs values that are not negative; missing values gave an unclear error from `quantile()`,
  or, with `t1` and `t2` given, were counted as below the threshold, and negative values counted a pixel as
  similar and dissimilar at once), and when the two objects share no pixel name (the similarity was `NaN`
  without a warning). `t1` and `t2` must be single numbers, `minQuantile` and `minPixels` between 0 and 1, and
  `foldChange` not negative. `verbose = TRUE` prints one message instead of one per gene.
* `linearRegression()`, `pixelClass()` and `savePlots()` use the assay that `spatialSimilarity()` classified
  when `assayName` is not given (they used the first assay, so the scatter plot could show other values than
  the ones classified).
* `savePlots()`:
  * without pixel geometries, panel 2 shows the second sample (it showed the first one again);
  * uses `assayName` for every panel (the expression panels of rasterized objects always showed the first assay);
  * no longer attaches ggplot2, gridExtra and patchwork to the search path; patchwork is now a suggested
    package, and gridExtra is not used;
  * the expression panels have the sample names as titles, and their colour bars the same size, in both
    branches.
* `plotCorrelationGeneExp()`:
  * draws pixels with negative values (the axes started at 0, so they were dropped);
  * shows the greater of the two p-values, or `NA` if either is `NA` (it showed `pValuePermuteX` when
    `pValuePermuteY` was `NA`), rounded to 3 significant digits;
  * also takes a `compareSpatial()` result, and shows its `padj`;
  * stops with a clear error for a gene that is not in the results, and no longer prints "Ignoring unknown
    labels".
* `pixelClass()` no longer prints "Coordinate system already present". It needs the sf package (now
  suggested) for rasterized objects.
* `spatialCorrelationGeneExpWithinSample()` stops with a clear error when the object has no row names.
* Documentation:
  * `spatialSimilarity()`: the proportions are called proportions, `numPixelInThresh` counts the pixels above
    the threshold in either object (not both), the `verbose` argument is documented, and the details explain
    the thresholds and the fold-change band (both ends included);
  * `plotCorrelationGeneExp()` and `linearRegression()` describe what they draw (the former had the
    description of `spatialCorrelationGeneExpWithinSample()`);
  * `spatialCorrelation()`: `X` and `Y` are vectors (a matrix with one row or one column is accepted);
  * `simRanPatternRasts`: each dataset keeps 1201 to 1381 of 5000 simulated cells, on 272 to 288 pixels, and
    its `colData` columns are listed; `speKidney` is ordered A, C, B;
  * return values are described with `\describe{}` lists (their item names were lost in the help pages), and
    typos are fixed;
  * the examples run in under 5 seconds each with 2 threads, without starting BiocParallel workers, and the
    `savePlots()` example plots the samples it compared;
  * `?STcompare` describes the package, and `citation("STcompare")` gives the Bioinformatics paper.

## Package size, namespace and dependencies

* The authors' precomputed results (`inst/extdata/*.RData`, 32 MB) are no longer part of the package. They are
  kept in the source repository, in `bench/published/`, as the reference of the validation script and the test
  fixtures; the byte-identical `brainCorrelation_1.RData` was removed. The source tarball is now 2.4 MB with
  the built vignettes (2.1 MB without; it was 50 MB) and the installed package 2.6 MB (34 MB).
  `system.file("extdata", ...)` no longer finds these files.
* The authors' analysis scripts (`inst/scripts/`) are no longer installed with the package: they needed
  GitHub-only packages and the authors' file paths. They are kept in the source repository with the results
  they wrote, in `bench/published/scripts/`.
* The NAMESPACE is generated by roxygen2. Only the documented user-facing functions are exported:
  `getGenePixelDF()` and `assignFill()` are internal now, and `threshold()` is removed (`spatialSimilarity()` no
  longer uses it). The `print()`, `summary()` and `as.data.frame()` methods of `compareSpatial()` results are
  registered as S3 methods instead of being exported as functions.
* DESCRIPTION: R 4.5 is required, as SEraster is in Bioconductor from release 3.21, which requires R 4.5; the
  title is spelled correctly; `biocViews` lets `BiocManager::install()` and
  `remotes::install_github()` find the Bioconductor dependencies; `URL` and `BugReports` point to GitHub;
  `SystemRequirements: C++17`; ggplot2 (>= 3.5.0) is required; patchwork and sf are suggested, class and
  reshape2 (no longer used) are not, and `Config/Needs/website` lists the packages that only the case-study
  articles of the website need.

## Documentation

* The tutorials run the analyses live with `compareSpatial()` instead of loading precomputed results.
  "Getting started" uses the built-in data, including a null-calibration check on 50 independent pairs of
  `simRanPatternRasts` that replaces `simRanPatternResults.RData`. The case studies "Acute kidney injury
  (10x Visium)" and "Comparison of MERFISH and Visium for mouse brain" download and cache their inputs and run
  in one to two minutes on 8 threads.
* New articles: "How STcompare works" (the tests step by step, with figures from the built-in data, and why
  rarely detected genes are not tested) and "Parameters, performance and reproducibility" (choosing the
  settings, run times, threads, seeds, and the legacy functions compared with `compareSpatial()`). Their
  formulas render offline (MathML), and the figures of all tutorials have alternative text.
* The case studies no longer need MERINGUE or scatterbar. The spatially variable genes of the published
  analyses ship in `inst/extdata/vignette-aki-svg-genes.txt` and `vignette-brain-svg-genes.txt`, and
  `data-raw/vignette_gene_lists.R` in the source repository regenerates them.
* The README has the installation with BiocManager (and remotes, which BiocManager needs to install from
  GitHub), the R and C++17 compiler requirements, a quick start, links to every tutorial, and corrected
  descriptions: alignment is done with STalign and rasterization with SEraster, and it is "Pearson's", not
  "Person", correlation.
* The website has Tutorials and Articles menus, a grouped function reference with descriptive titles, and a
  changelog. The case studies and the installation page are pkgdown articles, not package vignettes.
