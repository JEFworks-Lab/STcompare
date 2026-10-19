# Compare the spatial expression patterns of genes between two samples

For every gene measured in two samples on one shared grid of pixels,
`compareSpatial()` tests whether the gene's expression patterns are
correlated beyond what the spatial autocorrelation of each pattern alone
would produce (spatial correlation test), and how often the two samples'
values at the same pixel are within a fold change of each other (spatial
similarity). It returns one row per gene.

## Usage

``` r
compareSpatial(
  x,
  y = NULL,
  assay = 1,
  genes = NULL,
  tests = c("correlation", "similarity"),
  nPermutations = 10000,
  exceedances = 10,
  delta = c(0.01, 0.05, seq(0.1, 0.9, 0.1)),
  maxDistPrctile = 0.25,
  minDetected = NULL,
  surrogate = c("remap", "gaussian"),
  foldChange = 1,
  minQuantile = 0.05,
  minPixels = 0.1,
  adjustMethod = "BH",
  seed = 0L,
  nThreads = getOption("STcompare.nThreads", 1L),
  progress = interactive(),
  verbose = TRUE,
  keepNulls = FALSE
)
```

## Arguments

- x, y:

  Two `SpatialExperiment` objects with the same pixels, for example the
  elements of the list returned by one
  [`SEraster::rasterizeGeneExpression()`](https://rdrr.io/pkg/SEraster/man/rasterizeGeneExpression.html)
  call. Alternatively `x` is a list of two such objects and `y` is
  `NULL`, such as the list returned by
  `SEraster::rasterizeGeneExpression(list(a = ..., b = ...))`; the names
  of the list, if any, are used as the sample labels (otherwise `"x"`
  and `"y"`).

- assay:

  The assay to use: a name or an index, used for both objects, or two
  names or two indices, one per object. Default: the first assay.

- genes:

  The genes (row names) to compare. Default `NULL`: every gene present
  in both objects (a message says how many were dropped). Genes given
  here must be present in both objects.

- tests:

  Which tests to run: `"correlation"` (the spatial correlation test with
  permutations) and/or `"similarity"` (no permutations).

- nPermutations:

  The largest number of permutations per gene; the smallest possible
  p-value is `1 / (nPermutations + 1)`. Default 10000.

- exceedances:

  The number of null correlations at least as extreme as the observed
  one after which a gene stops early (\\h\\ of Besag and Clifford). The
  relative standard error of a p-value is about `1 / sqrt(exceedances)`.
  Default 10; `Inf` gives every gene exactly `nPermutations`
  permutations.

- delta:

  The candidate smoothing bandwidths, as shares of the pixels in the
  kernel's window (values in (0, 1\]). Values that leave fewer than 2
  pixels in the window (\\\lfloor N\delta \rfloor \< 2\\ for \\N\\
  shared pixels) are dropped with a message. The default adds 0.01 and
  0.05 to `seq(0.1, 0.9, 0.1)`, the grid the authors used for the kidney
  and MERFISH analyses: with the shorter grid, delta star is often its
  smallest value.

- maxDistPrctile:

  The quantile of the pairwise pixel distances up to which the
  variograms are computed. Default 0.25.

- minDetected:

  The smallest share of the shared pixels at which a gene must be
  detected, in each sample, to be tested for correlation (see *Rarely
  detected genes* in Details); `0` tests every gene that is not
  constant. Default `NULL`, which depends on `surrogate`: `0` (no
  filter) with `"remap"`, whose p-values are calibrated for rarely
  detected genes, and \\\sqrt{N}\\ of the \\N\\ shared pixels (18 of
  311, 47 of 2170) with `"gaussian"`, whose p-values are far too small
  for such genes. The number of pixels asked for is recorded in
  `attr(result, "params")$minDetectedPixels`.

- surrogate:

  How the surrogates of the correlation test take their values (see
  *Rank-remapped surrogates* in Details): `"remap"` (the default) gives
  each surrogate the gene's own values, in the surrogate's rank order,
  so that it has exactly the gene's distribution of values, which
  calibrates the p-values of sparse and skewed genes; `"gaussian"` uses
  them as the smoothing and the added noise leave them, the surrogates
  of the legacy functions and of the published analyses.

- foldChange:

  The similarity band: pixels with \\\|\log_2(y / x)\| \le\\
  `foldChange` are similar. Default 1 (within two-fold).

- minQuantile:

  The quantile of each gene's values in each sample used as its
  expression threshold for the similarity. Default 0.05.

- minPixels:

  The smallest share of the pixels that must pass the threshold for the
  similarity to be computed. Default 0.1.

- adjustMethod:

  The multiple-testing adjustment of `pValuePermute` across the tested
  genes, which gives `pValuePermuteAdj`; one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html).
  Default `"BH"`.

- seed:

  An integer seed for the permutations and the noise. Default 0.

- nThreads:

  The number of threads. Default `getOption("STcompare.nThreads", 1L)`,
  so it can be set once per session with
  `options(STcompare.nThreads = 8)`. The results do not depend on it.

- progress:

  Show a progress line (with an estimate of the remaining time) while
  the permutations run. Default
  [`interactive()`](https://rdrr.io/r/base/interactive.html). It is
  written as a message, so
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) also hides
  it.

- verbose:

  Show informative messages: genes or deltas that were dropped, genes
  whose delta star is at the edge of the grid, p-values limited by
  `nPermutations`, and the total time. Default `TRUE`.

- keepNulls:

  Keep the null correlations and delta stars of every permutation in
  `attr(result, "details")`. Default `FALSE`.

## Value

A data frame of class `"STcompareResult"` with one row per gene (row
names are the genes) and the columns below, in this order. The columns
of a test that was not run are absent. The column names are those of the
legacy functions,
[`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
and
[`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md),
where the two report the same quantity; the differences in meaning are
stated for each column. In the formulas, \\L\\ is the number of
permutations of the gene (the `nPermutations` column), and \\b_x\\ and
\\b_y\\ count its null correlations whose absolute value is at least the
observed \\\|r\|\\ when `x` and when `y` is permuted.

- `gene`:

  The gene (also the row name).

- `nPixels`:

  The number of shared pixels, the same for every gene.

- `correlationCoef`:

  Pearson's correlation \\r\\ of the gene's values in `x` and `y` over
  the shared pixels ([`cor()`](https://rdrr.io/r/stats/cor.html)).
  Reported for every gene with finite values, tested or not. The same as
  in
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md).

- `pValueNaive`:

  The analytical p-value of
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html), which assumes
  that the pixels are independent. For reference only: under spatial
  autocorrelation it is far too small. The same as in
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md).

- `pValuePermuteX`:

  The sequential permutation p-value of the direction that permutes `x`,
  \\(b_x + 1) / (L + 1)\\, at the common stopping point \\L\\ of the two
  directions. Not adjusted. Descriptive: when a gene stops early, the
  direction that did not reach `exceedances` has an inflated value. In
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
  this column holds \\(b_x + 1) / (B + 1)\\ for a fixed number \\B\\ of
  permutations, also unadjusted; its adjustment across genes is that
  function's `pValuePermuteXAdj`.

- `pValuePermuteY`:

  The same for the direction that permutes `y`, \\(b_y + 1) / (L + 1)\\.
  Not adjusted.

- `pValuePermute`:

  The permutation p-value of the gene, which combines the two
  directions: `exceedances / L` if the gene stopped early
  (`stop = "exceedances"`), otherwise \\(\max(b_x, b_y) + 1) / (L +
  1)\\. This is the larger of the two directions' sequential p-values,
  which is valid without assuming that the two directions are
  independent (see Details); it equals
  `max(pValuePermuteX, pValuePermuteY)` for a gene that reached
  `nPermutations` (`stop = "limit"`) and is slightly smaller,
  `exceedances / L` instead of `(exceedances + 1) / (L + 1)`, for a gene
  that stopped early. Not adjusted; its smallest possible value is
  `1 / (nPermutations + 1)`. The legacy functions have no such column:
  they adjust each direction on its own (`pValuePermuteXAdj`,
  `pValuePermuteYAdj`).

- `pValuePermuteAdj`:

  `pValuePermute` adjusted across the tested genes with `adjustMethod`
  ([`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); skipped
  and failed genes are not counted). This is the column to use to call a
  gene significant.

- `nPermutations`:

  The number of permutations \\L\\ the gene used in each direction:
  fewer than the `nPermutations` argument when the gene stopped early,
  equal to it otherwise. The sum of this column is the work of the call.

- `stop`:

  Why the permutations stopped: `"exceedances"` (early: a direction
  reached `exceedances` null correlations at least as extreme as the
  observed one), `"limit"` (the `nPermutations` argument was reached
  without that, so `pValuePermute` is \\(\max(b_x, b_y) + 1) / (L +
  1)\\), `"skipped"` (not tested) or `"failed"`.

- `deltaStarMedianX`, `deltaStarMedianY`:

  The median over the permutations of delta star, the bandwidth of
  `delta` whose surrogate variogram matched best, when `x` and when `y`
  is permuted. The same quantities as in
  [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md).

- `deltaGridEdge`:

  `TRUE` when, in either direction, more than half of the permutations
  chose the smallest or the largest delta of the grid used
  (`attr(result, "params")$deltaGrid`): the best bandwidth may lie
  outside the grid, so consider extending `delta`. `NA` with fewer than
  3 deltas.

- `percentSimilarity`:

  The proportion (between 0 and 1, despite the name) of the kept pixels
  with \\\|\log_2(y / x)\| \le\\ `foldChange`: the gene's spatial
  similarity. The same as in
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).

- `percentDissimilarityX`:

  The proportion of the kept pixels with \\\log_2(y / x) \<
  -\\`foldChange`: higher in `x`. The same as in
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).

- `percentDissimilarityY`:

  The proportion of the kept pixels with \\\log_2(y / x) \>\\
  `foldChange`: higher in `y`. The same as in
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).

- `numPixelInThresh`:

  The number of kept pixels: those whose value is above `t1` in `x` or
  above `t2` in `y`. The same as in
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).

- `t1`, `t2`:

  The gene's expression thresholds in `x` and in `y`: its `minQuantile`
  quantile in each sample. The same as in
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).

- `status`:

  `"ok"`, `"skipped"` (not tested for correlation: constant or rarely
  detected in a sample) or `"failed"` (see Details).

- `message`:

  Why a gene was skipped or failed, or why its similarity is `NA`
  (negative values); empty otherwise.

The result has four attributes:

- `params`:

  A list of the settings used: every argument after its default was
  applied (`tests`, `nPermutations`, `exceedances`, `delta`,
  `maxDistPrctile`, `minDetected`, `surrogate`, `foldChange`,
  `minQuantile`, `minPixels`, `adjustMethod`, `seed`, `nThreads`,
  `progress`, `verbose`, `keepNulls`), plus `samples` (the two sample
  labels), `assay` (the assay used in each sample), `nPixels` (the
  number of shared pixels), `genes`, `deltaGrid` (the deltas actually
  used, after the unusable ones were dropped) and `minDetectedPixels`
  (the number of detected pixels a gene needed in each sample to be
  tested; `NA` without the correlation test).

- `call`:

  The matched call.

- `runtime`:

  The run time, as
  [`proc.time()`](https://rdrr.io/r/base/proc.time.html) measures it.

- `details`:

  Only with `keepNulls = TRUE`: a list with one element per gene, each a
  list of the null correlations (`nullX`, `nullY`) and the delta stars
  (`deltaStarX`, `deltaStarY`) of its permutations, in permutation
  order, when `x` and when `y` is permuted (empty for genes that were
  not tested).

[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (which
returns a plain data frame) have methods for it
([`STcompareResult`](https://jef.works/STcompare/reference/STcompareResult.md));
a selection of rows keeps the class and the attributes, and a selection
of columns is a plain data frame.

## Details

**Input.** The two samples must share their pixels: rasterize them
together, in one call such as
`SEraster::rasterizeGeneExpression(list(x, y), ...)`, so that a pixel
name means the same location in both. The pixels used are those whose
names (column names) are in both objects, in the column order of `x`,
and their coordinates
([`SpatialExperiment::spatialCoords()`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html))
must agree.

**Spatial correlation test.** The observed statistic is Pearson's
correlation `r` of the gene's values in `x` and `y` over the shared
pixels. Its null distribution comes from spatially autocorrelated
permutations (Viladomat et al. 2014): the values of `x` are shuffled
across the pixels, smoothed with a Gaussian kernel whose window holds a
share \\\delta\\ of the pixels, and rescaled with added noise so that
their variogram matches the variogram of `x`; the \\\delta\\ of `delta`
whose variogram matches best (delta star) is kept, the values of this
surrogate are replaced by the gene's own values in the surrogate's rank
order (`surrogate = "remap"`, the default; see *Rank-remapped
surrogates*), and its correlation with `y` is one null correlation. The
same is done with `y` permuted and correlated with `x`. If there are
more than 1000 shared pixels, the variograms use a random subsample of
1000 of them, drawn from `seed` with R's default random number generator
whatever [`RNGkind()`](https://rdrr.io/r/base/Random.html) is set to.
The smoother and the variogram are those of the other correlation
functions of the package
([`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)).

**Adaptive p-values.** The permutations of a gene stop as soon as its
p-value is clearly not small (Besag and Clifford 1991). Both directions
(`x` permuted, `y` permuted) run in lockstep, and a gene stops at the
first permutation \\L\\ at which either direction has `exceedances` null
correlations whose absolute value is at least \\\|r\|\\, or at
`nPermutations`. Its p-value (`pValuePermute`) is `exceedances / L` if
it stopped early, and otherwise \\(\max(b_x, b_y) + 1) / (L + 1)\\,
where \\b_x\\ and \\b_y\\ count the exceedances of each direction. This
is the larger of the two directions' sequential p-values, each \\h / L\\
when the direction stopped the gene (\\b = h\\) and \\(b + 1) / (L +
1)\\ at the limit, so it needs no assumption about how the two
directions are related: if the observed correlation is distributed like
the null correlations of either direction, \\P(p \le \alpha) \le
\alpha\\ at every \\\alpha\\. The surrogates approximate that null
distribution; with `surrogate = "gaussian"` the approximation fails for
genes detected in few pixels (see *Rarely detected genes*). A gene that
is clearly not significant stops after a few dozen permutations, and a
significant gene runs to `nPermutations`, which sets the smallest
possible p-value, `1 / (nPermutations + 1)`. With `exceedances = Inf`
every gene gets exactly `nPermutations` permutations. The columns
`pValuePermuteX` and `pValuePermuteY` always hold \\(b + 1) / (L + 1)\\,
so `pValuePermute` equals their maximum for a gene that reached
`nPermutations` and is slightly smaller, \\h / L\\ instead of \\(h + 1)
/ (L + 1)\\, for a gene that stopped early. `pValuePermute` is adjusted
across the tested genes with `adjustMethod`, which gives
`pValuePermuteAdj`.

**Which p-value should I report?** `pValuePermuteAdj`: the permutation
p-value of the gene, adjusted for the number of genes tested. A gene
with `pValuePermuteAdj < 0.05` has spatial expression patterns that are
correlated between the samples beyond what the spatial autocorrelation
of each pattern explains, at a false discovery rate of 5 percent (with
the default `adjustMethod = "BH"`). The unadjusted `pValuePermute` is
the per-gene p-value to report next to it, and its smallest possible
value is `1 / (nPermutations + 1)`. `pValuePermuteX` and
`pValuePermuteY` are descriptive: they are not adjusted, and when a gene
stops early the direction that did not reach `exceedances` has an
inflated value. `pValueNaive` ignores the spatial autocorrelation and is
far too small; it is reported only to show how much the permutation test
matters. The legacy
[`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
returns the same unadjusted `pValuePermuteX` and `pValuePermuteY` (for a
fixed number of permutations), plus each of them adjusted across genes,
`pValuePermuteXAdj` and `pValuePermuteYAdj`; all four are per direction,
and there is no combined p-value: the published analyses called a gene
significant when both adjusted p-values were below 0.05.

**Rank-remapped surrogates.** With `surrogate = "remap"` (the default),
the values of each surrogate are replaced by the gene's own values,
placed in the rank order of the surrogate (the amplitude adjustment of
the AAFT surrogates of Theiler et al. 1992, in a single step): the
surrogate keeps the spatial arrangement of the smoothed, rescaled
permutation and has exactly the distribution of values of the gene,
zeros included. Delta star is chosen before the remapping, and for a
gene with a spatial pattern the remapping is a monotone distortion of
the surrogate's amplitudes, which changes its variogram a little. For a
gene without spatial structure the ranks of its surrogates are close to
uniformly random permutations, so its remapped surrogates are plain
random permutations of its values and the test is the exact permutation
test, whatever the gene's distribution (exactly so when the other gene
has no spatial structure either). Many rearrangements of a gene detected
in few pixels give exactly the observed correlation; as in the exact
test, such ties count as exceedances (a null correlation within a
relative 1e-9 of \\\|r\|\\ counts, because
[`cor()`](https://rdrr.io/r/stats/cor.html) rounds each rearrangement
differently). In simulations on the kidney and brain grids
(`bench/calibration-results.md` in the source repository), remapped
surrogates kept \\P(p \le \alpha) \le \alpha\\ for independent genes
detected in 1 to 100 percent of the pixels, had the same power as
gaussian surrogates on correlated Gaussian fields, recovered the
published kidney and brain genes at least as well, and cost about 7
percent more per permutation. `surrogate = "gaussian"` uses the
surrogates as the smoothing and the added noise leave them. They are the
surrogates of the other correlation functions of the package and of the
published analyses, and are kept for comparability with them: their
values are close to normally distributed, which gives p-values that are
far too small for genes detected in few pixels (next paragraph) and
conservative ones for skewed genes with a spatial pattern.

**Rarely detected genes.** Gaussian surrogates are smoothed and mixed
with Gaussian noise, so their values are close to normally distributed.
A gene detected in few pixels has a few high values among many zeros,
and the correlation of two such genes is dominated by the few pixels
where both are detected: chance coincidences give extreme correlations
more often than gaussian surrogates do, so small p-values are too small,
ten times or more for genes detected in one or two pixels. With
`surrogate = "gaussian"`, genes detected in fewer than a share
`minDetected` of the shared pixels of either sample are therefore not
tested for correlation (`status = "skipped"`). A gene is detected at a
pixel where its value is above its lowest value over the shared pixels
(for counts and normalized counts: where it is not zero). With gaussian
surrogates the default asks for \\\sqrt{N}\\ of the \\N\\ shared pixels
(18 of 311, 47 of 2170): the departure from the surrogates grows roughly
as \\N / (k_x k_y)\\ for genes detected in \\k_x\\ and \\k_y\\ pixels,
and this removes the genes for which it is largest; sparse or
zero-inflated genes above it can still get p-values that are somewhat
too small, so with gaussian surrogates test the genes with a spatial
pattern in both samples, or raise `minDetected`. Remapped surrogates
have the gene's own values, and the problem disappears: in the
simulations their p-values were calibrated down to genes detected in a
single pixel, so with `surrogate = "remap"` the default `minDetected` is
0 and every gene that is not constant in a sample is tested. Such genes
cost little (a gene without spatial structure stops after about ten
permutations), and a gene detected in a single pixel cannot get a
p-value below \\1 / N\\: the exact test is discrete.

**Spatial similarity.** The same computation as
[`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md):
for each gene, pixels are kept when their value is above the gene's
`minQuantile` quantile in `x` or in `y`; zeros among them are replaced
by 1e-4; `percentSimilarity` is the proportion of kept pixels with
\\\|\log_2(y / x)\| \le\\ `foldChange`. It is `NA` when fewer than a
share `minPixels` of the pixels are kept. A fold change needs values
that are not negative (counts or normalized expression on a linear
scale): for a gene with negative values the similarity is `NA`, and a
warning lists such genes.

**Reproducibility.** The permutations and the noise of each gene and
direction come from their own random stream, which depends only on
`seed`, the gene's name, the direction and the permutation number. The
results therefore do not depend on `nThreads`, on the order of the
genes, or on which other genes are in the call (except
`pValuePermuteAdj`, which is adjusted over the tested genes), and the
global random number generator (`.Random.seed`) is left unchanged. The
permutations are drawn over the shared pixels in the column order of `x`
(the order of `y` does not matter), and each direction has its own
streams: reordering the pixels of `x`, or swapping `x` and `y`, gives
other permutations and p-values that differ within Monte Carlo error.
Results are deterministic on a given platform; the normal deviates use
the C library's [`log()`](https://rdrr.io/r/base/Log.html) and can
differ in the last bits between platforms.

**Run time.** It is proportional to the total number of permutations
(the sum of the `nPermutations` column), and the work of all genes is
shared by `nThreads` threads. With the defaults, genes that are clearly
not significant use a few dozen permutations each and significant genes
10000, so the run time is dominated by the significant genes. A smaller
`nPermutations`, such as 1000, was 6 to 7 times faster on the published
data and still resolves p-values down to 1 / 1001.

**Skipped and failed genes.** A gene that is constant in a sample (for
example never detected) or detected in fewer pixels than `minDetected`
asks for (by default only with `surrogate = "gaussian"`) is not tested
for correlation: it gets `status = "skipped"`, `NA` p-values and the
reason in `message` (its similarity is still computed), and one message
counts such genes. A gene that cannot be tested gets
`status = "failed"`, `NA` in the statistics that could not be computed,
and the reason in `message`; one warning lists the failed genes, and the
other genes are not affected. This happens for genes with missing or
infinite values on the shared pixels (every statistic is `NA`) and when
the permutations of a gene fail (for example a degenerate variogram
fit).

## References

Viladomat J, Mazumder R, McInturff A, McCauley DJ, Hastie T (2014).
Assessing the significance of global and local correlations under
spatial autocorrelation: a nonparametric approach. *Biometrics*
70(2):409-418.
[doi:10.1111/biom.12139](https://doi.org/10.1111/biom.12139)

Besag J, Clifford P (1991). Sequential Monte Carlo p-values.
*Biometrika* 78(2):301-304.
[doi:10.1093/biomet/78.2.301](https://doi.org/10.1093/biomet/78.2.301)

Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD (1992). Testing
for nonlinearity in time series: the method of surrogate data. *Physica
D* 58(1-4):77-94.
[doi:10.1016/0167-2789(92)90102-S](https://doi.org/10.1016/0167-2789%2892%2990102-S)

## See also

[`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
and
[`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md),
which reproduce the published analyses exactly. The tutorials:
[`vignette("getting-started-with-STcompare", package = "STcompare")`](https://jef.works/STcompare/articles/getting-started-with-STcompare.md),
[`vignette("how-STcompare-works", package = "STcompare")`](https://jef.works/STcompare/articles/how-STcompare-works.md)
and
[`vignette("parameters-performance-reproducibility", package = "STcompare")`](https://jef.works/STcompare/articles/parameters-performance-reproducibility.md).

## Examples

``` r
data(speKidney)
# rasterize the three samples together, onto one grid of pixels
rast <- SEraster::rasterizeGeneExpression(speKidney, assay_name = "counts",
                                          resolution = 0.2, fun = "mean",
                                          square = FALSE)

# A and B have opposite spatial patterns, A and C similar ones
ab <- compareSpatial(list(A = rast$A, B = rast$B), nPermutations = 999,
                     nThreads = 2)
#> compareSpatial: 1 gene on 273 shared pixels in 0.2 s (999 permutations in total, 2 threads)
ab
#> STcompareResult: 1 gene on 273 shared pixels; A vs B (assay pixelval)
#> Correlation: adaptive permutation p-values (stop at 10 exceedances, at most 999 permutations); 11 deltas from 0.01 to 0.9; rank-remapped surrogates
#>   pValuePermuteAdj < 0.05: 1 gene (0 with correlationCoef > 0, 1 with correlationCoef < 0)
#>   stopped early: 0; reached nPermutations: 1; skipped: 0; failed: 0
#> Similarity (|log2(y / x)| <= 1): median percentSimilarity 0.535 (0 genes with NA)
#> 
#>  gene nPixels correlationCoef pValueNaive pValuePermuteX pValuePermuteY
#>  Gene     273          -0.947   5.65e-136          0.001          0.001
#>  pValuePermute pValuePermuteAdj nPermutations  stop deltaStarMedianX
#>          0.001            0.001           999 limit              0.2
#>  deltaStarMedianY deltaGridEdge percentSimilarity percentDissimilarityX
#>               0.2         FALSE             0.535                 0.267
#>  percentDissimilarityY numPixelInThresh   t1   t2 status message
#>                  0.198              273 9.87 7.67     ok        
ac <- compareSpatial(list(A = rast$A, C = rast$C), nPermutations = 999,
                     nThreads = 2)
#> compareSpatial: 1 gene on 277 shared pixels in 0.2 s (999 permutations in total, 2 threads)
summary(ac)
#> compareSpatial: 1 gene on 277 shared pixels (A vs C); 1 ok, 0 skipped, 0 failed
#> Correlation, pValuePermuteAdj < 0.05: 1 gene (1 with correlationCoef > 0, 0 with correlationCoef < 0)
#>   stopped early: 0; reached nPermutations: 1; skipped: 0; failed: 0
#>   nPermutations per gene: median 999, max 999; 999 in total
#>   settings: nPermutations = 999, exceedances = 10, 11 deltas (0.01 to 0.9), rank-remapped surrogates, BH adjustment, seed 0
#> Similarity, percentSimilarity:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>       0       0       0       0       0       0 
#> Run time: 0.2 s
as.data.frame(ac)[, c("gene", "correlationCoef", "pValuePermute",
                      "pValuePermuteAdj", "nPermutations", "percentSimilarity")]
#>      gene correlationCoef pValuePermute pValuePermuteAdj nPermutations
#> Gene Gene       0.9431531         0.001            0.001           999
#>      percentSimilarity
#> Gene                 0
```
