# Spatial correlation test with more permutations for promising genes

Function to calculate Pearson's correlation between assays from two
SpatialExperiment datasets using an iterative permutation strategy. It
first evaluates all genes with a smaller number of permutations, then
only reruns genes that remain potentially significant with progressively
larger permutation counts to refine their empirical p-values while
accounting for original degree of autocorrelation.

## Usage

``` r
spatialCorrelationGeneExpIterPermutations(
  input,
  alpha = 0.05,
  nPermutations = c(100, 1000),
  deltaX = NULL,
  deltaY = NULL,
  maxDistPrctile = 0.25,
  returnPermutations = FALSE,
  assayName = NULL,
  nThreads = 1,
  BPPARAM = NULL,
  verbose = TRUE,
  seed = 0,
  adjustMethod = "BH"
)
```

## Arguments

- input:

  `list` List of two SpatialExperiment objects with matched spatial
  locations. The first element corresponds to the first
  SpatialExperiment (`X`), and the second to the second
  SpatialExperiment (`Y`). The SpatialCoords of the two
  SpatialExperiment objects should be on the same coordinate framework
  and observations at the same coordinate location in both datasets
  should have the same row names. If the SpatialExperiment objects do
  not have shared locations, use
  [`SEraster::rasterizeGeneExpression()`](https://rdrr.io/pkg/SEraster/man/rasterizeGeneExpression.html)
  to generate SpatialExperiment objects with shared pixel locations. See
  `assayName` parameter if the SpatialExperiment objects have more than
  one assay.

- alpha:

  `numeric`: significance threshold used to decide which genes should be
  rerun at the next permutation level. After each iteration, a gene is
  carried forward only when both its unadjusted (raw) p-values
  `pValuePermuteX` and `pValuePermuteY` are less than
  `100 * alpha / nPermutations[k]` (the screening never uses the
  adjusted p-values, which are computed once, at the end), that is, less
  than `alpha` after 100 permutations and less than `alpha / 10` after
  1000 permutations. When `100 * alpha` is a whole number, as for the
  default, these are the genes with fewer than `100 * alpha` exceedances
  (null correlations at least as large as the observed correlation in
  absolute value) in each direction. Otherwise a gene with exactly
  `floor(100 * alpha)` exceedances may not be carried forward: for
  example, with `alpha = 0.053` and 100 permutations, a gene can have at
  most 4 exceedances, not 5. Default is `0.05`.

- nPermutations:

  `numeric vector`: numbers of permutations to use across iterative
  rounds. The vector is applied from smallest to largest. All genes are
  first tested with `nPermutations[1]`; genes passing the screening rule
  defined by `alpha` are rerun with `nPermutations[2]`, and so on. A
  gene's smallest possible unadjusted p-value is \\1 / (B + 1)\\, where
  \\B\\ is the number of permutations in the last round in which it was
  tested. Default is `c(100, 1000)`.

- deltaX:

  `list`: List of single numerics or list of numeric vectors to use for
  delta, the parameter controlling the degree of smoothing in
  permutations of X. The length of the list should be the same as the
  number of rows in the SpatialExperiment. Delta is a proportion
  calculated by dividing k neighbors by N total observations (columns)
  in X, where k is the number of neighbors in the permutation of X that
  should be within the radius smoothed by the Gaussian kernel to achieve
  the amount of autocorrelation present in the original X. If a single
  delta is not known, a sequence of deltas can be inputted and the best
  delta will be found such that it minimizes the sum of squares of the
  residuals between the variogram of the permutation generated from the
  delta and the variogram of the target. Default is `NULL`. If no value
  is supplied for `deltaX`, `seq(0.1,0.9,0.1)`, the sequence of every
  0.1 from 0.1 to 0.9, will be used to find the best delta for each row
  (gene) in X.

- deltaY:

  `list`: List of single numerics or list of numeric vectors to use for
  delta, the parameter controlling the degree of smoothing in
  permutations of Y. `deltaY` is like `deltaX` but for permuting data in
  Y instead of X. Default is `NULL`. If no value is supplied for
  `deltaY`, `seq(0.1,0.9,0.1)`, the sequence of every 0.1 from 0.1 to
  0.9, will be used to find the best delta for permutations for each row
  (gene) in Y.

- maxDistPrctile:

  `numeric`: percentile of distances between pixels to use as max
  distance when calculating variograms. Default = 0.25. At greater
  distances the variogram is less precise because there are fewer pairs
  of points with that distance between them. Therefore, since the goal
  is to minimize the difference between the variogram of X and those of
  its permutations, the variogram should be subsetted to the percentile
  that is more robust.

- returnPermutations:

  `logical`: indicate whether the dataframe returned as output will have
  columns with the values of the permutations used to calculate the null
  correlations and empirical p-values. Default is `FALSE`.

- assayName:

  `character` or `integer` A character string or numeric specifying the
  assay in the SpatialExperiment to use. Default is `NULL`. If no value
  is supplied for `assayName`, then the first assay is used as a
  default.

- nThreads:

  `integer`: Number of threads of the compiled code. Default = 1. The
  permutations of all genes are distributed over the threads, and the
  results do not depend on the number of threads. We recommend the
  number of cores available (`parallel::detectCores(logical = FALSE)`).

- BPPARAM:

  `BiocParallelParam`: Optional. If not `NULL`, its number of workers
  (`BiocParallel::bpnworkers(BPPARAM)`) is used as the number of threads
  instead of `nThreads`. No BiocParallel back-end is started (nothing is
  forked). Default is `NULL`.

- verbose:

  `logical`: if `TRUE` (default), print a message when the computation
  starts (genes, permutation rounds and threads) and when it ends
  (elapsed time and the number of genes tested in each round).

- seed:

  `integer`: Seed for the random number generator. Default `0`. The
  permutations are drawn after `set.seed(seed)` with the session's
  random number generator kinds
  ([`RNGkind()`](https://rdrr.io/r/base/Random.html)), and the noise of
  permutation \\b\\ after `set.seed(seed + b)` with `"L'Ecuyer-CMRG"`;
  every gene uses the same seed. The results depend on the seed only,
  and the global random number generator state (`.Random.seed`) is left
  unchanged.

- adjustMethod:

  `character`: multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), applied
  once, after the last round, across all genes to the final
  `pValuePermuteX` and `pValuePermuteY` columns separately; it gives the
  columns `pValuePermuteXAdj` and `pValuePermuteYAdj`, and the
  unadjusted columns are kept. Must be one of `p.adjust.methods`.
  Default is `"BH"`.

## Value

The output is returned as a `data.frame`. The rownames are the rownames
of the SpatialExperiments, and each row reflects the last permutation
round in which that gene was evaluated. The columns and their contents
are as follows:

- `correlationCoef`:

  Pearson's correlation coefficient.

- `pValueNaive`:

  the analytical p-value naively assuming independent observations

- `pValuePermuteX`:

  the empirical p-value from the null generated by permutations of
  observations in X, \\(b + 1) / (B + 1)\\ with the gene's final number
  of permutations \\B\\ (see
  [`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md)).
  Not adjusted for multiple testing

- `pValuePermuteY`:

  the empirical p-value from the null generated by permutations of
  observations in Y, computed like `pValuePermuteX`. Not adjusted

- `pValuePermuteXAdj`:

  `pValuePermuteX` adjusted across all genes with `adjustMethod`
  ([`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); genes
  with an `NA` p-value are not counted), after the last round

- `pValuePermuteYAdj`:

  `pValuePermuteY` adjusted across all genes with `adjustMethod`,
  separately from `pValuePermuteXAdj`. The published analyses called a
  gene significant when both adjusted p-values were below 0.05

- `deltaStarMedianX`:

  the median delta star (the delta which minimizes the difference
  between the variogram of the permutation and the variogram of
  observations) across permutations of X

- `deltaStarMedianY`:

  the median delta star across permutations of Y

- `deltaStarX`:

  list of delta star for all permutations of X

- `deltaStarY`:

  list of delta star for all permutations of Y

- `nullCorrelationsX`:

  list of B x 1 matrices: the correlation coefficients for Y and all
  permutations of X

- `nullCorrelationsY`:

  list of B x 1 matrices: the correlation coefficients for X and all
  permutations of Y

- `permutationsX`:

  (optional) an N x B matrix, where N is the length of X and B is the
  final `nPermutations` used for that gene. Each column is the resulting
  values of a permutation of X

- `permutationsY`:

  (optional) an N x B matrix, where N is the length of Y and B is the
  final `nPermutations` used for that gene. Each column is the resulting
  values of a permutation of Y

## Details

Every round tests the genes as
[`spatialCorrelationGeneExp`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
does, in one call of the compiled code on `nThreads` threads. A gene
carried forward from a round with \\B_k\\ permutations to a round with
\\B\_{k+1}\\ keeps its first \\B_k\\ permutations and adds permutations
\\B_k + 1, ..., B\_{k+1}\\. Permutation \\b\\ depends only on the seed
and on \\b\\, so the result is the same as a fresh run with \\B\_{k+1}\\
permutations, at a fraction of the cost.

Genes for which a null cannot be computed (see
[`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md))
get `NA` in every column computed from permutations, with one warning
per gene, and are not carried forward.

## Examples

``` r
data(speKidney)
rastKidney <- SEraster::rasterizeGeneExpression(
  speKidney,
  assay_name = "counts",
  resolution = 0.2,
  fun = "mean",
  square = FALSE
)

rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)

corr <- spatialCorrelationGeneExpIterPermutations(
  rastGexpListAB,
  nPermutations = c(100, 1000),
  nThreads = 2
)
#> spatialCorrelationGeneExpIterPermutations: 1 gene(s) x 2 directions on 273 shared pixels, permutation rounds 100, 1000, 2 thread(s)
#> spatialCorrelationGeneExpIterPermutations: done in 0.1 s; genes tested per round: 1 at B = 100, 1 at B = 1000 (0 gene(s) with NA permutation p-values)

corr[, 1:6]
#>      correlationCoef   pValueNaive pValuePermuteX pValuePermuteY
#> Gene      -0.9472813 5.652003e-136    0.000999001    0.000999001
#>      pValuePermuteXAdj pValuePermuteYAdj
#> Gene       0.000999001       0.000999001
# the number of permutations of each gene
lengths(corr$nullCorrelationsX)
#> [1] 1000
```
