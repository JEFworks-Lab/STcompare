# Spatial correlation test for every gene of two samples

Function to calculate Pearson's correlation between assays from two
SpatialExperiment datasets. To replace the analytical p-value which
results in a high false positive rate for autocorrelated spatial
patterns, it calculates empirical p-values from empirical null
distributions generated from permuting the datasets and then smoothing
to maintain the original degree of autocorrelation

## Usage

``` r
spatialCorrelationGeneExp(
  input,
  nPermutations = 100,
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

- nPermutations:

  `integer` or `double`: number of permutations to generate to build the
  empirical null distribution. This number will determine the precision
  of the p-value. Default is `100`, such that the smallest possible
  unadjusted p-value is \\1 / 101\\, about 0.0099

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
  a column with the values of the permutations used to calculate the
  null correlations and the empirical p-value. Default is `FALSE`

- assayName:

  `character` or `integer` A character string or numeric specifying the
  assay in the SpatialExperiment to use. Default is `NULL`. If no value
  is supplied for `assayName`, then the first assay is used as a default

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
  starts (genes, permutations and threads) and when it ends (elapsed
  time).

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
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). It is
  applied across all genes, separately to the `pValuePermuteX` column
  and to the `pValuePermuteY` column, and gives the columns
  `pValuePermuteXAdj` and `pValuePermuteYAdj`; the unadjusted columns
  are kept. Must be one of `p.adjust.methods`; with `"none"` the
  adjusted columns equal the unadjusted ones. Default is `"BH"`.

## Value

The output is returned as a `data.frame`. The rownames are the rownames
of the SpatialExperiments. The names of the columns and their contents
are as follows:

- `correlationCoef`:

  Pearson's correlation coefficient.

- `pValueNaive`:

  the analytical p-value naively assuming independent observations

- `pValuePermuteX`:

  the empirical p-value when creating a null from permutations of
  observations in X, \\(b + 1) / (B + 1)\\ (see
  [`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md)).
  Not adjusted for multiple testing

- `pValuePermuteY`:

  the empirical p-value when creating a null from permutations of
  observations in Y, computed like `pValuePermuteX`. Not adjusted

- `pValuePermuteXAdj`:

  `pValuePermuteX` adjusted across genes with `adjustMethod`
  ([`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); genes
  with an `NA` p-value are not counted)

- `pValuePermuteYAdj`:

  `pValuePermuteY` adjusted across genes with `adjustMethod`, separately
  from `pValuePermuteXAdj`. The published analyses called a gene
  significant when both adjusted p-values were below 0.05

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

  (optional) an N x B matrix, where N is the length of X and B is
  `nPermutations`. Each column is the resulting values of a permutation
  of X

- `permutationsY`:

  (optional) an N x B matrix, where N is the length of Y and B is
  `nPermutations`. Each column is the resulting values of a permutation
  of Y

## Details

For every gene (row), the correlation of its values in the two datasets
over their shared pixels is tested as in
[`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md):
with nulls from permutations of X and from permutations of Y. All genes
are computed in one call of the compiled code, which shares the
smoothing operators and the shuffles of the locations between genes and
distributes the permutations of all genes over `nThreads` threads. The
results are the same as for one gene at a time and do not depend on the
number of threads.

Genes for which a null cannot be computed (see
[`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md):
for example a gene that is constant or has missing values on the shared
pixels, or \\N \times delta \< 2\\) get `NA` in every column computed
from permutations, and one warning per such gene gives the reason.

## Examples

``` r
data(speKidney)

##### Rasterize to get pixels at matched spatial locations #####
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
               assay_name = 'counts', resolution = 0.2, fun = "mean",
               square = FALSE)

##### Use STcompare to calculate Pearson's correlation coefficient #####
rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
rastGexpListAC <- list(A = rastKidney$A, C = rastKidney$C)

negCorrelation <- spatialCorrelationGeneExp(rastGexpListAB, nThreads = 2)
#> spatialCorrelationGeneExp: 1 gene(s) x 2 directions on 273 shared pixels, 100 permutations, 2 thread(s)
#> spatialCorrelationGeneExp: done in 0.0 s (0 of 1 gene(s) with NA permutation p-values)
posCorrelation <- spatialCorrelationGeneExp(rastGexpListAC, nThreads = 2)
#> spatialCorrelationGeneExp: 1 gene(s) x 2 directions on 277 shared pixels, 100 permutations, 2 thread(s)
#> spatialCorrelationGeneExp: done in 0.0 s (0 of 1 gene(s) with NA permutation p-values)

negCorrelation
#>      correlationCoef   pValueNaive pValuePermuteX pValuePermuteY
#> Gene      -0.9472813 5.652003e-136     0.00990099     0.00990099
#>      pValuePermuteXAdj pValuePermuteYAdj deltaStarMedianX deltaStarMedianY
#> Gene        0.00990099        0.00990099              0.2              0.2
#>        deltaStarX   deltaStarY nullCorrelationsX nullCorrelationsY
#> Gene 0.2, 0.2.... 0.2, 0.2....      -0.16074....      -0.08458....
posCorrelation
#>      correlationCoef   pValueNaive pValuePermuteX pValuePermuteY
#> Gene       0.9431531 1.409195e-133     0.00990099     0.00990099
#>      pValuePermuteXAdj pValuePermuteYAdj deltaStarMedianX deltaStarMedianY
#> Gene        0.00990099        0.00990099              0.2              0.2
#>        deltaStarX   deltaStarY nullCorrelationsX nullCorrelationsY
#> Gene 0.2, 0.4.... 0.2, 0.2....      -0.00632....      0.052321....
```
