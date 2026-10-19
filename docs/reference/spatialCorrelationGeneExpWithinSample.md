# Spatial correlation test between pairs of genes of one sample

Function to calculate Pearson's correlation between rows from one
SpatialExperiment dataset. To replace the analytical p-value which
results in a high false positive rate for autocorrelated spatial
patterns, it calculates empirical p-values from empirical null
distributions generated from permuting the data and then smoothing to
maintain the original degree of autocorrelation

## Usage

``` r
spatialCorrelationGeneExpWithinSample(
  input,
  nPermutations = 100,
  delta = NULL,
  maxDistPrctile = 0.25,
  returnPermutations = FALSE,
  assayName = NULL,
  nThreads = 1,
  BPPARAM = NULL,
  verbose = TRUE,
  seed = 0
)
```

## Arguments

- input:

  `SpatialExperiment` A SpatialExperiment object. See `assayName`
  parameter if the SpatialExperiment object has more than one assay.

- nPermutations:

  `integer` or `double`: number of permutations to generate to build the
  empirical null distribution. This number will determine the precision
  of the p-value. Default is `100`, such that the smallest possible
  p-value is \\1 / 101\\, about 0.0099

- delta:

  `list`: List of single numerics or list of numeric vectors to use for
  delta, the parameter controlling the degree of smoothing in
  permutations of each row (gene). The length of the list should be the
  same as the number of rows in the SpatialExperiment. Delta is a
  proportion calculated by dividing k neighbors by N total observations
  (columns), where k is the number of neighbors in the permutation that
  should be within the radius smoothed by the Gaussian kernel to achieve
  the amount of autocorrelation present in the original data. If a
  single delta is not known, a sequence of deltas can be inputted and
  the best delta will be found such that it minimizes the sum of squares
  of the residuals between the variogram of the permutation generated
  from the delta and the variogram of the target. Default is `NULL`. If
  no value is supplied for `delta`, `seq(0.1,0.9,0.1)`, the sequence of
  every 0.1 from 0.1 to 0.9, will be used to find the best delta for
  each row (gene).

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
  starts (genes, pairs, permutations and threads) and when it ends
  (elapsed time).

- seed:

  `integer`: Seed for the random number generator. Default `0`. The
  permutations are drawn after `set.seed(seed)` with the session's
  random number generator kinds
  ([`RNGkind()`](https://rdrr.io/r/base/Random.html)), and the noise of
  permutation \\b\\ after `set.seed(seed + b)` with `"L'Ecuyer-CMRG"`;
  every gene uses the same seed. The results depend on the seed only,
  and the global random number generator state (`.Random.seed`) is left
  unchanged.

## Value

The output is returned as a `data.frame` with one row per pair of rows
(genes), in the order of `combn(rownames(input), 2)`. The rownames are
arbitrary (`"cor"`, `"cor1"`, ...). The names of the columns and their
contents are as follows:

- `correlationCoef`:

  Pearson's correlation coefficient.

- `pValueNaive`:

  the analytical p-value naively assuming independent observations

- `pValuePermuteX`:

  the p-value when creating an empirical null from permutations of the
  first gene of the pair, computed as \\(b + 1) / (B + 1)\\ (see
  [`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md))

- `pValuePermuteY`:

  the p-value when creating an empirical null from permutations of the
  second gene of the pair

- `deltaStarMedianX`:

  the median delta star (the delta which minimizes the difference
  between the variogram of the permutation and the variogram of
  observations) across permutations of the first gene

- `deltaStarMedianY`:

  the median delta star across permutations of the second gene

- `deltaStarX`:

  list of delta star for all permutations of the first gene

- `deltaStarY`:

  list of delta star for all permutations of the second gene

- `nullCorrelationsX`:

  list of B x 1 matrices: the correlation coefficients for the second
  gene and all permutations of the first

- `nullCorrelationsY`:

  list of B x 1 matrices: the correlation coefficients for the first
  gene and all permutations of the second

- `permutationsX`:

  (optional) an N x B matrix, where N is the number of pixels and B is
  `nPermutations`. Each column is the resulting values of a permutation
  of the first gene

- `permutationsY`:

  (optional) the same for the second gene

- `first`:

  the name of the first row in the pair

- `second`:

  the name of the second row in the pair

## Details

Every pair of rows (genes) is tested as in
[`spatialCorrelation`](https://jef.works/STcompare/reference/spatialCorrelation.md),
with the first gene of the pair as X and the second as Y. A gene's
permutations do not depend on its partner, so the compiled code computes
the permutations of every gene once and correlates them with every other
gene, on `nThreads` threads. The results are the same as for one pair at
a time and do not depend on the number of threads. Pairs with a gene for
which a null cannot be computed (for example a constant gene) get `NA`
in every column computed from permutations, and one warning per such
pair gives the reason. The empirical p-values are not adjusted for
multiple testing.

## Examples

``` r
data(speKidney)

##### Rasterize to get pixels at matched spatial locations #####
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
               assay_name = 'counts', resolution = 0.2, fun = "mean",
               square = FALSE)

# one SpatialExperiment whose three rows are the gene in samples A, C and
# B, on the pixels the three samples share
shared <- Reduce(intersect, lapply(rastKidney, colnames))
expr <- t(sapply(rastKidney, function(s) {
  SummarizedExperiment::assay(s)[1, shared]
}))
speACB <- SpatialExperiment::SpatialExperiment(
  assays = list(pixelval = expr),
  spatialCoords = SpatialExperiment::spatialCoords(rastKidney$A)[shared, ])

sc_within_sample <- spatialCorrelationGeneExpWithinSample(speACB,
                                                          nThreads = 2)
#> spatialCorrelationGeneExpWithinSample: 3 genes (3 pairs) on 270 pixels, 100 permutations, 2 thread(s)
#> spatialCorrelationGeneExpWithinSample: done in 0.0 s (0 of 3 pairs with NA permutation p-values)
sc_within_sample[, c("first", "second", "correlationCoef",
                     "pValuePermuteX", "pValuePermuteY")]
#>      first second correlationCoef pValuePermuteX pValuePermuteY
#> cor      A      C       0.9457538     0.00990099     0.00990099
#> cor1     A      B      -0.9467256     0.00990099     0.00990099
#> cor2     C      B      -0.9475692     0.00990099     0.00990099
```
