# Spatial correlation test of two vectors of values at the same locations

Function to calculate Pearson's correlation between two spatial
datasets. To replace the analytical p-value which results in a high
false positive rate for autocorrelated spatial patterns, it calculates
empirical p-values from empirical null distributions generated from
permuting the datasets and then smoothing to maintain the original
degree of autocorrelation

## Usage

``` r
spatialCorrelation(
  X,
  Y,
  pos,
  nPermutations = 100,
  deltaX = NULL,
  deltaY = NULL,
  maxDistPrctile = 0.25,
  returnPermutations = FALSE,
  nThreads = 1,
  BPPARAM = NULL,
  seed = 0
)
```

## Arguments

- X:

  `numeric`: a numeric vector with N observations (a matrix with one row
  or one column is also accepted)

- Y:

  `numeric`: a numeric vector with N observations, like `X`

- pos:

  `matrix`: an N x 2 matrix of the spatial x,y coordinates of
  observations

- nPermutations:

  `integer` or `double`: number of permutations to generate to build the
  empirical null distribution. This number will determine the precision
  of the p-value. Default is `100`, such that the smallest possible
  p-value is \\1 / 101\\, about 0.0099

- deltaX:

  `numeric`: A single numeric or a numeric vector for controlling the
  degree of smoothing in permutations of X. Delta is a proportion
  calculated by dividing k neighbors by N total observations in X, where
  k is the number of neighbors in the permutation of X that should be
  within the radius smoothed by the Gaussian kernel to achieve the
  amount of autocorrelation present in the original X. If a single delta
  is not known, a sequence of deltas can be inputted and the best delta
  will be found such that it minimizes the sum of squares of the
  residuals between the variogram of the permutation generated from the
  delta and the variogram of the target. Default is `NULL`. If no value
  is supplied for `deltaX`, `seq(0.1,0.9,0.1)`, the sequence of every
  0.1 from 0.1 to 0.9, will be used to find the best delta for X.

- deltaY:

  `numeric`: A single numeric or numeric vector for controlling the
  degree of smoothing in permutations of Y. `deltaY` is like `deltaX`
  but for observation in Y instead of X. Default is `NULL`. If no value
  is supplied for `deltaY`, `seq(0.1,0.9,0.1)`, the sequence of every
  0.1 from 0.1 to 0.9, will be used to find the best delta for
  permutations of Y.

- maxDistPrctile:

  `numeric`: percentile of distances between pixels to use as max
  distance when calculating variograms. Default = 0.25. At greater
  distances the variogram is less precise because there are fewer pairs
  of points with that distance between them. Therefore, since the goal
  is to minimize the difference between the variogram of X and those of
  its permutations, the variogram should be subsetted to the percentile
  that is more robust.

- returnPermutations:

  `logical`: `FALSE` (default) indicate whether the outputted dataframe
  will have a column with the values of the permutations used to
  calculate the null correlations and the empirical p-value.

- nThreads:

  `integer`: Number of threads of the compiled code. Default = 1. The
  permutations are distributed over the threads, and the results do not
  depend on the number of threads. We recommend the number of cores
  available (`parallel::detectCores(logical = FALSE)`).

- BPPARAM:

  `BiocParallelParam`: Optional. If not `NULL`, its number of workers
  (`BiocParallel::bpnworkers(BPPARAM)`) is used as the number of threads
  instead of `nThreads`. No BiocParallel back-end is started (nothing is
  forked). Default is `NULL`.

- seed:

  `integer`: Seed for the random number generator. Default `0`. The
  permutations are drawn after `set.seed(seed)` with the session's
  random number generator kinds
  ([`RNGkind()`](https://rdrr.io/r/base/Random.html)), and the noise of
  permutation \\b\\ after `set.seed(seed + b)` with `"L'Ecuyer-CMRG"`.
  The results depend on the seed only, and the global random number
  generator state (`.Random.seed`) is left unchanged.

## Value

The output is returned as a `data.frame` with one row (named `"cor"`)
containing the columns:

- `correlationCoef`:

  Pearson's correlation coefficient.

- `pValueNaive`:

  the analytical p-value naively assuming independent observations

- `pValuePermuteX`:

  the empirical p-value when creating a null from permutations of
  observations in X, computed as \\(b + 1) / (B + 1)\\ where \\b\\ is
  the number of null correlations whose absolute value is at least the
  absolute value of the observed correlation

- `pValuePermuteY`:

  the empirical p-value when creating a null from permutations of
  observations in Y, computed like `pValuePermuteX`

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

  list of a B x 1 matrix: the correlation coefficients for pairing Y and
  all permutations of X

- `nullCorrelationsY`:

  list of a B x 1 matrix: the correlation coefficients for pairing X and
  all permutations of Y

- `permutationsX`:

  (optional) an N x B matrix, where N is the length of X and B is
  `nPermutations`. Each column is the resulting values of a permutation
  of X

- `permutationsY`:

  (optional) an N x B matrix, where N is the length of Y and B is
  `nPermutations`. Each column is the resulting values of a permutation
  of Y

## Details

The null distribution of the correlation is built twice: by permuting X
while keeping Y fixed, and by permuting Y while keeping X fixed (see
[`viladomatCorrelation`](https://jef.works/STcompare/reference/viladomatCorrelation.md)
for how a permutation keeps the autocorrelation). Both directions use
the same seed, so they use the same shuffles of the locations. The
empirical p-value of a direction is \\(b + 1) / (B + 1)\\, where \\b\\
is the number of null correlations whose absolute value is at least the
absolute value of the observed correlation and \\B\\ is `nPermutations`;
it is never 0.

The computation is done by compiled code (C++) that reimplements the
`locfit` smoother and the `geoR` variogram of the original R
implementation and gives the same results: the same delta star for every
permutation, and the same null correlations up to floating-point
rounding (about 1e-14 on the published analyses). Results do not depend
on the number of threads, and the global random number generator state
is left unchanged.

If the null cannot be computed in a direction, the columns computed from
permutations are `NA` (an NA row) and a warning gives the reason. This
happens when X or Y is constant or has missing values, when \\N \times
delta \< 2\\ for a delta (fewer than 2 locations in the smoothing
neighbourhood), when `maxDistPrctile` is so small that the variogram has
fewer than 2 bins, or with exactly duplicated coordinates and a delta so
small that the neighbourhood holds no more than the copies of a
location. `correlationCoef` and `pValueNaive` are still computed by
[`cor.test()`](https://rdrr.io/r/stats/cor.test.html) where possible (on
the complete pairs if there are missing values).

## Examples

``` r
data(quakes)

#remove duplicated positions
quakes_data <- quakes[!duplicated(cbind(quakes$lat, quakes$long)),]

cor <- spatialCorrelation(X = quakes_data$depth,
                          Y = quakes_data$mag,
                          pos = cbind(quakes_data$lat, quakes_data$long),
                          nThreads = 2)
cor
#>     correlationCoef  pValueNaive pValuePermuteX pValuePermuteY deltaStarMedianX
#> cor       -0.230266 1.779041e-13     0.03960396      0.1386139              0.8
#>     deltaStarMedianY   deltaStarX   deltaStarY nullCorrelationsX
#> cor              0.5 0.9, 0.8.... 0.3, 0.7....      -0.15504....
#>     nullCorrelationsY
#> cor      -0.10144....

# plot the delta star (the delta which minimizes the difference between the
# variogram of the permutation and the variogram of observations) for all
# permutations to see if clear peak found in the range inputted
hist(cor$deltaStarX[[1]])

hist(cor$deltaStarY[[1]])


#plot null distribution of correlations to see if normally distributed
hist(cor$nullCorrelationsX[[1]])

hist(cor$nullCorrelationsY[[1]])


#example of inputting specific range for deltaX and deltaY
cor2 <- spatialCorrelation(X = quakes_data$depth,
                          Y = quakes_data$mag,
                          pos = cbind(quakes_data$lat, quakes_data$long),
                          deltaX = seq(0.05, 0.9, 0.05),
                          deltaY = seq(0.02, 0.5, 0.02),
                          nThreads = 2)

cor2
#>     correlationCoef  pValueNaive pValuePermuteX pValuePermuteY deltaStarMedianX
#> cor       -0.230266 1.779041e-13     0.03960396     0.02970297             0.75
#>     deltaStarMedianY   deltaStarX   deltaStarY nullCorrelationsX
#> cor             0.26 0.9, 0.8.... 0.5, 0.1....      -0.15141....
#>     nullCorrelationsY
#> cor      -0.10726....

hist(cor2$deltaStarX[[1]])

hist(cor2$deltaStarY[[1]])

hist(cor2$nullCorrelationsX[[1]])

hist(cor2$nullCorrelationsY[[1]])


#visualizations of the spatial data to verify negative correlation
library(ggplot2)
p1 <- ggplot2::ggplot(quakes_data,
                      ggplot2::aes(x = long, y = lat, color = depth)) +
ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_gradient(low = "lightblue", high = "blue") +
  ggplot2::labs(title = "Locations of Earthquakes off Fiji",
                x = "Longitude", y = "Latitude",
                color = "Depth (km)") +
  ggplot2::theme_minimal() +
  ggplot2::coord_quickmap()

p2 <- ggplot2::ggplot(quakes_data,
                      ggplot2::aes(x = long, y = lat, color = mag)) +
ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_gradient(low = "lightblue", high = "blue") +
  ggplot2::labs(title = "Locations of Earthquakes off Fiji",
                x = "Longitude", y = "Latitude",
                color = "Richter Magnitude") +
  ggplot2::theme_minimal() +
  ggplot2::coord_quickmap()

p1

p2
```
