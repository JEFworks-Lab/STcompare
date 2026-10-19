# Spatially autocorrelated surrogates and permutation p-value for one direction

Function to calculate Pearson's correlation between two spatial
datasets, X and Y. To replace the analytical p-value which results in a
high false positive rate for autocorrelated spatial patterns, it
calculates empirical p-values from empirical null distributions
generated from permuting dataset X by randomly shuffling the values and
then smoothing to maintain the original degree of autocorrelation of X

## Usage

``` r
viladomatCorrelation(
  data,
  delta,
  maxDistPrctile,
  nPermutations,
  nThreads = 1,
  BPPARAM = NULL,
  seed = 0
)
```

## Arguments

- data:

  `matrix` An N x 4 matrix with the first column as the values of X, the
  second column as the values of Y, the third column as the
  x-coordinates, and the fourth column as the y-coordinates.

- delta:

  `numeric vector` Given a point in dataset X, the percentage of
  neighbors which should be within the smoothing kernel. This can be a
  single numeric or a vector of candidate values; the best one is chosen
  for every permutation. Each delta must satisfy \\N \times delta \ge
  2\\.

- maxDistPrctile:

  `numeric`: Percentile of distances between pixels to use as max
  distance when calculating variograms. At greater distances the
  variogram is less precise because there are fewer pairs of points with
  that distance between them. Therefore, since the goal is to minimize
  the difference between the variogram of X and those of its
  permutations, the variogram should be subsetted to the percentile that
  is more robust.

- nPermutations:

  `integer`: Number of permutations to generate to build the empirical
  null distribution. This number will determine the precision of the
  p-value. The smallest possible p-value is \\1 / (nPermutations + 1)\\,
  for example about 0.0099 when `nPermutations <- 100`

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

The output is returned as a list.

- `deltaStarMedian`:

  numeric, the median of the deltas that minimize the residual sum of
  squares across each permutation

- `deltaStar`:

  numeric vector of length `nPermutations`: the delta that minimizes the
  residual sum of squares for each permutation

- `pValueGlobal`:

  numeric, empirical p-value for the Pearson's correlation of X and Y,
  computed as \\(b + 1) / (B + 1)\\ where \\b\\ is the number of null
  correlations whose absolute value is at least the absolute value of
  the observed correlation

- `nullCorGlobal`:

  a B x 1 matrix, where B is `nPermutations`. This matrix is the
  correlation coefficients between the permutations and Y that compose
  that null distribution used to calculate the empirical p-value

- `permutations`:

  an N x B matrix, where B is `nPermutations`. Each column is the
  resulting values of a permutation of X

## Details

Each permutation shuffles the values of X across the locations, smooths
the shuffled values with a Gaussian kernel whose nearest-neighbour
bandwidth covers a proportion delta of the locations, and rescales the
smoothed values and adds Gaussian noise so that their variogram matches
the variogram of X (Viladomat et al. 2014). This is repeated for every
delta, and the delta whose variogram matches best (delta star) gives the
permuted field. If X has more than 1000 values, the variograms use a
random subsample of 1000 locations.

The computation is done by compiled code (C++), which reimplements the
local-constant kernel smoother of `locfit` (its adaptive k-d tree with
vertex interpolation, `nn = delta`, `kern = "gauss"`) and the binned
variogram of
[`geoR::variog()`](https://rdrr.io/pkg/geoR/man/variog.html). It gives
the same results as the original R implementation: the same delta star
for every permutation and the same null correlations up to
floating-point rounding (about 1e-14 on the published analyses).

If the null cannot be computed (for example X is constant or has missing
values, \\N \times delta \< 2\\ for a delta, or the variogram has fewer
than 2 bins), every element of the result is `NA` and a warning gives
the reason. If Y has missing values or is constant, the permutations are
returned, and `nullCorGlobal` and `pValueGlobal` are `NA` with a
warning.

## References

adapted from: Viladomat, Júlia et al. “Assessing the significance of
global and local correlations under spatial autocorrelation: a
nonparametric approach.” Biometrics vol. 70,2 (2014): 409-18.
doi:10.1111/biom.12139

## Examples

``` r
data(quakes)

#remove duplicated positions
quakes_data <- quakes[!duplicated(cbind(quakes$lat, quakes$long)),]

data <- cbind(quakes_data$depth,
              quakes_data$mag,
              quakes_data$lat,
              quakes_data$long)

#sequence of deltas
delta <- seq(0.1,0.9,0.1)

# maximum distance for the variogram set at the 25% percentile of
# the distribution of pairs of distances:
maxDistPrctile <- 0.25

#number of permutations
nPermutations <- 10

resultsPermuteX <- viladomatCorrelation(data, delta, maxDistPrctile, nPermutations)
resultsPermuteX$pValueGlobal
#> [1] 0.09090909
resultsPermuteX$deltaStar
#>  [1] 0.9 0.8 0.4 0.8 0.9 0.5 0.8 0.7 0.9 0.8
```
