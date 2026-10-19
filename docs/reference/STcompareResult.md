# Methods for compareSpatial() results

[`print()`](https://rdrr.io/r/base/print.html) shows a header (genes,
pixels, samples and the test settings), the number of genes with
`pValuePermuteAdj < 0.05` by the sign of `correlationCoef`, how many
genes stopped early or were skipped, the median `percentSimilarity`, and
the first rows. [`summary()`](https://rdrr.io/r/base/summary.html)
returns these counts (and more: the distribution of `nPermutations`, the
settings, the genes with `deltaGridEdge`, `summary(percentSimilarity)`)
as an object with its own [`print()`](https://rdrr.io/r/base/print.html)
method; `alpha` sets the threshold on `pValuePermuteAdj`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the table as a plain data frame, without the class and the attributes.
Selecting rows with `[` keeps the class and the attributes (`details` is
subset to the selected genes); selecting columns gives a plain data
frame. The columns are described in
[`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md).

## Usage

``` r
# S3 method for class 'STcompareResult'
print(x, n = 6L, ...)

# S3 method for class 'STcompareResult'
summary(object, alpha = 0.05, ...)

# S3 method for class 'summary.STcompareResult'
print(x, ...)

# S3 method for class 'STcompareResult'
x[i, j, ..., drop]

# S3 method for class 'STcompareResult'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x, object:

  A result of
  [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md).

- n:

  The number of rows to print.

- alpha:

  The threshold on `pValuePermuteAdj` below which
  [`summary()`](https://rdrr.io/r/base/summary.html) counts a gene as
  significant.

- i, j, drop:

  Rows and columns to select, as for data frames.

- row.names, optional, ...:

  Passed on or ignored, as for the generics.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`"summary.STcompareResult"`;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) a data
frame.

## Examples

``` r
data(simRanPatternRasts)
# two independent simulated fields, rasterized onto one grid of pixels
res <- compareSpatial(simRanPatternRasts[[1]], simRanPatternRasts[[2]],
                      nPermutations = 199, nThreads = 2)
#> compareSpatial: 1 gene on 275 shared pixels in 0.1 s (11 permutations in total, 2 threads)
print(res)
#> STcompareResult: 1 gene on 275 shared pixels; x vs y (assay pixelval)
#> Correlation: adaptive permutation p-values (stop at 10 exceedances, at most 199 permutations); 11 deltas from 0.01 to 0.9; rank-remapped surrogates
#>   pValuePermuteAdj < 0.05: 0 genes (0 with correlationCoef > 0, 0 with correlationCoef < 0)
#>   stopped early: 1; reached nPermutations: 0; skipped: 0; failed: 0
#> Similarity (|log2(y / x)| <= 1): median percentSimilarity 1 (0 genes with NA)
#> 
#>  gene nPixels correlationCoef pValueNaive pValuePermuteX pValuePermuteY
#>     1     275         -0.0416       0.492          0.833          0.917
#>  pValuePermute pValuePermuteAdj nPermutations        stop deltaStarMedianX
#>          0.909            0.909            11 exceedances              0.2
#>  deltaStarMedianY deltaGridEdge percentSimilarity percentDissimilarityX
#>               0.2         FALSE                 1                     0
#>  percentDissimilarityY numPixelInThresh   t1   t2 status message
#>                      0              275 8.75 8.93     ok        
summary(res)
#> compareSpatial: 1 gene on 275 shared pixels (x vs y); 1 ok, 0 skipped, 0 failed
#> Correlation, pValuePermuteAdj < 0.05: 0 genes (0 with correlationCoef > 0, 0 with correlationCoef < 0)
#>   stopped early: 1; reached nPermutations: 0; skipped: 0; failed: 0
#>   nPermutations per gene: median 11, max 11; 11 in total
#>   settings: nPermutations = 199, exceedances = 10, 11 deltas (0.01 to 0.9), rank-remapped surrogates, BH adjustment, seed 0
#> Similarity, percentSimilarity:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>       1       1       1       1       1       1 
#> Run time: 0.1 s
head(as.data.frame(res))
#>   gene nPixels correlationCoef pValueNaive pValuePermuteX pValuePermuteY
#> 1    1     275     -0.04158115   0.4922708      0.8333333      0.9166667
#>   pValuePermute pValuePermuteAdj nPermutations        stop deltaStarMedianX
#> 1     0.9090909        0.9090909            11 exceedances              0.2
#>   deltaStarMedianY deltaGridEdge percentSimilarity percentDissimilarityX
#> 1              0.2         FALSE                 1                     0
#>   percentDissimilarityY numPixelInThresh       t1       t2 status message
#> 1                     0              275 8.748878 8.928013     ok        
```
