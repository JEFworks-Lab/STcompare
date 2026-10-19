# Computes spatial similarity of gene expression between two SpatialExperiments.

This function calculates the spatial similarity of gene expression
patterns between two SpatialExperiment objects based on thresholding and
fold-change criteria. Similarity is defined as the proportion of pixels
where the fold change in gene expression falls within a specified range.

## Usage

``` r
spatialSimilarity(
  input,
  t1 = NULL,
  t2 = NULL,
  minQuantile = 0.05,
  minPixels = 0.1,
  foldChange = 1,
  assayName = NULL,
  verbose = FALSE
)
```

## Arguments

- input:

  List of two SpatialExperiment objects (only the first two elements are
  used). The first element corresponds to the first spatial experiment
  (`x`), and the second to the second spatial experiment (`y`).
  Similarity is computed as the proportion of pixels satisfying
  `-b <= log2(y/x) <= b` for each gene, where `b = foldChange` (default
  1).

- t1:

  `numeric`: Gene expression threshold for the first spatial experiment
  (`x`), a single number used for every gene. Only pixels with values
  greater than `t1` in `x` or greater than `t2` in `y` are used to
  calculate the similarity score. Default is `NULL`: each gene's
  threshold is then its `minQuantile` quantile in `x`.

- t2:

  `numeric`: Gene expression threshold for the second spatial experiment
  (`y`), like `t1`. Default is `NULL`: each gene's threshold is then its
  `minQuantile` quantile in `y`.

- minQuantile:

  `numeric`: The quantile (between 0 and 1) of each gene's values in
  each experiment used as its threshold when `t1` or `t2` is not
  supplied. Default is `0.05`.

- minPixels:

  `numeric`: The smallest proportion (between 0 and 1) of the shared
  pixels that must pass the thresholds. If fewer pixels pass, a spatial
  similarity score is not calculated for the gene (`NA`). Default is
  `0.1`.

- foldChange:

  `numeric`: The similarity band on the log2 scale: pixels with
  \\\|\log_2(y / x)\| \le\\ `foldChange` are similar. Default is `1`,
  that is, within two-fold.

- assayName:

  A character string or numeric specifying the assay in the Spatial
  Experiment to use. Default is `NULL`. If no value is supplied for
  `assayName`, then the first assay is used as a default. The assay used
  is recorded in the result, and
  [`linearRegression()`](https://jef.works/STcompare/reference/linearRegression.md),
  [`pixelClass()`](https://jef.works/STcompare/reference/pixelClass.md)
  and
  [`savePlots()`](https://jef.works/STcompare/reference/savePlots.md)
  use it by default.

- verbose:

  `logical`: if `TRUE`, print a message with the number of genes and
  shared pixels compared. Default is `FALSE`.

## Value

A list containing:

- `similarityTable`:

  A data frame with one row per gene and the following columns:

  `gene`

  :   Gene name.

  `percentSimilarity`

  :   Proportion (between 0 and 1) of the kept pixels that are similar:
      the similarity score S.

  `percentDissimilarityX`

  :   Proportion of the kept pixels for which `log2(y/x) < -foldChange`
      (higher in the first experiment).

  `percentDissimilarityY`

  :   Proportion of the kept pixels for which `log2(y/x) > foldChange`
      (higher in the second experiment).

  `similarPixelID`

  :   List of pixel IDs classified as similar.

  `dissimilarPixelIDX`

  :   List of pixel IDs for `log2(y/x) < -foldChange`.

  `dissimilarPixelIDY`

  :   List of pixel IDs for `log2(y/x) > foldChange`.

  `numPixelInThresh`

  :   Number of pixels kept: above the threshold in either experiment.

  `pixelIDInThresh`

  :   List of the kept pixel IDs.

  `numPixelOutThresh`

  :   Number of pixels removed: at or below the threshold in both
      experiments.

  `pixelIDOutThresh`

  :   List of the removed pixel IDs.

  `t1`

  :   Threshold value for this gene in first spatial experiment.

  `t2`

  :   Threshold value for this gene in second spatial experiment.

  For a gene with fewer than `minPixels` kept pixels, the three
  proportions and the three lists of classified pixels are `NA`.

- `pixelLogTransformation`:

  A data frame with one row per gene that has a similarity score:
  `gene`, and `log`, a list with the `log2(y/x)` values of its kept
  pixels (in the order of `pixelIDInThresh`).

- `parameters`:

  A list with `foldChange`, `minPixels`, the `input` objects, and the
  `assayName` (the assay used), `minQuantile`, `t1` and `t2` arguments.

## Details

The genes compared are those in both objects (row names), and the pixels
compared are those whose names (column names) are in both objects:
rasterize the two samples in one call, such as
`SEraster::rasterizeGeneExpression(list(x, y), ...)`, so that a pixel
name means the same location in both.

For each gene, a pixel is kept when its value in the first object is
above `t1` or its value in the second object is above `t2` (strictly),
so the pixels removed are those at or below the threshold in both
objects. Zeros among the kept pixels are replaced by 1e-4. With \\l =
\log_2(y / x)\\ for a kept pixel, where \\x\\ and \\y\\ are its values
in the first and the second object, the pixel is similar when \\-b \le l
\le b\\ for \\b\\ = `foldChange` (the default 1 means within two-fold,
both ends included), higher in the first object when \\l \< -b\\, and
higher in the second object when \\l \> b\\. If fewer than a share
`minPixels` of the shared pixels are kept, the gene gets no similarity
(`NA`).

The values must be finite and not negative, such as counts or normalized
expression on a linear scale: a fold change is not defined for negative
values, and genes with negative, missing or infinite values on the
shared pixels stop the function with an error that names them. All genes
are computed at once by the helper that also computes the similarity
columns of
[`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md),
which gives the same numbers.

## See also

[`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md),
which computes the similarity together with the spatial correlation
test;
[`linearRegression()`](https://jef.works/STcompare/reference/linearRegression.md),
[`pixelClass()`](https://jef.works/STcompare/reference/pixelClass.md)
and [`savePlots()`](https://jef.works/STcompare/reference/savePlots.md)
to plot the result.

## Examples

``` r
data(speKidney)
##### Rasterize to get pixels at matched spatial locations #####
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
               assay_name = 'counts', resolution = 0.2, fun = "mean",
               square = FALSE)

s <- spatialSimilarity(list(rastKidney$A, rastKidney$C))
s$similarityTable[, c("gene", "percentSimilarity", "percentDissimilarityX",
                      "percentDissimilarityY", "numPixelInThresh")]
#>   gene percentSimilarity percentDissimilarityX percentDissimilarityY
#> 1 Gene                 0                     0                     1
#>   numPixelInThresh
#> 1              272
```
