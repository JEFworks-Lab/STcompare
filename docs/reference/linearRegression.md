# Plot a gene's values in two samples, coloured by similarity class

This function creates a scatter plot comparing gene expression levels
between two spatial experiments for a specified gene. It colors data
points based on similarity classification and overlays fold-change
threshold lines. No regression line is fitted: the solid line is \\y =
x\\. Both axes end at the larger of the two 95% quantiles of the gene's
values, so the pixels with the highest values are not drawn (with a
warning from ggplot2).

## Usage

``` r
linearRegression(input, gene, assayName = NULL)
```

## Arguments

- input:

  A list. Results from
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md).
  This includes the similarity table, log-transformed pixel data, and
  analysis parameters.

- gene:

  Character. The name of the gene to visualize.

- assayName:

  A character string or numeric specifying the assay in the Spatial
  Experiment to use. Default is `NULL`: the assay that
  [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  used (the first assay for results of STcompare 0.1.0 and earlier,
  which did not record it).

## Value

A ggplot2 scatter plot displaying gene expression values from two
spatial experiments. Data points are colored as follows:

- **blue**:

  Pixels classified as similar (within the fold-change threshold).

- **yellow**:

  Pixels with greater expression in dataset X than Y.

- **red**:

  Pixels with greater expression in dataset Y than X.

- **grey**:

  Pixels with gene expression below the threshold in both experiments.

The plot includes:

- **Solid line:**:

  y = x (equal expression).

- **Dashed lines:**:

  Fold-change similarity thresholds (upper and lower bounds).

## See also

[`plotGene()`](https://jef.works/STcompare/reference/plotGene.md), the
standardized figure of a gene in a
[`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
result, which includes this scatter plot.

## Examples

``` r
data(speKidney)
##### Rasterize to get pixels at matched spatial locations #####
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
               assay_name = 'counts', resolution = 0.2, fun = "mean",
               square = FALSE)

s <- spatialSimilarity(list(A = rastKidney$A, C = rastKidney$C))
linearRegression(s, "Gene")
#> Warning: Removed 14 rows containing missing values or values outside the scale range
#> (`geom_point()`).

```
