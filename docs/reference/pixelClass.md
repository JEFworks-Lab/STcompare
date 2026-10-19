# Map the similarity class of every pixel

This function visualizes the spatial distribution of gene expression
similarity by classifying pixels into four categories: below threshold,
similar, higher in the first dataset and higher in the second dataset.
The plot is drawn with the pixel geometries of rasterized
SpatialExperiment objects (`colData()` column `geometry`, which needs
the sf package), otherwise with points at the spatial coordinates.

## Usage

``` r
pixelClass(input, gene, assayName = NULL)
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
  used. The classification comes from `input`, so the plot does not
  depend on it.

## Value

A ggplot2 spatial plot displaying classified pixels, where:

- **blue**:

  Pixels classified as similar (within the fold-change threshold).

- **yellow**:

  Pixels with greater expression in dataset X than Y.

- **red**:

  Pixels with greater expression in dataset Y than X.

- **grey**:

  Pixels with gene expression below the threshold in both experiments.

The plot title includes the gene name and its similarity score.

## See also

[`plotGene()`](https://jef.works/STcompare/reference/plotGene.md), the
standardized figure of a gene in a
[`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
result, which includes this map.

## Examples

``` r
data(speKidney)
##### Rasterize to get pixels at matched spatial locations #####
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
                assay_name = 'counts', resolution = 0.2, fun = "mean",
                square = FALSE)
s <- spatialSimilarity(list(A = rastKidney$A, B = rastKidney$B))
pixelClass(s, "Gene")

```
