# Package index

## Compare two samples

The main function: the spatial correlation test and the spatial
similarity for every gene, and the methods of its result.

- [`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md)
  : Compare the spatial expression patterns of genes between two samples
- [`print(`*`<STcompareResult>`*`)`](https://jef.works/STcompare/reference/STcompareResult.md)
  [`summary(`*`<STcompareResult>`*`)`](https://jef.works/STcompare/reference/STcompareResult.md)
  [`print(`*`<summary.STcompareResult>`*`)`](https://jef.works/STcompare/reference/STcompareResult.md)
  [`` `[`( ``*`<STcompareResult>`*`)`](https://jef.works/STcompare/reference/STcompareResult.md)
  [`as.data.frame(`*`<STcompareResult>`*`)`](https://jef.works/STcompare/reference/STcompareResult.md)
  : Methods for compareSpatial() results

## Spatial correlation, legacy functions

The correlation functions of earlier versions. They reproduce the
published analyses exactly.

- [`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md)
  : Spatial correlation test for every gene of two samples
- [`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
  : Spatial correlation test with more permutations for promising genes
- [`spatialCorrelationGeneExpWithinSample()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpWithinSample.md)
  : Spatial correlation test between pairs of genes of one sample
- [`spatialCorrelation()`](https://jef.works/STcompare/reference/spatialCorrelation.md)
  : Spatial correlation test of two vectors of values at the same
  locations
- [`viladomatCorrelation()`](https://jef.works/STcompare/reference/viladomatCorrelation.md)
  : Spatially autocorrelated surrogates and permutation p-value for one
  direction

## Visualization

Standard figures of a compareSpatial() result: one gene with its maps,
similarity classes and statistics, and the overview of all genes. The
legacy functions plot single genes from the results of the legacy
functions.

- [`plotGene()`](https://jef.works/STcompare/reference/plotGene.md) :
  Standard figure of one gene in a compareSpatial() result
- [`plotOverview()`](https://jef.works/STcompare/reference/plotOverview.md)
  : Overview figure of a compareSpatial() result
- [`linearRegression()`](https://jef.works/STcompare/reference/linearRegression.md)
  : Plot a gene's values in two samples, coloured by similarity class
- [`pixelClass()`](https://jef.works/STcompare/reference/pixelClass.md)
  : Map the similarity class of every pixel
- [`plotCorrelationGeneExp()`](https://jef.works/STcompare/reference/plotCorrelationGeneExp.md)
  : Plot a gene's values in two samples at the shared pixels
- [`savePlots()`](https://jef.works/STcompare/reference/savePlots.md) :
  Maps, pixel classes and scatter plot of genes, optionally saved as PDF

## Spatial similarity, legacy function

The similarity with the class of every pixel.

- [`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md)
  : Computes spatial similarity of gene expression between two
  SpatialExperiments.

## Example data

- [`speKidney`](https://jef.works/STcompare/reference/speKidney.md) :
  Simulated Kidney SpatialExperiment Data
- [`simRanPatternRasts`](https://jef.works/STcompare/reference/simRanPatternRasts.md)
  : Simulated Random Spatially Patterned Kidney Datasets
