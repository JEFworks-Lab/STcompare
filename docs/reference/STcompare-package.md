# STcompare: Framework for Spatial Transcriptomic Data Comparison

Compares the spatial gene expression patterns of two spatial
transcriptomics datasets of structurally matched tissues, rasterized
onto a shared grid of pixels. For each gene, a spatial correlation test
gives a permutation p-value that accounts for the spatial
autocorrelation of both patterns (Viladomat et al. (2014)
[doi:10.1111/biom.12139](https://doi.org/10.1111/biom.12139) ), and a
spatial fold-change similarity measures how often the expression at
matched locations is within a fold change. The method is described in
Clifton et al. (2026)
[doi:10.1093/bioinformatics/btag644](https://doi.org/10.1093/bioinformatics/btag644)
.

## Details

The main function is
[`compareSpatial()`](https://jef.works/STcompare/reference/compareSpatial.md):
for every gene measured in two samples rasterized onto one grid of
pixels, it tests whether the spatial patterns are correlated beyond what
their spatial autocorrelation alone would give, and measures how often
the expression at matched pixels is within a fold change. The functions
of the published analyses
([`spatialCorrelationGeneExp()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExp.md),
[`spatialCorrelationGeneExpIterPermutations()`](https://jef.works/STcompare/reference/spatialCorrelationGeneExpIterPermutations.md)
and
[`spatialSimilarity()`](https://jef.works/STcompare/reference/spatialSimilarity.md))
reproduce them exactly.

## References

Clifton K, Jiang V, Peixoto RdS, Singh S, Matsuura R, Rabb H, Fan J
(2026). STcompare: comparative spatial transcriptomics data analysis of
structurally matched tissues to characterize differentially spatially
patterned genes. *Bioinformatics* 42(9):btag644.
[doi:10.1093/bioinformatics/btag644](https://doi.org/10.1093/bioinformatics/btag644)

## See also

Useful links:

- <https://jef.works/STcompare/>

- <https://github.com/JEFworks-Lab/STcompare>

- Report bugs at <https://github.com/JEFworks-Lab/STcompare/issues>

## Author

**Maintainer**: Kalen Clifton <kclifto2@jh.edu>
([ORCID](https://orcid.org/0000-0002-4095-8104))

Authors:

- Vivien Jiang <vjiang1@jhu.edu>
  ([ORCID](https://orcid.org/0009-0002-1417-4389))

- Jean Fan <jeanfan@jhu.edu>
  ([ORCID](https://orcid.org/0000-0002-0212-5451)) \[reviewer\]
