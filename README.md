
# STcompare

`STcompare` is an R package for comparing spatial gene expression patterns across spatial transcriptomics (ST) datasets, addressing cases where traditional non-spatial differential expression methods fail to capture spatially distinct variation. [jef.works/STcompare](jef.works/STcompare) is the `STcompare` R documentation website.  	

## Overview 


STcompare enables differential spatial comparison tests orthogonal to traditional differential gene (DGE) analysis. With DGE methods, simulated spatially distinct gene expression patterns A and B are identified as not differentially expressed because they have the same mean gene expression.  Conversely, highly spatially similar gene expression patterns A and C are differentially expressed because they have very different mean gene expression. 


<p align="center">
  <img src="https://github.com/JEFworks-Lab/STcompare/blob/main/images/overview_figure1.png?raw=true" width="700">
</p>

`STcompare` enables two differential spatial comparison tests: 

1.	**Spatial Correlation:** Person correlation assumes that each sample is independent and identically distributed random variables. However, in terms of gene expression, the gene expression of one cell influences the gene expression of the neighboring cells. Spatial Correlation computes the empirical p-value instead to consider the dependent relationship. 

2.	**Spatial Fold Change:** Computes a Similarity score for each gene in the comparison, compares the change in expression at matched locations. 

To make comparisons, `STcompare` relies on matched spatial locations in the ST datasets. As such, single cell resolution ST datasets first must be first aligned using [`SEraster`](https://github.com/JEFworks-Lab/SEraster) such that shared tissue structures are coincide in the same common coordinate framework and subsequently must be rasterized onto a shared pixel grid such that there is a one-to-one correspondence of units to be compared.  

## Installation 

```r
require(remotes)
remotes::install_github('JEFworks-Lab/STcompare')
```

## Tutorials 

1. [Getting started with `STcompare`](https://github.com/JEFworks-Lab/STcompare/blob/main/doc/getting-started-with-STcompare.md)

## Citation 

Clifton K, Jiang V, Singh S, Matsuura R, Rabb H, Fan J., STcompare: Comparative spatial transcriptomics data analysis to characterize differentially spatially patterned genes. bioRxiv 2025.11.21.689847. doi: https://doi.org/10.1101/2025.11.21.689847
