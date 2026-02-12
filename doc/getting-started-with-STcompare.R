## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  dpi = 72,
  dev = "png",
  fig.path = "getting-started-with-STcompare-figures/"

)

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
# require(remotes)
# remotes::install_github('JEFworks-Lab/STcompare')
# 

## ----import_packages, message=FALSE, warning=FALSE----------------------------
library(STcompare)
library(SpatialExperiment)
library(SEraster)
library(ggplot2)
library(patchwork)
library(viridis)



## ----load_data----------------------------------------------------------------

data("speKidney") 
head(speKidney)

## ----spatialExpression,  message=FALSE, warning=FALSE-------------------------

# This function takes the Spatial Experiment Object and creates a plot 
# of the spatial patterns of gene expression. 
plotSpatialPatterns <- function(spatialExperimentObj, name) {
  df <- cbind(
    as.data.frame(spatialCoords(spatialExperimentObj)),
    as.data.frame(colData(spatialExperimentObj))
  )
  
  df$Gene <- as.numeric(assay(spatialExperimentObj, "counts")["Gene", ])
  
  ggplot(df, aes(x = x, y = y, color = Gene)) +
    geom_point(size = 0.6) +
    coord_equal() +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    scale_color_viridis(option = "D") +
    labs(
      x = "x",
      y = "y",
      color = "Gene",
      title = paste("Simulated Kidney", name)
    )
}

pA <- plotSpatialPatterns(speKidney$A, "A")
pB <- plotSpatialPatterns(speKidney$B, "B")
pC <- plotSpatialPatterns(speKidney$C, "C")

pA + pB + pC


## ----box_plot,  message=FALSE, warning=FALSE----------------------------------

# Convert the counts for each kidney (A, B, C) into one dataframe.
# Each kidney’s expression counts are extracted from the "counts" assay and stored as a numeric vector.
df <- data.frame(
  value = c(
    as.numeric(assay(speKidney$A, "counts")), 
    as.numeric(assay(speKidney$B, "counts")), 
    as.numeric(assay(speKidney$C, "counts"))
  ),
  
  # Create a "sample" column showing which kidney each pixel belongs to.
  # The `rep()` call ensures that the correct number of labels (A/B/C) is assigned
  # based on the length of each counts vector.
  sample = factor(rep(
    c("A", "B", "C"),
    times = c(
      length(as.numeric(assay(speKidney$A, "counts"))), 
      length(as.numeric(assay(speKidney$B, "counts"))), 
      length(as.numeric(assay(speKidney$C, "counts")))
    )
  ))
)

# Plot a boxplot comparing the distribution of pixel intensities across the 3 kidneys.
# This represents a traditional magnitude-only comparison, 
# which cannot distinguish pattern similarity/differences.
compare_box_plot <- ggplot(df, aes(x = sample, y = value, fill = sample)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.7) +
  labs(x = "Sample", y = "Pixel Intensity", title = "Expression Magnitude Across Samples") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
compare_box_plot



## ----rasterization,  message=FALSE, warning=FALSE-----------------------------
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
                assay_name = 'counts', resolution = 0.2,
                square = FALSE)

# After rasterization, the output is a SpatialExperiment object, but the spatial units are now raster pixels rather than individual cells or spots. The assay has been converted to pixelval, and additional metadata (num_cell, cellID_list, geometry) records which original cells contributed to each pixel.
head(rastKidney)

# These are the plots to visualize what the kidneys looks like 
pA <- plotRaster(rastKidney$A, plotTitle = "Rasterized Kidney A")
pB <- plotRaster(rastKidney$B, plotTitle = "Rasterized Kidney B")
pC <- plotRaster(rastKidney$C, plotTitle = "Rasterized Kidney C")

pA + pB + pC


## ----spatialCorrelationAB, eval = TRUE, message=FALSE, warning=FALSE----------

# From the list of rasterized kidneys, rastKidney, we will take a subset of rastKidney 
# rastGexpListAB is a list of two kidneys, A and B 
rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)

# spatialCorrelationGeneExp takes input of a list of two SpatialExperiment objects -- rastGexpListAB
# nThreads, the default is 1, should be set to the number of cores available to allow for parallel computing 
scAB <- spatialCorrelationGeneExp(rastGexpListAB, nThreads = 1)

# correlationCoef is showing a negative linear relationship 
# Both the naive and the permuted p-values (pValuePermuteY and pValuePermuteX) are showing that the correlation is significant 
head(scAB)


## ----spatialCorrelationABplot, eval = TRUE, message=FALSE, warning=FALSE------
# visualization of the negative correlation 
# plotCorrelationGeneExp needs the list of rasterized SpatialExperiment objects, the table from spatialCorrelationGeneExp of the same objects, 
# and the gene name you are trying to plot. In the case, the gene name is "Gene". 
expAB <- plotCorrelationGeneExp(rastGexpListAB, scAB, "Gene")
expAB


## ----spatialCorrelationAC,  message=FALSE, warning=FALSE----------------------

# From the list of rasterized kidneys, rastKidney, we will take a subset of rastKidney 
# rastGexpListAC is a list of two kidneys, A and C
rastGexpListAC <- list(A = rastKidney$A, C = rastKidney$C)

# correlationCoef is showing a positive linear relationship 
# Both the naive and the permuted p-values (pValuePermuteY and pValuePermuteX) are showing that the correlation is significant 
scAC <- spatialCorrelationGeneExp(rastGexpListAC)
head(scAC)

# visualization of the positive correlation
# plotCorrelationGeneExp needs the list of rasterized SpatialExperiment objects, the table from spatialCorrelationGeneExp of the same objects, 
# and the gene name you are trying to plot. In the case, the gene name is "Gene".
expAC <- plotCorrelationGeneExp(rastGexpListAC, scAC, "Gene")
expAC


## ----show_ran_kid, eval=TRUE, echo=TRUE, message=FALSE------------------------

# using plotRaster from the SEraster package to visualize kidney 1 and 7
plotRaster(simRanPatternRasts[[1]], plotTitle = "simRanPatternRasts Kidney 1")
plotRaster(simRanPatternRasts[[7]], plotTitle = "simRanPatternRasts Kidney 7")

# taking a subset of simRanPatternRasts
# rastGexpList is a list of rasterize SpatialExperiment objects kidney 1 and kidney 7 
rastGexpList <- list(kidney_1 = simRanPatternRasts[[1]], kidney_7 = simRanPatternRasts[[7]])

# finding the correlation and the p-value of kidney 1 and kidney 7 
sc <- spatialCorrelationGeneExp(rastGexpList)

# the naive p-value is showing that the correlation is significant 
# while both of the permuted p-values are showing not significant p-values 
head(sc) 

# the plot further shows that there isn't a correlation between kidney 1 and kidney 7
plotCorrelationGeneExp(rastGexpList, sc, "1")



## ----load-results, eval=TRUE, echo=TRUE, message=FALSE------------------------

# Full script to generate the precomputed vignette results available at:
# system.file("scripts", "simRanPatternSpatialCorrelation.R", package = "STcompare")

# load data 
data_path <- system.file("extdata", "simRanPatternResults.RData", package = "STcompare")
load(data_path)

# Var1 and Var2 is every pairwise combination between the 100 kidneys 
# cors is the correlation coefficient of that pair 
# corspv is the naive p-value for that pair 
# corspv_correct is the permuted p-value for that pair (chosen to be the higher of pValuePermuteY and pValuePermuteX)
head(cors_df)

# remove the na values 
cors_df <- na.omit(cors_df)

naive_vs_correct_p_value <- ggplot(cors_df) +
  geom_point(aes(x =cors, y = -log10(corspv)), alpha = 0.1, size = 0.5, color = "blue") +
  geom_point(aes(x =cors, y = -log10(corspv_corrected)), alpha = 0.1, size = 0.5, color = "green") +
  theme_classic() +
  labs(x = "Correlation", y = "-log10(p-value)",
       title = "Naïve vs Spatially Corrected p-values") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "black")  +
  ylim(min(-log10(cors_df$corspv), na.rm = TRUE), max(-log10(cors_df$corspv), na.rm = TRUE)) +
  labs(x = "Correlation" , y = "-log10(p-value)")
naive_vs_correct_p_value



## ----spatial_similarity,  message=FALSE, warning=FALSE------------------------

# spatialSimilarity also takes in a list of 2 SpatialExperiment objects you want to compare between 
# t1 and t2 are the thresholding value for the first and second object in the list. 
# The threshold is used to remove pixels that have no or little expression for a given gene in the pixel.
# Failure to remove those those pixels will result in a trivial similarity score due to noise. 
# If t1 and t2 aren't provided, minQuantile is used to threshold based on the quantile of expression. 
# The default for minQuantile is 0.05, meaning by default the pixels at the bottom 5% are removed. 
# minPixels is the percentage of pixels left after thresholding. The default is 0.1. 
# When there is not enough pixels left, the similarity score will also be trivial.
# foldChange is the number of fold that are considered similar. The default is 1 fold. 

# compare kidney A and B 
sAB <- spatialSimilarity(
  list(A = rastKidney$A, B = rastKidney$B)
)

# compare kidney A and C 
sAC <- spatialSimilarity(
  list(A = rastKidney$A, C = rastKidney$C)
)


## ----lr_pc_AB,  message=FALSE, warning=FALSE----------------------------------

# To get the linear regression and the pixel classification plots, 
# the inputs are the spatialSimilarity output and the gene name 

# The Linear regression shows the areas in which gene expression 
# each pixels falls, the pixel is either: higher in B, similar, higher in A  
lrAB <- linearRegression(input=sAB, gene = "Gene")
lrAB

# Similar to the linear regression plot, the pixel classification plot takes 
# the plot of the rasterized kidney and classifies each pixel with the pixels 
# that are below the threshold are gray 
pcAB <- pixelClass(input=sAB, gene="Gene")
pcAB


## ----lr_pc_AC,  message=FALSE, warning=FALSE----------------------------------

# all the pixels have higher expression in C than in A, 
# and are outside of the similar expression boundaries

lrAC <- linearRegression(input=sAC, gene = "Gene")
lrAC

pcAC <- pixelClass(input=sAC, gene="Gene")
pcAC

