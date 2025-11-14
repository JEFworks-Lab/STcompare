## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 8, 
  fig.height = 6, 
  dpi = 72, 
  dev = "png"
)

## ----eval = FALSE, message=FALSE, warning=FALSE-------------------------------
# require(remotes)
# remotes::install_github('JEFworks-Lab/STcompare')
# 

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(STcompare)
library(SpatialExperiment)
library(SEraster)
library(ggplot2)


## ----load_data----------------------------------------------------------------

data("speKidney") 

## ----rasterization,  message=FALSE, warning=FALSE-----------------------------
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
                assay_name = 'counts', resolution = 0.2,
                square = FALSE)

pA <- plotRaster(rastKidney$A, plotTitle = "Simulated Kidney A")
pA

pB <- plotRaster(rastKidney$B, plotTitle = "Simulated Kidney B")
pB

pC <- plotRaster(rastKidney$C, plotTitle = "Simulated Kidney C")
pC



## ----box_plot,  message=FALSE, warning=FALSE----------------------------------

df <- data.frame(
  value = c(
    as.numeric(assay(rastKidney$A, "pixelval")), 
    as.numeric(assay(rastKidney$B, "pixelval")), 
    as.numeric(assay(rastKidney$C, "pixelval"))
  ),
  sample = factor(rep(c("A", "B", "C"),
                      times = c(
                        length(as.numeric(assay(rastKidney$A, "pixelval"))), 
                        length(as.numeric(assay(rastKidney$B, "pixelval"))), 
                        length(as.numeric(assay(rastKidney$C, "pixelval")))
                      )
                  ))
)

compare_box_plot <- ggplot(df, aes(x = sample, y = value, fill = sample)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.7) +
  labs(x = "Sample", y = "Pixel Intensity", title = "Expression Magnitude Across Samples") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
compare_box_plot



## ----spatial_similarity,  message=FALSE, warning=FALSE------------------------

# compare kidney A and B 
sAB <- spatialSimilarity(
  list(A = rastKidney$A, B = rastKidney$B)
)
  
# compare kidney A and C 
sAC <- spatialSimilarity(
  list(A = rastKidney$A, C = rastKidney$C)
)
  
# compare kidney B and C 
sBC <- spatialSimilarity(
  list(B = rastKidney$B, C = rastKidney$C)
)


## ----lr_pc_AB,  message=FALSE, warning=FALSE----------------------------------

# linear regression and pixel classification plots 

lrAB <- linearRegression(input=sAB, gene = "Gene")
lrAB

pcAB <- pixelClass(input=sAB, gene="Gene")
pcAB


## ----lr_pc_AC,  message=FALSE, warning=FALSE----------------------------------

lrAC <- linearRegression(input=sAC, gene = "Gene")
lrAC

pcAC <- pixelClass(input=sAC, gene="Gene")
pcAC

## ----lr_pc_BC,  message=FALSE, warning=FALSE----------------------------------

lrBC <- linearRegression(input=sBC, gene = "Gene")
lrBC

pcBC <- pixelClass(input=sBC, gene="Gene")
pcBC



## ----spatialCorrelationAB,  message=FALSE, warning=FALSE----------------------

rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
scAB <- spatialCorrelationGeneExp(rastGexpListAB)
scAB

rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
expAB <- plotCorrelationGeneExp(rastGexpListAB, scAB, "Gene")
expAB


## ----spatialCorrelationAC,  message=FALSE, warning=FALSE----------------------


rastGexpListAC <- list(A = rastKidney$A, C = rastKidney$C)
scAC <- spatialCorrelationGeneExp(rastGexpListAC)
scAC

expAC <- plotCorrelationGeneExp(rastGexpListAC, scAC, "Gene")
expAC


## ----spatialCorrelationBC,  message=FALSE, warning=FALSE----------------------

rastGexpListBC <- list(B = rastKidney$B, C = rastKidney$C)
scBC <- spatialCorrelationGeneExp(rastGexpListBC)
scBC


rastGexpListBC <- list(B = rastKidney$B, C = rastKidney$C)
expBC <- plotCorrelationGeneExp(rastGexpListBC, scBC, "Gene")
expBC



## ----show_ran_kid, eval=TRUE, echo=TRUE, message=FALSE------------------------

plotRaster(simRanPatternRasts[[1]], plotTitle = "simRanPatternRasts Kidney 1")
plotRaster(simRanPatternRasts[[7]], plotTitle = "simRanPatternRasts Kidney 7")

rastGexpList <- list(kidney_1 = simRanPatternRasts[[1]], kidney_7 = simRanPatternRasts[[7]])

sc <- spatialCorrelationGeneExp(rastGexpList)
sc

plotCorrelationGeneExp(rastGexpList, sc, "1")



## ----load-results, eval=TRUE, echo=TRUE, message=FALSE------------------------

# This code was run once to generate the precomputed vignette results.
# Full script available at:
# system.file("scripts", "simRanPatternSpatialCorrelation.R", package = "STcompare")

data_path <- system.file("extdata", "simRanPatternResults.RData", package = "STcompare")
load(data_path)
corrected


