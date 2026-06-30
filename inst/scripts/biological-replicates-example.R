# Introduction ####
# To demonstrate STcompare on real ST data, we compared two biological replicates of coronal slices of the adult mouse brain from different animals but the same location with respect to bregma assayed by the single-cell resolution ST technology MERFISH. 
# Because these are biological replicates, we hypothesized that genes should be significantly positively correlated with high similarity, particularly for SVGs which exhibit coordinated expression variation within tissue structures that likely have correspondence across replicates. 

# Load Libraries ####

library(STcompare)
library(SpatialExperiment)
library(SEraster)
library(viridis)
library(MERINGUE)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(BiocGenerics)
library(BiocParallel)
devtools::load_all()

# Load Data ####

# For `STcompare` to produce meaningful comparisons, the tissues must first be spatially aligned so that corresponding structures are matched across samples.
# We used the `STalign` package to align these two replicate MERFISH mouse brain tissues. [STalign Tutorial]("https://jef.works/STalign/notebooks/merfish-merfish-alignment.html")

s2r2_url <- "https://zenodo.org/records/10724029/files/STalign_S2R2.csv.gz?download=1"
s2r3_url <- "https://zenodo.org/records/10724029/files/STalign_S2R3_to_S2R2.csv.gz?download=1"

# read the data from the url 
read_csv_gz_url <- function(x, ...) {
  u  <- url(x, open = "rb")
  gz <- gzcon(u, text = TRUE)
  on.exit({
    try(close(gz), silent = TRUE)
    try(close(u),  silent = TRUE)
  }, add = TRUE)

  read.csv(gz, ...)
}

target <- read_csv_gz_url(s2r2_url)
source <- read_csv_gz_url(s2r3_url)

# Convert Data to SpatialExperiment Objects ####

## 1) target ####################################################################
dim(target)
target[1:5, 1:10]
# The first column is the cell id, the second and third columns are the x,y coordinates of the cell. The following columns are the counts of each gene per cell.


# extract the cell positions (x, y)
pos_target <- target[, c('x', 'y')]
rownames(pos_target) <- target$X
head(pos_target)

# after the third column, the remaining columns are the genes
# making this the counts matrix for spatial expression 
gene_target <- target[, 4:dim(target)[2]]
rownames(gene_target) <- target$X
gene_target[1:5, 1:5]

# for spatial experiments, convert from a dataframe to a dgCMatrix
class(gene_target)
target_sparse <- as(t(gene_target), "dgCMatrix")

# format into SpatialExperiment
spe_target <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = target_sparse),
  spatialCoords = as.matrix(pos_target)
)

## 2) source  ###################################################################

dim(source)
source[1:5, 1:10]
# The first column is the cell id, 
# the second and third columns are the unaligned x,y coordinates of the cell
# the fourth and fifth columns are the x,y coordinates aligned to the target with STalign
# the sixth and seventh columns are the x,y coordinates aligned to the target with an affine alignment

# extract the cell positions (x, y) after alignment 
pos_source <- source[, c('STalign_x', 'STalign_y')]
rownames(pos_source) <- source$X
colnames(pos_source) <- c("x", "y")
head(pos_source)

# the genes in source start at the 8th column 
gene_source <- source[, 8:dim(source)[2]]
rownames(gene_source) <- source$X
gene_source[1:5, 1:5]

# need to convert from data.frame to dgCMatrix
class(gene_source)
source_sparse <- as(t(gene_source), "dgCMatrix")

spe_source <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = as(t(gene_source), "dgCMatrix")),
  spatialCoords = as.matrix(pos_source)
)

identical(rownames(spe_target), rownames(spe_source))

# Rasterize ####

# To ensure matched spatial locations for comparison, we performed rasterization using SEraster to obtain matching 200µm pixel locations across the two samples.

spe_list <- list(target = spe_target, source = spe_source)

res = 200 
output <- SEraster::rasterizeGeneExpression(
  spe_list, 
  resolution = res, 
  BPPARAM = BiocParallel::SerialParam()
)

# Identify Spatially Variables Genes ####

# We identified 415 SVGs shared across both replicates by using MERINGUE to generate a list of SVGs for each replicate from the 483 genes assayed by MERFISH and then selecting the intersection of these lists. 
moransI <- function(SE_temp, filterDist = resolution){
  
  # Get neighbor-relationships
  w <- MERINGUE::getSpatialNeighbors(spatialCoords(SE_temp), filterDist = filterDist)
  par(mfrow=c(1,1))
  plotNetwork(spatialCoords(SE_temp), w)
  
  # Identify significantly spatially auto-correlated genes
  start_time <- Sys.time()
  I <- MERINGUE::getSpatialPatterns(as.matrix(assay(SE_temp)), w)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  results.filter <- MERINGUE::filterSpatialPatterns(mat = assay(SE_temp),
                                          I = I,
                                          w = w,
                                          adjustPv = TRUE,
                                          alpha = 0.05,
                                          minPercentCells = 0.01,
                                          verbose = TRUE, details = TRUE)
  
  return(results.filter)
}

start_time <- Sys.time()
moransI_merfish_target <- moransI(output$target, filterDist = res)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 46.27861 secs

start_time <- Sys.time()
moransI_merfish_source <- moransI(output$source, filterDist = res)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 44.94201 secs

svg_merfish_target <- moransI_merfish_target  %>% rownames()
svg_merfish_source <- moransI_merfish_source  %>% rownames()
save(svg_merfish_target, svg_merfish_source, file = "~/ST_compare/data/merfish_data/moransI_20260622.RData")
load(file = "~/ST_compare/data/merfish_data/moransI_20260622.RData")

svg_int <- intersect(svg_merfish_target, svg_merfish_source)
svg_uni <- union(svg_merfish_target, svg_merfish_source)
non_svg <- rownames(output$target)[!(rownames(output$target) %in% svg_int)]

length(svg_int) #[1] 489
length(non_svg) #[1] 160

# Using STcompare ####
## Spatial Correlation ####

# We applied STcompare’s spatial correlation test the genes, calculating a Pearson’s correlation coefficient r and an empirical p-value p_E for each
genes_only <- rownames(output$source)[!grepl("Blank", rownames(output$source))]
output <- list(target = output$target[genes_only,], source = output$source[genes_only,])

# To account for the smaller structures observed in the brain we added two additional lower values of 0.01 and 0.05 to the sequence of values tested to find the optimal delta for the Gaussian kernel to generate the permutations.
deltaList <- rep(list(c(0.01, 0.05, seq(0.1, 0.9, .1))), length(rownames(output$source)))

#shortList <- list(target = output$target[1:5,], source = output$source[1:5,])
start_time <- Sys.time()
merfishCorrelation <- spatialCorrelationGeneExpIterPermutations(
  output, 
  deltaX = deltaList, deltaY = deltaList,
  nThreads = 20,
  BPPARAM = BiocParallel::MulticoreParam()
)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 16.80652 hours

save(merfishCorrelation, file = "~/ST_compare/data/merfish_data/merfish_correlation_delta_0_01_0_9_BH_Iter1000_20260623.RData")
load(file = "~/ST_compare/data/merfish_data/merfish_correlation_delta_0_01_0_9_BH_Iter1000_20260623.RData")

## Results for SVGs ####

#number of genes
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>%
  dim()

#number of genes as svgs
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dim() #[1] 415  11

#check for NA
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(is.na(pValuePermuteX)) %>% 
  dim()

#replace NA with rows from previous run from 20260622
tmp_env <- new.env()
load("~/ST_compare/data/merfish_data/merfish_correlation_delta_0_01_0_9_BH_Iter1000_20260622.RData", envir = tmp_env)
merfishCorrelation20260622 <- tmp_env$merfishCorrelation

na_rows <- which(is.na(merfishCorrelation$pValuePermuteX))
row_ids <- rownames(merfishCorrelation)[na_rows]
merfishCorrelation[row_ids, ] <- merfishCorrelation20260622[row_ids, ]

#check for NA
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(is.na(pValuePermuteX)) %>% 
  dim()

save(merfishCorrelation, file = "~/github/STcompare/inst/extdata/merfishCorrelation.RData")

# Full script to generate the precomputed vignette results available at:
# system.file("scripts", "visiumKidneySpatialCorrelation.R", package = "STcompare")

# load data generated from the script 
load(system.file("extdata", "merfishCorrelation.RData", package = "STcompare"))
load(file = "~/github/STcompare/inst/extdata/merfishCorrelation.RData")

## Results for SVGs ####

#number of genes
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>%
  dim()

#number of genes as svgs
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dim() #[1] 415  11

# number of genes as svgs that are significantly positively correlated
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  dim()

# number of genes as svgs that are significantly negatively correlated
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef < 0) %>% 
  dim()

# number of genes as svgs that are not significantly  correlated
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  dim()

# names of genes as svgs that are significantly positively correlated
svgSigPos <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  rownames()

# names of genes as svgs that are not significantly  correlated
svgNotSig <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  rownames()

# In total, we identified 93% (384/415) of genes as significantly positively correlated and zero genes as significantly negative correlated with an adjusted empirical p < 0.05 (Figure 13d), suggesting that the most SVGs have not changed in their spatial patterning, as expected. 

## Results for non-SVGs ####

#number of genes classified as non svgs 
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dim()

# number of genes as non svgs that are significantly positively correlated
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  dim()

# number of genes as non svgs that are not significantly  correlated
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  dim()

# names of genes as non svgs that are significantly positively correlated
nonSVGsig <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  rownames()

# names of genes as non svgs that are not significantly correlated
nonSVGnotsig <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  rownames()

# When performing STcompare’s spatial correlation test on the remaining 68 genes not identified as SVG, 81% (55/68) were identified as not significantly positively correlated using Benjamini-Hochberg adjusted empirical p < 0.05. 
# The high percentage of non-SVGs identified as not significantly positively correlated is consistent with the expectation that genes not exhibiting autocorrelated spatial expression patterns are less likely to have correlated patterns across replicates. 

# Here we visualize the correlation coefficient and empirical p-value for each gene, colored by whether the gene is an SVG or not. The dashed line indicates the threshold for significance at p < 0.05.
plt <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::mutate(SVG = case_when(rownamesCol %in% non_svg ~ "not",
                                rownamesCol %in% svg_int ~ "svg",
                         .default = "Blank")) %>%
  dplyr::filter(SVG %in% c("svg", "not")) %>% 
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x= correlationCoef, y = -log10(pValueEmpiricalRound), color = SVG), alpha = 0.25, size = 3) +
  ggplot2::xlim(NA,1) +
  ggplot2::scale_color_manual(values = c("svg" = "green", "not" = "blue")) +
  ggplot2::theme_classic() +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "black")  + 
  ggplot2::labs(x = "Correlation Coefficient, r" , y = "-log10(p-value)", 
       title = "Rate of significance for SVGs versus non-SVGs"
       )

plt

##### Use STcompare to calculate fold-change similarity score #####
similarityMerfish <- spatialSimilarity(output, foldChange = 1)

df <- data.frame(gene = rownames(merfishCorrelation[svgSigPos,]),
                 correlation = merfishCorrelation[svgSigPos, "correlationCoef" ], 
                 similar = similarityMerfish$similarityTable$percentSimilarity[similarityMerfish$similarityTable$gene %in% svgSigPos])
head(df)

pltCorSim <- ggplot2::ggplot(df, ggplot2::aes(x = correlation, y = similar)) +
  ggplot2::geom_point(size = 1, alpha = 1, color = "black") +
  ggplot2::xlim(c(0,1)) +
  ggplot2::ylim(c(0,1)) +
  ggplot2::coord_fixed() +
  ggplot2::theme_classic() + 
  ggplot2::labs(x = "Correlation Coefficient, r" , y = "Fold-Change Similarity, S", 
       title = "Similarity for positively correlated SVGs"
       )

pltCorSim

# Visualize the spatial expression of genes with high and low similarity scores across the two replicates.
# Gene with high similarity score
gene <- "Gabbr1"

# create a color scale that is consistent across the two replicates for the same gene

sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(output[[1]])),
                            rownames(SpatialExperiment::spatialCoords(output[[2]])))

minList <- min(c(SummarizedExperiment::assay(output[[1]])[gene, sharedPixels], SummarizedExperiment::assay(output[[2]])[gene, sharedPixels]))
maxList <- max(c(SummarizedExperiment::assay(output[[1]])[gene, sharedPixels], SummarizedExperiment::assay(output[[2]])[gene, sharedPixels]))

sc <- ggplot2::scale_fill_gradientn(colors=viridis::viridis(20), 
                           limits = c(minList, maxList),
                           oob=scales::squish,
                           name = gene)

# plot the spatial expression of the gene across the two replicates
pltHighS2R2 <- SEraster::plotRaster(output[[1]][gene, sharedPixels], name = gene) + 
  sc + 
  ggplot2::theme_void()
pltHighS2R3 <- SEraster::plotRaster(output[[2]][gene, sharedPixels], name = gene) + 
  sc + 
  ggplot2::theme_void() 


# plot scatterplot of the gene's expression across the two replicates with points colored by similarity classification
simPlotHigh <- linearRegression(similarityMerfish, gene) + ggplot2::theme_classic() + ggplot2::coord_fixed()

# plot spatial pattern of similarity score
classPlotHigh <- pixelClass(similarityMerfish, gene)

# correlation coefficient and empirical p-value for the gene
merfishCorrelation[gene,]

# Gene with lower similarity score
gene <- "Adgrl1"

# create a color scale that is consistent across the two replicates for the same gene

minList <- min(c(SummarizedExperiment::assay(output[[1]])[gene, sharedPixels], SummarizedExperiment::assay(output[[2]])[gene, sharedPixels]))
maxList <- max(c(SummarizedExperiment::assay(output[[1]])[gene, sharedPixels], SummarizedExperiment::assay(output[[2]])[gene, sharedPixels]))

sc <- ggplot2::scale_fill_gradientn(colors=viridis::viridis(20), 
                           limits = c(minList, maxList),
                           oob=scales::squish,
                           name = gene)

# plot the spatial expression of the gene across the two replicates
pltLowS2R2 <- SEraster::plotRaster(output[[1]][gene, sharedPixels], name = gene) + 
  sc + 
  ggplot2::theme_void()
pltLowS2R3 <- SEraster::plotRaster(output[[2]][gene, sharedPixels], name = gene) + 
  sc + 
  ggplot2::theme_void() 

# plot scatterplot of the gene's expression across the two replicates with points colored by similarity classification
simPlotLow <- linearRegression(similarityMerfish, gene) + ggplot2::theme_classic() + ggplot2::coord_fixed()

# plot spatial pattern of similarity score
classPlotLow <- pixelClass(similarityMerfish, gene)

# correlation coefficient and empirical p-value for the gene
merfishCorrelation[gene,]

#grab repeating legends from plots
LegendHighS2R2 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltHighS2R2), "guide-box")
LegendHighS2R3 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltHighS2R3), "guide-box")
LegendLowS2R2 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltLowS2R2), "guide-box")
LegendLowS2R3 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltLowS2R3), "guide-box")

#multi-panel figure of the spatial expression of genes with high and low similarity scores across the two replicates, scatterplots of the gene's expression across the two replicates with points colored by similarity classification, and spatial pattern of similarity score.
gridExtra::grid.arrange(gridExtra::arrangeGrob(pltHighS2R2 + ggplot2::theme(legend.position="none"),
                                               pltHighS2R3 + ggplot2::theme(legend.position="none"), 
                                               LegendHighS2R2,
                                               simPlotHigh, classPlotHigh + ggplot2::theme(legend.position="none"),
                                               pltLowS2R2 + ggplot2::theme(legend.position="none"),
                                               pltLowS2R3 + ggplot2::theme(legend.position="none"),
                                               LegendLowS2R2,
                                               simPlotLow, classPlotLow + ggplot2::theme(legend.position="none"), 
                                               ncol=5, nrow = 2, widths = c(3,3,1,3,3))
                        )


##### Use STcompare to calculate fold-change similarity score #####
similarityMerfish <- spatialSimilarity(output, foldChange = 1)

df <- data.frame(gene = rownames(merfishCorrelation_affine[svgSigPos,]),
                 correlation = merfishCorrelation_affine[svgSigPos, "correlationCoef" ], 
                 similar = similarityMerfish$similarityTable$percentSimilarity[similarityMerfish$similarityTable$gene %in% svgSigPos])
head(df)

pltCorSim <- ggplot2::ggplot(df, ggplot2::aes(x = correlation, y = similar)) +
  ggplot2::geom_point(size = 1, alpha = 1, color = "black") +
  ggplot2::xlim(c(0,1)) +
  ggplot2::ylim(c(0,1)) +
  ggplot2::coord_fixed() +
  ggplot2::theme_classic() + 
  ggplot2::labs(x = "Correlation Coefficient, r" , y = "Fold-Change Similarity, S", 
       title = "Similarity for positively correlated SVGs"
       )

pltCorSim

#Repeat the above for the affine aligned source coordinates to demonstrate that the results are dependent on the alignment method used.

# extract the cell positions (x, y) after alignment 
pos_source_affine <- source[, c('affine_x', 'affine_y')]
rownames(pos_source_affine) <- source$X
colnames(pos_source_affine) <- c("x", "y")
head(pos_source_affine)

spe_source_affine <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = as(t(gene_source), "dgCMatrix")),
  spatialCoords = as.matrix(pos_source_affine)
)

df <- rbind(data.frame(pos = SpatialExperiment::spatialCoords(spe_target), sample = "target"),
            data.frame(pos = SpatialExperiment::spatialCoords(spe_source_affine), sample = "source_affine"))

pltTarget <- df %>%
  dplyr::filter(sample %in% c("target")) %>%
  ggplot2::ggplot(ggplot2::aes(x = pos.x, y = pos.y)) + 
  ggplot2::geom_point(size = 0.5, alpha = 0.25, color = "gray") +
  ggplot2::theme_void() +
  ggplot2::labs(title = "Target")

pltSourceAffine <- df %>%
  dplyr::filter(sample %in% c("source_affine")) %>%
  ggplot2::ggplot(ggplot2::aes(x = pos.x, y = pos.y)) + 
  ggplot2::geom_point(size = 0.5, alpha = 0.25, color = "gray") +
  ggplot2::theme_void() +
  ggplot2::labs(title = "Source (Affine)")

pltBoth <- df %>%
  dplyr::filter(sample %in% c("target", "source_affine")) %>%
  ggplot2::ggplot(ggplot2::aes(x = pos.x, y = pos.y, color = sample)) + 
  ggplot2::geom_point(size = 0.5, alpha = 0.25) +
  ggplot2::scale_color_manual(values = c("target" = "blue", "source_affine" = "green")) +
  ggplot2::theme_void() +
  ggplot2::labs(title = "Affine Alignment")

# Rasterize ####

# To ensure matched spatial locations for comparison, we performed rasterization using SEraster to obtain matching 200µm pixel locations across the two samples.

spe_list_affine <- list(target = spe_target, source = spe_source_affine)

res = 200 
output_affine <- SEraster::rasterizeGeneExpression(
  spe_list_affine, 
  resolution = res, 
  BPPARAM = BiocParallel::SerialParam()
)

sharedPixels_affine <- intersect(rownames(SpatialExperiment::spatialCoords(output_affine[[1]])),
                            rownames(SpatialExperiment::spatialCoords(output_affine[[2]])))

plt1 <- SEraster::plotRaster(output_affine[[1]][, sharedPixels_affine], name = "total counts") + ggplot2::theme_void()
plt2 <- SEraster::plotRaster(output_affine[[2]][, sharedPixels_affine], name = "total counts") + ggplot2::theme_void()

plt1 + plt2

# We applied STcompare’s spatial correlation test the genes, calculating a Pearson’s correlation coefficient r and an empirical p-value p_E for each
genes_only <- rownames(output_affine$source)[!grepl("Blank", rownames(output_affine$source))]
output_affine <- list(target = output_affine$target[genes_only,], source = output_affine$source[genes_only,])

# To account for the smaller structures observed in the brain we added two additional lower values of 0.01 and 0.05 to the sequence of values tested to find the optimal delta for the Gaussian kernel to generate the permutations.
deltaList <- rep(list(c(0.01, 0.05, seq(0.1, 0.9, .1))), length(rownames(output_affine$source)))

#shortList <- list(target = output_affine$target[1:5,], source = output_affine$source[1:5,])
start_time <- Sys.time()
merfishCorrelation_affine <- spatialCorrelationGeneExpIterPermutations(
  output_affine, 
  deltaX = deltaList, deltaY = deltaList,
  nThreads = 20,
  BPPARAM = BiocParallel::MulticoreParam()
)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 7.653494 hours

save(merfishCorrelation_affine, file = "~/ST_compare/data/merfish_data/merfish_correlation_affine_delta_0_01_0_9_BH_Iter1000_20260628.RData")
load(file = "~/ST_compare/data/merfish_data/merfish_correlation_affine_delta_0_01_0_9_BH_Iter1000_20260628.RData")
save(merfishCorrelation_affine, file = "~/github/STcompare/inst/extdata/merfishCorrelation_affine.RData")
load(file = "~/github/STcompare/inst/extdata/merfishCorrelation_affine.RData")


## Results for SVGs ####

#number of genes
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>%
  dim()

#number of genes as svgs
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dim() #[1] 415  11

# number of genes as svgs that are significantly positively correlated
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  dim()

# number of genes as svgs that are significantly negatively correlated
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef < 0) %>% 
  dim()

# number of genes as svgs that are not significantly  correlated
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  dim()

# names of genes as svgs that are significantly positively correlated
svgSigPos <- merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  rownames()

# names of genes as svgs that are not significantly  correlated
svgNotSig <- merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  rownames()

# In total, we identified 93% (384/415) of genes as significantly positively correlated and zero genes as significantly negative correlated with an adjusted empirical p < 0.05 (Figure 13d), suggesting that the most SVGs have not changed in their spatial patterning, as expected. 

## Results for non-SVGs ####

#number of genes classified as non svgs 
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dim()

# number of genes as non svgs that are significantly positively correlated
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  dim()

# number of genes as non svgs that are not significantly  correlated
merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  dim()

# names of genes as non svgs that are significantly positively correlated
nonSVGsig <- merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
  rownames()

# names of genes as non svgs that are not significantly correlated
nonSVGnotsig <- merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  rownames()

# When performing STcompare’s spatial correlation test on the remaining 68 genes not identified as SVG, 81% (55/68) were identified as not significantly positively correlated using Benjamini-Hochberg adjusted empirical p < 0.05. 
# The high percentage of non-SVGs identified as not significantly positively correlated is consistent with the expectation that genes not exhibiting autocorrelated spatial expression patterns are less likely to have correlated patterns across replicates. 

# Here we visualize the correlation coefficient and empirical p-value for each gene, colored by whether the gene is an SVG or not. The dashed line indicates the threshold for significance at p < 0.05.
plt <- merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::mutate(SVG = case_when(rownamesCol %in% non_svg ~ "not",
                                rownamesCol %in% svg_int ~ "svg",
                         .default = "Blank")) %>%
  dplyr::filter(SVG %in% c("svg", "not")) %>% 
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x= correlationCoef, y = -log10(pValueEmpiricalRound), color = SVG), alpha = 0.25, size = 3) +
  ggplot2::xlim(NA,1) +
  ggplot2::scale_color_manual(values = c("svg" = "green", "not" = "blue")) +
  ggplot2::theme_classic() +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "black")  + 
  ggplot2::labs(x = "Correlation Coefficient, r" , y = "-log10(p-value)", 
       title = "Rate of significance for SVGs versus non-SVGs"
       )

plt

##### Use STcompare to calculate fold-change similarity score #####
similarityMerfish_affine <- spatialSimilarity(output_affine, foldChange = 1)

# Compare results for the affine aligned source coordinates to demonstrate that the results are dependent on the alignment method used.
pValueAffine <- merfishCorrelation_affine %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation_affine)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::mutate(SVG = case_when(rownamesCol %in% non_svg ~ "not",
                                rownamesCol %in% svg_int ~ "svg",
                         .default = "Blank")) %>%
  dplyr::filter(SVG %in% c("svg", "not")) %>% 
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  dplyr::pull(pValueEmpiricalRound)

pValue <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::mutate(SVG = case_when(rownamesCol %in% non_svg ~ "not",
                                rownamesCol %in% svg_int ~ "svg",
                         .default = "Blank")) %>%
  dplyr::filter(SVG %in% c("svg", "not")) %>% 
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  dplyr::pull(pValueEmpiricalRound)


df <- data.frame(gene = rownames(merfishCorrelation_affine),
                 correlation = merfishCorrelation[, "correlationCoef" ],
                 correlation_affine = merfishCorrelation_affine[, "correlationCoef" ],
                 pValue = pValue,
                 pValueAffine = pValueAffine,
                 similar = similarityMerfish$similarityTable$percentSimilarity[similarityMerfish$similarityTable$gene %in% genes_only],
                 similarAffine = similarityMerfish_affine$similarityTable$percentSimilarity[similarityMerfish_affine$similarityTable$gene %in% genes_only]
                 )
head(df)

df <- df %>%
  dplyr::mutate(SVG = case_when(gene %in% non_svg ~ "not",
                                gene %in% svg_int ~ "svg",
                         .default = "Blank")) %>%
  dplyr::filter(SVG %in% c("svg", "not")) %>%
  dplyr::mutate(SVG = factor(SVG, levels = c("svg", "not"))) %>%
  dplyr::mutate(correlationDiff = correlation_affine - correlation) %>%
  dplyr::mutate(similarityDiff = similarAffine - similar) %>%
  dplyr::mutate(pValueOutcome = case_when(pValue < 0.05 & pValueAffine < 0.05 ~ "both",
                                          pValue < 0.05 & pValueAffine >= 0.05 ~ "pValue",
                                          pValue >= 0.05 & pValueAffine < 0.05 ~ "pValueAffine",
                                          .default = "neither")) %>%
  dplyr::mutate(pValueOutcome = factor(pValueOutcome, levels = c("both", "pValue", "pValueAffine", "neither")))
  

#Comparison of correlation coefficients for SVGs versus non-SVGs
pltCorrelation <- df %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x= correlation, y = correlation_affine, color = SVG), alpha = 0.25, size = 3) +
  ggplot2::xlim(0,1) +
  ggplot2::ylim(0,1) +
  ggplot2::coord_fixed() +
  ggplot2::scale_color_manual(values = c("svg" = "green", "not" = "blue")) +
  ggplot2::theme_classic() +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = "black")  +
  ggplot2::labs(x = "Correlation Coefficient, r" , y = "Correlation Coefficient, r (affine)", 
       title = "Comparison of correlation coefficients for SVGs versus non-SVGs"
       )


#correlation coefficients difference for SVGs versus non-SVGs with p-value outcomes
pltCorrelationDiff <- df %>%
  ggplot2::ggplot(ggplot2::aes(x= SVG, y = correlationDiff)) + 
  ggplot2::geom_violin(alpha = 0.25, size = 1) +
  ggplot2::geom_jitter(ggplot2::aes(x= SVG, y = correlationDiff, color = pValueOutcome), width = 0.2, alpha = 0.5, size = 3) +
  #scale_color_manual(values = c("both" = "black", "pValue" = "red", "pValueAffine" = "blue", "neither" = "gray")) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[1:4]) +
  #viridis::scale_color_viridis(discrete = TRUE, name = "p-value outcome", option = "rocket") +
  #ggthemes::scale_color_colorblind()
  ggplot2::theme_classic() +
  ggplot2::labs(x = "" , y = "correlation coefficient difference, r_affine - r", 
       title = "Differences in correlation coefficients for affine vs STalign for SVGs versus non-SVGs"
       )

LegendCorrDiff <- gtable::gtable_filter(ggplot2::ggplotGrob(pltCorrelationDiff), "guide-box")


#similarity difference for SVGs versus non-SVGs with p-value outcomes
df %>%
  ggplot2::ggplot(ggplot2::aes(x= SVG, y = similarityDiff)) + 
  ggplot2::geom_violin(alpha = 0.25, size = 1) +
  ggplot2::geom_jitter(ggplot2::aes(x= SVG, y = similarityDiff, color = pValueOutcome), width = 0.2, alpha = 0.25) +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "" , y = "similarity difference, S_affine - S", 
       title = "Differences in similarity scores for affine vs STalign for SVGs versus non-SVGs"
       )

# Comparison of similarity scores for SVGs versus non-SVGs
pltSimilarity <- df %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x= similar, y = similarAffine, color = SVG), alpha = 0.25, size = 3) +
  ggplot2::xlim(0,1) +
  ggplot2::ylim(0,1) +
  ggplot2::coord_fixed() +
  ggplot2::scale_color_manual(values = c("svg" = "green", "not" = "blue")) +
  ggplot2::theme_classic() +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = "black")  +
  ggplot2::labs(x = "Fold-Change Similarity, S" , y = "Fold-Change Similarity, S (affine)", 
       title = "Comparison of similarity scores for SVGs versus non-SVGs"
       )

#multi-panel figure of the spatial expression of genes with high and low similarity scores across the two replicates, scatterplots of the gene's expression across the two replicates with points colored by similarity classification, and spatial pattern of similarity score.
gridExtra::grid.arrange(gridExtra::arrangeGrob(pltTarget,
                                               pltSourceAffine,
                                               pltBoth,
                                               ncol=3, nrow = 1, widths = c(5,5,5)))

file_name = paste0("~/ST_compare/plots/affine_comparison_corr_sim_20260630.pdf")
pdf(file = file_name, width = 14, height = 6, onefile=FALSE)                                             
gridExtra::grid.arrange(gridExtra::arrangeGrob(pltCorrelation + ggplot2::theme(legend.position="none"),
                                               pltCorrelationDiff + ggplot2::theme(legend.position="none"),
                                               LegendCorrDiff,
                                               pltSimilarity + ggplot2::theme(legend.position="none"),
                                               ncol=4, nrow = 1, widths = c(5,8,1,5))
                        )
dev.off()


correlationDiffSVG <- df %>%
  dplyr::filter(SVG == "svg") %>%
  dplyr::pull(correlationDiff)

correlationDiffnonSVG <- df %>%
  dplyr::filter(SVG == "not") %>%
  dplyr::pull(correlationDiff)

t.test(correlationDiffSVG, correlationDiffnonSVG)


similarityDiffSVG <- df %>%
  dplyr::filter(SVG == "svg") %>%
  dplyr::pull(similarityDiff)

similarityDiffnonSVG <- df %>%
  dplyr::filter(SVG == "not") %>%
  dplyr::pull(similarityDiff)

t.test(similarityDiffSVG, similarityDiffnonSVG)
wilcox.test(similarityDiffSVG, similarityDiffnonSVG)

# return a gene with a large difference in correlation coefficient between the two alignment methods and a significant p-value for the original alignment method but not for the affine alignment method.
df %>%
  dplyr::filter(pValueOutcome == "pValue") %>%
  dplyr::filter(SVG == "svg") %>%
  dplyr::arrange(correlationDiff) %>%
  head()


gene <- "Drd4"
gene <- "Gabbr1"

# create a color scale that is consistent across the two replicates for the same gene

sharedPixels<- intersect(rownames(SpatialExperiment::spatialCoords(output[[1]])),
                            rownames(SpatialExperiment::spatialCoords(output[[2]])))


minList <- min(c(SummarizedExperiment::assay(output[[1]])[gene, sharedPixels], SummarizedExperiment::assay(output[[2]])[gene, sharedPixels]))
maxList <- max(c(SummarizedExperiment::assay(output[[1]])[gene, sharedPixels], SummarizedExperiment::assay(output[[2]])[gene, sharedPixels]))

sc <- ggplot2::scale_fill_gradientn(colors=viridis::viridis(20), 
                           limits = c(minList, maxList),
                           oob=scales::squish,
                           name = gene)

# plot the spatial expression of the gene across the two replicates
pltS2R2 <- SEraster::plotRaster(output[[1]][gene, sharedPixels], name = gene) + 
  sc + 
  ggplot2::theme_void()
pltS2R3 <- SEraster::plotRaster(output[[2]][gene, sharedPixels], name = gene) + 
  sc + 
  ggplot2::theme_void() 

#grab repeating legends from plots
LegendS2R2 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltS2R2), "guide-box")
LegendS2R3 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltS2R3), "guide-box")


sharedPixels_affine <- intersect(rownames(SpatialExperiment::spatialCoords(output_affine[[1]])),
                            rownames(SpatialExperiment::spatialCoords(output_affine[[2]])))


minList <- min(c(SummarizedExperiment::assay(output_affine[[1]])[gene, sharedPixels_affine], SummarizedExperiment::assay(output_affine[[2]])[gene, sharedPixels_affine]))
maxList <- max(c(SummarizedExperiment::assay(output_affine[[1]])[gene, sharedPixels_affine], SummarizedExperiment::assay(output_affine[[2]])[gene, sharedPixels_affine]))

sc <- ggplot2::scale_fill_gradientn(colors=viridis::viridis(20), 
                           limits = c(minList, maxList),
                           oob=scales::squish,
                           name = gene)

# plot the spatial expression of the gene across the two replicates
pltS2R2Affine <- SEraster::plotRaster(output_affine[[1]][gene, sharedPixels_affine], name = gene) + 
  sc + 
  ggplot2::theme_void()
pltS2R3Affine <- SEraster::plotRaster(output_affine[[2]][gene, sharedPixels_affine], name = gene) + 
  sc + 
  ggplot2::theme_void() 

#grab repeating legends from plots
LegendAffineS2R2 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltS2R2Affine), "guide-box")
LegendAffineS2R3 <- gtable::gtable_filter(ggplot2::ggplotGrob(pltS2R3Affine), "guide-box")


# plot scatterplot of the gene's expression across the two replicates with points colored by similarity classification
simPlot <- linearRegression(similarityMerfish, gene) + ggplot2::theme_classic() + ggplot2::coord_fixed()
simPlotAffine <- linearRegression(similarityMerfish_affine, gene) + ggplot2::theme_classic() + ggplot2::coord_fixed()


# plot spatial pattern of similarity score
classPlot <- pixelClass(similarityMerfish, gene)
classPlotAffine <- pixelClass(similarityMerfish_affine, gene)

#multi-panel figure of the spatial expression of genes with high and low similarity scores across the two replicates, scatterplots of the gene's expression across the two replicates with points colored by similarity classification, and spatial pattern of similarity score.
file_name = paste0("~/ST_compare/plots/affine_example_Gabbr1_20260630.pdf")
pdf(file = file_name, width = 14, height = 6, onefile=FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(pltS2R2 + ggplot2::theme(legend.position="none"),
                                               pltS2R3 + ggplot2::theme(legend.position="none"),
                                               LegendS2R2,
                                               simPlot, classPlot,
                                               pltS2R2Affine + ggplot2::theme(legend.position="none"),
                                               pltS2R3Affine + ggplot2::theme(legend.position="none"),
                                               LegendAffineS2R2,
                                               simPlotAffine, classPlotAffine, 
                                               ncol=5, nrow = 2, widths = c(3,3,1,3,3))
                        )
dev.off()

