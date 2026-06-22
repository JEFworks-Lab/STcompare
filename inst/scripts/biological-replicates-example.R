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

# To account for the smaller structures observed in the brain we added two additional lower values of 0.01 and 0.05 to the sequence of values tested to find the optimal delta for the Gaussian kernel to generate the permutations.
deltaList <- rep(list(c(0.01, 0.05, seq(0.1, 0.9, .1))), length(rownames(output$source)))

shortList <- list(target = output$target[1:5,], source = output$source[1:5,])
start_time <- Sys.time()
merfishCorrelation <- spatialCorrelationGeneExpIterPermutations(
  shortList, 
  deltaX = deltaList, deltaY = deltaList,
  nThreads = 5,
  BPPARAM = BiocParallel::MulticoreParam()
)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 17.40188 mins

## Results for SVGs ####

#number of genes as svgs
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dim()

# number of genes as svgs that are significantly positively correlated
merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::filter(rownamesCol %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%
  dplyr::filter(correlationCoef > 0) %>% 
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

# In total, we identified 86% (358/415) of genes as significantly positively correlated and zero genes as significantly negative correlated with an adjusted empirical p < 0.05 (Figure 13d), suggesting that the most SVGs have not changed in their spatial patterning, as expected. 


#The 14% (57/415) of SVGs without significant positive correlation may reflect low expression 

p1 <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::mutate(SVG = case_when(rownames(merfishCorrelation) %in% non_svg ~ "not",
                         rownames(merfishCorrelation) %in% svg_int ~ "svg",
                         .default = "")) %>%
  dplyr::mutate(Blank = case_when(!grepl("Blank", rownamesCol) ~ "Genes",
                                  grepl("Blank", rownamesCol) ~ "Blanks",
                                .default = "")) %>%
  dplyr::mutate(Blank = factor(Blank, levels = c("Genes", "Blanks"))) %>%
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::mutate(SIG = case_when(pValuePermuteX < 0.05 & pValuePermuteY < 0.05 ~ "sig",
                                pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05 ~ "not",
                                .default = "")) %>%
  dplyr::mutate(meanGexpS2R2 = Matrix::rowMeans(SummarizedExperiment::assay(output[[1]]))[names(Matrix::rowMeans(SummarizedExperiment::assay(output[[1]]))) %in% rownamesCol]) %>%
  dplyr::mutate(meanGexpS2R3 = Matrix::rowMeans(SummarizedExperiment::assay(output[[2]]))[names(Matrix::rowMeans(SummarizedExperiment::assay(output[[2]]))) %in% rownamesCol]) %>%
  dplyr::mutate(names = case_when(SIG %in% c("not") & Blank %in% c("Genes") ~ rownamesCol,
                                  SIG %in% c("sig") ~ "",
                                  .default = "")) %>%
  dplyr::mutate(names = case_when(rownamesCol %in% c("Gpr17", 
                                                     #"Agtr1a", 
                                                     #"Tas1r3", 
                                                     "Lhcgr", 
                                                     #"Sucnr1", 
                                                     #"Gnrhr", 
                                                     #"Rxfp4", 
                                                     "Ccr1", 
                                                     "Taar1") ~ rownamesCol,
                                  .default = "")) %>%
  dplyr::filter(SVG %in% c("svg")) %>% 
  ggplot2::ggplot(ggplot2::aes(x=log10(meanGexpS2R2) , y = correlationCoef, color = SIG, 
                      label = names,
                      #shape = Blank
                      )) + 
    ggplot2::geom_point(alpha = 1, size = 1) +
    ggplot2::scale_color_manual(values = c("sig" = "green", "not" = "blue")) +
    ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
    #viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "log10(mean gene expression in S2R2)"  , y =  "Correlation", 
                  title = "SVGs expression in S2R2"
    )

p2 <- merfishCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(merfishCorrelation)) %>% 
  dplyr::mutate(SVG = case_when(rownames(merfishCorrelation) %in% non_svg ~ "not",
                                rownames(merfishCorrelation) %in% svg_int ~ "svg",
                                .default = "")) %>%
  dplyr::mutate(Blank = case_when(!grepl("Blank", rownamesCol) ~ "Genes",
                                  grepl("Blank", rownamesCol) ~ "Blanks",
                                  .default = "")) %>%
  dplyr::mutate(Blank = factor(Blank, levels = c("Genes", "Blanks"))) %>%
  dplyr::filter(!grepl("Blank", rownamesCol)) %>% 
  dplyr::mutate(SIG = case_when(pValuePermuteX < 0.05 & pValuePermuteY < 0.05 ~ "sig",
                                pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05 ~ "not",
                                .default = "")) %>%
  dplyr::mutate(meanGexpS2R2 = Matrix::rowMeans(SummarizedExperiment::assay(output[[1]]))[names(Matrix::rowMeans(SummarizedExperiment::assay(output[[1]]))) %in% rownamesCol]) %>%
  dplyr::mutate(meanGexpS2R3 = Matrix::rowMeans(SummarizedExperiment::assay(output[[2]]))[names(Matrix::rowMeans(SummarizedExperiment::assay(output[[2]]))) %in% rownamesCol]) %>%
  dplyr::mutate(names = case_when(SIG %in% c("not") & Blank %in% c("Genes") ~ rownamesCol,
                                  SIG %in% c("sig") ~ "",
                                  .default = "")) %>%
  dplyr::mutate(names = case_when(rownamesCol %in% c("Gpr17", 
                                                     #"Agtr1a", 
                                                     #"Tas1r3", 
                                                     "Lhcgr", 
                                                     #"Sucnr1", 
                                                     #"Gnrhr", 
                                                     #"Rxfp4", 
                                                     "Ccr1", 
                                                     "Taar1") ~ rownamesCol,
                                  .default = "")) %>%
  dplyr::filter(SVG %in% c("svg")) %>% 
  ggplot2::ggplot(ggplot2::aes(x=log10(meanGexpS2R3) , y = correlationCoef, color = SIG, 
                      label = names,
                      #shape = Blank
                      )) + 
  ggplot2::geom_point(alpha = 1, size = 1) +
  ggplot2::scale_color_manual(values = c("sig" = "green", "not" = "blue")) +
  ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  #viridis::scale_color_viridis(discrete = TRUE) +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "log10(mean gene expression in S2R2)"  , y =  "Correlation", 
                title = "SVGs expression in S2R3"
  )

p1 + p2


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

# When performing STcompare’s spatial correlation test on the remaining 68 genes not identified as SVG, 94% (64/68) were identified as not significantly positively correlated using an Holm’s adjusted empirical p < 0.05. 
# The high percentage of non-SVGs identified as not significantly positively correlated is consistent with the expectation that genes not exhibiting autocorrelated spatial expression patterns are less likely to have correlated patterns across replicates. 