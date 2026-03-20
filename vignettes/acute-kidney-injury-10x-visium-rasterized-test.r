library(STcompare)
library(SpatialExperiment)
library(MERINGUE)
library(rhdf5)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)
library(BiocGenerics)

library(scatterbar)


# The Visium mouse kidney datasets can be found Zenodo: https://doi.org/10.5281/zenodo.17676991
# IL3 - ischemia acute kidney injury dataset 
# NL3 - control dataset 

# AKI counts ###################################################################
zenodo_url <- "https://zenodo.org/records/19074288/files/IL3_filtered_feature_bc_matrix.h5?download=1"
temp_h5 <- tempfile(fileext = ".h5")
download.file(zenodo_url, temp_h5, mode = "wb")

# this is a list of the internal groups and datasets stored in the HDF5 file 
# we are interested in the "matrix/barcodes" dataset 
rhdf5::h5ls(temp_h5)

# Read the CSC-encoded matrix components from the HDF5 file
file_barcodes <- as.character(rhdf5::h5read(temp_h5, "matrix/barcodes"))
barcodes_matrix <- rhdf5::h5read(temp_h5, "matrix") # stores the sparse matrix components 

# convert CSC encoding into an R dgCMatrix
aki_counts <- Matrix::sparseMatrix(
  dims = barcodes_matrix$shape,
  i = as.numeric(barcodes_matrix$indices),
  p = as.numeric(barcodes_matrix$indptr),
  x = as.numeric(barcodes_matrix$data),
  index1 = FALSE
)
colnames(aki_counts) <- file_barcodes
rownames(aki_counts) <- barcodes_matrix[["features"]]$name

unlink(temp_h5)
head(aki_counts)

# Control counts ###############################################################

zenodo_url <- "https://zenodo.org/records/19074288/files/NL3_filtered_feature_bc_matrix.h5?download=1"
temp_h5 <- tempfile(fileext = ".h5")
download.file(zenodo_url, temp_h5, mode = "wb")

# Read the CSC-encoded matrix components from the HDF5 file
file_barcodes <- as.character(rhdf5::h5read(temp_h5, "matrix/barcodes"))
barcodes_matrix <- rhdf5::h5read(temp_h5, "matrix") # stores the sparse matrix components 

# convert CSC encoding into an R dgCMatrix
control_counts <- Matrix::sparseMatrix(
  dims = barcodes_matrix$shape,
  i = as.numeric(barcodes_matrix$indices),
  p = as.numeric(barcodes_matrix$indptr),
  x = as.numeric(barcodes_matrix$data),
  index1 = FALSE
)
colnames(control_counts) <- file_barcodes
rownames(control_counts) <- barcodes_matrix[["features"]]$name

unlink(temp_h5)
head(control_counts)

# positions ####################################################################
zenodo_url <- "https://zenodo.org/records/19074288/files/IL3_tissue_positions.csv?download=1"
aki_pos <- read.csv(zenodo_url, header = TRUE, stringsAsFactors = FALSE)

zenodo_url <- "https://zenodo.org/records/19074288/files/NL3_tissue_positions.csv?download=1"
ctrl_pos <- read.csv(zenodo_url, header = TRUE, stringsAsFactors = FALSE)

head(ctrl_pos)
head(aki_pos)

# number of columns in counts is not the same number of rows in pos
dim(aki_counts)
dim(aki_pos)

dim(control_counts)
dim(ctrl_pos)

# limit the counts to the spots with positions and make sure they are in the same order
aki_counts <- aki_counts[, colnames(aki_counts) %in% aki_pos$barcode]
control_counts <- control_counts[, colnames(control_counts) %in% ctrl_pos$barcode]

# number of columns in counts is the same number of rows in pos
dim(aki_counts)
dim(aki_pos)

dim(control_counts)
dim(ctrl_pos)

# plot positions ###############################################################

# set the barcode column as the rownames
rownames(ctrl_pos) <- ctrl_pos$barcode
ctrl_pos <- ctrl_pos[,-1]
colnames(ctrl_pos) <- c("x", "y")

rownames(aki_pos) <- aki_pos$barcode
aki_pos <- aki_pos[,-1]
colnames(aki_pos) <- c("x", "y")

ctrl_pos$group <- "Control"
aki_pos$group  <- "AKI"

# visualize both control and AKI positions on the same plot
df <- rbind(ctrl_pos, aki_pos)
pos_plot <- ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point() +
  scale_color_manual(values = c("Control" = "blue", "AKI" = "red")) +
  coord_fixed() +
  labs(x = "x", y = "y", color = "Group") +
  theme_classic()
pos_plot

# rotate the tissue by 90-degrees for consistency with the paper 
ctrl_pos_rot <- transform(ctrl_pos, x = y, y = -x + max(ctrl_pos$x))
aki_pos_rot <- transform(aki_pos, x = y, y = -x + max(aki_pos$x))

df_rot <- rbind(ctrl_pos_rot, aki_pos_rot)
pos_rot_plot <- ggplot(df_rot, aes(x = x, y = y, color = group)) +
  geom_point(size=0.5) +
  scale_color_manual(values = c("Control" = "blue", "AKI" = "red")) +
  coord_fixed() +
  labs(x = "x", y = "y", color = "Group") +
  theme_classic()
pos_rot_plot

# center align X coordinates
ctrl_pos_aligned <- ctrl_pos_rot
ctrl_pos_aligned[,1] <- ctrl_pos_aligned[,1] - mean(ctrl_pos_aligned[,1])
aki_pos_aligned <- aki_pos_rot
aki_pos_aligned[,1] <- aki_pos_aligned[,1] - mean(aki_pos_aligned[,1])
# simply stretch Y positions to align
aki_pos_aligned[,2] <- aki_pos_aligned[,2]*max(ctrl_pos_aligned[,2])/max(aki_pos_aligned[,2])

# visualize to assess alignment
df_aligned <- rbind(ctrl_pos_aligned, aki_pos_aligned)
pos_aligned_plot <- ggplot(df_aligned, aes(x = x, y = y, color = group)) +
  geom_point(size=0.5) +
  scale_color_manual(values = c("Control" = "blue", "AKI" = "red")) +
  coord_fixed() +
  labs(x = "x", y = "y", color = "Group") +
  theme_classic()
pos_aligned_plot

AKI_ctrl_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = control_counts),
  spatialCoords = as.matrix(ctrl_pos_aligned[,1:2]),
)

AKI_aki_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = aki_counts),
  spatialCoords = as.matrix(aki_pos_aligned[,1:2])
)

# the region labels were generated in the paper using a clustering approach and are used here to ensure that we are comparing shared regions across the two samples
load(file = "~/ST_compare/data/kidney_data/region_labels_20241106.RData")

# rename names of region labels so can used shared spot locations
ctrl_region <- new.com[grepl("CTRL", names(new.com))]
names(ctrl_region) <- sub("^[^_]*_", "", names(ctrl_region))

aki_region <- new.com[grepl("IR", names(new.com))]
names(aki_region) <- sub("^[^_]*_", "", names(aki_region))

# save the cluster labels in spatial experiments
AKI_ctrl_SE$region <- ctrl_region[colnames(AKI_ctrl_SE)]
AKI_aki_SE$region <- aki_region[colnames(AKI_aki_SE)]

#save the cluster labels as one-hot encoded matrices
ctrl_region_onehot <- model.matrix(~ region - 1, data = colData(AKI_ctrl_SE))
aki_region_onehot <- model.matrix(~ region - 1, data = colData(AKI_aki_SE))
colnames(ctrl_region_onehot) <- gsub("region", "", colnames(ctrl_region_onehot))
colnames(aki_region_onehot) <- gsub("region", "", colnames(aki_region_onehot))

#append spatial coordinates to the one-hot encoded matrices to align
ctrl_region_onehot <- cbind(ctrl_pos, ctrl_region_onehot)
aki_region_onehot <- cbind(aki_pos, aki_region_onehot)

#save one-hot encoded matrices as csv files to use as input for the region-based alignment methods
write.csv(ctrl_region_onehot, file = "~/ST_compare/data/kidney_data/ctrl_region_onehot.csv", row.names = TRUE)
write.csv(aki_region_onehot, file = "~/ST_compare/data/kidney_data/aki_region_onehot.csv", row.names = TRUE)


load(file = "~/ST_compare/data/kidney_data/harmonized_clusters_20240626.RData")

# rename names of region labels so can used shared spot locations
ctrl_emb <- emb.harmony[grepl("CTRL", rownames(emb.harmony)),]
rownames(ctrl_emb) <- sub("^[^_]*_", "", rownames(ctrl_emb))
colnames(ctrl_emb) <- c("tSNE1", "tSNE2")

aki_emb <- emb.harmony[grepl("IR", rownames(emb.harmony)),]
rownames(aki_emb) <- sub("^[^_]*_", "", rownames(aki_emb))
colnames(aki_emb) <- c("tSNE1", "tSNE2")

df1 <- data.frame(SpatialExperiment::spatialCoords(AKI_ctrl_SE), dataset = "ctrl", region = AKI_ctrl_SE$region, ctrl_emb)
df2 <- data.frame(SpatialExperiment::spatialCoords(AKI_aki_SE), dataset = "aki", region = AKI_aki_SE$region, aki_emb)
df5 <- rbind(df1, df2)

# spatial plot of control with region cluster labels
pctrl <- ggplot2::ggplot(df1, ggplot2::aes(x = x, y = y, color = region)) +
  ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void()

# spatial plot of ischemia with region cluster labels
paki <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y, color = region)) +
  ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void()

# tSNE plot with region cluster labels
ptSNEregion <- ggplot2::ggplot(df5, ggplot2::aes(x = tSNE1, y = tSNE2, color = region)) +
  ggplot2::geom_point(size = 0.4) +
  ggplot2::coord_fixed() +
  ggplot2::theme_classic()
ptSNEregion

# tSNE plot with dataset labels to show harmonization
ptSNEdataset <- ggplot2::ggplot(df5, ggplot2::aes(x = tSNE1, y = tSNE2, color = dataset)) +
  ggplot2::geom_point(size = 0.4) +
  ggplot2::coord_fixed() +
  ggplot2::theme_classic()
ptSNEdataset

# rasterize the gene expression for both samples using SEraster
input <- list(AKI_ctrl = AKI_ctrl_SE, AKI_aki = AKI_aki_SE)
rast <- SEraster::rasterizeGeneExpression(
  input,
  resolution = 5,
  fun        = "sum",
  square     = FALSE,
  assay_name = 'counts'
)

# add a CPM normalization assay to the rasterized SpatialExperiment
assay(rast$AKI_ctrl, "CPM") <- Matrix::t(Matrix::t(assay(rast$AKI_ctrl))/Matrix::colSums(assay(rast$AKI_ctrl)))*1e6
assay(rast$AKI_aki, "CPM") <- Matrix::t(Matrix::t(assay(rast$AKI_aki))/Matrix::colSums(assay(rast$AKI_aki)))*1e6

# plot total expression for new pixels
p1 <- SEraster::plotRaster(rast$AKI_ctrl, assay_name = "CPM", name = "total expression")
p2 <- SEraster::plotRaster(rast$AKI_aki, assay_name = "CPM", name = "total expression")
p1 + p2

# rasterize the region labels for both samples using rasterizeCellType
rast_clusters <- SEraster::rasterizeCellType(input, 
  col_name = "region", 
  resolution = 5, 
  fun = "mean", 
  square = FALSE,
  BPPARAM = BiocParallel::MulticoreParam())

#check that plots look correct and that same coordinate system created from rasterizeGeneExpression and rasterizeCellType
# plot spatial patterns of regions proportion for both samples
sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(rast_clusters$AKI_ctrl)),
                          rownames(SpatialExperiment::spatialCoords(rast_clusters$AKI_aki)))
p3 <- SEraster::plotRaster(rast_clusters$AKI_ctrl[,sharedPixels], feature_name = "cortex", name = "region")
p4 <- SEraster::plotRaster(rast_clusters$AKI_aki[,sharedPixels], feature_name = "cortex", name = "region")
p3 + p4

p5 <- SEraster::plotRaster(rast_clusters$AKI_ctrl[,sharedPixels], feature_name = "interface", name = "region")
p6 <- SEraster::plotRaster(rast_clusters$AKI_aki[,sharedPixels], feature_name = "interface", name = "region")
p5 + p6

p7 <- SEraster::plotRaster(rast_clusters$AKI_ctrl[,sharedPixels], feature_name = "medulla", name = "region")
p8 <- SEraster::plotRaster(rast_clusters$AKI_aki[,sharedPixels], feature_name = "medulla", name = "region")
p7 + p8

#use scatterbar to visualize the proportion of each region in each sample
dfproportions_ctrl <- data.frame(t(as.matrix(SummarizedExperiment::assay(rast_clusters$AKI_ctrl)[,sharedPixels])))
dfcoords_ctrl <- data.frame(SpatialExperiment::spatialCoords(rast_clusters$AKI_ctrl)[sharedPixels,])

p9 <- scatterbar::scatterbar(data = dfproportions_ctrl, 
                              xy = dfcoords_ctrl, 
                              size_x = 4, size_y = 4, 
                              padding_x = 0.01, padding_y = 0.01) + coord_fixed()

dfproportions_aki <- data.frame(t(as.matrix(SummarizedExperiment::assay(rast_clusters$AKI_aki)[,sharedPixels])))
dfcoords_aki <- data.frame(SpatialExperiment::spatialCoords(rast_clusters$AKI_aki)[sharedPixels,])

p10 <- scatterbar::scatterbar(data = dfproportions_aki, xy = dfcoords_aki, size_x = 4, size_y = 4, padding_x = 0.01, padding_y = 0.01) + coord_fixed()
p9 + p10

# Returns a dataframe with these columns: 
# observed: Observed Moran's I statistic measuring spatial autocorrelation
# expected: Expected Moran's I under null hypothesis of random distribution
# sd: Standard deviation of Moran's I under the null hypothesis
# p.value: Statistical significance of the spatial pattern
# p.adj: Adjusted p-value (FDR-corrected) for multiple testing correction
moransI <- function(SE_temp, filterDist, assayName=1){
  
  # use MERINGUE to get the neighbor-relationships 
  w <- MERINGUE::getSpatialNeighbors(SpatialExperiment::spatialCoords(SE_temp), filterDist = filterDist)
  par(mfrow=c(1,1))
  MERINGUE::plotNetwork(SpatialExperiment::spatialCoords(SE_temp), w)
  
  # Identify significantly spatially auto-correlated genes
  I <- MERINGUE::getSpatialPatterns(SummarizedExperiment::assays(SE_temp)[[assayName]], w)
  
  return(I)
}

moransI_ctrl <- moransI(rast$AKI_ctrl, filterDist = 10, assayName = "CPM")
head(moransI_ctrl)

moransI_aki <- moransI(rast$AKI_aki, filterDist = 10, assayName = "CPM")
head(moransI_aki)

# Filter for the genes that are statistically significantly autocorrelated 
svg_ctrl <- moransI_ctrl %>% filter(p.adj == 0) %>% rownames()
svg_aki <- moransI_aki %>% filter(p.adj == 0) %>% rownames()
svg_int <- intersect(svg_ctrl, svg_aki) # genes that are SVGs in both conditions 
svg_uni <- union(svg_ctrl, svg_aki) # genes that are SVG in one condition but not the other 

# Out of 13638 genes: 
# there are 1498 genes that are SVGs in the control dataset 
length(svg_ctrl)

# there are 2308 genes that are SVGs in the aki dataset  
length(svg_aki)

# there are 943 genes that are SVGs in both datasets   
length(svg_int)

# there are 2863 genes that are SVGs in either datasets 
length(svg_uni)

#return names in `svg_int` that are not in `rownames(rast$AKI_ctrl)`
setdiff(svg_int, rownames(rast$AKI_ctrl))

# remove leading X only if followed by a number
svg_int <- sub("^X(?=[0-9])", "", svg_int, perl = TRUE)
# replace . with -
svg_int <- gsub("\\.", "-", svg_int)

#check if any differences remain
setdiff(svg_int, rownames(rast$AKI_ctrl))

# there are 32285 genes in both samples 
# which rows in the matrix that the gene is expressed in greater than 5% of all spots 
# this has to be true in both datasets to be considered a "good gene" 
good.genes <- names(which(
  (Matrix::rowSums(SummarizedExperiment::assay(rast$AKI_ctrl, 'CPM') > 0) / ncol(SummarizedExperiment::assay(rast$AKI_ctrl, 'CPM'))*100 > 5) &
    (Matrix::rowSums(SummarizedExperiment::assay(rast$AKI_aki, 'CPM') > 0) / ncol(SummarizedExperiment::assay(rast$AKI_aki, 'CPM'))*100 > 5)))
length(good.genes)

# choosing the genes that are SVGs in both conditions 
genes_chosen <- svg_int

## alternative: choosing the genes that are expressed in greater than 5% of all spots
# genes_chosen <- good.genes

# The input is the spatial experiment for the genes that are SVGs in both conditions 
input <- list('AKI_ctrl'= rast$AKI_ctrl[genes_chosen,],
              'AKI_aki'= rast$AKI_aki[genes_chosen,])

# running spatialCorrelationGeneExp 
set.seed(0)
start_time <- Sys.time()
deltaList <- replicate(length(genes_chosen), c(0.01, 0.05, seq(0.1, 0.9, .1)), simplify = FALSE)
kidneyCorrelation <- STcompare::spatialCorrelationGeneExp(
  input, 
  #nPermutations = c(100,1000), # increased the number of permutations to account for stochasticity 
  deltaX = deltaList, 
  deltaY = deltaList, 
  returnPermutations = FALSE, 
  assayName = "CPM", 
  nThreads = 5,
  verbose = TRUE
)
end_time <- Sys.time()
print(end_time - start_time) 

#save the correlation results as an RData file
save(kidneyCorrelation, file = "~/ST_compare/data/kidney_data/kidneyCorrelation_20260320.RData")
#load the correlation results as an RData file
load("~/ST_compare/data/kidney_data/kidneyCorrelation_20260320.RData")

# look at results
head(kidneyCorrelation)

# saving the significantly positively correlated svg genes 
svgSigPos <- kidneyCorrelation %>%
  dplyr::filter(rownames(kidneyCorrelation) %in% svg_int) %>%
  dplyr::filter(p.adjust(pValuePermuteX) < 0.05 & p.adjust(pValuePermuteY) < 0.05) %>%  
  dplyr::filter(correlationCoef > 0) %>% 
  rownames()
# 360 genes are SVGs that are significantly positively correlated 
length(svgSigPos)

# saving the significantly negatively correlated svg genes 
svgSigNeg <- kidneyCorrelation %>%
  dplyr::filter(rownames(kidneyCorrelation) %in% svg_int) %>%
  dplyr::filter(p.adjust(pValuePermuteX) < 0.05 & p.adjust(pValuePermuteY) < 0.05) %>%  
  dplyr::filter(correlationCoef < 0) %>% 
  rownames()
# 5 genes are SVGs that are significantly negatively correlated 
length(svgSigNeg)

# visualize the significantly positively correlate and significantly negatively correlated svg genes  
# Figure 2J - Rate of Significance for SVGs 
fig_2j <- kidneyCorrelation %>%
  dplyr::mutate(Sig = dplyr::case_when(rownames(kidneyCorrelation) %in% svgSigPos ~ "SigPos",
                                       rownames(kidneyCorrelation) %in% svgSigNeg ~ "SigNeg",
                                       .default = "")) %>%
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  ggplot2::ggplot(ggplot2::aes(x= correlationCoef, 
                               y = -log10(pValueEmpiricalRound), 
                               color = Sig)) + 
  ggplot2::geom_point(alpha = 0.25, size = 1) +
  ggplot2::xlim(NA,1) +
  ggplot2::scale_color_manual(values = c("SigPos" = "green", "SigNeg" = "blue")) +
  ggplot2::theme_classic() +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "black")  + 
  ggplot2::labs(x = "Correlation" , y = "-log10(p-value)", title = "Rate of Significance for SVGs")
fig_2j

# Adding a column for the empirical p-value (pValueEmpirical) 
# defined as the larger of the two x and y permuted p-values (greater of pValuePermuteX and pValuePermuteY)
kidneyCorrelationEmpirical <- kidneyCorrelation %>%
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX))


# kidneyCorrPos <- kidneyCorrelation[svgSigPos,]
# head(kidneyCorrPos[order(kidneyCorrPos$correlationCoef, decreasing=TRUE),])
# gene1 <- 'Slc5a3'
# kidneyCorrNeg <- kidneyCorrelation[svgSigNeg,]
# head(kidneyCorrNeg[order(kidneyCorrNeg$correlationCoef, decreasing=FALSE),])
# gene2 <- 'Ech1'

#filter the correlation dataframe to only the genes that are significantly positively correlated
kidneyCorrPos <- kidneyCorrelation[svgSigPos,]
#arrange this resulting dataframe in descending order for the correlation coefficient 
kidneyCorrPos %>% arrange(desc(correlationCoef))
#save the name of the significantly positively correlated gene with the highest correlation coefficient
gene1 <- kidneyCorrPos %>% arrange(desc(correlationCoef)) %>% rownames() %>% head(1)

#filter the correlation dataframe to only the genes that are significantly negatively correlated
kidneyCorrNeg <- kidneyCorrelation[svgSigNeg,]
#arrange this resulting dataframe in ascending order for the correlation coefficient 
kidneyCorrNeg %>% arrange(correlationCoef)
#save the name of the significantly negatively correlated gene with the lowest correlation coefficient
gene2 <- kidneyCorrNeg %>% arrange(correlationCoef) %>% rownames() %>% head(1)


# plot spatial patterns of example genes
sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(rast$AKI_ctrl)),
                          rownames(SpatialExperiment::spatialCoords(rast$AKI_aki)))

gene1_plt_ctrl  <- SEraster::plotRaster(rast$AKI_ctrl[gene1, sharedPixels], assay_name = "CPM", name = "CPM expression")
gene1_plt_aki <- SEraster::plotRaster(rast$AKI_aki[gene1, sharedPixels], assay_name = "CPM", name = "CPM expression")

gene1_plt_ctrl + gene1_plt_aki

gene2_plt_ctrl  <- SEraster::plotRaster(rast$AKI_ctrl[gene2, sharedPixels], assay_name = "CPM", name = "CPM expression")
gene2_plt_aki <- SEraster::plotRaster(rast$AKI_aki[gene2, sharedPixels], assay_name = "CPM", name = "CPM expression")

gene2_plt_ctrl + gene2_plt_aki

# Figure 2h - example of correlation visualization 
# use the plotCorrelationGeneExp function to generate the correlation plots 

gene1_corr <- plotCorrelationGeneExp(
  speList = input, 
  spatialCorrelation = kidneyCorrelation, 
  geneName = gene1, 
  assayName = "CPM"
)
gene2_corr <- plotCorrelationGeneExp(
  speList = input, 
  spatialCorrelation = kidneyCorrelation, 
  geneName = gene2, 
  assayName = "CPM"
)

gene1_corr

gene1_plt_ctrl | gene1_plt_aki | gene1_corr
gene2_plt_ctrl | gene2_plt_aki | gene2_corr

# Finding the spatial similarity
# If t1 and t2 are null, then the default threshold is the 0.05 quantile of gene expression. 
start_time <- Sys.time()
ss <- spatialSimilarity(input, foldChange = 1, assayName = "CPM", t1 = NULL, t2 = NULL)
end_time <- Sys.time()
print(end_time-start_time)


# Figure 2I - fold change linear regression and spot classification 
gene1_lr <- linearRegression(ss, gene=gene1, assayName = "CPM") + ggplot2::theme_classic() + ggplot2::coord_fixed()
gene2_lr <- linearRegression(ss, gene=gene2, assayName = "CPM") + ggplot2::theme_classic() + ggplot2::coord_fixed()

gene1_pc <- pixelClass(ss, gene=gene1, assayName = "CPM")
gene2_pc <- pixelClass(ss, gene=gene2, assayName = "CPM")

gene1_lr | gene1_pc
gene2_lr | gene2_pc

# 2k - Comparing Correlation Metric to Similiarty metric 

# create a df of correlation and similarity values for each gene 
sim_corr <- data.frame(gene = ss$similarityTable$gene[ss$similarityTable$gene %in% svg_int],
                       correlation = kidneyCorrelation$correlationCoef[rownames(kidneyCorrelation) %in% svg_int], 
                       similar = ss$similarityTable$percentSimilarity[ss$similarityTable$gene %in% svg_int],
                       percentDissimilarityX =ss$similarityTable$percentDissimilarityX[ss$similarityTable$gene %in% svg_int],
                       percentDissimilarityY = ss$similarityTable$percentDissimilarityY[ss$similarityTable$gene %in% svg_int])
sim_corr

sim_corr_plt <- sim_corr %>%
  dplyr::mutate(Sig = case_when(gene %in% svgSigPos ~ "SigPos",
                                gene %in% svgSigNeg ~ "SigNeg",
                                .default = "NotSig")) %>%
  dplyr::filter(Sig %in% c("SigPos", "SigNeg")) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x = correlation, y = similar, color = Sig), alpha = 0.25, size = 1) +
  ggplot2::xlim(NA,1) +
  ggplot2::coord_fixed() +
  ggplot2::scale_color_manual(values = c("SigPos" = "green", "SigNeg" = "blue")) +
  ggplot2::theme_classic() + 
  ggplot2::labs(x = "Correlation" , y = "Similarity", 
                title = "Similarity for significantly positively and negatively correlated SVGs")

sim_corr_plt



# paper figures

# visualize the significantly positively correlate and significantly negatively correlated svg genes  
# Figure 2J - Rate of Significance for SVGs 
fig_2j <- kidneyCorrelation %>%
  dplyr::mutate(Sig = dplyr::case_when(rownames(kidneyCorrelation) %in% svgSigPos ~ "SigPos",
                                       rownames(kidneyCorrelation) %in% svgSigNeg ~ "SigNeg",
                                       .default = "")) %>%
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  ggplot2::ggplot(ggplot2::aes(x= correlationCoef, 
                               y = -log10(pValueEmpiricalRound), 
                               color = Sig)) + 
  ggplot2::geom_point(alpha = 0.25, size = 1) +
  ggplot2::xlim(NA,1) +
  ggplot2::scale_color_manual(values = c("SigPos" = "green", "SigNeg" = "blue")) +
  ggplot2::theme_classic() +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "black")  + 
  ggplot2::labs(x = "Correlation" , y = "-log10(p-value)", title = "Rate of Significance for SVGs")
fig_2j

sim_corr_plt <- sim_corr %>%
  dplyr::mutate(Sig = case_when(gene %in% svgSigPos ~ "SigPos",
                                gene %in% svgSigNeg ~ "SigNeg",
                                .default = "NotSig")) %>%
  dplyr::filter(Sig %in% c("SigPos", "SigNeg")) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_point(ggplot2::aes(x = correlation, y = similar, color = Sig), alpha = 0.25, size = 1) +
  ggplot2::xlim(NA,1) +
  ggplot2::coord_fixed() +
  ggplot2::scale_color_manual(values = c("SigPos" = "green", "SigNeg" = "blue")) +
  ggplot2::theme_classic() + 
  ggplot2::labs(x = "Correlation" , y = "Similarity", 
                title = "Similarity for significantly positively and negatively correlated SVGs")

sim_corr_plt

gene1 <- sim_corr %>%
  dplyr::mutate(Sig = case_when(gene %in% svgSigPos ~ "SigPos",
                                gene %in% svgSigNeg ~ "SigNeg",
                                .default = "NotSig")) %>%
  dplyr::filter(Sig %in% c("SigPos")) %>%
  dplyr::arrange(desc(correlation)) %>%
  dplyr::arrange(similar) %>%
  dplyr::pull(gene) %>%
  head(1)

gene2 <- sim_corr %>%
  dplyr::mutate(Sig = case_when(gene %in% svgSigPos ~ "SigPos",
                                gene %in% svgSigNeg ~ "SigNeg",
                                .default = "NotSig")) %>%
  dplyr::filter(Sig %in% c("SigNeg")) %>%
  dplyr::mutate(diff = abs(percentDissimilarityX - percentDissimilarityY)) %>%
  dplyr::arrange(diff) %>%
  dplyr::pull(gene) %>%
  head(1)

# plot spatial patterns of example genes
sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(rast$AKI_ctrl)),
                          rownames(SpatialExperiment::spatialCoords(rast$AKI_aki)))

minList1 <- quantile(c(SummarizedExperiment::assay(rast$AKI_ctrl, "CPM")[gene1, sharedPixels], SummarizedExperiment::assay(rast$AKI_aki, "CPM")[gene1, sharedPixels]), 0.05)
maxList1 <- quantile(c(SummarizedExperiment::assay(rast$AKI_ctrl, "CPM")[gene1, sharedPixels], SummarizedExperiment::assay(rast$AKI_aki, "CPM")[gene1, sharedPixels]), 0.95)

sc1 <- ggplot2::scale_fill_gradientn(colors=viridis::viridis(20), 
                           limits = c(minList1, maxList1),
                           oob=scales::squish,
                           name = gene1)

minList2 <- quantile(c(SummarizedExperiment::assay(rast$AKI_ctrl, "CPM")[gene2, sharedPixels], SummarizedExperiment::assay(rast$AKI_aki, "CPM")[gene2, sharedPixels]), 0.05)
maxList2 <- quantile(c(SummarizedExperiment::assay(rast$AKI_ctrl, "CPM")[gene2, sharedPixels], SummarizedExperiment::assay(rast$AKI_aki, "CPM")[gene2, sharedPixels]), 0.95)

sc2 <- ggplot2::scale_fill_gradientn(colors=viridis::viridis(20), 
                           limits = c(minList2, maxList2),
                           oob=scales::squish,
                           name = gene2)

gene1_plt_ctrl  <- SEraster::plotRaster(rast$AKI_ctrl[gene1, sharedPixels], assay_name = "CPM", name = "CPM expression") + sc1
gene1_plt_aki <- SEraster::plotRaster(rast$AKI_aki[gene1, sharedPixels], assay_name = "CPM", name = "CPM expression") + sc1

gene1_plt_ctrl + gene1_plt_aki

gene2_plt_ctrl  <- SEraster::plotRaster(rast$AKI_ctrl[gene2, sharedPixels], assay_name = "CPM", name = "CPM expression") + sc2
gene2_plt_aki <- SEraster::plotRaster(rast$AKI_aki[gene2, sharedPixels], assay_name = "CPM", name = "CPM expression") + sc2

gene2_plt_ctrl + gene2_plt_aki

gene1_corr <- plotCorrelationGeneExp(
  speList = input, 
  spatialCorrelation = kidneyCorrelation, 
  geneName = gene1, 
  assayName = "CPM"
)
gene2_corr <- plotCorrelationGeneExp(
  speList = input, 
  spatialCorrelation = kidneyCorrelation, 
  geneName = gene2, 
  assayName = "CPM"
)


gene1_plt_ctrl | gene1_plt_aki | gene1_corr
gene2_plt_ctrl | gene2_plt_aki | gene2_corr

# Figure 2I - fold change linear regression and spot classification 
gene1_lr <- linearRegression(ss, gene=gene1, assayName = "CPM") + ggplot2::theme_classic() + ggplot2::coord_fixed()
gene2_lr <- linearRegression(ss, gene=gene2, assayName = "CPM") + ggplot2::theme_classic() + ggplot2::coord_fixed()

gene1_pc <- pixelClass(ss, gene=gene1, assayName = "CPM")
gene2_pc <- pixelClass(ss, gene=gene2, assayName = "CPM")

gene1_lr | gene1_pc
gene2_lr | gene2_pc

# figure 2f 
file_name = paste0("~/ST_compare/plots/figure2_kidney_spots_regions_20260320.pdf")
pdf(file = file_name, width = 4, height = 4, onefile=FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(p9, ncol=1, nrow=1),
                        gridExtra::arrangeGrob(p10, ncol=1, nrow=1))
dev.off()

# supplemental figure 
file_name = paste0("~/ST_compare/plots/figure2_kidney_tSNE_20260320.pdf")
pdf(file = file_name, width = 8, height = 8, onefile=FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(pctrl, paki, ncol=2, nrow=1),
                        gridExtra::arrangeGrob(ptSNEregion, ptSNEdataset, ncol=2, nrow=1))
dev.off()

#figure 2ghi


file_name = paste0("~/ST_compare/plots/kidney_plots_20251121_Vnn1_Nbl1.pdf")
pdf(file = file_name, width = 17, height = 4, onefile=FALSE)

#grab  legends from plots
# LegendHighS2R2 <- gtable::gtable_filter(ggplot2::ggplotGrob(gene1_plt_ctrl), "guide-box")
# LegendHighS2R3 <- gtable::gtable_filter(ggplot2::ggplotGrob(gene1_plt_aki), "guide-box")
# LegendLowS2R2 <- gtable::gtable_filter(ggplot2::ggplotGrob(gene2_plt_ctrl), "guide-box")
# LegendLowS2R3 <- gtable::gtable_filter(ggplot2::ggplotGrob(gene2_plt_aki), "guide-box")
# gridExtra::grid.arrange(gridExtra::arrangeGrob(gene1_plt_ctrl + ggplot2::theme(legend.position="none"),
#                                                gene1_plt_aki + ggplot2::theme(legend.position="none"),
#                                                gene1_corr,
#                                                gene1_lr, gene1_pc + ggplot2::theme(legend.position="none"),
#                                                LegendHighS2R2, LegendHighS2R3,
#                                                gene2_plt_ctrl + ggplot2::theme(legend.position="none"),
#                                                gene2_plt_aki + ggplot2::theme(legend.position="none"),
#                                                gene2_corr,
#                                                gene2_lr, gene2_pc + ggplot2::theme(legend.position="none"),
#                                                LegendLowS2R2, LegendLowS2R3,
#                                                ncol=7, nrow = 2, widths = c(3,3,3,3,3,1,1))
# )

#grab  legends from plots
Legend1 <- gtable::gtable_filter(ggplot2::ggplotGrob(gene1_plt_ctrl), "guide-box")
Legend2 <- gtable::gtable_filter(ggplot2::ggplotGrob(gene2_plt_ctrl), "guide-box")
gridExtra::grid.arrange(gridExtra::arrangeGrob(gene1_plt_ctrl + ggplot2::theme(legend.position="none"),
                                               gene1_plt_aki + ggplot2::theme(legend.position="none"),
                                               gene1_corr,
                                               gene1_lr, gene1_pc + ggplot2::theme(legend.position="none"),
                                               Legend1,
                                               gene2_plt_ctrl + ggplot2::theme(legend.position="none"),
                                               gene2_plt_aki + ggplot2::theme(legend.position="none"),
                                               gene2_corr,
                                               gene2_lr, gene2_pc + ggplot2::theme(legend.position="none"),
                                               Legend2,
                                               ncol=6, nrow = 2, widths = c(3,3,3,3,3,1))
)
dev.off()

