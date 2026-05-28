library(STcompare)
library(rhdf5)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(BiocGenerics)
library(SEraster)
library(rjson)
library(scatterbar)
devtools::load_all()


# The MERFISH data can be found at Zenodo: https://zenodo.org/records/10724029
# S2R3 - Slice 2 Replicate 3

# MERFISH counts #################################################
zenodo_url <- "https://zenodo.org/records/10724029/files/STalign_S2R3_to_Visium.csv.gz?download=1"

temp <- tempfile(fileext = ".csv.gz")
download.file(zenodo_url, destfile = temp, mode = "wb")
df <- read.csv(gzfile(temp), row.names = 1)

head(df)

MERFISH.gexp <- df[, -(1:5)]
MERFISH.gexp <- as.matrix(t(MERFISH.gexp))
head(MERFISH.gexp)

# MERFISH positions #################################################

#subset area most informative (had high weight) in mapping with Visium
pos.aligned <- df[, c("STalign_y", "STalign_x")]
rownames(pos.aligned) <- rownames(df)
annot <- df[, c("Pmatch")]
names(annot) <- names(df)
vi <- annot > 0.95
MERFISH.pos.good <- pos.aligned[vi,]
colnames(MERFISH.pos.good) <- c("y", "x")

#rotate 90 clockwise to match the conventional orientation
MERFISH.pos.good <- cbind(MERFISH.pos.good[,"x"], -MERFISH.pos.good[,"y"])
rownames(MERFISH.pos.good) <- rownames(pos.aligned[vi,])

#switch names of x and y to match the convention of scatterbar
colnames(MERFISH.pos.good) <- c("x", "y")



# The Visium can be found on 10X Genomics:https://www.10xgenomics.com/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0 

# Visium counts  ###############################################################

url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"

temp_tar <- tempfile(fileext = ".tar.gz")
download.file(url, destfile = temp_tar, mode = "wb")

temp_dir <- tempdir()
untar(temp_tar, exdir = temp_dir)

matrix_path <- file.path(temp_dir, "filtered_feature_bc_matrix", "matrix.mtx.gz")
barcodes_path <- file.path(temp_dir, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
features_path <- file.path(temp_dir, "filtered_feature_bc_matrix", "features.tsv.gz")

Visium.gexp.raw <- Matrix::readMM(matrix_path)
Visium.barcodes <- read.csv(barcodes_path, header=FALSE)
Visium.genenames <- read.csv(features_path, sep='\t', header=FALSE)


dim(Visium.gexp.raw)
length(Visium.barcodes[,1])
length(Visium.genenames[,2])
rownames(Visium.gexp.raw) <- Visium.genenames[,2]
colnames(Visium.gexp.raw) <- Visium.barcodes[,1]
head(Visium.gexp.raw)



url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/Visium_FFPE_Mouse_Brain_spatial.tar.gz"

temp_tar <- tempfile(fileext = ".tar.gz")
download.file(url, destfile = temp_tar, mode = "wb")

temp_dir <- tempdir()
untar(temp_tar, exdir = temp_dir)

pos_path <- file.path(temp_dir, "spatial", "tissue_positions_list.csv")
scale_path <- file.path(temp_dir, "spatial", "scalefactors_json.json")

Visium.pos.info <- read.csv(pos_path, header = FALSE, row.names = 1)
head(Visium.pos.info)
Visium.pos <- Visium.pos.info[, 4:5]
Visium.pos.good <- Visium.pos[colnames(Visium.gexp.raw),]

convert <- fromJSON(file = scale_path)
range(Visium.pos.good*convert$tissue_hires_scalef)
Visium.pos.good <- Visium.pos.good*convert$tissue_hires_scalef
colnames(Visium.pos.good) <- c("y", "x")

#rotate 90 clockwise to match the conventional orientation
Visium.pos.good  <- cbind(Visium.pos.good[,"x"], -Visium.pos.good[,"y"])
rownames(Visium.pos.good) <- colnames(Visium.gexp.raw)

#switch names of x and y to match the convention of scatterbar
colnames(Visium.pos.good) <- c("x", "y")

genes.have <- intersect(rownames(Visium.gexp.raw), rownames(MERFISH.gexp))



MERFISH_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = MERFISH.gexp[genes.have, rownames(MERFISH.pos.good)],
                libnorm = MERINGUE::normalizeCounts(MERFISH.gexp[genes.have, rownames(MERFISH.pos.good)], log = FALSE)),
  spatialCoords = as.matrix(MERFISH.pos.good),
  
)

Visium_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = Visium.gexp.raw[genes.have,],
                libnorm = MERINGUE::normalizeCounts(Visium.gexp.raw[genes.have,], log=FALSE)),
  spatialCoords = as.matrix(Visium.pos.good)

)

# Create a list of SpatialExperiment objects to rasterize
spe_list <- list(MERFISH = MERFISH_SE, Visium = Visium_SE)
# Choose a rasterization resolution that is appropriate for the data. The resolution should be large enough to capture the spatial patterns of interest but not too large that it smooths out important details. 
resolution = 20
# Rasterize the gene expression data for each sample onto the same coordinate system. This will create a new list of SpatialExperiment objects that contains the rasterized gene expression values.
rast <- SEraster::rasterizeGeneExpression(spe_list, 
                                            resolution = resolution,
                                            assay_name = 'libnorm',
                                            fun = "mean",
                                            square = FALSE,
                                            BPPARAM = BiocParallel::MulticoreParam())

# add an assay that computes log transformation after mean-rasterization
SummarizedExperiment::assay(rast[["MERFISH"]], "lognorm") <- log10(SummarizedExperiment::assay(rast[["MERFISH"]]) + 1)
SummarizedExperiment::assay(rast[["Visium"]], "lognorm") <- log10(SummarizedExperiment::assay(rast[["Visium"]]) + 1)



# gene is expressed in greater than 1% of all spots 
# this has to be true in both datasets to be considered a "good gene" 
good.genes <- names(which(
  (Matrix::rowSums(SummarizedExperiment::assay(rast[["MERFISH"]], "lognorm") > 0) / dim(rast[["MERFISH"]])[2]*100 > 1) &
  (Matrix::rowSums(SummarizedExperiment::assay(rast[["Visium"]], "lognorm") > 0) / dim(rast[["Visium"]])[2]*100 > 1)))
length(good.genes)

# Reduced the number of genes from 466 --> 325 genes 
rast[["MERFISH"]] <- rast[["MERFISH"]][good.genes,]
rast[["Visium"]] <- rast[["Visium"]][good.genes,]



moransI <- function(SE_temp, filterDist, assayName=1){
  
  # use MERINGUE to get the neighbor-relationships 
  w <- MERINGUE::getSpatialNeighbors(SpatialExperiment::spatialCoords(SE_temp), filterDist = filterDist)
  par(mfrow=c(1,1))
  MERINGUE::plotNetwork(SpatialExperiment::spatialCoords(SE_temp), w)
  
  # Identify significantly spatially auto-correlated genes
  I <- MERINGUE::getSpatialPatterns(SummarizedExperiment::assays(SE_temp)[[assayName]], w)
  
  # filter for significant spatially auto-correlated genes with multiple testing correction and a minimum percentage of cells driving the pattern
  results.filter <- MERINGUE::filterSpatialPatterns(mat = SummarizedExperiment::assays(SE_temp)[[assayName]],
                                          I = I,
                                          w = w,
                                          adjustPv = TRUE,
                                          alpha = 0.05,
                                          minPercentCells = 0.01,
                                          verbose = TRUE, details = TRUE)

  return(results.filter)
}


moransI_MERFISH <- moransI(rast[["MERFISH"]], filterDist = resolution + 1, assayName = "lognorm")

moransI_Visium <- moransI(Visium_SE, filterDist = 25, assayName = "libnorm")


# Running the function `moransI` above already filtered the dataframe to only the genes that are statistically significantly autocorrelated so we can just take the rownames of the resulting dataframes to get the list of spatially variable genes in each dataset.
svg_MERFISH <- moransI_MERFISH %>% rownames()
svg_Visium <- moransI_Visium %>% rownames()

# genes that are SVGs in both conditions 
svg_int <- intersect(svg_MERFISH, svg_Visium) 

# genes that are not SVGs in both conditions
non_svg <- good.genes[!(good.genes %in% svg_int)] 

# Out of 325 genes: 
# there are 324 genes that are SVGs in the MERFISH dataset 
print("SVG MERFISH")
length(svg_MERFISH)

# there are 230 genes that are SVGs in the Visium dataset  
print("SVG Visium")
length(svg_Visium)

# there are 230 genes that are SVGs in both datasets   
print("SVG Intersect")
length(svg_int)

# there are 95 genes that are not SVGs in both datasets 
print("Not SVG")
length(non_svg)


# running spatialCorrelationGeneExpIterPermutations 

print("Brain Correlation")
start_time <- Sys.time()
brainCorrelation <- STcompare::spatialCorrelationGeneExpIterPermutations(
  rast,
  assayName = "lognorm", 
  returnPermutations = FALSE,
  nThreads = 22,
  seed = 0
)
end_time <- Sys.time()
print(end_time - start_time) 

save(brainCorrelation, file = file.path("inst", "extdata", "brain-MERFISH-10x-visium", "brainCorrelation.RData"))


# saving the significantly positively correlated svg genes 
svgSigPos <- brainCorrelation %>%
  dplyr::filter(rownames(brainCorrelation) %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%  
  dplyr::filter(correlationCoef > 0) %>%
  dplyr::arrange(desc(correlationCoef)) %>% 
  rownames()


# saving the significantly negatively correlated svg genes 
svgSigNeg <- brainCorrelation %>%
  dplyr::filter(rownames(brainCorrelation) %in% svg_int) %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%  
  dplyr::filter(correlationCoef < 0) %>% 
  rownames()


# saving the not significantly correlated non-svg genes
nonsvgNotSig <- brainCorrelation %>%
  dplyr::mutate(rownamesCol = rownames(brainCorrelation)) %>% 
  dplyr::filter(rownamesCol %in% non_svg) %>%
  dplyr::filter(pValuePermuteX >= 0.05 | pValuePermuteY >= 0.05) %>% 
  dplyr::arrange(correlationCoef) %>% 
  rownames()

geneSig <- svgSigPos[1]
geneNotSig <- nonsvgNotSig[1]

#get the shared pixels between the two samples to show only the pixels that are used for the one-to-one comparisons
sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(rast$MERFISH)),
                           rownames(SpatialExperiment::spatialCoords(rast$Visium)))



print("Spatial Similarity")
start_time <- Sys.time()
ss <- spatialSimilarity(rast, assayName = "lognorm", foldChange = 1, t1 = NULL, t2 = NULL)
end_time <- Sys.time()
print(end_time-start_time)




zenodo_url <- "https://zenodo.org/records/19582556/files/STalign_cell_type_transcriptional_correlations.csv.gz?download=1"
temp <- tempfile(fileext = ".csv.gz")
download.file(zenodo_url, destfile = temp, mode = "wb")

# cmat <- read.csv(gzfile(temp), row.names = 1)
cmat <- read.csv(gzfile(temp), row.names = 1, check.names = FALSE)

colnames(cmat) <- rownames(cmat)

zenodo_url <- "https://zenodo.org/records/19582556/files/STalign_S2R3_cell_type_annotations.csv.gz?download=1"
temp <- tempfile(fileext = ".csv.gz")
download.file(zenodo_url, destfile = temp, mode = "wb")

ctMERFISH <- read.csv(gzfile(temp), row.names = 1, stringsAsFactors = FALSE)
tempNames <- rownames(ctMERFISH)
levelNames <- rownames(cmat)
ctMERFISH <- factor(ctMERFISH$x, levels = levelNames)
names(ctMERFISH) <- tempNames

zenodo_url <- "https://zenodo.org/records/19582556/files/STalign_Visium_cell_type_annotations.csv.gz?download=1"
temp <- tempfile(fileext = ".csv.gz")
download.file(zenodo_url, destfile = temp, mode = "wb")

ctVisium <- read.csv(gzfile(temp), row.names = 1, stringsAsFactors = FALSE)
colnames(ctVisium) <- rownames(cmat)

cmat <- as.matrix(cmat)

# make list of the most transcriptionally similar cell types across the two technologies
best.cts <- unlist(lapply(1:nrow(cmat), function(i) {
  if (Matrix::diag(cmat)[i] == max(cmat[i,]) & Matrix::diag(cmat)[i] == max(cmat[,i])) {
    return(rownames(cmat)[i])
  }
}))


# assign colors to each cell type for visualization, prioritizing showing the best transcriptionally aligned cell types in color and the rest of the cell types in grey
set.seed(111111)
ccols <- rep('#eeeeee', nrow(cmat))
names(ccols) <- rownames(cmat)
ccols[best.cts] <- sample(rainbow(length(best.cts)), length(best.cts))
ccols
#names(ccols) <- make.names(names(ccols))

#plot the deconvolved cell type proportions for Visium
dfproportions <- data.frame(ctVisium)
colnames(dfproportions) <- names(ccols)
dfcoords <- data.frame(SpatialExperiment::spatialCoords(Visium_SE))


# make new spatial experiment for Visium with cluster labels as the assay
Visium_SE_ct <- SpatialExperiment::SpatialExperiment(
  assays = list(celltypes = t(ctVisium)),
  spatialCoords = as.matrix(SpatialExperiment::spatialCoords(Visium_SE))

)
# make new spatial experiment for MERFISH with one-hot encoded cluster labels as the assay
one_hot <- model.matrix(~ x - 1, data = data.frame( x = ctMERFISH[intersect(colnames(MERFISH_SE), names(ctMERFISH))]))
colnames(one_hot) <- gsub("x", "", colnames(one_hot))

MERFISH_SE_ct <- SpatialExperiment::SpatialExperiment(
  assays = list(celltypes = t(one_hot)),
  spatialCoords = as.matrix(SpatialExperiment::spatialCoords(MERFISH_SE)[intersect(colnames(MERFISH_SE), names(ctMERFISH)), ])
)


# rasterize the cell type labels for Visium using rasterizeGeneExpression since we have the deconvolved proportions for Visium rather than discrete cluster labels

input <- list(Visium = Visium_SE_ct, MERFISH = MERFISH_SE_ct)
rast_ct <- SEraster::rasterizeGeneExpression(input, 
  assay_name = "celltypes", 
  resolution = resolution, 
  fun = "mean", 
  square = FALSE,
  BPPARAM = BiocParallel::MulticoreParam())



sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(rast_ct$Visium)),
                           rownames(SpatialExperiment::spatialCoords(rast_ct$MERFISH)))

#For Visium, use scatterbar to visualize the proportion of each region in each sample
dfproportions_Visium <- data.frame(t(as.matrix(SummarizedExperiment::assay(rast_ct$Visium)[, sharedPixels])))
colnames(dfproportions_Visium) <- names(ccols)
dfcoords_Visium <- data.frame(SpatialExperiment::spatialCoords(rast_ct$Visium)[sharedPixels, ])


#For MERFISH, use scatterbar to visualize the proportion of each region in each sample
dfproportions_MERFISH <- data.frame(t(as.matrix(SummarizedExperiment::assay(rast_ct$MERFISH)[, sharedPixels])))
colnames(dfproportions_MERFISH) <- names(ccols)
dfcoords_MERFISH <- data.frame(SpatialExperiment::spatialCoords(rast_ct$MERFISH)[sharedPixels, ])


# running spatialCorrelationGeneExp 

print("ctCorrelation spatialCorrelationGeneExpIterPermutations")
start_time <- Sys.time()
ctCorrelation <- STcompare::spatialCorrelationGeneExpIterPermutations(
  rast_ct, 
  returnPermutations = FALSE,
  nThreads = 22
)
end_time <- Sys.time()
print(end_time - start_time) 

save(ctCorrelation, file = file.path("inst", "extdata", "brain-MERFISH-10x-visium", "ctCorrelation.RData"))



