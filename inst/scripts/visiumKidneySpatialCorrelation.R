
library(STcompare)
library(SpatialExperiment)
library(MERINGUE)
library(rhdf5)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)
library(BiocGenerics)


# The initlization is the same as the tutorial

# The Visium mouse kidney datasets can be found Zenodo: https://doi.org/10.5281/zenodo.17676991
# IL3 - ischemia acute kidney injury dataset
# NL3 - control dataset

# AKI counts ###################################################################
zenodo_url <- "https://zenodo.org/records/19074288/files/IL3_filtered_feature_bc_matrix.h5?download=1"
temp_h5 <- tempfile(fileext = ".h5")
download.file(zenodo_url, temp_h5, mode = "wb")


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

# positions ####################################################################

zenodo_url <- "https://zenodo.org/records/19074288/files/IL3_tissue_positions.csv?download=1"
aki_pos <- read.csv(zenodo_url, header = TRUE, stringsAsFactors = FALSE)

zenodo_url <- "https://zenodo.org/records/19074288/files/NL3_tissue_positions.csv?download=1"
ctrl_pos <- read.csv(zenodo_url, header = TRUE, stringsAsFactors = FALSE)

# limit the counts to the spots with positions and make sure they are in the same order
aki_counts <- aki_counts[, colnames(aki_counts) %in% aki_pos$barcode]
control_counts <- control_counts[, colnames(control_counts) %in% ctrl_pos$barcode]

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

zenodo_url <- "https://zenodo.org/records/19486091/files/aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz?download=1"
aki_STalign <- read.csv(gzcon(url(zenodo_url, open = "rb"), text = TRUE), header = TRUE, row.names =1)

# store the aligned positions
aki_pos_aligned <- aki_STalign[,c("aligned_x", "aligned_y")]

#rename columns and add group for consistency with ctrl_pos_rot
colnames(aki_pos_aligned) <- c("x", "y")
aki_pos_aligned$group  <- "AKI"
# rotate the tissue by 90-degrees for consistency with the paper
aki_pos_rot <- transform(aki_pos_aligned, x = y, y = -x + max(aki_pos_aligned$x))
aki_pos_aligned <- aki_pos_rot

ctrl_pos_rot <- transform(ctrl_pos, x = y, y = -x + max(ctrl_pos$x))


AKI_ctrl_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = control_counts),
  spatialCoords = as.matrix(ctrl_pos_rot[,1:2]),
)

AKI_aki_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = aki_counts),
  spatialCoords = as.matrix(aki_pos_aligned[,1:2])
)

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



w <- MERINGUE::getSpatialNeighbors(SpatialExperiment::spatialCoords(rast$AKI_ctrl), filterDist = 10)

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
moransI_aki <- moransI(rast$AKI_aki, filterDist = 10, assayName = "CPM")

svg_ctrl <- moransI_ctrl %>% filter(p.adj == 0) %>% rownames()
svg_aki <- moransI_aki %>% filter(p.adj == 0) %>% rownames()

svg_int <- intersect(svg_ctrl, svg_aki) # genes that are SVGs in both conditions

# choosing the genes that are SVGs in both conditions
genes_chosen <- svg_int

# alternative: choosing the genes that are expressed in greater than 5% of all spots
# genes_chosen <- good.genes

rownames(rast$AKI_ctrl) <- make.names(rownames(rast$AKI_ctrl), unique = TRUE)
rownames(rast$AKI_aki)  <- make.names(rownames(rast$AKI_aki),  unique = TRUE)

# The input is the spatial experiments subset to the chosen genes of interest
input <- list('AKI_ctrl'= rast$AKI_ctrl[genes_chosen,],
             'AKI_aki'= rast$AKI_aki[genes_chosen,])

# `deltaX` and `deltaY` are the parameters that controlling the degree of smoothing in permutations of datasets X and Y, respectively
deltaList <- replicate(length(genes_chosen), c(0.01, 0.05, seq(0.1, 0.9, .1)), simplify = FALSE)

# creates 100 permutations and then if a gene has an empirical p-value is less than 0.05, creates 1000 permutations for more accurate p-value estimation
nPermutations <- c(100, 1000)

# running spatialCorrelationGeneExp
start_time <- Sys.time()
deltaList <- replicate(length(genes_chosen), c(0.01, 0.05, seq(0.1, 0.9, .1)), simplify = FALSE)

kidneyCorrelation <- spatialCorrelationGeneExpIterPermutations(
  input,
  nPermutations = nPermutations,
  deltaX = deltaList,
  deltaY = deltaList,
  assayName = "CPM", # use the CPM normalized assay for this analysis
  nThreads = 22, # parallelize genes across threads
  BPPARAM = BiocParallel::MulticoreParam(), #setting up parallel-execution back-end
  verbose = TRUE, # prints name of gene as progress updates, which is helpful for long runs
  seed = 0) # set seed for reproducibility of the permutation results
end_time <- Sys.time()
print(end_time - start_time)

save(kidneyCorrelation, file = file.path("inst", "extdata", "kidneyCorrelation.RData"))
