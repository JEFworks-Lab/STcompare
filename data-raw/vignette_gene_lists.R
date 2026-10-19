#!/usr/bin/env Rscript
# data-raw/vignette_gene_lists.R
#
# Regenerate the gene lists that the case-study articles read from inst/extdata/:
#   * vignette-aki-svg-genes.txt: the 1046 genes of the published AKI analysis (the row names of
#     bench/published/kidneyCorrelation.RData), in the gene symbols of the 10x feature files. The published
#     analysis named the genes with make.names(unique = TRUE), which changed 19 symbols (mt.Nd1 for mt-Nd1);
#     all 1046 genes are first occurrences of their symbol, so the articles' make.unique() keeps their names.
#     They are the genes spatially variable in both sections by MERINGUE's Moran's I (adjusted p-value of 0
#     in both); MERINGUE 1.0 on the rebuilt input gives the same list up to 4 genes (the p.adj == 0 criterion
#     depends on floating-point underflow), which this script reports;
#   * vignette-brain-svg-genes.txt: the 230 of the 325 analysed brain genes that are spatially variable in
#     both technologies, recomputed with MERINGUE as in the published analysis.
#
# Needs the GitHub package MERINGUE (remotes::install_github("JEFworks-Lab/MERINGUE")) and rhdf5, and the
# rasterized inputs: run data-raw/build_inputs_aki.R and data-raw/build_inputs_brain.R first. Run from the
# repository root:
#   Rscript data-raw/vignette_gene_lists.R              # writes inst/extdata/vignette-*-svg-genes.txt
#   Rscript data-raw/vignette_gene_lists.R --out=DIR    # writes them to DIR instead (to compare)
# Runtime about 4 minutes. The lists it writes are identical to the shipped ones (checked 2026-10-05).
source("data-raw/download_data.R")
suppressPackageStartupMessages({ library(SpatialExperiment); library(SummarizedExperiment); library(Matrix) })
if (!requireNamespace("MERINGUE", quietly = TRUE)) {
  stop("MERINGUE is needed: remotes::install_github(\"JEFworks-Lab/MERINGUE\")", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(hit <- grep("^--out=", args, value = TRUE))) sub("^--out=", "", hit[1]) else file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

moran <- function(mat, coords, filterDist, filter = FALSE) {
  w <- MERINGUE::getSpatialNeighbors(coords, filterDist = filterDist)
  I <- MERINGUE::getSpatialPatterns(mat, w)
  if (!filter) return(I)
  MERINGUE::filterSpatialPatterns(mat = mat, I = I, w = w, adjustPv = TRUE, alpha = 0.05, minPercentCells = 0.01,
                                  details = TRUE)
}
write_list <- function(genes, file, header) {
  path <- file.path(out_dir, file)
  writeLines(c(paste("#", header), genes), path)
  cat(sprintf("wrote %s (%d genes)\n", path, length(genes)))
}

# --- AKI: the published genes, in the 10x symbols ------------------------------------------------------
aki <- readRDS(stc_input_file("aki_rast.rds"))
f <- stc_download(group = "aki", quiet = TRUE)
symbols <- as.character(rhdf5::h5read(f[["NL3_filtered_feature_bc_matrix.h5"]], "matrix/features")$name)
idx <- match(aki$genes, make.names(symbols, unique = TRUE))
stopifnot(!anyNA(idx), identical(make.unique(symbols)[idx], symbols[idx]))
akiGenes <- symbols[idx]
iC <- moran(assay(aki$rast$AKI_ctrl, "CPM"), spatialCoords(aki$rast$AKI_ctrl), 10)
iA <- moran(assay(aki$rast$AKI_aki, "CPM"), spatialCoords(aki$rast$AKI_aki), 10)
svgMeringue <- intersect(rownames(iC)[which(iC$p.adj == 0)], rownames(iA)[which(iA$p.adj == 0)])
cat("AKI: MERINGUE on the rebuilt input: published only:", paste(setdiff(aki$genes, svgMeringue), collapse = ", "),
    "; MERINGUE only:", paste(setdiff(svgMeringue, aki$genes), collapse = ", "), "\n")
write_list(akiGenes, "vignette-aki-svg-genes.txt", c(
  "Spatially variable genes used in the article \"Acute kidney injury (10x Visium)\".",
  "These are the 1046 genes of the published analysis (Clifton et al. 2026, Figure 2f-j): genes spatially",
  "variable in both the sham control (NL3) and the AKI (IL3) kidney section by Moran's I",
  "(MERINGUE::getSpatialPatterns() on the CPM of the hexagonal pixels of resolution 5, neighbours within",
  "a distance of 10; adjusted p-value of 0 in both sections). Gene symbols as in the 10x feature files;",
  "one gene per line, in the order of the feature files."))

# --- brain: Moran's I of the MERFISH pixels and of the Visium spots (before rasterization) --------------
brain <- readRDS(stc_input_file("brain_merfish_visium_rast.rds"))
f <- stc_download(group = "brain", quiet = TRUE)
td <- stc_cache_dir("extracted", "visium_brain")
untar(f[["Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"]], exdir = td)
untar(f[["Visium_FFPE_Mouse_Brain_spatial.tar.gz"]], exdir = td)
counts <- Matrix::readMM(file.path(td, "filtered_feature_bc_matrix", "matrix.mtx.gz"))
barcodes <- read.csv(file.path(td, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), header = FALSE)[, 1]
features <- read.csv(file.path(td, "filtered_feature_bc_matrix", "features.tsv.gz"), sep = "\t", header = FALSE)[, 2]
dimnames(counts) <- list(features, barcodes)
positions <- read.csv(file.path(td, "spatial", "tissue_positions_list.csv"), header = FALSE, row.names = 1)
scale <- jsonlite::fromJSON(file.path(td, "spatial", "scalefactors_json.json"))$tissue_hires_scalef
visiumXY <- cbind(x = positions[barcodes, 5], y = -positions[barcodes, 4]) * scale
rownames(visiumXY) <- barcodes
# the gene names of the MERFISH table (after its row names and 5 columns of cell data), as read.csv() names them
merfishGenes <- names(utils::read.csv(gzfile(stc_download(files = "STalign_S2R3_to_Visium.csv.gz", quiet = TRUE)),
                                      nrows = 1))[-(1:6)]
genesHave <- intersect(features, merfishGenes)
libnorm <- Matrix::t(Matrix::t(counts[genesHave, ]) / Matrix::colSums(counts[genesHave, ])) * 1e6
iM <- moran(assay(brain$rast$MERFISH, "lognorm"), spatialCoords(brain$rast$MERFISH), 21, filter = TRUE)
iV <- moran(libnorm, visiumXY, 25, filter = TRUE)
brainGenes <- intersect(rownames(iM), rownames(iV))
cat("brain:", length(brainGenes), "of", nrow(brain$rast$MERFISH), "analysed genes are spatially variable in both\n")
write_list(brainGenes, "vignette-brain-svg-genes.txt", c(
  "Spatially variable genes used in the article \"Comparison of MERFISH and Visium for mouse brain\".",
  "These are the 230 of the 325 analysed genes that are spatially variable in both technologies by Moran's I",
  "(MERINGUE::getSpatialPatterns() and filterSpatialPatterns(adjustPv = TRUE, alpha = 0.05,",
  "minPercentCells = 0.01)), computed as in the published analysis (Clifton et al. 2026): on log10(x + 1) of the",
  "MERFISH pixels of resolution 20 (neighbours within 21) and on the library-normalized Visium spots (neighbours",
  "within 25). One gene per line."))
