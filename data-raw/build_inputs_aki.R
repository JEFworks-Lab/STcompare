#!/usr/bin/env Rscript
# data-raw/build_inputs_aki.R
#
# Rebuild the rasterized AKI kidney 10x Visium input (IL3 = ischemic AKI, NL3 = sham control) exactly as
# inst/scripts/visiumKidneySpatialCorrelation.R and vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd
# do, but from the md5-verified cached downloads (see data-raw/download_data.R).
#
# Run from the repository root:
#   Rscript data-raw/build_inputs_aki.R
# Output (cache, never the repository):
#   <cache>/data-raw/inputs/aki_rast.rds
#     list(rast = list(AKI_ctrl, AKI_aki) rasterized SpatialExperiments with assays counts and CPM,
#          shared = shared pixel IDs (311), genes = the 1046 genes of inst/extdata/kidneyCorrelation.RData,
#          meta = provenance)
# Check: Pearson r of every published gene must reproduce inst/extdata/kidneyCorrelation.RData to 1e-12.
# Needs rhdf5 (Bioconductor). Runtime about 1 min.
source("data-raw/download_data.R")
suppressPackageStartupMessages({ library(SpatialExperiment); library(SummarizedExperiment); library(Matrix) })
if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Package 'rhdf5' is required: BiocManager::install('rhdf5')")
t_start <- Sys.time()
f <- stc_download(group = "aki")

read_h5 <- function(h5) {
  bc <- as.character(rhdf5::h5read(h5, "matrix/barcodes"))
  bm <- rhdf5::h5read(h5, "matrix")
  m <- Matrix::sparseMatrix(dims = bm$shape, i = as.numeric(bm$indices), p = as.numeric(bm$indptr),
                            x = as.numeric(bm$data), index1 = FALSE)
  colnames(m) <- bc
  rownames(m) <- bm[["features"]]$name
  m
}
aki_counts <- read_h5(f[["IL3_filtered_feature_bc_matrix.h5"]])
control_counts <- read_h5(f[["NL3_filtered_feature_bc_matrix.h5"]])
cat("h5 dims: IL3 (AKI)", dim(aki_counts), "; NL3 (control)", dim(control_counts), "\n")

aki_pos <- read.csv(f[["IL3_tissue_positions.csv"]], header = TRUE, stringsAsFactors = FALSE)
ctrl_pos <- read.csv(f[["NL3_tissue_positions.csv"]], header = TRUE, stringsAsFactors = FALSE)
aki_counts <- aki_counts[, colnames(aki_counts) %in% aki_pos$barcode]
control_counts <- control_counts[, colnames(control_counts) %in% ctrl_pos$barcode]
cat("in-tissue dims: IL3", dim(aki_counts), "; NL3", dim(control_counts), "\n")

rownames(ctrl_pos) <- ctrl_pos$barcode; ctrl_pos <- ctrl_pos[, -1]; colnames(ctrl_pos) <- c("x", "y")
rownames(aki_pos) <- aki_pos$barcode; aki_pos <- aki_pos[, -1]; colnames(aki_pos) <- c("x", "y")
ctrl_pos$group <- "Control"; aki_pos$group <- "AKI"

aki_STalign <- read.csv(gzfile(f[["aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz"]]),
                        header = TRUE, row.names = 1)
aki_pos_aligned <- aki_STalign[, c("aligned_x", "aligned_y")]
colnames(aki_pos_aligned) <- c("x", "y")
aki_pos_aligned$group <- "AKI"
# rotate both tissues by 90 degrees, as in the published script
aki_pos_aligned <- transform(aki_pos_aligned, x = y, y = -x + max(aki_pos_aligned$x))
ctrl_pos_rot <- transform(ctrl_pos, x = y, y = -x + max(ctrl_pos$x))
# the published script passes the coordinates without reordering them to the count columns
stopifnot(identical(colnames(control_counts), rownames(ctrl_pos_rot)),
          identical(colnames(aki_counts), rownames(aki_pos_aligned)))

AKI_ctrl_SE <- SpatialExperiment(assays = list(counts = control_counts), spatialCoords = as.matrix(ctrl_pos_rot[, 1:2]))
AKI_aki_SE <- SpatialExperiment(assays = list(counts = aki_counts), spatialCoords = as.matrix(aki_pos_aligned[, 1:2]))
t0 <- Sys.time()
rast <- SEraster::rasterizeGeneExpression(list(AKI_ctrl = AKI_ctrl_SE, AKI_aki = AKI_aki_SE),
                                          resolution = 5, fun = "sum", square = FALSE, assay_name = "counts",
                                          BPPARAM = BiocParallel::SerialParam())
cat(sprintf("rasterized in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
assay(rast$AKI_ctrl, "CPM") <- Matrix::t(Matrix::t(assay(rast$AKI_ctrl)) / Matrix::colSums(assay(rast$AKI_ctrl))) * 1e6
assay(rast$AKI_aki, "CPM") <- Matrix::t(Matrix::t(assay(rast$AKI_aki)) / Matrix::colSums(assay(rast$AKI_aki))) * 1e6
rownames(rast$AKI_ctrl) <- make.names(rownames(rast$AKI_ctrl), unique = TRUE)
rownames(rast$AKI_aki) <- make.names(rownames(rast$AKI_aki), unique = TRUE)
shared <- intersect(rownames(spatialCoords(rast$AKI_ctrl)), rownames(spatialCoords(rast$AKI_aki)))
for (nm in names(rast)) cat(sprintf("raster %s: %d genes x %d pixels\n", nm, nrow(rast[[nm]]), ncol(rast[[nm]])))
cat("shared pixels:", length(shared), "\n")

# reproduction check against the published results (computed by the package authors on the original inputs)
e <- new.env(); load(file.path("inst", "extdata", "kidneyCorrelation.RData"), envir = e)
kc <- e$kidneyCorrelation
genes <- rownames(kc)
stopifnot(all(genes %in% rownames(rast$AKI_ctrl)))
X <- as.matrix(assay(rast$AKI_ctrl, "CPM")[genes, shared])
Y <- as.matrix(assay(rast$AKI_aki, "CPM")[genes, shared])
r <- vapply(seq_along(genes), function(i) cor(X[i, ], Y[i, ]), 0)
dr <- max(abs(r - kc$correlationCoef))
cat(sprintf("reproduction of kidneyCorrelation$correlationCoef: max|dr| = %.3g (%d/%d bit-identical)\n",
            dr, sum(r == kc$correlationCoef), length(r)))
if (!(dr <= 1e-12)) stop("Rebuilt AKI input does not reproduce the published correlations (max |dr| = ", dr, ")")

meta <- stc_meta("data-raw/build_inputs_aki.R",
                 inputs = stc_manifest[stc_manifest$group == "aki", c("file", "md5", "url")])
stc_save_rds(list(rast = rast, shared = shared, genes = genes, meta = meta), stc_input_file("aki_rast.rds"))
cat(sprintf("done in %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
