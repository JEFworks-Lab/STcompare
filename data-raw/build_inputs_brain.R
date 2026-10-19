#!/usr/bin/env Rscript
# data-raw/build_inputs_brain.R
#
# Rebuild the rasterized MERFISH-vs-Visium mouse brain inputs exactly as
# bench/published/scripts/brain-MERFISH-10x-visium.R and vignettes/articles/brain-MERFISH-10x-visium.Rmd do,
# from the md5-verified cached downloads:
#   * gene expression: MERFISH slice 2 replicate 3 aligned to Visium with STalign (cells with Pmatch > 0.95)
#     vs 10x Visium FFPE adult mouse brain; libnorm = counts / library size * 1e6 on the shared genes;
#     SEraster resolution 20 (hires pixels), fun = "mean", hexagons; lognorm = log10(x + 1); genes detected
#     in > 1% of pixels in both datasets (325 genes, 2170 shared pixels);
#   * cell types (input of bench/published/brain-MERFISH-10x-visium/ctCorrelation.RData): Visium deconvolved
#     proportions vs MERFISH one-hot cell-type labels, same rasterization (16 cell types, 2174 shared pixels).
# MERINGUE::normalizeCounts(log = FALSE) is re-implemented below (verbatim from MERINGUE/R/process.R, GPL-3,
# JEFworks-Lab) so that MERINGUE is not needed.
#
# Run from the repository root:
#   Rscript data-raw/build_inputs_brain.R
# Outputs (cache, never the repository):
#   <cache>/data-raw/inputs/brain_merfish_visium_rast.rds  list(rast = list(MERFISH, Visium), shared, meta)
#   <cache>/data-raw/inputs/brain_celltype_rast.rds        list(rast = list(Visium, MERFISH), shared, meta)
# Checks: Pearson r reproduces brainCorrelation.RData and ctCorrelation.RData to 1e-12.
# Runtime about 1.5 min.
source("data-raw/download_data.R")
suppressPackageStartupMessages({ library(SpatialExperiment); library(SummarizedExperiment); library(Matrix) })
t_start <- Sys.time()
f <- stc_download(group = c("brain", "celltype"))

normalizeCounts <- function(counts, normFactor = NULL, depthScale = 1e+06, pseudo = 1, log = TRUE) {
  if (!any(class(counts) %in% c("dgCMatrix", "dgTMatrix"))) counts <- Matrix::Matrix(counts, sparse = TRUE)
  if (is.null(normFactor)) normFactor <- Matrix::colSums(counts)
  counts <- Matrix::t(Matrix::t(counts) / normFactor)
  counts <- counts * depthScale
  if (log) counts <- log10(counts + pseudo)
  counts
}
read_scalefactors <- function(path) {
  if (requireNamespace("rjson", quietly = TRUE)) rjson::fromJSON(file = path) else jsonlite::fromJSON(path)
}

# --- MERFISH S2R3 aligned to Visium -------------------------------------------------------------------
df <- read.csv(gzfile(f[["STalign_S2R3_to_Visium.csv.gz"]]), row.names = 1)
MERFISH.gexp <- as.matrix(t(df[, -(1:5)]))
pos.aligned <- df[, c("STalign_y", "STalign_x")]; rownames(pos.aligned) <- rownames(df)
annot <- df[, c("Pmatch")]; names(annot) <- names(df)
vi <- annot > 0.95
cat("MERFISH cells:", nrow(df), "; Pmatch > 0.95:", sum(vi), "; genes:", nrow(MERFISH.gexp), "\n")
MERFISH.pos.good <- pos.aligned[vi, ]; colnames(MERFISH.pos.good) <- c("y", "x")
MERFISH.pos.good <- cbind(MERFISH.pos.good[, "x"], -MERFISH.pos.good[, "y"])   # rotate 90 degrees clockwise
rownames(MERFISH.pos.good) <- rownames(pos.aligned[vi, ]); colnames(MERFISH.pos.good) <- c("x", "y")

# --- Visium FFPE adult mouse brain ---------------------------------------------------------------------
td <- stc_cache_dir("extracted", "visium_brain")
untar(f[["Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"]], exdir = td)
untar(f[["Visium_FFPE_Mouse_Brain_spatial.tar.gz"]], exdir = td)
Visium.gexp.raw <- Matrix::readMM(file.path(td, "filtered_feature_bc_matrix", "matrix.mtx.gz"))
Visium.barcodes <- read.csv(file.path(td, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), header = FALSE)
Visium.genenames <- read.csv(file.path(td, "filtered_feature_bc_matrix", "features.tsv.gz"), sep = "\t", header = FALSE)
rownames(Visium.gexp.raw) <- Visium.genenames[, 2]; colnames(Visium.gexp.raw) <- Visium.barcodes[, 1]
cat("Visium filtered matrix:", dim(Visium.gexp.raw), "\n")
Visium.pos.info <- read.csv(file.path(td, "spatial", "tissue_positions_list.csv"), header = FALSE, row.names = 1)
Visium.pos.good <- Visium.pos.info[, 4:5][colnames(Visium.gexp.raw), ]
convert <- read_scalefactors(file.path(td, "spatial", "scalefactors_json.json"))
Visium.pos.good <- Visium.pos.good * convert$tissue_hires_scalef
colnames(Visium.pos.good) <- c("y", "x")
Visium.pos.good <- cbind(Visium.pos.good[, "x"], -Visium.pos.good[, "y"])     # rotate 90 degrees clockwise
rownames(Visium.pos.good) <- colnames(Visium.gexp.raw); colnames(Visium.pos.good) <- c("x", "y")
genes.have <- intersect(rownames(Visium.gexp.raw), rownames(MERFISH.gexp))
cat("genes shared by MERFISH and Visium:", length(genes.have), "\n")

MERFISH_SE <- SpatialExperiment(
  assays = list(counts = MERFISH.gexp[genes.have, rownames(MERFISH.pos.good)],
                libnorm = normalizeCounts(MERFISH.gexp[genes.have, rownames(MERFISH.pos.good)], log = FALSE)),
  spatialCoords = as.matrix(MERFISH.pos.good))
Visium_SE <- SpatialExperiment(
  assays = list(counts = Visium.gexp.raw[genes.have, ],
                libnorm = normalizeCounts(Visium.gexp.raw[genes.have, ], log = FALSE)),
  spatialCoords = as.matrix(Visium.pos.good))
t0 <- Sys.time()
rast <- SEraster::rasterizeGeneExpression(list(MERFISH = MERFISH_SE, Visium = Visium_SE), resolution = 20,
                                          assay_name = "libnorm", fun = "mean", square = FALSE,
                                          BPPARAM = BiocParallel::SerialParam())
cat(sprintf("rasterized in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
assay(rast[["MERFISH"]], "lognorm") <- log10(assay(rast[["MERFISH"]]) + 1)
assay(rast[["Visium"]], "lognorm") <- log10(assay(rast[["Visium"]]) + 1)
good.genes <- names(which(
  (Matrix::rowSums(assay(rast[["MERFISH"]], "lognorm") > 0) / dim(rast[["MERFISH"]])[2] * 100 > 1) &
  (Matrix::rowSums(assay(rast[["Visium"]], "lognorm") > 0) / dim(rast[["Visium"]])[2] * 100 > 1)))
rast[["MERFISH"]] <- rast[["MERFISH"]][good.genes, ]
rast[["Visium"]] <- rast[["Visium"]][good.genes, ]
shared <- intersect(rownames(spatialCoords(rast$MERFISH)), rownames(spatialCoords(rast$Visium)))
cat(sprintf("raster: %d genes; MERFISH %d px, Visium %d px, shared %d\n",
            length(good.genes), ncol(rast$MERFISH), ncol(rast$Visium), length(shared)))

e <- new.env(); load(stc_published_file("brain-MERFISH-10x-visium", "brainCorrelation.RData"), envir = e)
ref <- e$brainCorrelation
stopifnot(identical(rownames(ref), good.genes))
X <- as.matrix(assay(rast$MERFISH, "lognorm")[rownames(ref), shared])
Y <- as.matrix(assay(rast$Visium, "lognorm")[rownames(ref), shared])
r <- vapply(seq_len(nrow(ref)), function(i) cor(X[i, ], Y[i, ]), 0)
dr <- max(abs(r - ref$correlationCoef))
cat(sprintf("reproduction of brainCorrelation$correlationCoef: max|dr| = %.3g (%d/%d bit-identical)\n",
            dr, sum(r == ref$correlationCoef), nrow(ref)))
if (!(dr <= 1e-12)) stop("Rebuilt brain input does not reproduce the published correlations (max |dr| = ", dr, ")")
meta <- stc_meta("data-raw/build_inputs_brain.R",
                 inputs = stc_manifest[stc_manifest$group == "brain", c("file", "md5", "url")])
stc_save_rds(list(rast = rast, shared = shared, meta = meta), stc_input_file("brain_merfish_visium_rast.rds"))

# --- cell types (ctCorrelation input; X = Visium, Y = MERFISH) ------------------------------------------
cmat <- read.csv(gzfile(f[["STalign_cell_type_transcriptional_correlations.csv.gz"]]), row.names = 1, check.names = FALSE)
colnames(cmat) <- rownames(cmat)
ctM <- read.csv(gzfile(f[["STalign_S2R3_cell_type_annotations.csv.gz"]]), row.names = 1, stringsAsFactors = FALSE)
tn <- rownames(ctM); ctM <- factor(ctM$x, levels = rownames(cmat)); names(ctM) <- tn
ctV <- read.csv(gzfile(f[["STalign_Visium_cell_type_annotations.csv.gz"]]), row.names = 1, stringsAsFactors = FALSE)
colnames(ctV) <- rownames(cmat)
stopifnot(identical(rownames(ctV), colnames(Visium.gexp.raw)))
Vse <- SpatialExperiment(assays = list(celltypes = t(ctV)), spatialCoords = Visium.pos.good)
keep <- intersect(rownames(MERFISH.pos.good), names(ctM))
oh <- model.matrix(~ x - 1, data = data.frame(x = ctM[keep])); colnames(oh) <- gsub("x", "", colnames(oh))
Mse <- SpatialExperiment(assays = list(celltypes = t(oh)), spatialCoords = MERFISH.pos.good[keep, ])
rc <- SEraster::rasterizeGeneExpression(list(Visium = Vse, MERFISH = Mse), assay_name = "celltypes", resolution = 20,
                                        fun = "mean", square = FALSE, BPPARAM = BiocParallel::SerialParam())
sh_ct <- intersect(colnames(rc$Visium), colnames(rc$MERFISH))
e <- new.env(); load(stc_published_file("brain-MERFISH-10x-visium", "ctCorrelation.RData"), envir = e)
ref_ct <- e$ctCorrelation
r_ct <- vapply(rownames(ref_ct), function(g) cor(as.numeric(assay(rc$Visium)[g, sh_ct]), as.numeric(assay(rc$MERFISH)[g, sh_ct])), 0)
dr_ct <- max(abs(r_ct - ref_ct$correlationCoef))
cat(sprintf("cell types: %d annotated cells, %d types, shared %d px; reproduction of ctCorrelation: max|dr| = %.3g\n",
            length(keep), nrow(rc$Visium), length(sh_ct), dr_ct))
if (!(dr_ct <= 1e-12)) stop("Rebuilt cell-type input does not reproduce ctCorrelation (max |dr| = ", dr_ct, ")")
meta_ct <- stc_meta("data-raw/build_inputs_brain.R",
                    inputs = stc_manifest[stc_manifest$group %in% c("brain", "celltype"), c("file", "md5", "url")])
stc_save_rds(list(rast = rc, shared = sh_ct, meta = meta_ct), stc_input_file("brain_celltype_rast.rds"))
cat(sprintf("done in %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
