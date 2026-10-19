#!/usr/bin/env Rscript
# data-raw/build_inputs_merfish_replicates.R
#
# Rebuild the rasterized MERFISH biological-replicate inputs exactly as
# bench/published/scripts/biological-replicates-example.R does, from the md5-verified cached downloads: target = slice 2 replicate 2 (S2R2), source = S2R3 aligned to S2R2
# with STalign (coordinates STalign_x/y) and, separately, with the affine-only alignment (affine_x/y);
# SEraster resolution 200 with its defaults (first assay = counts, fun = "mean", square pixels); "Blank-*"
# control probes dropped (483 genes; 1371 shared pixels for STalign, 1299 for affine).
# These inputs are not used by the committed fixtures; they are the tier-3 (full-scale benchmark) data.
#
# Run from the repository root:
#   Rscript data-raw/build_inputs_merfish_replicates.R
# Output (cache, never the repository):
#   <cache>/data-raw/inputs/merfish_replicates_rast.rds
#     list(rast = list(target, source) [STalign], rast_affine = list(target, source) [affine], meta)
# Checks: Pearson r reproduces merfishCorrelation.RData and merfishCorrelation_affine.RData to 1e-12.
# Runtime about 2 min.
source("data-raw/download_data.R")
suppressPackageStartupMessages({ library(SpatialExperiment); library(SummarizedExperiment); library(Matrix) })
t_start <- Sys.time()
f <- stc_download(group = "merfish")

target_tab <- read.csv(gzfile(f[["STalign_S2R2.csv.gz"]]))
source_tab <- read.csv(gzfile(f[["STalign_S2R3_to_S2R2.csv.gz"]]))
cat("S2R2 (target):", dim(target_tab), "; S2R3 -> S2R2 (source):", dim(source_tab), "\n")
# cell IDs are parsed as doubles (e.g. 1.00442548580637e+38), exactly as in the published script
pos_target <- target_tab[, c("x", "y")]; rownames(pos_target) <- target_tab$X
gene_target <- target_tab[, 4:ncol(target_tab)]; rownames(gene_target) <- target_tab$X
spe_target <- SpatialExperiment(assays = list(counts = as(t(gene_target), "dgCMatrix")),
                                spatialCoords = as.matrix(pos_target))
pos_source <- source_tab[, c("STalign_x", "STalign_y")]; rownames(pos_source) <- source_tab$X; colnames(pos_source) <- c("x", "y")
gene_source <- source_tab[, 8:ncol(source_tab)]; rownames(gene_source) <- source_tab$X
spe_source <- SpatialExperiment(assays = list(counts = as(t(gene_source), "dgCMatrix")),
                                spatialCoords = as.matrix(pos_source))
pos_affine <- source_tab[, c("affine_x", "affine_y")]; rownames(pos_affine) <- source_tab$X; colnames(pos_affine) <- c("x", "y")
spe_affine <- SpatialExperiment(assays = list(counts = as(t(gene_source), "dgCMatrix")),
                                spatialCoords = as.matrix(pos_affine))
stopifnot(identical(rownames(spe_target), rownames(spe_source)))

rasterize200 <- function(lst) SEraster::rasterizeGeneExpression(lst, resolution = 200, BPPARAM = BiocParallel::SerialParam())
t0 <- Sys.time(); out <- rasterize200(list(target = spe_target, source = spe_source))
cat(sprintf("rasterized (STalign) in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
t0 <- Sys.time(); out_aff <- rasterize200(list(target = spe_target, source = spe_affine))
cat(sprintf("rasterized (affine) in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
genes_only <- rownames(out$source)[!grepl("Blank", rownames(out$source))]
out <- list(target = out$target[genes_only, ], source = out$source[genes_only, ])
out_aff <- list(target = out_aff$target[genes_only, ], source = out_aff$source[genes_only, ])

check <- function(o, f, label) {
  sh <- intersect(rownames(spatialCoords(o[[1]])), rownames(spatialCoords(o[[2]])))
  e <- new.env(); n <- load(stc_published_file(f), envir = e); ref <- get(n, envir = e)
  X <- as.matrix(assay(o[[1]])[rownames(ref), sh]); Y <- as.matrix(assay(o[[2]])[rownames(ref), sh])
  r <- vapply(seq_len(nrow(ref)), function(i) cor(X[i, ], Y[i, ]), 0)
  dr <- max(abs(r - ref$correlationCoef))
  cat(sprintf("%s: %d genes, target %d px, source %d px, shared %d; reproduction of %s: max|dr| = %.3g (%d/%d bit-identical)\n",
              label, length(genes_only), ncol(o[[1]]), ncol(o[[2]]), length(sh), f, dr, sum(r == ref$correlationCoef), nrow(ref)))
  if (!(dr <= 1e-12)) stop("Rebuilt MERFISH input (", label, ") does not reproduce ", f, " (max |dr| = ", dr, ")")
  invisible(sh)
}
check(out, "merfishCorrelation.RData", "STalign")
check(out_aff, "merfishCorrelation_affine.RData", "affine")
meta <- stc_meta("data-raw/build_inputs_merfish_replicates.R",
                 inputs = stc_manifest[stc_manifest$group == "merfish", c("file", "md5", "url")])
stc_save_rds(list(rast = out, rast_affine = out_aff, meta = meta), stc_input_file("merfish_replicates_rast.rds"))
cat(sprintf("done in %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
