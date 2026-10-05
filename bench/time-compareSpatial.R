#!/usr/bin/env Rscript
# bench/time-compareSpatial.R
#
# Run time of compareSpatial() with its defaults (adaptive p-values: exceedances = 10, nPermutations = 10000;
# the extended delta grid; rank-remapped surrogates, surrogate = "remap", with no detection filter) on the
# realistic test genes and on the full inputs of the published AKI and brain analyses, and the number of
# permutations each gene used.
#
# Usage, from the repository root, with the package installed (R CMD INSTALL compiles the engine with R's
# optimising flags; devtools::load_all() compiles it with -O0 by default, about 10 times slower):
#
#   R CMD INSTALL .
#   Rscript bench/time-compareSpatial.R                  # 16 threads, every dataset available
#   Rscript bench/time-compareSpatial.R --threads=8 --datasets=aki_fixture,brain_fixture
#   Rscript bench/time-compareSpatial.R --datasets=aki,brain --nPermutations=1000
#
# Datasets: aki_fixture and brain_fixture (the 35 AKI genes on 311 pixels and the 30 brain genes on 2170
# pixels of tests/testthat/fixtures/realistic_fixture.rds), aki (the 1046 genes of the published AKI analysis,
# control vs AKI, assay CPM) and brain (the 325 genes of the published brain analysis, MERFISH vs Visium, assay
# lognorm). aki and brain read the inputs that data-raw/build_inputs_aki.R and data-raw/build_inputs_brain.R
# cache in <cache>/data-raw/inputs/ (<cache>: $STCOMPARE_DATA_CACHE or tools::R_user_dir("STcompare",
# "cache")); they are skipped when the cache is missing. --out=DIR saves the results there as RDS (never
# inside the repository).

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[length(hit)]) else default
}
threads <- as.integer(opt("threads", 16L))
n_perm <- as.numeric(opt("nPermutations", 10000))
datasets <- strsplit(opt("datasets", "aki_fixture,brain_fixture,aki,brain"), ",", fixed = TRUE)[[1]]
out_dir <- opt("out", NA_character_)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(script)) dirname(dirname(normalizePath(script))) else "."
cache <- Sys.getenv("STCOMPARE_DATA_CACHE", tools::R_user_dir("STcompare", which = "cache"))
inputs <- file.path(cache, "data-raw", "inputs")

suppressPackageStartupMessages({
  library(STcompare)
  library(SpatialExperiment)
})

# a pair of SpatialExperiment objects from a pair of the realistic fixture
fixture_pair <- function(P) {
  mk <- function(m) {
    dimnames(m) <- list(P$genes, P$pixel)
    SpatialExperiment(assays = list(counts = m),
                      spatialCoords = matrix(P$pos, ncol = 2, dimnames = list(P$pixel, c("x", "y"))))
  }
  list(x = mk(P$X), y = mk(P$Y), assay = "counts")
}
load_dataset <- function(id) {
  rf <- file.path(root, "tests", "testthat", "fixtures", "realistic_fixture.rds")
  switch(id,
    aki_fixture = fixture_pair(readRDS(rf)$pairs$aki),
    brain_fixture = fixture_pair(readRDS(rf)$pairs$brain),
    aki = {
      f <- file.path(inputs, "aki_rast.rds")
      if (!file.exists(f)) return(NULL)
      inp <- readRDS(f)
      list(x = inp$rast$AKI_ctrl[inp$genes, ], y = inp$rast$AKI_aki[inp$genes, ], assay = "CPM")
    },
    brain = {
      f <- file.path(inputs, "brain_merfish_visium_rast.rds")
      if (!file.exists(f)) return(NULL)
      inp <- readRDS(f)
      list(x = inp$rast$MERFISH, y = inp$rast$Visium, assay = "lognorm")
    },
    stop("unknown dataset ", id))
}

rows <- list()
for (id in datasets) {
  d <- load_dataset(id)
  if (is.null(d)) {
    message(id, ": cached input not found in ", inputs, "; skipped")
    next
  }
  t0 <- proc.time()[["elapsed"]]
  res <- compareSpatial(d$x, d$y, assay = d$assay, nPermutations = n_perm, nThreads = threads, progress = FALSE,
                        verbose = FALSE)
  wall <- proc.time()[["elapsed"]] - t0
  if (!is.na(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, file.path(out_dir, sprintf("compareSpatial-%s-B%d.rds", id, as.integer(n_perm))))
  }
  L <- res$nPermutations
  q <- stats::quantile(L, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  rows[[id]] <- data.frame(
    dataset = id, genes = nrow(res), pixels = res$nPixels[1], deltas = length(attr(res, "params")$deltaGrid),
    wall_s = round(wall, 1), s_per_gene = signif(wall / nrow(res), 3),
    perms_total = sum(L, na.rm = TRUE), perms_min = min(L, na.rm = TRUE), perms_q1 = q[1], perms_median = q[2],
    perms_q3 = q[3], perms_max = max(L, na.rm = TRUE), early = sum(res$stop == "exceedances"),
    limit = sum(res$stop == "limit"), skipped = sum(res$status == "skipped"), failed = sum(res$status == "failed"),
    padj_05 = sum(res$padj < 0.05, na.rm = TRUE),
    us_per_perm_thread = signif(1e6 * wall * threads / (2 * sum(L, na.rm = TRUE)), 3))
  message(sprintf("%s: %d genes in %.1f s", id, nrow(res), wall))
  mode <- attr(res, "params")$surrogate
}
tab <- do.call(rbind, rows)
cat(sprintf("\ncompareSpatial() with exceedances = 10, nPermutations = %d, %s surrogates, %d threads; %s, %s\n\n",
            as.integer(n_perm), mode, threads, R.version.string, Sys.info()[["machine"]]))
cat(paste0("| ", paste(names(tab), collapse = " | "), " |"), sep = "\n")
cat(paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"), sep = "\n")
for (i in seq_len(nrow(tab))) cat(paste0("| ", paste(format(unlist(tab[i, ]), scientific = FALSE), collapse = " | "), " |"), sep = "\n")
cat("\nperms_*: permutations per gene (the nPermutations column; both directions run that many). early: stopped by\n",
    "exceedances; limit: reached nPermutations; skipped: not tested (constant in a sample, or detected in fewer\n",
    "pixels than minDetected asks for).\n",
    "us_per_perm_thread: wall time x threads per permutation and direction (kept permutations only; the batches\n",
    "also compute a few that are discarded).\n", sep = "")
