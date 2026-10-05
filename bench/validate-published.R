#!/usr/bin/env Rscript
# bench/validate-published.R
#
# Main acceptance test of the compiled engine (dev/engine-spec.md, section 5): re-run every published
# STcompare analysis whose results are kept in bench/published (they shipped in inst/extdata until version
# 0.1.0), from inputs rebuilt from the public downloads (data-raw/build_inputs_*.R), with the authors'
# parameters (bench/published/scripts/), and compare gene by gene.
#
# Usage, from the repository root. By default the script installs the package it sits in (a copy of R/,
# src/, inst/, ...) into a temporary library with R CMD INSTALL, so that the engine is compiled with R's
# optimising flags, as users get it. (devtools::load_all() compiles with -O0 by default: about 10 times slower.)
#
#   Rscript bench/validate-published.R                       # mode "exported", 16 threads, all analyses
#   Rscript bench/validate-published.R --mode=internal       # call .stc_engine_correlate() directly
#   Rscript bench/validate-published.R --threads=8 --analyses=celltypes,brain
#
# Options (each also settable as an environment variable, shown in brackets, or, when the script is
# source()d, as a variable of the same name in the global environment, e.g. validate_mode <- "internal"):
#   --mode=exported|internal  [STCOMPARE_VALIDATE_MODE, validate_mode]
#       exported  call spatialCorrelationGeneExpIterPermutations() / spatialCorrelationGeneExp() exactly
#                 as the authors did (assayName, delta lists, nPermutations, seed, adjustMethod = "BH").
#       internal  call the engine (.stc_engine_correlate()) directly: round 1 for every gene, then the
#                 genes whose unadjusted pValuePermuteX and pValuePermuteY are both below
#                 (alpha / nPermutations[k]) * 100 continue the same session to nPermutations[k + 1]
#                 (prefix property); BH across genes at the end.
#   --threads=N               [STCOMPARE_VALIDATE_THREADS, validate_threads]   default 16
#   --analyses=a,b,...        [STCOMPARE_VALIDATE_ANALYSES, validate_analyses] default: all six, from
#                             aki_iter, aki_fixed, merfish_affine, merfish_stalign, brain, celltypes
#   --pkg=DIR|installed       package to test (default: the directory above bench/); "installed" uses
#                             library(STcompare)
#   --published=DIR           [STCOMPARE_PUBLISHED_DIR, validate_published] the published results (default:
#                             published/ beside this script, else <pkg>/bench/published)
#   --build=install|load_all  [STCOMPARE_VALIDATE_BUILD] how a source directory is loaded: "install" (default)
#                             R CMD INSTALL into a temporary library; "load_all" devtools::load_all() with
#                             debug = FALSE and recompile = TRUE (writes object files into src/)
#   --out=DIR                 per-gene CSVs and our results (RDS); default <cache>/bench/validate-published,
#                             outside the repository
#   --results=FILE            markdown summary (default <pkg>/bench/validation-results.md). The section of
#                             this mode is replaced; it is written only when all six analyses ran, or when
#                             --results is given explicitly
#   --label=TEXT              note for the section heading (for example "pre-integration")
#   --allow-legacy            run "exported" even when the exported functions are still the legacy R code
#                             (refused by default: the legacy code needs about 200 CPU-hours here)
#   --reuse                   do not recompute: compare and report the results that an earlier run of the
#                             same mode saved in --out (for example after changing a tolerance)
#
# Inputs: <cache>/data-raw/inputs/{aki_rast,merfish_replicates_rast,brain_merfish_visium_rast,
# brain_celltype_rast}.rds, where <cache> is $STCOMPARE_DATA_CACHE or tools::R_user_dir("STcompare",
# "cache"). Build them with data-raw/download_data.R and data-raw/build_inputs_*.R (see data-raw/README.md).
#
# Comparisons per gene (published value = reference):
#   correlationCoef     |difference| <= 1e-12 (r is bounded by 1, so this is relative to its scale; the
#                       data-raw/ input checks use the same rule). The per-gene |difference| / |r| is
#                       reported as well.
#   pValueNaive         |ours / published - 1| <= 1e-10 (both 0 counts as equal)
#   permutations        the same number per gene (same screening decisions)
#   deltaStarX/Y        identical for every stored permutation
#   nullCorrelationsX/Y max |difference| <= 1e-9 x max |published| per gene and direction (every stored null)
#   exceedances         identical counts of |null| > |r| (the legacy count) and |null| >= |r| (the current
#                       count), our nulls against our r, published nulls against the published r
#   adjusted p          ours equal p.adjust((b + 1) / (B + 1), "BH") recomputed from the PUBLISHED nulls,
#                       b = #(|null| >= |r|), across the genes of the analysis (relative tolerance 1e-12;
#                       any count difference moves a p-value by at least 1e-4 relative)
# The stored p-values themselves were computed with the legacy definition b / B (strict >) and, in the
# iterative analyses, BH-adjusted across genes (in kidneyCorrelationNoIter.RData never adjusted: the legacy
# spatialCorrelationGeneExp() adjusted one gene at a time). They are checked for consistency with the
# published nulls only, to recognise rows that came from another run (merfishCorrelation.RData).
#
# Every exported call is also checked to leave .Random.seed and RNGkind() unchanged.

options(warn = 1, stringsAsFactors = FALSE)

# ---- options -------------------------------------------------------------------------------------

.args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, env, var, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit)) return(sub(paste0("^--", name, "="), "", hit[length(hit)]))
  if (!is.null(var) && exists(var, envir = globalenv(), inherits = FALSE)) return(get(var, envir = globalenv()))
  if (!is.null(env) && nzchar(Sys.getenv(env))) return(Sys.getenv(env))
  default
}
flag <- function(name) any(.args == paste0("--", name))
unknown <- setdiff(sub("=.*$", "", .args),
                   c("--mode", "--threads", "--analyses", "--pkg", "--build", "--out", "--results", "--label",
                     "--published", "--allow-legacy", "--reuse"))
if (length(unknown)) stop("unknown option(s): ", paste(unknown, collapse = ", "), " (see the header of this script)")

script_path <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) normalizePath(f[1], mustWork = FALSE) else NA_character_
})
mode <- match.arg(as.character(opt("mode", "STCOMPARE_VALIDATE_MODE", "validate_mode", "exported")),
                  c("exported", "internal"))
threads <- as.integer(opt("threads", "STCOMPARE_VALIDATE_THREADS", "validate_threads", 16L))
if (is.na(threads) || threads < 1L) stop("--threads must be a positive integer")
all_analyses <- c("aki_iter", "aki_fixed", "merfish_affine", "merfish_stalign", "brain", "celltypes")
analyses <- opt("analyses", "STCOMPARE_VALIDATE_ANALYSES", "validate_analyses", all_analyses)
if (length(analyses) == 1L) analyses <- strsplit(analyses, ",", fixed = TRUE)[[1]]
analyses <- trimws(analyses)
if (!length(analyses) || !all(analyses %in% all_analyses)) {
  stop("--analyses must be a comma-separated subset of: ", paste(all_analyses, collapse = ", "))
}
analyses <- all_analyses[all_analyses %in% analyses]
pkg <- opt("pkg", NULL, "validate_pkg",
           if (!is.na(script_path)) dirname(dirname(script_path)) else ".")
label <- opt("label", NULL, "validate_label", "")
build <- match.arg(as.character(opt("build", "STCOMPARE_VALIDATE_BUILD", "validate_build", "install")),
                   c("install", "load_all"))
allow_legacy <- flag("allow-legacy")
reuse <- flag("reuse")

cache_root <- local({
  root <- Sys.getenv("STCOMPARE_DATA_CACHE", unset = "")
  if (!nzchar(root)) root <- tools::R_user_dir("STcompare", which = "cache")
  path.expand(root)
})
inputs_dir <- file.path(cache_root, "data-raw", "inputs")
out_dir <- opt("out", NULL, "validate_out", file.path(cache_root, "bench", "validate-published"))
results_given <- any(grepl("^--results=", .args)) || exists("validate_results", envir = globalenv(), inherits = FALSE)

# ---- package -------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (identical(pkg, "installed")) {
    library(STcompare)
    pkg_desc <- sprintf("installed STcompare %s (%s)", utils::packageVersion("STcompare"),
                        dirname(system.file(package = "STcompare")))
    pkg_root <- NA_character_
  } else {
    pkg_root <- normalizePath(pkg, mustWork = TRUE)
    if (!file.exists(file.path(pkg_root, "DESCRIPTION")) ||
        !any(grepl("^Package:[[:space:]]*STcompare[[:space:]]*$", readLines(file.path(pkg_root, "DESCRIPTION"))))) {
      stop("--pkg must be the STcompare source directory (or \"installed\"): ", pkg_root)
    }
    version <- read.dcf(file.path(pkg_root, "DESCRIPTION"), "Version")[1]
    if (build == "install") {
      # a copy of the sources (no object files), installed with R's default (optimising) flags
      tmp <- tempfile("validate-published-")
      src <- file.path(tmp, "STcompare")
      lib <- file.path(tmp, "lib")
      dir.create(src, recursive = TRUE)
      dir.create(lib)
      parts <- intersect(c("DESCRIPTION", "NAMESPACE", "R", "src", "inst", "data"), list.files(pkg_root))
      file.copy(file.path(pkg_root, parts), src, recursive = TRUE)
      unlink(list.files(file.path(src, "src"), pattern = "[.](o|so|dll)$", full.names = TRUE))
      install_log <- file.path(tmp, "install.log")
      cat(sprintf("installing %s into a temporary library (R CMD INSTALL) ...\n", pkg_root))
      st <- system2(file.path(R.home("bin"), "R"),
                    c("CMD", "INSTALL", "--no-docs", "--no-test-load", "--no-multiarch", "-l", shQuote(lib), shQuote(src)),
                    stdout = install_log, stderr = install_log)
      if (!identical(st, 0L)) {
        cat(tail(readLines(install_log), 30), sep = "\n")
        stop("R CMD INSTALL failed (log: ", install_log, ")", call. = FALSE)
      }
      library(STcompare, lib.loc = lib)
      cc <- grep("stc_engine[.]cpp", readLines(install_log), value = TRUE)
      flags <- if (length(cc)) unique(regmatches(cc[1], gregexpr("(^| )-(O[0-9a-z]*|g|march=[^ ]+|mcpu=[^ ]+)(?= |$)",
                                                                   cc[1], perl = TRUE))[[1]]) else character(0)
      pkg_desc <- sprintf("STcompare %s from %s, installed with R CMD INSTALL (engine compiled with %s)", version,
                          pkg_root, if (length(flags)) paste(trimws(flags), collapse = " ") else "R's default flags")
    } else {
      devtools::load_all(pkg_root, quiet = TRUE, export_all = FALSE, debug = FALSE, recompile = TRUE)
      pkg_desc <- sprintf("STcompare %s loaded with devtools::load_all(\"%s\", debug = FALSE, recompile = TRUE)",
                          version, pkg_root)
    }
  }
})
# The cached inputs are SpatialExperiment objects: load the namespaces with their S4 methods before the first
# readRDS() (otherwise the first rownames() of a freshly unserialised object can return NULL).
for (p in c("SummarizedExperiment", "SpatialExperiment")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("package ", p, " is required")
}
results_md <- opt("results", NULL, "validate_results",
                  if (!is.na(pkg_root)) file.path(pkg_root, "bench", "validation-results.md") else NA_character_)
# The published results (.RData) are not part of the package: they are kept in bench/published.
published_dir <- opt("published", "STCOMPARE_PUBLISHED_DIR", "validate_published",
                     if (!is.na(script_path)) file.path(dirname(script_path), "published")
                     else file.path(if (is.na(pkg_root)) "." else pkg_root, "bench", "published"))
if (!dir.exists(published_dir)) stop("--published: directory of the published results not found: ", published_dir, call. = FALSE)
ns <- asNamespace("STcompare")
engine_correlate <- get(".stc_engine_correlate", envir = ns)

# A fingerprint of the engine sources, so that the results file says which build was tested.
engine_fingerprint <- local({
  if (is.na(pkg_root)) return(NA_character_)
  f <- sort(c(Sys.glob(file.path(pkg_root, "src", "stc_*")), Sys.glob(file.path(pkg_root, "R", "*.R"))))
  tf <- tempfile()
  on.exit(unlink(tf))
  writeLines(paste(basename(f), unname(tools::md5sum(f))), tf)
  substr(unname(tools::md5sum(tf)), 1, 12)
})
git_commit <- local({
  if (is.na(pkg_root)) return(NA_character_)
  out <- suppressWarnings(tryCatch(system2("git", c("-C", shQuote(pkg_root), "rev-parse", "--short", "HEAD"),
                                           stdout = TRUE, stderr = FALSE), error = function(e) character(0)))
  if (length(out) && is.null(attr(out, "status"))) out[1] else NA_character_
})

# The exported functions still running the legacy R code would take about 200 CPU-hours: refuse.
legacy_exported <- exists("matchingVariograms", envir = ns, inherits = FALSE) ||
  any(grepl("matchingVariograms|geoR::variog", deparse(body(get("viladomatCorrelation", envir = ns)))))
if (mode == "exported" && legacy_exported && !allow_legacy && !reuse) {
  stop("The exported functions of this build still run the legacy R implementation (matchingVariograms() ",
       "exists), which would need about 200 CPU-hours for these analyses. Use --mode=internal, or a build in ",
       "which the exported functions call the engine (or pass --allow-legacy).", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!is.na(pkg_root) && startsWith(paste0(normalizePath(out_dir), "/"), paste0(pkg_root, "/"))) {
  stop("--out must be outside the package source directory: ", out_dir, call. = FALSE)
}

# ---- the published analyses ----------------------------------------------------------------------
# Parameters from bench/published/scripts/ (dev/investigation/04-datasets-and-test-tiers.md, section 1.3). X is the
# first element of the authors' input list, Y the second. pair() takes the cached input (data-raw/).

grid11 <- c(0.01, 0.05, seq(0.1, 0.9, .1))  # the authors' expression (kidney and MERFISH scripts)
specs <- list(
  aki_iter = list(
    title = "AKI kidney Visium, iterative", ref = "kidneyCorrelation.RData", input = "aki_rast.rds",
    pair = function(inp) list(AKI_ctrl = inp$rast$AKI_ctrl, AKI_aki = inp$rast$AKI_aki),
    xy = "X = control (NL3), Y = AKI (IL3)", assay = "CPM", delta = grid11, nPermutations = c(100, 1000),
    iterative = TRUE, stored = "BH", script = "bench/published/scripts/visiumKidneySpatialCorrelation.R",
    authors_time = NA_character_, authors_hours = NA_real_, authors_threads = "22"),
  aki_fixed = list(
    title = "AKI kidney Visium, fixed B = 100", ref = "kidneyCorrelationNoIter.RData", input = "aki_rast.rds",
    pair = function(inp) list(AKI_ctrl = inp$rast$AKI_ctrl, AKI_aki = inp$rast$AKI_aki),
    xy = "X = control (NL3), Y = AKI (IL3)", assay = "CPM", delta = grid11, nPermutations = 100,
    iterative = FALSE, stored = "raw", script = "bench/published/scripts/KidneyNoIter.R",
    authors_time = NA_character_, authors_hours = NA_real_, authors_threads = "22"),
  merfish_affine = list(
    title = "MERFISH replicates, affine", ref = "merfishCorrelation_affine.RData",
    input = "merfish_replicates_rast.rds",
    pair = function(inp) list(target = inp$rast_affine$target, source = inp$rast_affine$source),
    xy = "X = S2R2 (target), Y = S2R3 (source, affine)", assay = NULL, delta = grid11,
    nPermutations = c(100, 1000), iterative = TRUE, stored = "BH",
    script = "bench/published/scripts/biological-replicates-example.R",
    authors_time = "7.65 h", authors_hours = 7.653494, authors_threads = "20 (MulticoreParam())"),
  merfish_stalign = list(
    title = "MERFISH replicates, STalign", ref = "merfishCorrelation.RData",
    input = "merfish_replicates_rast.rds",
    pair = function(inp) list(target = inp$rast$target, source = inp$rast$source),
    xy = "X = S2R2 (target), Y = S2R3 (source, STalign)", assay = NULL, delta = grid11,
    nPermutations = c(100, 1000), iterative = TRUE, stored = "BH", composite = TRUE,
    script = "bench/published/scripts/biological-replicates-example.R",
    authors_time = "16.8 h", authors_hours = 16.80652, authors_threads = "20 (MulticoreParam())"),
  brain = list(
    title = "Brain MERFISH vs Visium", ref = file.path("brain-MERFISH-10x-visium", "brainCorrelation.RData"),
    input = "brain_merfish_visium_rast.rds",
    pair = function(inp) list(MERFISH = inp$rast$MERFISH, Visium = inp$rast$Visium),
    xy = "X = MERFISH, Y = Visium", assay = "lognorm", delta = NULL, nPermutations = c(100, 1000),
    iterative = TRUE, stored = "BH", script = "bench/published/scripts/brain-MERFISH-10x-visium.R",
    authors_time = "1.78 h", authors_hours = 1.78, authors_threads = "22"),
  celltypes = list(
    title = "Brain cell types", ref = file.path("brain-MERFISH-10x-visium", "ctCorrelation.RData"),
    input = "brain_celltype_rast.rds",
    pair = function(inp) list(Visium = inp$rast$Visium, MERFISH = inp$rast$MERFISH),
    xy = "X = Visium, Y = MERFISH", assay = NULL, delta = NULL, nPermutations = c(100, 1000),
    iterative = TRUE, stored = "BH", script = "bench/published/scripts/brain-MERFISH-10x-visium.R",
    authors_time = "6.95 min", authors_hours = 6.95 / 60, authors_threads = "22")
)
SEED <- 0
ALPHA <- 0.05
MAXDIST <- 0.25
TOL_R <- 1e-12
TOL_PNAIVE <- 1e-10
TOL_NULL <- 1e-9
TOL_PADJ <- 1e-12

load_ref <- function(f) {
  path <- file.path(published_dir, f)
  if (!file.exists(path)) stop("published result not found: ", path)
  e <- new.env()
  n <- load(path, envir = e)
  get(n[1], envir = e)
}
input_cache <- new.env()
load_input <- function(f) {
  if (is.null(input_cache[[f]])) {
    path <- file.path(inputs_dir, f)
    if (!file.exists(path)) {
      stop("cached input not found: ", path, "\nBuild it from the repository root with Rscript data-raw/",
           c(aki_rast.rds = "build_inputs_aki.R", merfish_replicates_rast.rds = "build_inputs_merfish_replicates.R",
             brain_merfish_visium_rast.rds = "build_inputs_brain.R", brain_celltype_rast.rds = "build_inputs_brain.R")[[f]],
           " (see data-raw/README.md), or set STCOMPARE_DATA_CACHE.", call. = FALSE)
    }
    input_cache[[f]] <- readRDS(path)
  }
  input_cache[[f]]
}
as_num <- function(x) as.numeric(unlist(x, use.names = FALSE))
is_na_row <- function(x) length(x) == 0L || (length(x) == 1L && is.na(x))
rng_state <- function() list(kind = RNGkind(), seed = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
fmt_time <- function(s) {
  if (is.na(s)) return("")
  if (s < 120) sprintf("%.1f s", s) else if (s < 7200) sprintf("%.1f min", s / 60) else sprintf("%.2f h", s / 3600)
}

# ---- running one analysis ------------------------------------------------------------------------

# Our results in a common form: one element per gene of the published table, in its order.
run_exported <- function(spec, input, genes) {
  G <- length(genes)
  delta <- if (is.null(spec$delta)) NULL else rep(list(spec$delta), G)
  before <- rng_state()
  t0 <- Sys.time()
  if (spec$iterative) {
    res <- spatialCorrelationGeneExpIterPermutations(
      input, alpha = ALPHA, nPermutations = spec$nPermutations, deltaX = delta, deltaY = delta,
      maxDistPrctile = MAXDIST, returnPermutations = FALSE, assayName = spec$assay, nThreads = threads,
      BPPARAM = NULL, verbose = FALSE, seed = SEED, adjustMethod = "BH")
  } else {
    res <- spatialCorrelationGeneExp(
      input, nPermutations = spec$nPermutations, deltaX = delta, deltaY = delta, maxDistPrctile = MAXDIST,
      returnPermutations = FALSE, assayName = spec$assay, nThreads = threads, BPPARAM = NULL,
      verbose = FALSE, seed = SEED, adjustMethod = "BH")
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rng_ok <- identical(rng_state(), before)
  if (!identical(rownames(res), genes)) {
    if (!setequal(rownames(res), genes)) stop("the exported function returned other genes than the published table")
    res <- res[genes, , drop = FALSE]
  }
  nullX <- lapply(res$nullCorrelationsX, as_num)
  nullY <- lapply(res$nullCorrelationsY, as_num)
  B <- vapply(nullX, function(v) if (is_na_row(v)) NA_integer_ else length(v), 0L)
  list(r = unname(as.numeric(res$correlationCoef)), rY = unname(as.numeric(res$correlationCoef)),
       pnaive = as.numeric(res$pValueNaive), B = B, nullX = nullX, nullY = nullY,
       dsX = lapply(res$deltaStarX, as_num), dsY = lapply(res$deltaStarY, as_num),
       adjX = as.numeric(res$pValuePermuteX), adjY = as.numeric(res$pValuePermuteY),
       bX_engine = rep(NA_integer_, G), bY_engine = rep(NA_integer_, G),
       status = ifelse(is.na(B), "failed", "ok"), message = rep("", G),
       rounds = NULL, elapsed = elapsed, rng_ok = rng_ok)
}

run_internal <- function(spec, input, genes) {
  sx <- input[[1]]
  sy <- input[[2]]
  shared <- intersect(rownames(SpatialExperiment::spatialCoords(sx)), rownames(SpatialExperiment::spatialCoords(sy)))
  pos <- SpatialExperiment::spatialCoords(sx)[shared, , drop = FALSE]
  assay <- if (is.null(spec$assay)) 1L else spec$assay  # assayName = NULL: the first assay
  X <- t(as.matrix(SummarizedExperiment::assay(sx, assay)[genes, shared, drop = FALSE]))
  Y <- t(as.matrix(SummarizedExperiment::assay(sy, assay)[genes, shared, drop = FALSE]))
  G <- length(genes)
  delta <- if (is.null(spec$delta)) NULL else rep(list(spec$delta), G)
  before <- rng_state()
  t0 <- Sys.time()
  # correlationCoef and pValueNaive as the exported functions compute them
  ct <- lapply(seq_len(G), function(g) tryCatch(suppressWarnings(stats::cor.test(X[, g], Y[, g])),
                                                error = function(e) list(estimate = NA_real_, p.value = NA_real_)))
  r <- vapply(ct, function(o) unname(as.numeric(o$estimate)), 0)
  pnaive <- vapply(ct, function(o) as.numeric(o$p.value), 0)
  rY <- suppressWarnings(vapply(seq_len(G), function(g) stats::cor(Y[, g], X[, g]), 0))
  B <- spec$nPermutations
  res <- engine_correlate(X, Y, pos, deltaX = delta, deltaY = delta, nPermutations = B[1], seed = SEED,
                          maxDistPrctile = MAXDIST, nThreads = threads)
  state <- attr(res, "state")
  cur <- list(L = res$L, bX = res$bX, bY = res$bY, nullX = unclass(res$nullX), nullY = unclass(res$nullY),
              dsX = unclass(res$deltaStarX), dsY = unclass(res$deltaStarY), status = res$status,
              message = res$message)
  screen <- function(bX, bY, L, status, k) {
    t <- (ALPHA / B[k]) * 100  # the legacy screening threshold, as R/iterativePermutations.R computes it
    status == "ok" & !is.na(bX) & !is.na(bY) & (bX + 1) / (L + 1) < t & (bY + 1) / (L + 1) < t
  }
  rounds <- data.frame(round = 1L, nPermutations = B[1], genes = G, failed = sum(res$status != "ok"),
                       promoted = NA_integer_)
  promote <- if (spec$iterative && length(B) > 1L) which(screen(res$bX, res$bY, res$L, res$status, 1L)) else integer(0)
  if (spec$iterative && length(B) > 1L) {
    rounds$promoted[1] <- length(promote)
    for (k in seq_along(B)[-1]) {
      if (!length(promote)) break
      rk <- engine_correlate(state = state, units = promote, nPermutations = B[k], nThreads = threads)
      cur$L[promote] <- rk$L
      cur$bX[promote] <- rk$bX
      cur$bY[promote] <- rk$bY
      cur$nullX[promote] <- unclass(rk$nullX)
      cur$nullY[promote] <- unclass(rk$nullY)
      cur$dsX[promote] <- unclass(rk$deltaStarX)
      cur$dsY[promote] <- unclass(rk$deltaStarY)
      cur$status[promote] <- rk$status
      cur$message[promote] <- rk$message
      nxt <- if (k < length(B)) promote[screen(rk$bX, rk$bY, rk$L, rk$status, k)] else integer(0)
      rounds <- rbind(rounds, data.frame(round = k, nPermutations = B[k], genes = length(promote),
                                         failed = sum(rk$status != "ok"),
                                         promoted = if (k < length(B)) length(nxt) else NA_integer_))
      promote <- nxt
    }
  }
  pX <- ifelse(cur$status == "ok", (cur$bX + 1) / (cur$L + 1), NA_real_)
  pY <- ifelse(cur$status == "ok", (cur$bY + 1) / (cur$L + 1), NA_real_)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rng_ok <- identical(rng_state(), before)
  list(r = r, rY = rY, pnaive = pnaive, B = ifelse(cur$status == "ok", cur$L, NA_integer_),
       nullX = lapply(cur$nullX, as_num), nullY = lapply(cur$nullY, as_num),
       dsX = lapply(cur$dsX, as_num), dsY = lapply(cur$dsY, as_num),
       adjX = stats::p.adjust(pX, "BH"), adjY = stats::p.adjust(pY, "BH"),
       bX_engine = cur$bX, bY_engine = cur$bY, status = cur$status, message = cur$message,
       rounds = rounds, elapsed = elapsed, rng_ok = rng_ok, N = length(shared))
}

# ---- comparing with the published table -----------------------------------------------------------

first_diff <- function(a, b) {
  n <- min(length(a), length(b))
  w <- which(a[seq_len(n)] != b[seq_len(n)])
  if (length(w)) w[1] else NA_integer_
}

compare <- function(spec, ref, ours) {
  genes <- rownames(ref)
  G <- length(genes)
  r_pub <- unname(as.numeric(ref$correlationCoef))
  rows <- vector("list", G)
  for (i in seq_len(G)) {
    pnx <- as_num(ref$nullCorrelationsX[[i]])
    pny <- as_num(ref$nullCorrelationsY[[i]])
    pdx <- as_num(ref$deltaStarX[[i]])
    pdy <- as_num(ref$deltaStarY[[i]])
    onx <- ours$nullX[[i]]
    ony <- ours$nullY[[i]]
    odx <- ours$dsX[[i]]
    ody <- ours$dsY[[i]]
    pub_na <- is_na_row(pnx)
    our_na <- is_na_row(onx) || ours$status[i] != "ok"
    B_pub <- if (pub_na) NA_integer_ else length(pnx)
    B_our <- if (our_na) NA_integer_ else length(onx)
    row <- list(gene = genes[i], status = ours$status[i], message = ours$message[i],
                na_pub = pub_na, na_ours = our_na,
                r_pub = r_pub[i], r_ours = ours$r[i], r_absdiff = abs(ours$r[i] - r_pub[i]),
                r_reldiff = abs(ours$r[i] - r_pub[i]) / abs(r_pub[i]),
                pnaive_pub = ref$pValueNaive[i], pnaive_ours = ours$pnaive[i],
                pnaive_reldiff = if (isTRUE(ref$pValueNaive[i] == 0 && ours$pnaive[i] == 0)) 0 else
                  abs(ours$pnaive[i] / ref$pValueNaive[i] - 1),
                B_pub = B_pub, B_ours = B_our)
    if (!pub_na && !our_na) {
      n <- min(B_pub, B_our)
      s <- seq_len(n)
      row$dsX_first_diff <- if (length(odx) == B_our && length(pdx) == B_pub) first_diff(odx, pdx) else 0L
      row$dsY_first_diff <- if (length(ody) == B_our && length(pdy) == B_pub) first_diff(ody, pdy) else 0L
      row$dsX_ndiff <- sum(odx[s] != pdx[s])
      row$dsY_ndiff <- sum(ody[s] != pdy[s])
      row$nullX_absdiff <- max(abs(onx[s] - pnx[s]))
      row$nullY_absdiff <- max(abs(ony[s] - pny[s]))
      row$nullX_reldiff <- row$nullX_absdiff / max(abs(pnx[s]))
      row$nullY_reldiff <- row$nullY_absdiff / max(abs(pny[s]))
      row$nullX_worst <- which.max(abs(onx[s] - pnx[s]))
      row$nullY_worst <- which.max(abs(ony[s] - pny[s]))
      row$bgtX_pub <- sum(abs(pnx) > abs(r_pub[i]))
      row$bgtY_pub <- sum(abs(pny) > abs(r_pub[i]))
      row$bgeX_pub <- sum(abs(pnx) >= abs(r_pub[i]))
      row$bgeY_pub <- sum(abs(pny) >= abs(r_pub[i]))
      row$bgtX_ours <- sum(abs(onx) > abs(ours$r[i]))
      row$bgtY_ours <- sum(abs(ony) > abs(ours$rY[i]))
      row$bgeX_ours <- sum(abs(onx) >= abs(ours$r[i]))
      row$bgeY_ours <- sum(abs(ony) >= abs(ours$rY[i]))
      row$bX_engine <- ours$bX_engine[i]
      row$bY_engine <- ours$bY_engine[i]
    } else {
      for (nm in c("dsX_first_diff", "dsY_first_diff", "dsX_ndiff", "dsY_ndiff", "nullX_worst", "nullY_worst",
                   "bgtX_pub", "bgtY_pub", "bgeX_pub", "bgeY_pub", "bgtX_ours", "bgtY_ours", "bgeX_ours",
                   "bgeY_ours")) row[[nm]] <- NA_integer_
      for (nm in c("nullX_absdiff", "nullY_absdiff", "nullX_reldiff", "nullY_reldiff")) row[[nm]] <- NA_real_
      row$bX_engine <- ours$bX_engine[i]
      row$bY_engine <- ours$bY_engine[i]
      if (!pub_na) {  # published counts even when ours failed
        row$bgtX_pub <- sum(abs(pnx) > abs(r_pub[i])); row$bgtY_pub <- sum(abs(pny) > abs(r_pub[i]))
        row$bgeX_pub <- sum(abs(pnx) >= abs(r_pub[i])); row$bgeY_pub <- sum(abs(pny) >= abs(r_pub[i]))
      }
    }
    rows[[i]] <- row
  }
  d <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  # p-values: the current definition from the published nulls, BH across genes; the stored legacy values
  d$pX_raw_pub <- (d$bgeX_pub + 1) / (d$B_pub + 1)
  d$pY_raw_pub <- (d$bgeY_pub + 1) / (d$B_pub + 1)
  d$pX_adj_pub <- stats::p.adjust(d$pX_raw_pub, "BH")
  d$pY_adj_pub <- stats::p.adjust(d$pY_raw_pub, "BH")
  d$pX_adj_ours <- ours$adjX
  d$pY_adj_ours <- ours$adjY
  d$pX_stored <- as.numeric(ref$pValuePermuteX)
  d$pY_stored <- as.numeric(ref$pValuePermuteY)
  legX <- d$bgtX_pub / d$B_pub
  legY <- d$bgtY_pub / d$B_pub
  if (spec$stored == "BH") {
    legX <- stats::p.adjust(legX, "BH")
    legY <- stats::p.adjust(legY, "BH")
  }
  d$stored_consistent <- (abs(d$pX_stored - legX) <= 1e-12 | (is.na(d$pX_stored) & is.na(legX))) &
    (abs(d$pY_stored - legY) <= 1e-12 | (is.na(d$pY_stored) & is.na(legY)))
  d$stored_consistent_X <- abs(d$pX_stored - legX) <= 1e-12
  d$composite <- isTRUE(spec$composite) & !d$stored_consistent

  # per-gene verdicts
  both <- !d$na_pub & !d$na_ours
  rel_ok <- function(a, b, tol) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol * abs(b))
  d$ok_na <- d$na_pub == d$na_ours
  d$ok_r <- (is.na(d$r_pub) & is.na(d$r_ours)) | (!is.na(d$r_absdiff) & d$r_absdiff <= TOL_R)
  d$ok_pnaive <- (is.na(d$pnaive_pub) & is.na(d$pnaive_ours)) | (!is.na(d$pnaive_reldiff) & d$pnaive_reldiff <= TOL_PNAIVE)
  d$ok_B <- (is.na(d$B_pub) & is.na(d$B_ours)) | (!is.na(d$B_pub) & !is.na(d$B_ours) & d$B_pub == d$B_ours)
  d$ok_deltaStar <- !both | (d$ok_B & d$dsX_ndiff == 0 & d$dsY_ndiff == 0 & d$dsX_first_diff %in% NA &
                               d$dsY_first_diff %in% NA)
  d$ok_nulls <- !both | (d$ok_B & d$nullX_reldiff <= TOL_NULL & d$nullY_reldiff <= TOL_NULL)
  d$ok_count_gt <- !both | (d$bgtX_ours == d$bgtX_pub & d$bgtY_ours == d$bgtY_pub)
  d$ok_count_ge <- !both | (d$bgeX_ours == d$bgeX_pub & d$bgeY_ours == d$bgeY_pub &
                              (is.na(d$bX_engine) | (d$bX_engine == d$bgeX_ours & d$bY_engine == d$bgeY_ours)))
  d$ok_padj <- rel_ok(d$pX_adj_ours, d$pX_adj_pub, TOL_PADJ) & rel_ok(d$pY_adj_ours, d$pY_adj_pub, TOL_PADJ)
  d$padj_identical <- (d$pX_adj_ours == d$pX_adj_pub | (is.na(d$pX_adj_ours) & is.na(d$pX_adj_pub))) &
    (d$pY_adj_ours == d$pY_adj_pub | (is.na(d$pY_adj_ours) & is.na(d$pY_adj_pub)))
  checks <- c("ok_na", "ok_r", "ok_pnaive", "ok_B", "ok_deltaStar", "ok_nulls", "ok_count_gt", "ok_count_ge", "ok_padj")
  for (cc in checks) d[[cc]][is.na(d[[cc]])] <- FALSE
  d$ok <- Reduce(`&`, d[checks])
  d
}

# ---- report --------------------------------------------------------------------------------------

summarise <- function(id, spec, d, ours, N) {
  main <- !d$composite
  cnt <- function(col, sel = main) sum(!d[[col]][sel])
  mx <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  data.frame(
    analysis = id, title = spec$title, reference = basename(spec$ref), genes = nrow(d),
    compared = sum(main), composite = sum(d$composite), N = N,
    deltas = if (is.null(spec$delta)) 9L else length(spec$delta),
    schedule = paste(spec$nPermutations, collapse = "/"),
    pub_B = paste(sapply(spec$nPermutations, function(b) sum(d$B_pub == b, na.rm = TRUE)), collapse = "/"),
    our_B = paste(sapply(spec$nPermutations, function(b) sum(d$B_ours == b, na.rm = TRUE)), collapse = "/"),
    mm_na = cnt("ok_na"), mm_r = cnt("ok_r"), mm_pnaive = cnt("ok_pnaive"), mm_B = cnt("ok_B"),
    mm_deltaStar = cnt("ok_deltaStar"), mm_nulls = cnt("ok_nulls"), mm_count_gt = cnt("ok_count_gt"),
    mm_count_ge = cnt("ok_count_ge"), mm_padj = cnt("ok_padj"), mm_any = cnt("ok"),
    padj_identical = sum(d$padj_identical[main], na.rm = TRUE),
    max_null_rel = mx(pmax(d$nullX_reldiff[main], d$nullY_reldiff[main])),
    max_null_abs = mx(pmax(d$nullX_absdiff[main], d$nullY_absdiff[main])),
    max_r_abs = mx(d$r_absdiff[main]), max_r_rel = mx(d$r_reldiff[main]),
    max_pnaive_rel = mx(d$pnaive_reldiff[main]),
    failed = sum(d$na_ours), rng_ok = ours$rng_ok,
    seconds = ours$elapsed, threads = ours$threads, s_per_gene = ours$elapsed / nrow(d),
    authors_time = spec$authors_time, authors_threads = spec$authors_threads,
    speedup = if (is.na(spec$authors_hours)) NA_real_ else spec$authors_hours * 3600 / ours$elapsed,
    stringsAsFactors = FALSE)
}

details_lines <- function(id, d, max_genes = 10L) {
  bad <- which(!d$ok & !d$composite)
  if (!length(bad)) return(character(0))
  checks <- c(ok_na = "NA status", ok_r = "r", ok_pnaive = "naive p", ok_B = "permutations",
              ok_deltaStar = "deltaStar", ok_nulls = "nulls", ok_count_gt = "count >", ok_count_ge = "count >=",
              ok_padj = "BH p")
  out <- sprintf("- `%s`: %d gene(s) with a mismatch%s", id, length(bad),
                 if (length(bad) > max_genes) sprintf(" (first %d listed)", max_genes) else "")
  for (i in head(bad, max_genes)) {
    what <- names(checks)[!unlist(d[i, names(checks)])]
    out <- c(out, sprintf("  - %s: %s; B pub/ours %s/%s; deltaStar first diff X %s, Y %s (%s/%s differ); null rel diff X %.3g (perm %s), Y %.3g (perm %s); counts >= X %s/%s, Y %s/%s; r %.17g vs %.17g%s",
                          d$gene[i], paste(checks[what], collapse = ", "), d$B_pub[i], d$B_ours[i],
                          d$dsX_first_diff[i], d$dsY_first_diff[i], d$dsX_ndiff[i], d$dsY_ndiff[i],
                          d$nullX_reldiff[i], d$nullX_worst[i], d$nullY_reldiff[i], d$nullY_worst[i],
                          d$bgeX_pub[i], d$bgeX_ours[i], d$bgeY_pub[i], d$bgeY_ours[i], d$r_pub[i], d$r_ours[i],
                          if (nzchar(d$message[i])) paste0("; ", d$message[i]) else ""))
  }
  out
}

cat(sprintf("validate-published: mode %s, %d threads, analyses %s\n  package: %s\n  inputs: %s\n  output: %s\n",
            mode, threads, paste(analyses, collapse = ", "), pkg_desc, inputs_dir, out_dir))
set.seed(20261004)  # a defined .Random.seed, so that the RNG check is meaningful
summaries <- list()
runs <- list()
details <- character(0)
composite_notes <- character(0)
rounds_notes <- character(0)
for (id in analyses) {
  spec <- specs[[id]]
  cat(sprintf("\n== %s: %s (%s)\n", id, spec$title, spec$ref))
  ref <- load_ref(spec$ref)
  genes <- rownames(ref)
  inp <- load_input(spec$input)
  pair <- spec$pair(inp)
  if (!all(genes %in% rownames(pair[[1]])) || !all(genes %in% rownames(pair[[2]]))) {
    stop(id, ": the cached input lacks genes of the published table; rebuild it with data-raw/build_inputs_*.R")
  }
  input <- lapply(pair, function(se) se[genes, ])
  N <- length(intersect(rownames(SpatialExperiment::spatialCoords(input[[1]])),
                        rownames(SpatialExperiment::spatialCoords(input[[2]]))))
  cat(sprintf("   %d genes, %d shared pixels, %s, assay %s, deltas %s, nPermutations %s\n", length(genes), N, spec$xy,
              if (is.null(spec$assay)) "NULL (the first)" else spec$assay, if (is.null(spec$delta)) "default (0.1..0.9)" else paste(spec$delta, collapse = ","),
              paste(spec$nPermutations, collapse = ", ")))
  saved <- file.path(out_dir, sprintf("%s_%s_results.rds", id, mode))
  if (reuse) {
    if (!file.exists(saved)) stop(id, ": --reuse, but there are no saved results: ", saved, call. = FALSE)
    sv <- readRDS(saved)
    if (!identical(sv$genes, genes)) stop(id, ": the saved results are for other genes: ", saved, call. = FALSE)
    ours <- sv$ours
    runs[[id]] <- sv$run
    cat(sprintf("   reusing the results computed %s (%s; engine fingerprint %s)\n", format(sv$run$time, "%Y-%m-%d %H:%M"),
                sv$run$pkg, sv$run$fingerprint))
  } else {
    ours <- if (mode == "exported") run_exported(spec, input, genes) else run_internal(spec, input, genes)
    ours$threads <- threads
    runs[[id]] <- list(time = Sys.time(), pkg = pkg_desc, fingerprint = engine_fingerprint, threads = threads)
    saveRDS(list(spec = spec[setdiff(names(spec), "pair")], genes = genes, ours = ours, mode = mode, run = runs[[id]]),
            saved)
  }
  cat(sprintf("   %s %s on %d threads (%.3f s per gene); global RNG state unchanged: %s\n",
              if (reuse) "had run in" else "ran in", fmt_time(ours$elapsed), ours$threads,
              ours$elapsed / length(genes), ours$rng_ok))
  if (!is.null(ours$rounds)) {
    print(ours$rounds, row.names = FALSE)
    rounds_notes <- c(rounds_notes, sprintf("- `%s`: %s", id, paste(sprintf(
      "round %d (B = %d): %d genes%s", ours$rounds$round, ours$rounds$nPermutations, ours$rounds$genes,
      ifelse(is.na(ours$rounds$promoted), "", sprintf(", %d promoted", ours$rounds$promoted))), collapse = "; ")))
  }
  d <- compare(spec, ref, ours)
  csv <- file.path(out_dir, sprintf("%s_%s.csv", id, mode))
  utils::write.csv(d, csv, row.names = FALSE)
  s <- summarise(id, spec, d, ours, N)
  summaries[[id]] <- s
  details <- c(details, details_lines(id, d))
  if (any(d$composite)) {
    cm <- d[d$composite, ]
    zero <- d$bgtX_pub %in% 0 & d$bgtY_pub %in% 0
    composite_notes <- c(composite_notes, sprintf(paste(
      "- `%s`: the stored p-values of %d rows are not BH(b / B) of the published nulls: %d through pValuePermuteX",
      "(the number in dev/investigation/04) and %d more through pValuePermuteY only (their raw pX is 0). They are",
      "%d of the %d rows whose stored p-values can show this at all: the other %d rows have no exceedance in either",
      "direction, so their stored p is 0 under any adjustment%s. The stored values therefore cannot tell which rows",
      "were copied from the earlier run (`bench/published/scripts/biological-replicates-example.R`, lines 210-217); they show",
      "that the stored BH adjustment was computed over other raw p-values than the published nulls give. These rows",
      "are compared like the others but not counted in the table: %d match on every check (r, naive p, permutations,",
      "deltaStar, nulls, both counts, BH p), %d differ%s; max relative null difference %.2g."),
      id, nrow(cm), sum(!cm$stored_consistent_X, na.rm = TRUE), sum(cm$stored_consistent_X, na.rm = TRUE),
      nrow(cm), sum(!zero), sum(zero),
      if (sum(!zero) > nrow(cm)) sprintf(", and %s %s the largest raw p-values, which BH leaves unchanged",
                                         paste(d$gene[!zero & !d$composite], collapse = ", "),
                                         if (sum(!zero & !d$composite) == 1L) "has" else "have") else "",
      sum(cm$ok), sum(!cm$ok),
      if (any(!cm$ok)) paste0(" (", paste(head(cm$gene[!cm$ok], 10), collapse = ", "), ")") else "",
      max(pmax(cm$nullX_reldiff, cm$nullY_reldiff), na.rm = TRUE)))
  }
  cat(sprintf("   compared %d genes%s: mismatches NA %d, r %d, naive p %d, permutations %d, deltaStar %d, nulls %d, count> %d, count>= %d, BH p %d (any %d); max rel null diff %.2g; BH p bit-identical %d/%d\n",
              s$compared, if (s$composite) sprintf(" (+%d composite rows reported separately)", s$composite) else "",
              s$mm_na, s$mm_r, s$mm_pnaive, s$mm_B, s$mm_deltaStar, s$mm_nulls, s$mm_count_gt, s$mm_count_ge,
              s$mm_padj, s$mm_any, s$max_null_rel, s$padj_identical, s$compared))
  cat(sprintf("   per-gene table: %s\n", csv))
}
S <- do.call(rbind, summaries)
utils::write.csv(S, file.path(out_dir, sprintf("summary_%s.csv", mode)), row.names = FALSE)

# ---- markdown summary -----------------------------------------------------------------------------

md_table <- function(df) {
  num <- vapply(df, is.numeric, NA)
  df[] <- lapply(df, function(v) ifelse(is.na(v), "", as.character(v)))
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(ifelse(num, "---:", "---"), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
f2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "g", digits = 2))
tab1 <- data.frame(
  Analysis = sprintf("%s (`%s`)", S$title, S$reference),
  Genes = ifelse(S$composite > 0, sprintf("%d (+%d)", S$compared, S$composite), as.character(S$compared)),
  N = S$N,
  `B 100/1000 pub` = ifelse(grepl("/", S$schedule), S$pub_B, paste0(S$pub_B, "/-")),
  `B 100/1000 ours` = ifelse(grepl("/", S$schedule), S$our_B, paste0(S$our_B, "/-")),
  r = S$mm_r, `naive p` = S$mm_pnaive, B = S$mm_B, `δ*` = S$mm_deltaStar, nulls = S$mm_nulls,
  `count >` = S$mm_count_gt, `count ≥` = S$mm_count_ge, `BH p` = S$mm_padj, NA. = S$mm_na,
  `max rel null diff` = f2(S$max_null_rel), `max abs null diff` = f2(S$max_null_abs),
  `max abs Δr` = f2(S$max_r_abs), `max Δr / r` = f2(S$max_r_rel),
  check.names = FALSE, stringsAsFactors = FALSE)
names(tab1)[names(tab1) == "NA."] <- "NA rows"
tab2 <- data.frame(
  Analysis = S$title, Genes = S$genes, `Wall time` = vapply(S$seconds, fmt_time, ""), Threads = S$threads,
  `s / gene` = sprintf("%.3f", S$s_per_gene), `Authors' time` = ifelse(is.na(S$authors_time), "not recorded", S$authors_time),
  `Authors' workers` = S$authors_threads,
  Speedup = ifelse(is.na(S$speedup), "", sprintf("%.0f×", S$speedup)),
  `RNG state unchanged` = ifelse(S$rng_ok, "yes", "NO"),
  check.names = FALSE, stringsAsFactors = FALSE)

section_key <- paste0("validate-published:", mode)
heading <- sprintf("## Mode `%s`%s", mode, if (nzchar(label)) paste0(" — ", label) else "")
sys <- Sys.info()
section <- c(
  sprintf("<!-- BEGIN %s -->", section_key),
  heading, "",
  if (reuse) {
    times <- do.call(c, lapply(runs, `[[`, "time"))
    sprintf("Re-reported %s (`--reuse`) from results computed %s to %s on %s (%s %s, %s), %s; %s threads.",
            format(Sys.time(), "%Y-%m-%d %H:%M %Z"), format(min(times), "%Y-%m-%d %H:%M"),
            format(max(times), "%H:%M %Z"),
            sys[["nodename"]], sys[["sysname"]], sys[["release"]], R.version$platform, R.version.string,
            paste(unique(vapply(runs, function(r) as.character(r$threads), "")), collapse = ", "))
  } else {
    sprintf("Run %s on %s (%s %s, %s), %s; %d threads.", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), sys[["nodename"]],
            sys[["sysname"]], sys[["release"]], R.version$platform, R.version.string, threads)
  },
  sprintf("Package: %s; engine/R source fingerprint `%s`%s.",
          paste(unique(vapply(runs, `[[`, "", "pkg")), collapse = "; "),
          paste(unique(vapply(runs, `[[`, "", "fingerprint")), collapse = ", "),
          if (!is.na(git_commit)) sprintf(", git HEAD `%s` (uncommitted changes are not part of this id)", git_commit) else ""),
  sprintf("Call: `Rscript bench/validate-published.R --mode=%s --threads=%s%s`. Per-gene tables: `%s` (outside the repository).",
          mode, paste(unique(vapply(runs, function(r) as.character(r$threads), "")), collapse = ","),
          if (nzchar(label)) sprintf(" --label=\"%s\"", label) else "", out_dir),
  "",
  "Mismatching genes by check (0 everywhere = the published analysis is reproduced). Genes = genes compared;",
  "\"(+k)\": rows of `merfishCorrelation.RData` whose stored p-values do not follow from its stored nulls",
  "(the table was patched with rows of an earlier run; see below). They are compared but counted separately.",
  "",
  md_table(tab1), "",
  md_table(tab2), "",
  "Speedup = the authors' reported wall time / ours (their 20-22 workers on an unrecorded machine, our threads here).",
  "",
  if (length(rounds_notes)) c("Screening rounds (ours):", "", rounds_notes, "") else NULL,
  if (length(composite_notes)) c("Stored p-values that do not follow from the stored nulls:", "", composite_notes, "") else NULL,
  if (length(details)) c("Mismatches:", "", details, "") else c("No mismatches.", ""),
  sprintf("<!-- END %s -->", section_key))

header <- c(
  "# Validation against the published analyses",
  "",
  "Written by `bench/validate-published.R` (see `bench/README.md`). Each analysis whose results are in",
  "`bench/published` is re-run from inputs rebuilt from the public downloads (`data-raw/`) with the authors'",
  "parameters (`bench/published/scripts/`), and compared gene by gene with the stored table:",
  "",
  "- `r`: |Δ correlationCoef| ≤ 1e-12; `naive p`: relative difference of pValueNaive ≤ 1e-10;",
  "- `B`: the same number of permutations per gene (the same screening decisions);",
  "- `δ*`: deltaStarX and deltaStarY identical for every stored permutation;",
  "- `nulls`: every stored null correlation within 1e-9 × the largest published |null| of that gene and direction;",
  "- `count >` / `count ≥`: identical numbers of |null| > |r| (the legacy count) and |null| ≥ |r| (the current one);",
  "- `BH p`: our final pValuePermuteX/Y equal p.adjust((b + 1) / (B + 1), \"BH\") recomputed from the published",
  "  nulls (relative tolerance 1e-12; the stored values used the legacy b / B and are not compared directly).",
  "",
  "On macOS arm64 the MERFISH correlationCoef differ from the stored ones by up to 2.6e-14 (10-11 of 483 are",
  "bit-identical) because R's cor() accumulates in long double, which is plain double there: the MERFISH analyses",
  "were evidently computed with extended precision (on Linux arm64, 467 and 469 of 483 are bit-identical, max",
  "difference 1.1e-16), while the AKI, brain and cell-type values reproduce best on macOS arm64. The engine is",
  "not involved in r or the naive p-value.",
  "",
  "All six analyses reproduce on macOS arm64, where the AKI, brain and cell-type results were computed. On Linux",
  "arm64 (Docker, Bioconductor 3.22, GCC 13) the MERFISH and cell-type analyses reproduce as well, but AKI and brain",
  "do not: glibc's hypot() puts 4 (AKI) and 16 (brain) lattice pairs that lie within a few ulp of the last",
  "variogram bin edge on the other side of it. geoR::variog() bins them the same way there, so the legacy R",
  "code does not reproduce these two published analyses on Linux either.",
  "")
if (!is.na(results_md) && (results_given || identical(analyses, all_analyses))) {
  old <- if (file.exists(results_md)) readLines(results_md, warn = FALSE) else header
  b <- grep(sprintf("<!-- BEGIN %s -->", section_key), old, fixed = TRUE)
  e <- grep(sprintf("<!-- END %s -->", section_key), old, fixed = TRUE)
  new <- if (length(b) == 1L && length(e) == 1L && e > b) {
    c(old[seq_len(b - 1L)], section, old[-seq_len(e)])
  } else {
    c(old, if (length(old) && nzchar(old[length(old)])) "" else NULL, section)
  }
  dir.create(dirname(results_md), recursive = TRUE, showWarnings = FALSE)
  writeLines(new, results_md)
  cat(sprintf("\nsummary written to %s\n", results_md))
} else {
  cat("\n(summary file not written: a subset of the analyses ran; pass --results=FILE to write it anyway)\n")
}
cat("\n", paste(md_table(tab1), collapse = "\n"), "\n\n", paste(md_table(tab2), collapse = "\n"), "\n", sep = "")
if (length(details)) cat("\nMismatches:\n", paste(details, collapse = "\n"), "\n", sep = "")
if (length(composite_notes)) cat("\n", paste(composite_notes, collapse = "\n"), "\n", sep = "")
invisible(S)
