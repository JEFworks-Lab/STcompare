#!/usr/bin/env Rscript
# data-raw/build_test_fixtures.R
#
# Builds the committed test fixtures from the CURRENT (legacy, pure R) implementation of STcompare:
#
#   tier0  tests/testthat/fixtures/kernel_fixture.rds       seven small cases with the step-by-step
#          intermediates of viladomatCorrelation()/matchingVariograms() and the expected nulls, deltaStar
#          and raw p-values in both directions; spatialSimilarity() outputs for speKidney; edge-case notes
#   tier1  tests/testthat/fixtures/calibration_fixture.rds  simRanPatternRasts as plain matrices, all 4950
#          independent null pairs, the shipped reference rates, the mixing recipe for power, and the
#          reference results (p-values, tail counts, deltaStar) of the current code for the 60 jobs of
#          tests/testthat/test-calibration.R
#   tier2  tests/testthat/fixtures/realistic_fixture.rds    AKI kidney Visium (35 genes x 311 pixels) and
#          MERFISH-vs-Visium brain (30 genes x 2170 pixels) subsets with golden nulls/deltaStar copied
#          from inst/extdata (published results)
#   verify re-runs the current code on every realistic gene (B = $STCOMPARE_VERIFY_B, default 100,
#          in parallel) and on the kernel cases, and stops unless everything matches
#
# Usage (from the repository root, on a Unix-alike; on Windows set STCOMPARE_BUILD_WORKERS=1):
#   Rscript data-raw/build_test_fixtures.R [tier0] [tier1] [tier2] [verify] [compare | force]
# Default: all four tiers. Each new fixture is compared with the committed file, ignoring `meta` (provenance,
# creation time, runtimes); it is written only if its content changed, or always with "force". With "compare"
# nothing is written: the differences are listed and the script stops with an error if there are any.
# tier0 and tier2 read rasterized inputs from the cache (<cache>/data-raw/inputs/, see download_data.R), so run
# data-raw/build_inputs_brain.R (tier0, tier2) and data-raw/build_inputs_aki.R (tier2) first; the cached inputs
# are checked against the published correlation coefficients before use. Parallel parts use
# $STCOMPARE_BUILD_WORKERS processes (default: detectCores() - 2).
#
# The fixtures must be built from the legacy R implementation: every intermediate is checked for
# bit-identity against the package functions (stopifnot), and content hashes of the code are recorded in meta.
# They encode the arithmetic of the build machine (see tests/testthat/helper-platform.R), whose signature is
# stored in meta$platform_signature.
rng_default <- c("Mersenne-Twister", "Inversion", "Rejection")
rng_at_startup <- RNGkind()
args <- commandArgs(trailingOnly = TRUE)
tiers_all <- c("tier0", "tier1", "tier2", "verify")
if (length(setdiff(args, c(tiers_all, "compare", "force"))))
  stop("Unknown argument(s): ", paste(setdiff(args, c(tiers_all, "compare", "force")), collapse = ", "),
       "\nUsage: Rscript data-raw/build_test_fixtures.R [tier0] [tier1] [tier2] [verify] [compare | force]")
tiers <- intersect(args, tiers_all)
if (!length(tiers)) tiers <- tiers_all
compare_only <- "compare" %in% args
force_write <- "force" %in% args
if (compare_only && force_write) stop("'compare' and 'force' are mutually exclusive")
source("data-raw/download_data.R")
source(file.path("tests", "testthat", "helper-platform.R"))      # stc_platform_signature(), shared with the tests
suppressPackageStartupMessages({ library(SpatialExperiment); library(SummarizedExperiment) })
devtools::load_all(".", quiet = TRUE, helpers = FALSE)

# The permutation stream of the fixtures is defined for R's default generators. Pin them whatever the session
# (for example ~/.Rprofile) chose, and again before every set.seed() below.
use_default_rng <- function() RNGkind(rng_default[1], rng_default[2], rng_default[3])
use_default_rng()
if (!identical(rng_at_startup, rng_default))
  message("Note: the session started with RNGkind() = ", paste(rng_at_startup, collapse = " / "),
          "; the build uses R's defaults (", paste(rng_default, collapse = " / "), ").")

fixdir <- file.path("tests", "testthat", "fixtures")
workers <- as.integer(Sys.getenv("STCOMPARE_BUILD_WORKERS", max(1L, parallel::detectCores() - 2L)))
if (.Platform$OS.type == "windows" && workers > 1) { message("Windows: running serially"); workers <- 1L }
`%||%` <- function(a, b) if (is.null(a)) b else a
secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))
compare_failed <- FALSE

# mclapply on Unix-alikes; lapply (errors captured the same way) with one worker or on Windows.
par_lapply <- function(X, FUN) {
  if (workers > 1) return(parallel::mclapply(X, FUN, mc.cores = workers, mc.preschedule = FALSE))
  lapply(X, function(x) tryCatch(FUN(x), error = function(e) structure(conditionMessage(e), class = "try-error")))
}

rng_semantics <- paste(
  "Permutation order: Mersenne-Twister (Inversion, Rejection sampling) from set.seed(seed) in the calling session;",
  "if N > 1000 the variogram subsample ids <- sample(N, 1000) is drawn first, then permutation i is",
  "X[sample.int(N, N)] for i = 1..B in order (sample(X, length(X))). Both directions use the same seed and",
  "therefore the same index permutations. Noise: permutation i is processed inside a BiocParallel task,",
  "where RNGkind is L'Ecuyer-CMRG (Inversion); matchingVariograms() calls set.seed(seed + i) and draws",
  "rnorm(N) once per delta, in delta order. Results do not depend on nThreads or on the BiocParallel backend.",
  "Prefix property: permutation i depends only on (seed, i), so the first k nulls of a B-permutation run",
  "equal a k-permutation run.")

# Differences between two fixture objects: added, removed or changed leaves (named lists are walked by name,
# unnamed lists by position; data frames, matrices and vectors are compared as a whole).
fixture_diff <- function(a, b, path = "") {
  node <- function(x) is.list(x) && !is.data.frame(x)
  key <- function(x) if (is.null(names(x)) || anyDuplicated(names(x)) || any(!nzchar(names(x)))) NULL else names(x)
  if (node(a) && node(b) && identical(is.null(key(a)), is.null(key(b))) && (!is.null(key(a)) || length(a) == length(b))) {
    ks <- if (is.null(key(a))) seq_along(a) else union(names(a), names(b))
    out <- character(0)
    for (k in ks) {
      p <- if (is.numeric(k)) sprintf("%s[[%d]]", path, k) else if (nzchar(path)) paste0(path, "$", k) else k
      if (is.character(k) && !k %in% names(b)) out <- c(out, paste("removed:", p))
      else if (is.character(k) && !k %in% names(a)) out <- c(out, paste("added:  ", p))
      else out <- c(out, fixture_diff(a[[k]], b[[k]], p))
    }
    return(out)
  }
  if (identical(a, b)) character(0) else paste("changed:", path)
}

save_fixture <- function(x, f) {
  p <- file.path(fixdir, f)
  old <- if (file.exists(p)) readRDS(p) else NULL
  content <- function(z) z[setdiff(names(z), "meta")]
  d <- if (is.null(old)) "added:   (no committed file)" else fixture_diff(content(old), content(x))
  cat(sprintf("  %s vs committed file, ignoring meta: %s\n", f,
              if (length(d)) sprintf("%d difference(s)", length(d)) else "identical"))
  if (length(d)) cat(paste0("    ", utils::head(d, 60), "\n"), if (length(d) > 60) sprintf("    ... and %d more\n", length(d) - 60), sep = "")
  # Changed input data (not added cases) usually means a stale or modified cache: the correlation checks
  # against inst/extdata cannot see, e.g., a rescaled gene, but the committed fixture can.
  changed_inputs <- grep("^changed: (.*\\$input\\$(X|Y|pos)|pairs\\$[^$]+\\$(X|Y|pos)|coords|fields)$", d, value = TRUE)
  if (length(changed_inputs) && !force_write && !compare_only)
    stop("The input data differ from those stored in the committed ", f, ":\n  ", paste(changed_inputs, collapse = "\n  "),
         "\nIf the cached inputs are stale, rebuild them with data-raw/build_inputs_*.R; if the change is intended, ",
         "rerun with 'force'. Nothing was written.", call. = FALSE)
  if (compare_only) {
    if (length(d)) compare_failed <<- TRUE
    cat("  compare mode: nothing written\n")
    return(invisible(!length(d)))
  }
  if (!length(d) && !force_write) {
    cat(sprintf("  %s unchanged; not rewritten ('force' rewrites it, e.g. to refresh meta)\n", p))
    return(invisible(TRUE))
  }
  dir.create(fixdir, recursive = TRUE, showWarnings = FALSE)
  con <- xzfile(p, "wb", compression = 9)
  base::saveRDS(x, con)                  # (SummarizedExperiment masks saveRDS)
  close(con)
  cat(sprintf("wrote %s (%.1f KB)\n", p, file.size(p) / 1024))
  invisible(p)
}

# Cached inputs: must have been written by `script`; callers also check them against the published results.
need_input <- function(name, script) {
  p <- stc_input_file(name)
  if (!file.exists(p)) stop("Missing cached input ", p, "\nRun first:  Rscript ", script, call. = FALSE)
  x <- readRDS(p)
  if (!identical(x$meta$script, script))
    stop("Cached input ", p, " was not written by ", script, " (meta$script = ", deparse(x$meta$script),
         "). Rebuild it:  Rscript ", script, call. = FALSE)
  x
}
load_extdata <- function(f) { e <- new.env(); n <- load(file.path("inst", "extdata", f), envir = e); get(n, envir = e) }
# Stop unless the cached input reproduces the published correlation coefficients and naive p-values of every
# published gene: catches stale inputs (older builder, other SEraster/sf versions) and modified files.
check_input <- function(label, X, Y, ref) {
  genes <- rownames(ref)
  r <- vapply(genes, function(g) stats::cor(X[g, ], Y[g, ]), 0)
  p <- vapply(genes, function(g) stats::cor.test(X[g, ], Y[g, ])$p.value, 0)
  dr <- max(abs(r - ref$correlationCoef))
  dp <- max(abs(p / ref$pValueNaive - 1))
  cat(sprintf("  %s: cached input vs published (%d genes): max|dr| = %.2g, max rel d(naive p) = %.2g\n",
              label, length(genes), dr, dp))
  if (!(dr <= 1e-12 && dp <= 1e-10))
    stop(label, ": the cached input does not reproduce the published correlations (max|dr| = ", dr,
         "); rebuild it with data-raw/build_inputs_*.R", call. = FALSE)
}

## =============================================================================================
## Tier 0: step-by-step capture of the legacy algorithm (checked bit-for-bit against the package)
## =============================================================================================

# One direction of viladomatCorrelation(): X is permuted, Y is not used here.
capture_direction <- function(X, pos, delta, maxDistPrctile, B, seed, detail_vectors) {
  N <- length(X); lat <- pos[, 1]; long <- pos[, 2]        # package: lat <- data[, 3], long <- data[, 4]
  use_default_rng()
  set.seed(seed)
  ids <- if (N > 1000) sample(N, 1000) else seq_len(N)
  prctile <- quantile(dist(cbind(lat[ids], long[ids])), probs = maxDistPrctile)
  target <- geoR::variog(data = X[ids], coords = cbind(long[ids], lat[ids]), max.dist = prctile,
                         option = "bin", messages = FALSE)
  perm_index <- vapply(seq_len(B), function(i) sample.int(N, N), integer(N))   # == the sample(X, N) stream
  on.exit(use_default_rng(), add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  # every permutation through the exported kernel, exactly as viladomatCorrelation() calls it
  mv <- lapply(seq_len(B), function(i) suppressWarnings(
    matchingVariograms(X[perm_index[, i]], long, lat, delta, target, prctile, ids, i, seed = seed + i)))
  # permutation 1 replayed step by step
  xr <- X[perm_index[, 1]]
  set.seed(seed + 1)
  per <- vector("list", length(delta))
  for (k in seq_along(delta)) {
    fit <- suppressWarnings(locfit::locfit(xr ~ locfit::lp(long, lat, nn = delta[k], deg = 0), kern = "gauss", maxk = 300))
    xd <- fitted(fit)
    v1 <- geoR::variog(data = xd[ids], coords = cbind(long[ids], lat[ids]), option = "bin", max.dist = prctile, messages = FALSE)
    bet <- as.numeric(lm(target$v ~ 1 + v1$v)$coefficients)
    noise <- rnorm(length(xd))
    hat <- xd * sqrt(abs(bet[2])) + noise * sqrt(abs(bet[1]))
    v2 <- geoR::variog(data = hat[ids], coords = cbind(long[ids], lat[ids]), option = "bin", max.dist = prctile, messages = FALSE)
    stopifnot(identical(v1$n, target$n), identical(v2$n, target$n), identical(v1$u, target$u))
    per[[k]] <- list(delta = delta[k], variog_fitted_v = v1$v, lm_coef = bet, variog_hat_v = v2$v,
                     rss = sum((v2$v - target$v)^2))
    if (detail_vectors) {
      # reference for a smoother that is not bit-compatible with locfit's default adaptive tree:
      # the exact local-constant Gaussian-kernel fit at every data point (locfit ev = dat())
      xd_exact <- fitted(suppressWarnings(locfit::locfit(xr ~ locfit::lp(long, lat, nn = delta[k], deg = 0),
                                                         kern = "gauss", maxk = 300, ev = locfit::dat())))
      per[[k]] <- c(per[[k]], list(fitted = as.numeric(xd), fitted_exact_evdat = as.numeric(xd_exact),
                                   noise = noise, hat = as.numeric(hat)))
      per[[k]]$hat_raw <- hat                                  # (self-check only, removed below)
    }
  }
  rss1 <- vapply(per, function(p) p$rss, 0)
  ks <- which.min(rss1)
  checks <- c(replay_rss_identical_to_matchingVariograms = identical(rss1, mv[[1]]$residus),
              replay_argmin_identical = identical(ks, mv[[1]]$delta.star.id),
              replay_hat_identical = if (detail_vectors) identical(per[[ks]]$hat_raw, mv[[1]]$hat.X.delta.star) else NA)
  hat_star <- if (detail_vectors) as.numeric(per[[ks]]$hat) else NULL
  if (detail_vectors) per <- lapply(per, function(p) { p$hat_raw <- NULL; p })
  list(N = N, ids = ids, prctile = unname(prctile), maxDistPrctile = maxDistPrctile,
       target_variog = list(u = target$u, v = target$v, n = target$n, bins.lim = target$bins.lim),
       perm_index = perm_index,
       residus = vapply(mv, function(m) m$residus, numeric(length(delta))),     # length(delta) x B
       delta_star_id = vapply(mv, function(m) m$delta.star.id, 1L),
       detail = list(i = 1L, noise_seed = seed + 1, per_delta = per, residus = rss1, delta_star_id = ks,
                     hat_star = hat_star),
       .permutations = do.call(cbind, lapply(mv, function(m) m$hat.X.delta.star)),
       .checks = checks)
}

perm_fingerprint <- function(P) rbind(mean = colMeans(P), sd = apply(P, 2, sd), min = apply(P, 2, min), max = apply(P, 2, max))

make_case <- function(name, X, Y, pos, pixel = NULL, deltaX = seq(0.1, 0.9, 0.1), deltaY = deltaX, B = 10, seed = 0,
                      maxDistPrctile = 0.25, detail_vectors = TRUE, notes = "", extra_input = list()) {
  t0 <- Sys.time()
  X <- unname(as.numeric(X)); Y <- unname(as.numeric(Y)); pos <- unname(as.matrix(pos)); N <- length(X)
  fwd <- capture_direction(X, pos, deltaX, maxDistPrctile, B, seed, detail_vectors)
  rev <- capture_direction(Y, pos, deltaY, maxDistPrctile, B, seed, FALSE)
  use_default_rng()
  pf <- suppressWarnings(viladomatCorrelation(data.frame(X = X, Y = Y, x = pos[, 1], y = pos[, 2]), deltaX,
                                              maxDistPrctile, B, nThreads = 1, seed = seed))
  use_default_rng()
  pr <- suppressWarnings(viladomatCorrelation(data.frame(X = Y, Y = X, x = pos[, 1], y = pos[, 2]), deltaY,
                                              maxDistPrctile, B, nThreads = 1, seed = seed))
  use_default_rng()
  sc <- suppressWarnings(spatialCorrelation(X, Y, pos, nPermutations = B, deltaX = deltaX, deltaY = deltaY,
                                            maxDistPrctile = maxDistPrctile, returnPermutations = TRUE, seed = seed))
  el <- secs(t0)
  use_default_rng()
  # self-checks: the captured intermediates reproduce the package bit for bit
  chk <- c(
    perm_stream = identical(X[fwd$perm_index[, 1]],
                            { set.seed(seed); if (N > 1000) invisible(sample(N, 1000)); sample(X, N) }),
    perm_index_same_in_both_directions = identical(fwd$perm_index, rev$perm_index),
    setNames(fwd$.checks, paste0("fwd_", names(fwd$.checks))),
    setNames(rev$.checks, paste0("rev_", names(rev$.checks))),
    fwd_all_permutations_identical_to_viladomat = identical(fwd$.permutations, unname(pf$permutations)),
    rev_all_permutations_identical_to_viladomat = identical(rev$.permutations, unname(pr$permutations)),
    fwd_deltaStar_identical = identical(deltaX[fwd$delta_star_id], as.numeric(pf$deltaStar)),
    rev_deltaStar_identical = identical(deltaY[rev$delta_star_id], as.numeric(pr$deltaStar)),
    spatialCorrelation_nullX_identical = identical(as.numeric(sc$nullCorrelationsX[[1]]), as.numeric(pf$nullCorGlobal)),
    spatialCorrelation_nullY_identical = identical(as.numeric(sc$nullCorrelationsY[[1]]), as.numeric(pr$nullCorGlobal)),
    spatialCorrelation_permX_identical = identical(unname(sc$permutationsX[[1]]), unname(pf$permutations)),
    spatialCorrelation_permY_identical = identical(unname(sc$permutationsY[[1]]), unname(pr$permutations)),
    spatialCorrelation_deltaStar_identical = identical(sc$deltaStarX[[1]], pf$deltaStar) && identical(sc$deltaStarY[[1]], pr$deltaStar),
    spatialCorrelation_p_identical = identical(sc$pValuePermuteX, pf$pValueGlobal) && identical(sc$pValuePermuteY, pr$pValueGlobal))
  cat(sprintf("  case %-24s N=%4d B=%2d deltas=%d/%d seed=%g maxDist=%.2f  r=%+.4f pX=%.2f pY=%.2f  (%.1fs)  self-checks: %d/%d TRUE\n",
              name, N, B, length(deltaX), length(deltaY), seed, maxDistPrctile, sc$correlationCoef, sc$pValuePermuteX,
              sc$pValuePermuteY, el, sum(chk, na.rm = TRUE), sum(!is.na(chk))))
  if (!all(chk, na.rm = TRUE)) { print(chk[!chk %in% TRUE]); stop("self-check failed for case ", name) }
  strip <- function(cap) { cap$.permutations <- NULL; cap$.checks <- NULL; cap }
  expected_dir <- function(p, xv, yv, keep_perm1) {
    r_obs <- as.vector(cor(xv, yv))
    list(r_obs = r_obs, deltaStar = as.numeric(p$deltaStar), deltaStarMedian = p$deltaStarMedian,
         nullCor = as.numeric(p$nullCorGlobal), nExtreme = sum(abs(p$nullCorGlobal) > abs(r_obs)),
         # legacy definition b / B, independent of the package's current p-value formula
         pValue = sum(abs(p$nullCorGlobal) > abs(r_obs)) / length(p$nullCorGlobal),
         perm_fingerprint = perm_fingerprint(p$permutations),
         perm1 = if (keep_perm1) as.numeric(p$permutations[, 1]) else NULL)
  }
  list(name = name, notes = notes,
       input = c(list(X = X, Y = Y, pos = pos, pixel = pixel, delta = if (identical(deltaX, deltaY)) deltaX else NULL,
                      maxDistPrctile = maxDistPrctile, nPermutations = B, seed = seed, deltaX = deltaX, deltaY = deltaY),
                 extra_input),
       intermediate = list(forward = strip(fwd), reverse = strip(rev)),
       expected = list(
         forward = expected_dir(pf, X, Y, N <= 1000),
         reverse = expected_dir(pr, Y, X, N <= 1000),
         spatialCorrelation = list(correlationCoef = unname(sc$correlationCoef), pValueNaive = sc$pValueNaive,
                                   pValuePermuteX = sc$pValuePermuteX, pValuePermuteY = sc$pValuePermuteY,
                                   deltaStarMedianX = sc$deltaStarMedianX, deltaStarMedianY = sc$deltaStarMedianY)),
       self_checks = chk, .runtime = el)
}

if ("tier0" %in% tiers) {
  cat("== Tier 0: kernel fixture ==\n")
  t_tier <- Sys.time()
  data(speKidney)
  rk <- SEraster::rasterizeGeneExpression(speKidney, assay_name = "counts", resolution = 0.2, fun = "mean", square = FALSE)
  shAB <- intersect(rownames(spatialCoords(rk$A)), rownames(spatialCoords(rk$B)))
  shAC <- intersect(rownames(spatialCoords(rk$A)), rownames(spatialCoords(rk$C)))
  getv <- function(s, px) as.numeric(assay(s)[1, px])
  b <- need_input("brain_merfish_visium_rast.rds", "data-raw/build_inputs_brain.R")
  bc <- load_extdata("brain-MERFISH-10x-visium/brainCorrelation.RData")
  px <- b$shared
  bX <- as.matrix(assay(b$rast$MERFISH, "lognorm")[rownames(bc), px])
  bY <- as.matrix(assay(b$rast$Visium, "lognorm")[rownames(bc), px])
  check_input("brain", bX, bY, bc)
  cases <- list()
  cases$kidney_AB <- make_case("kidney_AB", getv(rk$A, shAB), getv(rk$B, shAB), spatialCoords(rk$A)[shAB, ], pixel = shAB,
    notes = paste("data(speKidney) A vs B rasterized with SEraster (resolution 0.2, mean, hexagons), shared pixels;",
                  "negative control r ~ -0.95. Hexagonal lattice: many pixel pairs are tied at max.dist (75 pairs at",
                  "exactly the 25% distance quantile; 448 within 4 ulp of geoR's umax), so binning depends on",
                  "ulp-level distance arithmetic. Complete per-delta vectors are stored for permutation 1 of the",
                  "forward direction."))
  cases$kidney_AC <- make_case("kidney_AC", getv(rk$A, shAC), getv(rk$C, shAC), spatialCoords(rk$A)[shAC, ], pixel = shAC,
    detail_vectors = FALSE,
    notes = "data(speKidney) A vs C (same rasterization), positive control r ~ +0.94; summary intermediates only.")
  use_default_rng()
  set.seed(7)
  jit <- matrix(runif(2 * length(shAB), -0.02, 0.02), ncol = 2)
  cases$kidney_AB_jitter <- make_case("kidney_AB_jitter", getv(rk$A, shAB), getv(rk$B, shAB),
                                      spatialCoords(rk$A)[shAB, ] + jit, pixel = shAB,
    extra_input = list(jitter_recipe = "pos = kidney_AB pos + matrix(runif(2 * N, -0.02, 0.02), ncol = 2) after set.seed(7)"),
    notes = paste("kidney_AB with coordinates jittered by U(-0.02, 0.02): no tied distances. As in every data set,",
                  "the pair that defines geoR's umax still lies exactly on the last bin edge, so its binning depends",
                  "on the libm hypot() (see tests/testthat/helper-platform.R). Complete per-delta vectors for",
                  "permutation 1 (forward)."))
  q <- datasets::quakes
  q <- q[!duplicated(cbind(q$lat, q$long)), ][1:300, ]
  cases$quakes_irregular <- make_case("quakes_irregular", q$depth, q$mag, cbind(q$lat, q$long),
    extra_input = list(recipe = "q <- quakes[!duplicated(cbind(quakes$lat, quakes$long)), ][1:300, ]; X = q$depth, Y = q$mag, pos = cbind(q$lat, q$long)"),
    notes = paste("base R datasets::quakes (first 300 unique positions; the roxygen example data): irregularly",
                  "spaced points, but on a 0.01-degree grid, so distances can tie; no download. Complete per-delta",
                  "vectors for permutation 1 (forward)."))
  cases$brain_Oprk1_subsample <- make_case("brain_Oprk1_subsample", bX["Oprk1", ], bY["Oprk1", ],
    spatialCoords(b$rast$MERFISH)[px, ], pixel = px, B = 5, detail_vectors = FALSE,
    notes = paste("MERFISH vs Visium brain, gene Oprk1 (lognorm), N = 2170 > 1000: exercises the sample(N, 1000)",
                  "variogram subsample, drawn before the permutations; Visium side zero-inflated. Summary",
                  "intermediates only (no per-delta vectors, no permutation fields)."))
  cases$kidney_AB_jitter_seed17 <- make_case("kidney_AB_jitter_seed17", getv(rk$A, shAB), getv(rk$B, shAB),
    spatialCoords(rk$A)[shAB, ] + jit, pixel = shAB, deltaX = seq(0.1, 0.9, 0.1), deltaY = c(0.2, 0.5, 0.8), B = 3,
    seed = 17, maxDistPrctile = 0.3, detail_vectors = FALSE,
    extra_input = list(recipe = "inputs of kidney_AB_jitter; seed = 17, maxDistPrctile = 0.3, deltaX = seq(0.1, 0.9, 0.1), deltaY = c(0.2, 0.5, 0.8)"),
    notes = paste("Argument plumbing: seed != 0 (permutations from set.seed(17), noise from set.seed(17 + i)),",
                  "maxDistPrctile != 0.25 and different delta grids for the two directions (deltaX permutes X,",
                  "deltaY permutes Y). Summary intermediates only."))
  idx1k <- seq_len(1000)
  cases$brain_Oprk1_N1000 <- make_case("brain_Oprk1_N1000", bX["Oprk1", idx1k], bY["Oprk1", idx1k],
    spatialCoords(b$rast$MERFISH)[px[idx1k], ], pixel = px[idx1k], deltaX = c(0.1, 0.3), B = 2, detail_vectors = FALSE,
    extra_input = list(recipe = "first 1000 pixels of brain_Oprk1_subsample; deltas c(0.1, 0.3); B = 2"),
    notes = paste("N = 1000 exactly: the variogram subsample is drawn only when N > 1000, so here all points are",
                  "used and the permutations start right after set.seed(seed). Summary intermediates only."))

  sim <- function(a, b, ...) {
    s <- spatialSimilarity(list(a, b), ...)
    st <- s$similarityTable
    list(args = list(...),
         table = st[, c("gene", "percentSimilarity", "percentDissimilarityX", "percentDissimilarityY",
                        "numPixelInThresh", "numPixelOutThresh", "t1", "t2")],
         similarPixelID = st$similarPixelID[[1]], dissimilarPixelIDX = st$dissimilarPixelIDX[[1]],
         dissimilarPixelIDY = st$dissimilarPixelIDY[[1]], pixelIDInThresh = st$pixelIDInThresh[[1]],
         pixelIDOutThresh = st$pixelIDOutThresh[[1]], log2ratio = s$pixelLogTransformation$log[[1]])
  }
  similarity <- list(
    recipe = "SEraster::rasterizeGeneExpression(speKidney, assay_name = 'counts', resolution = 0.2, fun = 'mean', square = FALSE); spatialSimilarity(list(rast$A, rast$B), ...)",
    kidney_AB = sim(rk$A, rk$B), kidney_AC = sim(rk$A, rk$C), kidney_AC_foldChange2 = sim(rk$A, rk$C, foldChange = 2))
  cat(sprintf("  spatialSimilarity: A-B S=%.4f, A-C S=%.4f, A-C (foldChange 2) S=%.4f\n",
              similarity$kidney_AB$table$percentSimilarity, similarity$kidney_AC$table$percentSimilarity,
              similarity$kidney_AC_foldChange2$table$percentSimilarity))
  edge_cases <- data.frame(
    case = c("constant X", "all-zero X", "single non-zero pixel", "one NA in X",
             "1 <= delta*N < 2, e.g. N = 12 with delta 0.1, or N < 200 with delta 0.01",
             "delta*N < 1, e.g. delta = 0.001 with N = 273", "duplicated coordinates", "N = 30 points"),
    legacy_behaviour = c("error caught by spatialCorrelation(): r, naive p and empirical p all NA (cor.test warns sd is zero)",
                         "same as constant X", "runs; finite p-values",
                         "r computed with the NA dropped; pX = pY = NA",
                         "locfit segfaults ('C stack overflow') and kills the R session: never run in tests",
                         "error 'procv: no points with non-zero weight' caught -> NA p-values",
                         "runs", "runs with locfit warnings 'Estimated rdf < 1.0'"),
    tested = c("test-reference-kernel.R (error paths)", "no", "no", "test-reference-kernel.R (error paths)", "never",
               "test-reference-kernel.R (error paths)", "no", "no"),
    stringsAsFactors = FALSE)
  runtime <- vapply(cases, function(cs) cs$.runtime, 0)
  cases <- lapply(cases, function(cs) { cs$.runtime <- NULL; cs })
  fixture <- list(rng = rng_semantics, cases = cases, similarity = similarity, edge_cases = edge_cases,
                  notes = paste("Legacy R behaviour, pinned for a future C++ backend. Expected p-values use the legacy",
                                "definition (b / B with strict '>', can be 0, no multiple-testing adjustment; the package",
                                "now computes (b + 1) / (B + 1)); nulls and deltaStar are the",
                                "primary references. Tolerances: integers, bin counts and deltaStar exact; RNG draws",
                                "1e-14; floating point 1e-12. The values are exact only where geoR bins each case's",
                                "pixel pairs as on the build machine (meta$platform_signature, compared per case by",
                                "the tests)."))
  fixture$meta <- stc_meta("data-raw/build_test_fixtures.R tier0",
                           inputs = list(brain_merfish_visium_rast.rds = b$meta[c("created", "script", "inputs")]),
                           sources = stc_sources(b$meta$inputs$file),
                           goldens = file.path("inst", "extdata", "brain-MERFISH-10x-visium", "brainCorrelation.RData"),
                           extra = list(rng_kind_at_startup = rng_at_startup, runtime_sec_current_R = runtime,
                                        licence_notice = "tests/testthat/fixtures/README.md"))
  fixture$meta$platform_signature <- stc_platform_signature(stc_signature_sets("kernel_fixture.rds", fixture))
  fixture <- fixture[c("meta", setdiff(names(fixture), "meta"))]
  save_fixture(fixture, "kernel_fixture.rds")
  cat(sprintf("  tier0 done in %.1f s\n", secs(t_tier)))
}

## =============================================================================================
## Tier 1: calibration
## =============================================================================================

# mixing recipe for positive controls: same covariance model, population correlation rho with f_i
# (tests/testthat/helper-fixtures.R defines the same function as stc_mix())
mix_mu <- 10
mix <- function(fi, fj, rho) rho * (fi - mix_mu) + sqrt(1 - rho^2) * (fj - mix_mu) + mix_mu

calibration_job <- function(fields, coords, i, j, rho, B, seed = 0) {
  sh <- which(!is.na(fields[i, ]) & !is.na(fields[j, ]))
  X <- unname(fields[i, sh]); Yj <- unname(fields[j, sh])
  Y <- if (rho == 0) Yj else mix(X, Yj, rho)
  use_default_rng()
  o <- suppressWarnings(spatialCorrelation(X, Y, unname(coords[sh, , drop = FALSE]), nPermutations = B,
                                           BPPARAM = BiocParallel::SerialParam(), seed = seed))
  list(row = data.frame(i = i, j = j, rho = rho, N = length(sh), r = unname(o$correlationCoef), pNaive = o$pValueNaive,
                        pX = o$pValuePermuteX, pY = o$pValuePermuteY,
                        nExtremeX = sum(abs(o$nullCorrelationsX[[1]]) > abs(o$correlationCoef)),
                        nExtremeY = sum(abs(o$nullCorrelationsY[[1]]) > abs(o$correlationCoef)), B = B,
                        deltaStarMedianX = o$deltaStarMedianX, deltaStarMedianY = o$deltaStarMedianY),
       deltaStarX = as.numeric(o$deltaStarX[[1]]), deltaStarY = as.numeric(o$deltaStarY[[1]]))
}

if ("tier1" %in% tiers) {
  cat("== Tier 1: calibration fixture ==\n")
  t_tier <- Sys.time()
  data(simRanPatternRasts)
  ids <- lapply(simRanPatternRasts, function(s) rownames(spatialCoords(s)))
  all_px <- unique(unlist(ids)); all_px <- all_px[order(as.integer(sub("pixel", "", all_px)))]
  coords <- matrix(NA_real_, length(all_px), 2, dimnames = list(all_px, c("x", "y")))
  fields <- matrix(NA_real_, length(simRanPatternRasts), length(all_px), dimnames = list(NULL, all_px))
  ncell <- matrix(NA_integer_, length(simRanPatternRasts), length(all_px), dimnames = list(NULL, all_px))
  for (i in seq_along(simRanPatternRasts)) {
    s <- simRanPatternRasts[[i]]; p <- rownames(spatialCoords(s))
    coords[p, ] <- spatialCoords(s)[p, ]
    fields[i, p] <- as.numeric(assay(s, "pixelval")[1, p])
    ncell[i, colnames(s)] <- as.integer(colData(s)$num_cell)
  }
  pr <- t(combn(length(simRanPatternRasts), 2))
  naive <- t(apply(pr, 1, function(ij) {
    sh <- which(!is.na(fields[ij[1], ]) & !is.na(fields[ij[2], ]))
    ct <- cor.test(fields[ij[1], sh], fields[ij[2], sh]); c(length(sh), ct$estimate, ct$p.value)
  }))
  pairs <- data.frame(i = pr[, 1], j = pr[, 2], n_shared = as.integer(naive[, 1]), r = naive[, 2], p_naive = naive[, 3])
  ref <- load_extdata("simRanPatternResults.RData")
  cat(sprintf("  fields %d x %d pixels (NA where absent); pairs %d; naive p<0.05: %.1f%%; BH(naive)<0.05: %.1f%%; shipped empirical p<0.05: %.2f%%\n",
              nrow(fields), ncol(fields), nrow(pairs), 100 * mean(pairs$p_naive < 0.05),
              100 * mean(p.adjust(pairs$p_naive, "BH") < 0.05), 100 * mean(ref$corspv_corrected < 0.05)))
  # pairs used by tests/testthat/test-calibration.R: 40 DISJOINT null pairs (80 distinct fields, so the 40
  # p-values are independent and a binomial bound applies) and the first 20 of them mixed to rho = 0.6
  use_default_rng()
  set.seed(20261003)
  o <- sample(nrow(fields))
  null_pairs <- data.frame(i = pmin(o[seq(1, 79, 2)], o[seq(2, 80, 2)]), j = pmax(o[seq(1, 79, 2)], o[seq(2, 80, 2)]))
  test_jobs <- rbind(data.frame(null_pairs, rho = 0), data.frame(null_pairs[1:20, ], rho = 0.6))
  B_test <- 100
  t0 <- Sys.time()
  out <- par_lapply(seq_len(nrow(test_jobs)), function(k)
    calibration_job(fields, coords, test_jobs$i[k], test_jobs$j[k], test_jobs$rho[k], B_test))
  if (any(vapply(out, inherits, NA, "try-error"))) { print(out[vapply(out, inherits, NA, "try-error")]); stop("calibration job failed") }
  res <- do.call(rbind, lapply(out, `[[`, "row"))
  dsX <- sapply(out, `[[`, "deltaStarX")
  dsY <- sapply(out, `[[`, "deltaStarY")
  cat(sprintf("  reference run of %d test jobs (B = %d) on %d workers: %.1f s wall\n", nrow(res), B_test, workers, secs(t0)))
  for (rho in unique(res$rho)) {
    d <- res[res$rho == rho, ]
    cat(sprintf("    rho=%.1f (n=%d): mean r=%.3f; naive p<0.05: %.2f; pX<0.05: %.2f; pY<0.05: %.2f; max(pX,pY)<0.05: %.2f\n",
                rho, nrow(d), mean(d$r), mean(d$pNaive < 0.05), mean(d$pX < 0.05), mean(d$pY < 0.05), mean(pmax(d$pX, d$pY) < 0.05)))
  }
  cat(sprintf("    deltaStar share at 0.1 / interior / 0.9: X %.2f / %.2f / %.2f; Y %.2f / %.2f / %.2f\n",
              mean(dsX == 0.1), mean(dsX > 0.1 & dsX < 0.9), mean(dsX == 0.9),
              mean(dsY == 0.1), mean(dsY > 0.1 & dsY < 0.9), mean(dsY == 0.9)))
  mix_check <- list(fi = c(9, 10, 11.5, 7.25), fj = c(10.2, 8, 12, 13.5), rho = 0.6)
  mix_check$value <- mix(mix_check$fi, mix_check$fj, mix_check$rho)
  fixture <- list(rng = rng_semantics, coords = coords, fields = fields, num_cell = ncell, pairs = pairs,
                  reference_shipped = list(
                    source = paste("inst/extdata/simRanPatternResults.RData: all 9900 ordered pairs, pValuePermuteX only,",
                                   "0 replaced by 0.01; generated 2026-04-24 with spatialCorrelationGeneExp_test, i.e.",
                                   "before the 2026-04-28 screening-threshold fix, so it is a statistical reference only"),
                    cors_df = ref, empirical_rate_p_lt_0.05 = mean(ref$corspv_corrected < 0.05)),
                  mix_mu = mix_mu,
                  mix_recipe = "stc_mix(f_i, f_j, rho, mu = 10) = rho * (f_i - mu) + sqrt(1 - rho^2) * (f_j - mu) + mu",
                  mix_check = mix_check,
                  test_jobs = test_jobs,
                  reference_test_jobs = list(
                    B = B_test, seed = 0, results = res, deltaStarX = dsX, deltaStarY = dsY,
                    call = paste("spatialCorrelation(X = fields[i, sh], Y = (rho == 0 ? fields[j, sh] : stc_mix(fields[i, sh], fields[j, sh], rho)),",
                                 "pos = coords[sh, ], nPermutations = 100, seed = 0, BPPARAM = SerialParam()) with sh the pixels",
                                 "present in both fields and default deltas / maxDistPrctile; RNGkind() at R defaults"),
                    notes = paste("raw p-values as computed by the package at build time, tail counts nExtremeX/Y",
                                  "(strict '>'; tests convert them with the current p-value definition), deltaStar medians and the B x 60",
                                  "matrices of deltaStar (deltaStarX, deltaStarY) of the legacy R implementation")),
                  notes = paste("Each field: GRF (exponential / Matern nu = 0.5, range 0.1 in [0,1]^2) + N(0, 0.3) noise + 10,",
                                "cells kept inside a kidney-shaped region (~1240 of 5000 cells), rasterized on a shared hex grid",
                                "(spacing 0.2). All 4950 unordered pairs are independent nulls. stc_mix(f_i, f_j, rho) gives a field",
                                "with the same covariance model whose population correlation with f_i is rho (power)."))
  fixture$meta <- stc_meta("data-raw/build_test_fixtures.R tier1",
                           goldens = file.path("inst", "extdata", "simRanPatternResults.RData"),
                           extra = list(rng_kind_at_startup = rng_at_startup, runtime_sec_reference_jobs = secs(t0),
                                        licence_notice = "tests/testthat/fixtures/README.md"))
  fixture$meta$platform_signature <- stc_platform_signature(stc_signature_sets("calibration_fixture.rds", fixture))
  fixture <- fixture[c("meta", setdiff(names(fixture), "meta"))]
  save_fixture(fixture, "calibration_fixture.rds")
  cat(sprintf("  tier1 done in %.1f s\n", secs(t_tier)))
}

## =============================================================================================
## Tier 2: realistic regression subsets with golden values from inst/extdata
## =============================================================================================
legacy_p <- function(nc, r) sum(abs(nc) > abs(r)) / length(nc)

if ("tier2" %in% tiers) {
  cat("== Tier 2: realistic fixture ==\n")
  t_tier <- Sys.time()
  pick <- function(ref, X, Y) {
    data.frame(gene = rownames(ref), r = ref$correlationCoef,
               p100X = mapply(function(nc, r) legacy_p(nc[1:100], r), ref$nullCorrelationsX, ref$correlationCoef),
               p100Y = mapply(function(nc, r) legacy_p(nc[1:100], r), ref$nullCorrelationsY, ref$correlationCoef),
               nperm = lengths(ref$nullCorrelationsX), pX_final = ref$pValuePermuteX, pY_final = ref$pValuePermuteY,
               fz_X = rowMeans(X[rownames(ref), ] == 0), fz_Y = rowMeans(Y[rownames(ref), ] == 0),
               stringsAsFactors = FALSE, row.names = NULL)
  }
  pack <- function(label, rl, assayName, ref, ref_file, genes, classes, deltas, notes, B_keep = 100) {
    sh <- intersect(rownames(spatialCoords(rl[[1]])), rownames(spatialCoords(rl[[2]])))
    X <- as.matrix(assay(rl[[1]], assayName)[genes, sh]); Y <- as.matrix(assay(rl[[2]], assayName)[genes, sh])
    gold <- lapply(genes, function(g) {
      nx <- as.numeric(ref[g, "nullCorrelationsX"][[1]]); ny <- as.numeric(ref[g, "nullCorrelationsY"][[1]])
      list(correlationCoef = ref[g, "correlationCoef"], pValueNaive = ref[g, "pValueNaive"],
           pValuePermuteX_published = ref[g, "pValuePermuteX"], pValuePermuteY_published = ref[g, "pValuePermuteY"],
           nPermutations_published = length(nx),
           nullX = nx[seq_len(min(B_keep, length(nx)))], nullY = ny[seq_len(min(B_keep, length(ny)))],
           deltaStarX = as.numeric(ref[g, "deltaStarX"][[1]])[seq_len(min(B_keep, length(nx)))],
           deltaStarY = as.numeric(ref[g, "deltaStarY"][[1]])[seq_len(min(B_keep, length(ny)))],
           pRawX_first100 = legacy_p(nx[1:100], ref[g, "correlationCoef"]),
           pRawY_first100 = legacy_p(ny[1:100], ref[g, "correlationCoef"]))
    })
    names(gold) <- genes
    check_input(paste(label, "(selected genes)"), X, Y, ref[genes, ])
    list(label = label, source = ref_file, genes = genes, X = X, Y = Y, pos = unname(spatialCoords(rl[[1]])[sh, ]),
         pixel = sh, gene_class = setNames(classes, genes),
         params = list(delta = deltas, maxDistPrctile = 0.25, seed = 0, nPermutations_published = c(100, 1000),
                       alpha = 0.05, adjustMethod = "BH",
                       published_call = "spatialCorrelationGeneExpIterPermutations(seed = 0) (BH across genes after the 1000-permutation round)",
                       golden_nulls_kept = B_keep),
         golden = gold, notes = notes)
  }
  # --- AKI kidney (X = control NL3, Y = AKI IL3), resolution 5, CPM
  aki <- need_input("aki_rast.rds", "data-raw/build_inputs_aki.R")
  kc <- load_extdata("kidneyCorrelation.RData")
  shA <- aki$shared
  XA <- as.matrix(assay(aki$rast$AKI_ctrl, "CPM")[rownames(kc), shA]); YA <- as.matrix(assay(aki$rast$AKI_aki, "CPM")[rownames(kc), shA])
  check_input("AKI (all published genes)", XA, YA, kc)
  tabA <- pick(kc, XA, YA)
  sigA <- tabA$pX_final < 0.05 & tabA$pY_final < 0.05
  posA <- with(tabA[sigA & tabA$r > 0, ], gene[order(-r)])[c(1:5, seq(20, 400, length.out = 5))]   # strongest + spread
  negA <- with(tabA[sigA & tabA$r < 0, ], gene[order(r)])[1:10]                                     # most negative
  nulA <- with(tabA[tabA$p100X > 0.5 & tabA$p100Y > 0.5, ], gene[order(abs(r))])[1:10]              # smallest |r|
  bordA <- with(tabA[!sigA & tabA$nperm == 1000, ], gene[order(pmax(pX_final, pY_final))])[1:5]     # screened, BH p ~ 0.05
  genesA <- c(posA, negA, nulA, bordA)
  clsA <- rep(c("positive", "negative", "null", "borderline"), c(length(posA), length(negA), length(nulA), length(bordA)))
  akiFix <- pack("AKI_NL3_vs_IL3_res5_CPM", aki$rast, "CPM", kc, "inst/extdata/kidneyCorrelation.RData", genesA, clsA,
                 c(0.01, 0.05, seq(0.1, 0.9, .1)),
                 paste("Visium mouse kidney, sham control (NL3, X) vs ischemic AKI (IL3, Y); AKI spots aligned to the control",
                       "with STalign (affine, region one-hots); both rotated 90 degrees; SEraster resolution 5 (array-index",
                       "units, anisotropic), fun = sum, hexagons; CPM; genes from the 1046 shared SVGs."))
  akiFix$selection <- tabA[match(genesA, tabA$gene), ]
  # --- brain MERFISH (S2R3 -> Visium) vs Visium FFPE, resolution 20, mean of libnorm then log10(x + 1)
  br <- need_input("brain_merfish_visium_rast.rds", "data-raw/build_inputs_brain.R")
  bc <- load_extdata("brain-MERFISH-10x-visium/brainCorrelation.RData")
  XB <- as.matrix(assay(br$rast$MERFISH, "lognorm")[rownames(bc), br$shared]); YB <- as.matrix(assay(br$rast$Visium, "lognorm")[rownames(bc), br$shared])
  check_input("brain (all published genes)", XB, YB, bc)
  tabB <- pick(bc, XB, YB)
  sigB <- tabB$pX_final < 0.05 & tabB$pY_final < 0.05
  posB <- with(tabB[sigB & tabB$r > 0, ], gene[order(-r)])[c(1:6, seq(20, 120, length.out = 4))]
  nulB <- with(tabB[tabB$p100X > 0.5 & tabB$p100Y > 0.5, ], gene[order(abs(r))])[1:10]
  sparseB <- with(tabB[!(tabB$gene %in% c(posB, nulB)) & tabB$fz_Y > 0.95, ], gene[order(-fz_Y)])[1:5]
  bordB <- with(tabB[!sigB & tabB$nperm == 1000 & !(tabB$gene %in% c(posB, nulB, sparseB)), ], gene[order(pmax(pX_final, pY_final))])[1:5]
  genesB <- c(posB, nulB, sparseB, bordB)
  clsB <- rep(c("positive", "null", "sparse_Y", "borderline"), c(length(posB), length(nulB), length(sparseB), length(bordB)))
  brFix <- pack("brain_MERFISH_vs_Visium_res20_lognorm", br$rast, "lognorm", bc,
                "inst/extdata/brain-MERFISH-10x-visium/brainCorrelation.RData", genesB, clsB, seq(0.1, 0.9, 0.1),
                paste("MERFISH S2R3 (Pmatch > 0.95, STalign to Visium; X) vs Visium FFPE adult mouse brain (Y); libnorm = CPM",
                      "on shared genes; SEraster resolution 20 (hires pixels), fun = mean, hexagons; lognorm = log10(x + 1);",
                      "N = 2170 > 1000, so the variogram uses a 1000-pixel subsample."))
  brFix$selection <- tabB[match(genesB, tabB$gene), ]
  brFix$engineered_negatives <- list(
    genes = posB[1:5],
    recipe = "Yflip = max(Y[g, ]) - Y[g, ] (computed in the test; not stored)",
    expected = paste("correlationCoef(X, Yflip) == -correlationCoef(X, Y) and nullCorrelationsX(X, Yflip) == -nullCorrelationsX(X, Y)",
                     "to ~1e-15 with identical deltaStarX, so pValuePermuteX is unchanged; pValuePermuteY is only",
                     "statistically equal (the permuted field changes)."))
  sel_rules <- paste(
    "AKI: positive = significant (published pX, pY < 0.05) with r > 0, the 5 strongest plus 5 spread over ranks 20..400;",
    "negative = the 10 most negative significant genes; null = raw p100X, p100Y > 0.5 with the smallest |r|;",
    "borderline = not significant but screened into the 1000-permutation round, smallest max(pX, pY).",
    "Brain: positive = 6 strongest + 4 spread over ranks 20..120; null as AKI; sparse_Y = other genes with > 95% zeros on",
    "the Visium side (most zeros first); borderline as AKI. p100X/p100Y = legacy raw p over the first 100 stored nulls.")
  cat(sprintf("  AKI: %d genes x %d px (%s)\n", nrow(akiFix$X), ncol(akiFix$X), paste(names(table(clsA)), table(clsA), sep = ":", collapse = " ")))
  cat(sprintf("  brain: %d genes x %d px (%s) + %d engineered negatives\n", nrow(brFix$X), ncol(brFix$X),
              paste(names(table(clsB)), table(clsB), sep = ":", collapse = " "), length(brFix$engineered_negatives$genes)))
  op <- options(width = 200, max.print = 1e5)
  print(akiFix$selection, digits = 3, row.names = FALSE)
  print(brFix$selection, digits = 3, row.names = FALSE)
  options(op)
  fixture <- list(rng = rng_semantics, pairs = list(aki = akiFix, brain = brFix), selection_rules = sel_rules,
                  notes = paste("Golden values are copied from the published results in inst/extdata, computed by the",
                                "authors with the current code (20-22 workers; the machine is not recorded). They",
                                "reproduce bit for bit on macOS arm64 but not on Linux x86-64 or arm64 (glibc 2.39),",
                                "where geoR bins the AKI and brain pixel pairs differently, so they were most likely",
                                "computed on macOS arm64. Only the first 100 nulls/deltaStar are kept: by the prefix",
                                "property they equal a 100-permutation run. pValuePermuteX/Y_published are the final",
                                "BH-adjusted values of the iterative protocol, kept for reference only; tests use raw",
                                "nulls, deltaStar, r and the naive p-value."))
  fixture$meta <- stc_meta("data-raw/build_test_fixtures.R tier2",
                           inputs = list(aki_rast.rds = aki$meta[c("created", "script", "inputs")],
                                         brain_merfish_visium_rast.rds = br$meta[c("created", "script", "inputs")]),
                           sources = stc_sources(c(aki$meta$inputs$file, br$meta$inputs$file)),
                           goldens = file.path("inst", "extdata", c("kidneyCorrelation.RData",
                                                                    "brain-MERFISH-10x-visium/brainCorrelation.RData")),
                           extra = list(rng_kind_at_startup = rng_at_startup, licence_notice = "tests/testthat/fixtures/README.md"))
  fixture$meta$platform_signature <- stc_platform_signature(stc_signature_sets("realistic_fixture.rds", fixture))
  fixture <- fixture[c("meta", setdiff(names(fixture), "meta"))]
  save_fixture(fixture, "realistic_fixture.rds")
  cat(sprintf("  tier2 done in %.1f s\n", secs(t_tier)))
}

## =============================================================================================
## Verify: re-run the current code against the fixtures (same tolerances as the tests)
## =============================================================================================
if ("verify" %in% tiers) {
  Bv <- as.integer(Sys.getenv("STCOMPARE_VERIFY_B", "100"))
  t_tier <- Sys.time()
  ok <- TRUE
  kf <- readRDS(file.path(fixdir, "kernel_fixture.rds"))
  rf <- readRDS(file.path(fixdir, "realistic_fixture.rds"))
  # does this machine bin each coordinate set as the build machine did? (tests/testthat/helper-platform.R)
  platform_note <- function(fx, name) {
    ref <- fx$meta$platform_signature
    cur <- stc_platform_signature(stc_signature_sets(name, fx))
    bad <- Filter(length, lapply(setNames(names(ref$sets), names(ref$sets)), function(s) stc_set_diff(ref, cur, s)))
    if (length(bad)) cat(sprintf(paste0("  NOTE: geoR bins %s differently here than on the build machine (%s): exact ",
                                        "agreement is not expected for %s\n"),
                                 paste(names(bad), collapse = ", "), ref$info$platform, if (length(bad) == 1) "it" else "them"))
    invisible(names(bad))
  }
  B_real <- min(Bv, 100L)
  B_neg <- min(Bv, 20L)
  if (Bv > 100L) message("STCOMPARE_VERIFY_B = ", Bv, " exceeds the 100 stored nulls per gene; using B = 100")
  cat(sprintf("== Verify: current R code vs fixtures (kernel cases at their stored B; realistic genes at B = %d, engineered negatives at B = %d; %d workers) ==\n",
              B_real, B_neg, workers))
  platform_note(kf, "kernel_fixture.rds")
  for (cs in kf$cases) {
    use_default_rng()
    o <- suppressWarnings(spatialCorrelation(cs$input$X, cs$input$Y, cs$input$pos, nPermutations = cs$input$nPermutations,
                                             deltaX = cs$input$deltaX, deltaY = cs$input$deltaY,
                                             maxDistPrctile = cs$input$maxDistPrctile, seed = cs$input$seed))
    dn <- max(abs(as.numeric(o$nullCorrelationsX[[1]]) - cs$expected$forward$nullCor),
              abs(as.numeric(o$nullCorrelationsY[[1]]) - cs$expected$reverse$nullCor))
    ds <- identical(as.numeric(o$deltaStarX[[1]]), cs$expected$forward$deltaStar) &&
      identical(as.numeric(o$deltaStarY[[1]]), cs$expected$reverse$deltaStar)
    good <- dn <= 1e-12 && ds
    cat(sprintf("  kernel %-24s B=%2d  max|dnull| = %.2g, deltaStar identical: %s  -> %s\n", cs$name,
                cs$input$nPermutations, dn, ds, if (good) "ok" else "MISMATCH"))
    ok <- ok && good
  }
  platform_note(rf, "realistic_fixture.rds")
  jobs <- do.call(rbind, lapply(names(rf$pairs), function(nm) data.frame(pair = nm, gene = rf$pairs[[nm]]$genes, flip = FALSE)))
  jobs <- rbind(jobs, data.frame(pair = "brain", gene = rf$pairs$brain$engineered_negatives$genes, flip = TRUE))
  run_job <- function(k) {
    P <- rf$pairs[[jobs$pair[k]]]; g <- jobs$gene[k]; G <- P$golden[[g]]
    Y <- if (jobs$flip[k]) max(P$Y[g, ]) - P$Y[g, ] else P$Y[g, ]
    B <- if (jobs$flip[k]) B_neg else B_real
    t0 <- Sys.time()
    use_default_rng()
    o <- suppressWarnings(spatialCorrelation(P$X[g, ], Y, P$pos, nPermutations = B, deltaX = P$params$delta,
                                             deltaY = P$params$delta, maxDistPrctile = P$params$maxDistPrctile,
                                             BPPARAM = BiocParallel::SerialParam(), seed = P$params$seed))
    sgn <- if (jobs$flip[k]) -1 else 1
    nx <- as.numeric(o$nullCorrelationsX[[1]]); ny <- as.numeric(o$nullCorrelationsY[[1]])
    data.frame(pair = jobs$pair[k], gene = g, class = if (jobs$flip[k]) "engineered_negative" else P$gene_class[[g]], B = B,
               sec = secs(t0), dr = abs(o$correlationCoef - sgn * G$correlationCoef),
               rel_dp_naive = if (jobs$flip[k]) NA else abs(o$pValueNaive / G$pValueNaive - 1),
               max_dnullX = max(abs(nx - sgn * G$nullX[1:B])),
               max_dnullY = if (jobs$flip[k]) NA else max(abs(ny - G$nullY[1:B])),
               deltaStarX_identical = identical(as.numeric(o$deltaStarX[[1]]), G$deltaStarX[1:B]),
               deltaStarY_identical = if (jobs$flip[k]) NA else identical(as.numeric(o$deltaStarY[[1]]), G$deltaStarY[1:B]))
  }
  t0 <- Sys.time()
  res <- par_lapply(seq_len(nrow(jobs)), run_job)
  err <- vapply(res, inherits, NA, "try-error")
  if (any(err)) { print(res[err]); stop("verify jobs failed") }
  res <- do.call(rbind, res)
  op <- options(width = 250, max.print = 1e5)
  print(res, digits = 3, row.names = FALSE)
  options(op)
  pass <- with(res, dr <= 1e-12 & (is.na(rel_dp_naive) | rel_dp_naive <= 1e-10) & max_dnullX <= 1e-12 &
                 (is.na(max_dnullY) | max_dnullY <= 1e-12) & deltaStarX_identical & (is.na(deltaStarY_identical) | deltaStarY_identical))
  cat(sprintf("  realistic: %d/%d jobs reproduce the golden values (max|dnull| = %.2g; %.1f CPU-s in %.1f s wall)\n",
              sum(pass), length(pass), max(c(res$max_dnullX, res$max_dnullY), na.rm = TRUE), sum(res$sec), secs(t0)))
  if (!all(pass)) { print(res[!pass, ]); ok <- FALSE }
  cat(sprintf("  verify done in %.1f s\n", secs(t_tier)))
  if (!ok) stop("verification failed")
}

if (compare_failed) stop("compare: the rebuilt content differs from the committed fixtures (see above)")
