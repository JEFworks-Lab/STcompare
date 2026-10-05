#!/usr/bin/env Rscript
# bench/calibrate-surrogates.R
#
# Calibration, power and cost of the two surrogate modes of compareSpatial(): surrogate = "gaussian" (the
# Viladomat surrogates as the smoothing and the noise leave them, the surrogates of the legacy functions) and
# surrogate = "remap" (rank-remapped onto the gene's own values), per dev/surrogate-remap-spec.md, section 4.
# Every run uses the default settings of compareSpatial() (adaptive p-values with exceedances = 10 and
# nPermutations = 10000, the extended delta grid) unless stated, except that the surrogate mode and minDetected
# are always passed explicitly: the study was designed when "gaussian" with the sqrt(N) detection filter was
# the default of compareSpatial(); on its recommendation "remap" with no filter is the default now (minDetected
# = NULL means 0 pixels with "remap" and sqrt(N) pixels with "gaussian"), so compare() passes the sqrt(N) filter
# (sqrt_share()) where the study relied on the old default, and minDetected = 0 where stated, so that both modes
# test the same genes and the run reproduces.
#
#   dense   Null, dense fields: all 4950 pairs of the 100 independent simulated fields of data(simRanPatternRasts)
#           (as stored in tests/testthat/fixtures/calibration_fixture.rds). P(p <= alpha) per mode.
#   sparse  Null, sparse genes, with minDetected = 0:
#           (a) the setting of the acceptance review: 1800 independent genes detected in 1, 2, 3, 5, 10 or 20 of
#               the 311 AKI pixels (values round(rexp(k) * 20) + 1 at k random pixels, no spatial structure);
#               P(p <= alpha) per number of detected pixels and the BH false discoveries;
#           (b) zero-inflated spatial fields: Poisson counts of a Gaussian random field (exponential covariance)
#               on the AKI (311 pixels) and brain (2170 pixels) coordinates, with the mean chosen so that a share
#               of 1%, 3%, 10% or 30% of the pixels is detected, plus the dense field itself (100%);
#               P(p <= alpha) per detection fraction.
#   power   Mixed fields with population correlation rho in {0.2, 0.4, 0.6} (the stc_mix recipe of the calibration
#           fixture) on random pairs of the simulated fields, at the same settings and seeds for both modes.
#   real    The published AKI (1046 genes, CPM) and brain (325 genes, lognorm) inputs from the cache: significant
#           genes (padj < 0.05) per mode, their overlap, the overlap with the published results (both BH-adjusted
#           direction p-values below 0.05), and the wall time of each mode (the cost). Then the whole AKI raster
#           (all genes, minDetected = 0, nPermutations = 1000): tested and significant genes by number of
#           detected pixels, per mode.
#
# Usage, from the repository root, with the package installed (R CMD INSTALL; devtools::load_all() compiles the
# engine without optimisation, about 10 times slower):
#   Rscript bench/calibrate-surrogates.R                      # 16 threads, every part
#   Rscript bench/calibrate-surrogates.R --threads=8 --parts=dense,sparse
#   Rscript bench/calibrate-surrogates.R --quick              # a small run of every part, to check the script
# Options: --threads=N (default 16); --parts=a,b,... (dense, sparse, power, real; default all); --quick (small
# sizes); --out=DIR (per-gene results as RDS; default <cache>/bench/calibrate-surrogates, never the repository;
# <cache> is $STCOMPARE_DATA_CACHE or tools::R_user_dir("STcompare", "cache")); --results=FILE (the markdown
# tables; default bench/calibration-results.md, where only the part between the markers
# "<!-- calibrate-surrogates:begin -->" and "<!-- calibrate-surrogates:end -->" is replaced, so that the
# discussion written around them is kept; --quick writes no results file unless --results is given).
# The real-data part needs <cache>/data-raw/inputs/{aki_rast,brain_merfish_visium_rast}.rds (built by
# data-raw/build_inputs_aki.R and data-raw/build_inputs_brain.R) and the published results in bench/published;
# it is skipped when the cache is missing.

options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[length(hit)]) else default
}
flag <- function(name) any(args == paste0("--", name))
threads <- as.integer(opt("threads", 16L))
if (is.na(threads) || threads < 1L) stop("--threads must be a positive integer")
all_parts <- c("dense", "sparse", "power", "real")
parts <- strsplit(opt("parts", paste(all_parts, collapse = ",")), ",", fixed = TRUE)[[1]]
if (!all(parts %in% all_parts)) stop("--parts must be a subset of: ", paste(all_parts, collapse = ", "))
quick <- flag("quick")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(script)) dirname(dirname(normalizePath(script))) else "."
cache <- Sys.getenv("STCOMPARE_DATA_CACHE", tools::R_user_dir("STcompare", which = "cache"))
inputs <- file.path(cache, "data-raw", "inputs")
out_dir <- opt("out", file.path(cache, "bench", "calibrate-surrogates"))
results_md <- opt("results", if (quick) NA_character_ else file.path(root, "bench", "calibration-results.md"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (startsWith(paste0(normalizePath(out_dir), "/"), paste0(normalizePath(root), "/"))) {
  stop("--out must be outside the repository: ", out_dir)
}

suppressPackageStartupMessages({
  library(STcompare)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

modes <- c("gaussian", "remap")
alphas <- c(0.1, 0.05, 0.01, 0.001)
# sizes (--quick: a smoke run)
n_dense <- if (quick) 60L else Inf           # pairs of simulated fields (Inf: all 4950)
sparse_per_k <- if (quick) 10L else 300L     # genes per number of detected pixels (the review's setting)
sparse_per_bin <- if (quick) c(aki = 10L, brain = 6L) else c(aki = 400L, brain = 200L)
n_power <- if (quick) 20L else 300L          # pairs per rho
n_perm_akiall <- 1000L

# ---- helpers ------------------------------------------------------------------------------------

# a SpatialExperiment from a genes x pixels matrix and pixels x 2 coordinates (row names = pixels)
spe <- function(M, pos) SpatialExperiment(assays = list(v = M), spatialCoords = pos)
load_avg <- function() {
  up <- tryCatch(system("uptime", intern = TRUE), error = function(e) "")
  m <- regmatches(up, regexpr("load average[s]?: .*$", up))
  if (length(m)) sub("load average[s]?: ", "", m) else NA_character_
}
machine <- function() {
  cpu <- if (Sys.info()[["sysname"]] == "Darwin") {
    tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) "")
  } else if (file.exists("/proc/cpuinfo")) {
    sub("^.*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])
  } else {
    ""
  }
  sprintf("%s (%s, %d logical cores), %s", cpu, Sys.info()[["machine"]], parallel::detectCores(), R.version.string)
}
wall <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- expr
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}
fmt_s <- function(s) if (s < 120) sprintf("%.1f s", s) else sprintf("%.1f min", s / 60)
# "rate (se)" of P(p <= alpha) with the binomial standard error
rate <- function(p, alpha) {
  p <- p[!is.na(p)]
  r <- mean(p <= alpha)
  sprintf("%.4f (%.4f)", r, sqrt(r * (1 - r) / length(p)))
}
rates <- function(p) vapply(alphas, function(a) rate(p, a), "")
rate_names <- sprintf("P(p <= %g)", alphas)
md_table <- function(df) {
  fmt <- function(v) if (is.numeric(v)) format(v, scientific = FALSE, trim = TRUE) else as.character(v)
  cells <- vapply(df, fmt, character(nrow(df)))
  if (nrow(df) == 1L) cells <- matrix(cells, nrow = 1L)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(cells, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
stc_mix <- function(fi, fj, rho, mu = 10) rho * (fi - mu) + sqrt(1 - rho^2) * (fj - mu) + mu
# a Gaussian random field with exponential covariance exp(-d / range) at coordinates pos: G x N
grf <- function(pos, G, range, seed) {
  C <- exp(-as.matrix(stats::dist(pos)) / range)
  diag(C) <- diag(C) + 1e-8
  R <- chol(C)
  set.seed(seed)
  t(crossprod(R, matrix(stats::rnorm(nrow(pos) * G), nrow(pos), G)))
}
# Poisson counts of a field f with the mean exp(mu + f), mu chosen so that a share `detected` of the pixels is
# expected to be nonzero
poisson_field <- function(f, detected) {
  mu <- stats::uniroot(function(m) mean(1 - exp(-exp(m + f))) - detected, c(-25, 10))$root
  stats::rpois(length(f), exp(mu + f))
}
# The sqrt(N) detection filter of the N shared pixels (ceiling(sqrt(N)) pixels, as compareSpatial() resolves it): the
# default of both modes when this study was designed, now the default of "gaussian" only (see the header).
sqrt_share <- function(x, y) {
  N <- length(intersect(colnames(x), colnames(y)))
  ceiling(sqrt(N)) / N
}
compare <- function(x, y, mode, minDetected = sqrt_share(x, y), ...) {
  compareSpatial(x, y, tests = "correlation", surrogate = mode, minDetected = minDetected, nThreads = threads,
                 progress = FALSE, verbose = FALSE, ...)
}

cat(sprintf("calibrate-surrogates: %d threads, parts %s%s\n  %s\n  load average at start: %s\n", threads,
            paste(parts, collapse = ", "), if (quick) " (quick)" else "", machine(), load_avg()))
md <- character(0)
section <- function(...) md <<- c(md, "", ...)
timing <- data.frame(part = character(0), wall = character(0), load_after = character(0))
note_time <- function(part, seconds) {
  timing <<- rbind(timing, data.frame(part = part, wall = fmt_s(seconds), load_after = load_avg()))
  message(sprintf("%s: %s (load average now %s)", part, fmt_s(seconds), load_avg()))
}

rf <- readRDS(file.path(root, "tests", "testthat", "fixtures", "realistic_fixture.rds"))
cf <- readRDS(file.path(root, "tests", "testthat", "fixtures", "calibration_fixture.rds"))
coords_of <- function(P) matrix(P$pos, ncol = 2L, dimnames = list(P$pixel, c("x", "y")))

# ---- simulated field pairs (dense null and power) --------------------------------------------------

fields <- cf$fields  # 100 x 305, NA where a field has no pixel
all_pairs <- t(utils::combn(nrow(fields), 2L))
# one pair, both modes, one gene per call (the pairs have different shared pixels), with the pair's own seed
run_pair <- function(i, j, rho, seed) {
  sh <- which(!is.na(fields[i, ]) & !is.na(fields[j, ]))
  X <- unname(fields[i, sh])
  Y <- unname(fields[j, sh])
  if (rho != 0) Y <- stc_mix(X, Y, rho, cf$mix_mu)
  px <- colnames(fields)[sh]
  co <- cf$coords[sh, , drop = FALSE]
  x <- spe(matrix(X, 1L, dimnames = list("g", px)), co)
  y <- spe(matrix(Y, 1L, dimnames = list("g", px)), co)
  do.call(rbind, lapply(modes, function(m) {
    r <- compareSpatial(x, y, tests = "correlation", surrogate = m, minDetected = sqrt_share(x, y), nThreads = 1L,
                        progress = FALSE, verbose = FALSE, seed = seed)
    data.frame(i = i, j = j, rho = rho, mode = m, N = r$nPixels, r = r$r, p = r$p, L = r$nPermutations, stop = r$stop,
               status = r$status)
  }))
}
run_pairs <- function(pairs, rho) {
  res <- parallel::mclapply(seq_len(nrow(pairs)), function(k) {
    tryCatch(run_pair(pairs[k, 1L], pairs[k, 2L], rho, seed = k), error = function(e) {
      data.frame(i = pairs[k, 1L], j = pairs[k, 2L], rho = rho, mode = modes, N = NA_integer_, r = NA_real_, p = NA_real_,
                 L = NA_integer_, stop = "error", status = conditionMessage(e))
    })
  }, mc.cores = threads)
  do.call(rbind, res)
}

if ("dense" %in% parts) {
  pairs <- if (is.finite(n_dense)) all_pairs[seq_len(n_dense), , drop = FALSE] else all_pairs
  w <- wall(run_pairs(pairs, 0))
  d <- w$value
  saveRDS(d, file.path(out_dir, "dense_null.rds"))
  note_time("dense", w$seconds)
  tab <- do.call(rbind, lapply(modes, function(m) {
    s <- d[d$mode == m, ]
    data.frame(mode = m, pairs = nrow(s), failed = sum(s$status == "failed"), `median L` = stats::median(s$L, na.rm = TRUE),
               `BH < 0.05` = sum(stats::p.adjust(s$p, "BH") < 0.05, na.rm = TRUE),
               t(stats::setNames(rates(s$p), rate_names)), check.names = FALSE)
  }))
  section("## Null, dense fields", "",
          sprintf(paste("All %d pairs of the 100 independent simulated fields of `data(simRanPatternRasts)` (one gene per pair,",
                        "%d to %d shared pixels), default settings, both modes on the same pairs and seeds. P(p <= alpha)",
                        "with its binomial standard error; the pairs share fields, so the standard errors understate the",
                        "uncertainty somewhat. Wall time %s on %d processes."),
                  nrow(pairs), min(d$N, na.rm = TRUE), max(d$N, na.rm = TRUE), fmt_s(w$seconds), threads), "",
          md_table(tab))
  print(tab, row.names = FALSE)
}

# ---- sparse null --------------------------------------------------------------------------------

if ("sparse" %in% parts) {
  aki <- rf$pairs$aki
  pos_aki <- coords_of(aki)
  N <- nrow(pos_aki)
  # (a) the review's setting (scratch scripts a7_sparse_null.R and a7b_sparse_default.R of the acceptance review)
  ks <- c(1L, 2L, 3L, 5L, 10L, 20L)
  set.seed(123)
  mkgene <- function(k) {
    v <- numeric(N)
    v[sample(N, k)] <- round(stats::rexp(k) * 20) + 1
    v
  }
  kk <- rep(ks, each = sparse_per_k)
  X <- t(sapply(kk, mkgene))
  Y <- t(sapply(kk, mkgene))
  genes <- sprintf("s%04d", seq_along(kk))
  dimnames(X) <- dimnames(Y) <- list(genes, rownames(pos_aki))
  res_a <- lapply(modes, function(m) wall(compare(spe(X, pos_aki), spe(Y, pos_aki), m, minDetected = 0, seed = 4)))
  names(res_a) <- modes
  saveRDS(lapply(res_a, function(w) as.data.frame(w$value)), file.path(out_dir, "sparse_review.rds"))
  for (m in modes) note_time(paste0("sparse (a) ", m), res_a[[m]]$seconds)
  tab_a <- do.call(rbind, lapply(ks, function(k) {
    do.call(rbind, lapply(modes, function(m) {
      d <- as.data.frame(res_a[[m]]$value)[kk == k, ]
      data.frame(`detected pixels` = k, mode = m, genes = nrow(d), failed = sum(d$status == "failed"),
                 `median L` = stats::median(d$nPermutations, na.rm = TRUE), t(stats::setNames(rates(d$p), rate_names)),
                 check.names = FALSE)
    }))
  }))
  disc <- vapply(modes, function(m) {
    d <- as.data.frame(res_a[[m]]$value)
    sprintf("%s: %d BH false discoveries (padj < 0.05) among %d genes, %d p-values at the floor 1 / 10001, %s permutations in %s",
            m, sum(d$padj < 0.05, na.rm = TRUE), nrow(d), sum(d$p == 1 / 10001, na.rm = TRUE),
            format(sum(d$nPermutations, na.rm = TRUE), big.mark = ","), fmt_s(res_a[[m]]$seconds))
  }, "")
  section("## Null, sparse genes without spatial structure (the acceptance review's setting)", "",
          sprintf(paste("%d independent genes on the AKI coordinates (%d pixels): for each number k of detected pixels, %d",
                        "genes in each sample with values round(rexp(k) * 20) + 1 at k random pixels (set.seed(123), as",
                        "in the review). Default settings, minDetected = 0, seed 4, %d threads. The sqrt(N) filter (the",
                        "default minDetected of gaussian surrogates) would skip every gene detected in fewer than %d pixels."),
                  length(genes), N, sparse_per_k, threads, ceiling(sqrt(N))), "",
          md_table(tab_a), "", paste0("- ", disc))
  print(tab_a, row.names = FALSE)
  cat(disc, sep = "\n")

  # (b) zero-inflated spatial fields by detection fraction, on the AKI and brain coordinates
  fractions <- c(0.01, 0.03, 0.1, 0.3, 1)
  tab_b <- NULL
  res_b <- list()
  for (set in c("aki", "brain")) {
    P <- rf$pairs[[set]]
    pos <- coords_of(P)
    N <- nrow(pos)
    per <- sparse_per_bin[[set]]
    range <- 0.1 * max(apply(pos, 2L, function(v) diff(range(v))))
    G <- per * length(fractions)
    FX <- grf(pos, G, range, seed = 100L + nchar(set))
    FY <- grf(pos, G, range, seed = 200L + nchar(set))
    frac <- rep(fractions, each = per)
    set.seed(7)
    X <- t(vapply(seq_len(G), function(g) if (frac[g] < 1) poisson_field(FX[g, ], frac[g]) else FX[g, ], numeric(N)))
    Y <- t(vapply(seq_len(G), function(g) if (frac[g] < 1) poisson_field(FY[g, ], frac[g]) else FY[g, ], numeric(N)))
    genes <- sprintf("z%04d", seq_len(G))
    dimnames(X) <- dimnames(Y) <- list(genes, rownames(pos))
    detected <- pmin(rowSums(X > 0), rowSums(Y > 0))
    res <- lapply(modes, function(m) wall(compare(spe(X, pos), spe(Y, pos), m, minDetected = 0, seed = 5)))
    names(res) <- modes
    res_b[[set]] <- list(frac = frac, detected = detected, results = lapply(res, function(w) as.data.frame(w$value)))
    for (m in modes) note_time(sprintf("sparse (b) %s %s", set, m), res[[m]]$seconds)
    tab_b <- rbind(tab_b, do.call(rbind, lapply(fractions, function(f) {
      do.call(rbind, lapply(modes, function(m) {
        d <- as.data.frame(res[[m]]$value)[frac == f, ]
        data.frame(coordinates = sprintf("%s (N = %d)", set, N), `detection fraction` = f,
                   `detected pixels (min of x, y), median` = stats::median(detected[frac == f]), mode = m, genes = nrow(d),
                   `skipped (constant)` = sum(d$status == "skipped"), failed = sum(d$status == "failed"),
                   `median L` = stats::median(d$nPermutations, na.rm = TRUE),
                   `BH < 0.05` = sum(d$padj < 0.05, na.rm = TRUE), t(stats::setNames(rates(d$p), rate_names)),
                   check.names = FALSE)
      }))
    })))
  }
  saveRDS(res_b, file.path(out_dir, "sparse_fields.rds"))
  section("## Null, zero-inflated spatial fields by detection fraction", "",
          sprintf(paste("Independent pairs of Gaussian random fields (exponential covariance, range 10%% of the coordinate",
                        "span) turned into Poisson counts with the mean chosen so that the given share of pixels is",
                        "detected (nonzero); 100%% is the Gaussian field itself. %d genes per fraction on the AKI",
                        "coordinates and %d on the brain coordinates, both modes on the same genes, default settings,",
                        "minDetected = 0, seed 5, %d threads. padj is adjusted within each run (all fractions together)."),
                  sparse_per_bin[["aki"]], sparse_per_bin[["brain"]], threads), "",
          md_table(tab_b))
  print(tab_b, row.names = FALSE)
}

# ---- power --------------------------------------------------------------------------------------

if ("power" %in% parts) {
  set.seed(2026)
  rhos <- c(0.2, 0.4, 0.6)
  pw <- NULL
  res_p <- list()
  for (rho in rhos) {
    pairs <- all_pairs[sample(nrow(all_pairs), n_power), , drop = FALSE]
    w <- wall(run_pairs(pairs, rho))
    res_p[[as.character(rho)]] <- w$value
    note_time(sprintf("power rho = %g", rho), w$seconds)
    pw <- rbind(pw, do.call(rbind, lapply(modes, function(m) {
      s <- w$value[w$value$mode == m, ]
      data.frame(rho = rho, mode = m, pairs = nrow(s), `mean r` = round(mean(s$r, na.rm = TRUE), 3), failed = sum(s$status == "failed"),
                 `median L` = stats::median(s$L, na.rm = TRUE), `P(p <= 0.05)` = rate(s$p, 0.05), `P(p <= 0.01)` = rate(s$p, 0.01),
                 `P(p <= 0.001)` = rate(s$p, 0.001), check.names = FALSE)
    })))
  }
  saveRDS(res_p, file.path(out_dir, "power.rds"))
  section("## Power, mixed simulated fields", "",
          sprintf(paste("%d random pairs (i, j) of the simulated fields per rho, Y = stc_mix(f_i, f_j, rho) (the recipe of the",
                        "calibration fixture: a field with the same covariance model and population correlation rho with",
                        "X = f_i), default settings, the same pairs and seeds for both modes."), n_power), "",
          md_table(pw))
  print(pw, row.names = FALSE)
}

# ---- the published inputs and the whole AKI raster --------------------------------------------------

if ("real" %in% parts) {
  f_aki <- file.path(inputs, "aki_rast.rds")
  f_brain <- file.path(inputs, "brain_merfish_visium_rast.rds")
  if (!file.exists(f_aki) || !file.exists(f_brain)) {
    message("real: cached inputs not found in ", inputs, "; skipped")
  } else {
    load_ref <- function(f) {
      e <- new.env()
      nm <- load(file.path(root, "bench", "published", f), envir = e)
      get(nm[1], envir = e)
    }
    aki <- readRDS(f_aki)
    brain <- readRDS(f_brain)
    if (quick) {  # a few genes only
      aki$genes <- utils::head(aki$genes, 30L)
      brain$rast <- lapply(brain$rast, function(s) s[seq_len(20L), ])
    }
    sets <- list(
      aki = list(x = aki$rast$AKI_ctrl[aki$genes, ], y = aki$rast$AKI_aki[aki$genes, ], assay = "CPM",
                 ref = load_ref("kidneyCorrelation.RData"), title = "AKI kidney, control vs AKI (1046 published genes, 311 pixels)"),
      brain = list(x = brain$rast$MERFISH, y = brain$rast$Visium, assay = "lognorm",
                   ref = load_ref(file.path("brain-MERFISH-10x-visium", "brainCorrelation.RData")),
                   title = "Brain, MERFISH vs Visium (325 published genes, 2170 pixels)"))
    real <- NULL
    res_r <- list()
    for (id in names(sets)) {
      s <- sets[[id]]
      res <- lapply(modes, function(m) wall(compare(s$x, s$y, m, assay = s$assay, seed = 0)))
      names(res) <- modes
      for (m in modes) note_time(sprintf("real %s %s", id, m), res[[m]]$seconds)
      d <- lapply(res, function(w) as.data.frame(w$value))
      res_r[[id]] <- d
      pub <- rownames(s$ref)[s$ref$pValuePermuteX < 0.05 & s$ref$pValuePermuteY < 0.05]
      sig <- lapply(d, function(t) t$gene[!is.na(t$padj) & t$padj < 0.05])
      both <- intersect(sig$gaussian, sig$remap)
      real <- rbind(real, do.call(rbind, lapply(modes, function(m) {
        t <- d[[m]]
        data.frame(dataset = id, mode = m, genes = nrow(t), tested = sum(t$status == "ok"), skipped = sum(t$status == "skipped"),
                   failed = sum(t$status == "failed"), `padj < 0.05` = length(sig[[m]]),
                   `in both modes` = length(both), `also published` = length(intersect(sig[[m]], pub)),
                   `published only` = length(setdiff(pub, sig[[m]])), `published total` = length(pub),
                   `median L` = stats::median(t$nPermutations, na.rm = TRUE),
                   permutations = format(sum(t$nPermutations, na.rm = TRUE), big.mark = ","),
                   wall = fmt_s(res[[m]]$seconds),
                   `us per task-permutation x threads` = signif(1e6 * res[[m]]$seconds * threads / (2 * sum(t$nPermutations, na.rm = TRUE)), 3),
                   check.names = FALSE)
      })))
    }
    saveRDS(res_r, file.path(out_dir, "real.rds"))
    section("## The published inputs: significant genes and cost", "",
            sprintf(paste("compareSpatial() with the default settings on the inputs of the published analyses (seed 0, %d",
                          "threads). Significant: padj < 0.05. Published: both BH-adjusted direction p-values below 0.05 in",
                          "bench/published (100 then 1000 permutations). The last column is the wall time per permutation",
                          "and direction, times the threads, over the permutations kept."), threads), "",
            md_table(real))
    print(real, row.names = FALSE)

    # the whole AKI raster with minDetected = 0: how many sparse genes each mode calls significant
    x_all <- aki$rast$AKI_ctrl
    y_all <- aki$rast$AKI_aki
    if (quick) {
      set.seed(1)
      keep <- sort(sample(nrow(x_all), 400L))
      x_all <- x_all[keep, ]
      y_all <- y_all[keep, ]
    }
    sh <- aki$shared
    cpmX <- assay(x_all, "CPM")[, sh]
    cpmY <- assay(y_all, "CPM")[, sh]
    nz <- pmin(Matrix::rowSums(cpmX > 0), Matrix::rowSums(cpmY > 0))
    res <- lapply(modes, function(m) wall(compare(x_all, y_all, m, assay = "CPM", minDetected = 0, nPermutations = n_perm_akiall, seed = 0)))
    names(res) <- modes
    for (m in modes) note_time(sprintf("real AKI all genes %s", m), res[[m]]$seconds)
    d_all <- lapply(res, function(w) as.data.frame(w$value))
    saveRDS(d_all, file.path(out_dir, "aki_all.rds"))
    bins <- cut(nz[d_all$gaussian$gene], c(-1, 0, 2, 5, 17, 50, 100, 400), labels = c("0", "1-2", "3-5", "6-17", "18-50", "51-100", "> 100"))
    tab_all <- do.call(rbind, lapply(levels(bins), function(b) {
      do.call(rbind, lapply(modes, function(m) {
        t <- d_all[[m]][bins == b, ]
        sig <- !is.na(t$padj) & t$padj < 0.05
        data.frame(`detected pixels (min of x, y)` = b, mode = m, genes = nrow(t), tested = sum(t$status == "ok"),
                   `padj < 0.05` = sum(sig), `of them published genes` = sum(sig & t$gene %in% aki$genes),
                   `P(p <= 0.001)` = if (any(t$status == "ok")) rate(t$p, 0.001) else "",
                   `median L` = if (any(t$status == "ok")) stats::median(t$nPermutations, na.rm = TRUE) else NA,
                   check.names = FALSE)
      }))
    }))
    section("## The whole AKI raster with minDetected = 0", "",
            sprintf(paste("All %d genes of the AKI raster (CPM, 311 shared pixels), minDetected = 0, nPermutations = %d,",
                          "exceedances = 10, seed 0, %d threads: wall %s (gaussian) and %s (remap). Genes by the number of",
                          "pixels where both samples detect them; the sqrt(N) filter (the default minDetected of gaussian",
                          "surrogates) tests only the genes detected in at least 18 pixels of each sample. The 1046 published",
                          "genes are spatially variable genes."),
                    nrow(d_all$gaussian), n_perm_akiall, threads, fmt_s(res$gaussian$seconds), fmt_s(res$remap$seconds)), "",
            md_table(tab_all))
    print(tab_all, row.names = FALSE)
  }
}

# ---- results file -------------------------------------------------------------------------------

header <- c(sprintf("Run of `bench/calibrate-surrogates.R` on %s with %d threads%s.", format(Sys.time(), "%Y-%m-%d %H:%M"),
                    threads, if (quick) " (--quick)" else ""),
            sprintf("Machine: %s. Package: STcompare %s.", machine(), utils::packageVersion("STcompare")), "",
            "Wall time of each part and the load averages (1, 5, 15 min) after it:", "", md_table(timing))
md <- c("<!-- calibrate-surrogates:begin -->", header, md, "", "<!-- calibrate-surrogates:end -->")
if (!is.na(results_md)) {
  old <- if (file.exists(results_md)) readLines(results_md) else character(0)
  b <- grep("<!-- calibrate-surrogates:begin -->", old, fixed = TRUE)
  e <- grep("<!-- calibrate-surrogates:end -->", old, fixed = TRUE)
  new <- if (length(b) == 1L && length(e) == 1L && e > b) {
    c(old[seq_len(b - 1L)], md, old[seq.int(e + 1L, length(old))[seq.int(e + 1L, length(old)) <= length(old)]])
  } else {
    c(old, if (length(old)) "", md)
  }
  writeLines(new, results_md)
  cat("results written to", results_md, "\n")
} else {
  cat(md, sep = "\n")
}
