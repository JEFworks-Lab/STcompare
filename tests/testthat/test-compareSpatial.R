# Tests of compareSpatial() (R/compareSpatial.R; dev/compare-spatial-spec.md, section 8) and of the engine's
# independent streams (src/stc_rng.h, RNG_STREAMS in src/stc_engine.h): invariance (threads, batch schedule,
# gene order and subsets, the session's RNG kinds), adaptive stopping against an offline Besag-Clifford
# computation, the stream draws, the null calibration (slow tier), positive and negative controls, input
# errors and failed genes, the similarity helper against spatialSimilarity(), and the methods. Labels are
# explained in helper-fixtures.R. No legacy R implementation is run.

rf <- fx_read("realistic_fixture.rds")

# --- helpers ------------------------------------------------------------------------------------

# AKI genes of the realistic fixture as a pair of SpatialExperiment objects (311 shared pixels; with the
# default delta grid all 11 deltas are usable). With the shorter grid cs_delta (cheaper; delta star is then
# often its smallest value), exceedances = 3 and 60 permutations, and the default rank-remapped surrogates,
# Upk2, Agt, Calm1 and Fxyd3 stop early (after 3, 40, 47 and 59 permutations), and Gpx1 (no exceedance) and
# Slc6a20b (2 exceedances) reach the limit. (With gaussian surrogates Fxyd3 reaches the limit with 2
# exceedances; the other genes stop at the same permutations: test-surrogate-remap.R.)
cs_genes <- c("Gpx1", "Upk2", "Agt", "Calm1", "Fxyd3", "Slc6a20b")
cs_delta <- c(0.05, 0.1, 0.3, 0.6, 0.9)
cs_input <- function(genes = cs_genes, P = rf$pairs$aki) fx_spe_pair(P, genes)

# compareSpatial() without progress output or messages.
cs_run <- function(input, ..., nThreads = 1L) {
  compareSpatial(input, ..., nThreads = nThreads, progress = FALSE, verbose = FALSE)
}

# The table of a result, without the class and the attributes that describe the call.
cs_table <- function(res) as.data.frame(res)

# The value of an expression and its messages (muffled).
cs_messages <- function(expr) {
  m <- character(0)
  value <- withCallingHandlers(expr, message = function(cnd) {
    m <<- c(m, conditionMessage(cnd))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = m)
}

# The result with a batch schedule (option STcompare.compare_schedule; see .stc_compare_schedule()).
cs_with_schedule <- function(schedule, expr) {
  old <- options(STcompare.compare_schedule = schedule)
  on.exit(options(old))
  expr
}

# The exceedances of the default surrogate mode, "remap": |null| >= |r| within a relative 1e-9 (the ties of the
# permutation distribution, which cor() rounds differently for each rearrangement; kRemapTieRel in
# src/stc_engine.h). Gaussian surrogates count |null| >= |r| exactly (test-surrogate-remap.R).
cs_exceedances <- function(null, r) sum(abs(null) >= abs(r) * (1 - 1e-9))

# Besag-Clifford computed offline from complete sequences of nulls (h exceedances or n_max permutations,
# both directions in lockstep), in compareSpatial()'s terms and with the exceedance rule of the default mode.
cs_bc_offline <- function(nullX, nullY, r, h, n_max) {
  cx <- cumsum(abs(nullX) >= abs(r) * (1 - 1e-9))
  cy <- cumsum(abs(nullY) >= abs(r) * (1 - 1e-9))
  hit <- which(cx >= h | cy >= h)
  early <- length(hit) > 0L
  L <- if (early) hit[1] else n_max
  bX <- cx[L]
  bY <- cy[L]
  list(L = as.integer(L), stop = if (early) "exceedances" else "limit",
       p = if (early) max(bX, bY) / L else (max(bX, bY) + 1) / (L + 1),
       pX = (bX + 1) / (L + 1), pY = (bY + 1) / (L + 1))
}

# --- the result ---------------------------------------------------------------------------------------

test_that("portable: one row per gene with the specified columns, attributes and statistics", {
  genes <- c("Gpx1", "Upk2", "Agt")
  input <- cs_input(genes)
  res <- cs_run(input, nPermutations = 60, exceedances = 3, keepNulls = TRUE)
  expect_s3_class(res, c("STcompareResult", "data.frame"), exact = TRUE)
  expect_identical(names(res), c("gene", "nPixels", "r", "pNaive", "p", "padj", "pX", "pY", "nPermutations", "stop",
                                 "deltaStarMedianX", "deltaStarMedianY", "deltaGridEdge", "similarity",
                                 "dissimilarityX", "dissimilarityY", "nPixelsSimilarity", "thresholdX", "thresholdY",
                                 "status", "message"))
  expect_identical(res$gene, genes)
  expect_identical(rownames(res), genes)
  expect_false(any(vapply(res, is.list, NA)))
  expect_identical(res$nPixels, rep(311L, length(genes)))
  expect_true(all(res$status == "ok") && all(res$message == ""))
  expect_setequal(unique(res$stop), c("exceedances", "limit"))
  # r and the naive p-value are cor() and cor.test()'s
  P <- rf$pairs$aki
  for (g in genes) {
    ct <- stats::cor.test(P$X[g, ], P$Y[g, ])
    expect_identical(res[g, "r"], unname(ct$estimate), info = g)
    expect_identical(res[g, "pNaive"], ct$p.value, info = g)
  }
  expect_identical(res$padj, stats::p.adjust(res$p, "BH"))
  # the details hold the permutations behind every summary column
  det <- attr(res, "details")
  expect_identical(names(det), genes)
  grid <- attr(res, "params")$deltaGrid
  expect_identical(grid, c(0.01, 0.05, seq(0.1, 0.9, 0.1)))
  for (g in genes) {
    d <- det[[g]]
    L <- res[g, "nPermutations"]
    expect_length(d$nullX, L)
    expect_length(d$deltaStarY, L)
    bX <- cs_exceedances(d$nullX, res[g, "r"])
    bY <- cs_exceedances(d$nullY, res[g, "r"])
    expect_identical(res[g, "pX"], (bX + 1) / (L + 1), info = g)
    expect_identical(res[g, "pY"], (bY + 1) / (L + 1), info = g)
    expect_identical(res[g, "p"], if (res[g, "stop"] == "exceedances") max(bX, bY) / L else (max(bX, bY) + 1) / (L + 1),
                     info = g)
    expect_identical(res[g, "deltaStarMedianX"], stats::median(d$deltaStarX), info = g)
    expect_identical(res[g, "deltaStarMedianY"], stats::median(d$deltaStarY), info = g)
    at_edge <- function(ds) 2 * sum(ds %in% range(grid)) > length(ds)
    expect_identical(res[g, "deltaGridEdge"], at_edge(d$deltaStarX) || at_edge(d$deltaStarY), info = g)
  }
  # without the details, the table is the same (the engine then keeps only the current batch's nulls)
  expect_identical(cs_table(cs_run(input, nPermutations = 60, exceedances = 3)), cs_table(res))
  prm <- attr(res, "params")
  expect_identical(prm$samples, c("x", "y"))
  expect_identical(prm$nPixels, 311L)
  expect_identical(prm$exceedances, 3)
  expect_identical(prm$surrogate, "remap")      # the default mode, with no detection filter
  expect_identical(prm$minDetectedPixels, 0L)
  expect_true(is.call(attr(res, "call")))
  expect_s3_class(attr(res, "runtime"), "proc_time")
  # only the tests asked for
  expect_identical(names(cs_run(input, tests = "similarity")),
                   c("gene", "nPixels", "similarity", "dissimilarityX", "dissimilarityY", "nPixelsSimilarity",
                     "thresholdX", "thresholdY", "status", "message"))
  cor_only <- cs_run(input[[1]], input[[2]], genes = "Agt", tests = "correlation", nPermutations = 9)
  expect_identical(names(cor_only)[c(1:3, 13:15)], c("gene", "nPixels", "r", "deltaGridEdge", "status", "message"))
})

# --- invariance -----------------------------------------------------------------------------------

test_that("portable: results do not depend on threads, the batch schedule, gene order or the other genes", {
  input <- cs_input()
  run <- function(exceedances = 3, ...) {
    cs_table(cs_run(input, tests = "correlation", nPermutations = 60, exceedances = exceedances, delta = cs_delta, ...))
  }
  ref <- run()
  expect_identical(ref$stop, c("limit", "exceedances", "exceedances", "exceedances", "exceedances", "limit"))
  expect_identical(ref$nPermutations, c(60L, 3L, 40L, 47L, 59L, 60L))
  expect_true(any(ref$deltaGridEdge) && !all(ref$deltaGridEdge))
  expect_identical(run(nThreads = 2L), ref)
  for (s in list(list(first_batch = 1, growth = 1, max_batch = 7, chunk = 3),
                 list(first_batch = 64, max_batch = 64, chunk = 64),
                 list(first_batch = 5, growth = 3, max_batch = 1000, chunk = 1))) {
    expect_identical(cs_with_schedule(s, run(nThreads = 2L)), ref, info = paste(names(s), unlist(s), collapse = " "))
  }
  # fixed B: one batch, or many
  fixed <- run(exceedances = Inf)
  expect_identical(cs_with_schedule(list(max_batch = 13, chunk = 4), run(exceedances = Inf, nThreads = 2L)), fixed)
  # the genes in another order: the same rows (the adjustment does not depend on the order either)
  rev_genes <- rev(cs_genes)
  expect_identical(run(genes = rev_genes)[cs_genes, ], ref)
  rev_input <- lapply(input, function(s) s[rev_genes, ])
  expect_identical(cs_table(cs_run(rev_input, tests = "correlation", nPermutations = 60, exceedances = 3,
                                   delta = cs_delta))[cs_genes, ], ref)
  # a subset of the genes: the same rows except padj (adjusted over fewer genes)
  sub <- run(genes = c("Upk2", "Gpx1"))
  cols <- setdiff(names(ref), "padj")
  expect_identical(sub[, cols], ref[c("Upk2", "Gpx1"), cols])
  expect_identical(sub$padj, stats::p.adjust(sub$p, "BH"))
  # the pixels of y in another order (they are matched to the pixels of x by name): the same results
  shuffled <- list(input[[1]], input[[2]][, c(2:ncol(input[[2]]), 1L)])
  expect_identical(cs_table(cs_run(shuffled, tests = "correlation", nPermutations = 60, exceedances = 3,
                                   delta = cs_delta)), ref)
})

test_that("portable: compareSpatial() leaves the global RNG state alone and does not depend on the session's RNG kinds", {
  input <- cs_input(c("Gpx1", "Upk2"))
  run <- function(seed = 0L) cs_run(input, nPermutations = 20, exceedances = 3, delta = cs_delta, seed = seed, keepNulls = TRUE)
  same <- function(a, b) identical(cs_table(a), cs_table(b)) && identical(attr(a, "details"), attr(b, "details"))
  local_default_rng()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  ref <- run()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(42)
  before <- .Random.seed
  expect_true(same(run(), ref))
  expect_identical(.Random.seed, before)
  suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rounding"))
  set.seed(7)
  before <- .Random.seed
  expect_true(same(run(), ref))
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), c("L'Ecuyer-CMRG", "Box-Muller", "Rounding"))
  # another seed gives other permutations, hence other null correlations (with 20 permutations of these two
  # genes the summary table can coincide: Gpx1 has no exceedance and Upk2 stops at the third permutation)
  other <- run(1L)
  expect_false(identical(attr(other, "details")$Gpx1$nullX, attr(ref, "details")$Gpx1$nullX))
  expect_false(identical(attr(other, "details")$Upk2$nullY, attr(ref, "details")$Upk2$nullY))
})

test_that("portable: more than 1000 pixels: the variogram subsample depends on the seed only (default RNG kinds)", {
  # .stc_subsample() draws the legacy plan's subsample under R's default kinds whatever the session's kinds
  local_default_rng()
  set.seed(3)
  ids <- sample(2170, 1000)
  suppressWarnings(RNGkind("Knuth-TAOCP-2002", "Ahrens-Dieter", "Rounding"))
  set.seed(11)
  before <- .Random.seed
  expect_identical(.stc_subsample(2170L, 3L), ids)
  expect_identical(.Random.seed, before)
  expect_identical(.stc_subsample(311L, 3L), seq_len(311L))
  # through compareSpatial(): brain genes on 2170 pixels give the same results under other RNG kinds
  input <- cs_input(c("Slc17a7", "Efemp1"), P = rf$pairs$brain)
  run <- function() cs_table(cs_run(input, tests = "correlation", nPermutations = 4, exceedances = Inf, nThreads = 2L,
                                    delta = c(0.1, 0.3, 0.6)))
  a <- run()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  expect_identical(run(), a)
  expect_identical(a$nPixels, c(2170L, 2170L))
  expect_true(all(a$status == "ok"))
})

# --- adaptive stopping ------------------------------------------------------------------------------

test_that("portable: adaptive stopping equals Besag-Clifford computed offline on a fixed-length run from the same streams", {
  input <- cs_input()
  h <- 3
  n_max <- 60L
  fixed <- cs_run(input, tests = "correlation", nPermutations = n_max, exceedances = Inf, keepNulls = TRUE,
                  delta = cs_delta, nThreads = 2L)
  ad <- cs_run(input, tests = "correlation", nPermutations = n_max, exceedances = h, keepNulls = TRUE, delta = cs_delta)
  # exceedances = Inf is fixed B: every gene at n_max, p = max(pX, pY)
  expect_identical(fixed$stop, rep("limit", length(cs_genes)))
  expect_identical(fixed$nPermutations, rep(n_max, length(cs_genes)))
  expect_identical(fixed$p, pmax(fixed$pX, fixed$pY))
  expect_true(all(c("exceedances", "limit") %in% ad$stop))
  df <- attr(fixed, "details")
  da <- attr(ad, "details")
  for (g in cs_genes) {
    o <- cs_bc_offline(df[[g]]$nullX, df[[g]]$nullY, fixed[g, "r"], h, n_max)
    L <- o$L
    expect_identical(ad[g, "nPermutations"], L, info = g)
    expect_identical(ad[g, "stop"], o$stop, info = g)
    expect_identical(ad[g, "p"], o$p, info = g)
    expect_identical(c(ad[g, "pX"], ad[g, "pY"]), c(o$pX, o$pY), info = g)
    # the adaptive run's permutations are the first L of the fixed run's
    for (k in names(da[[g]])) expect_identical(da[[g]][[k]], df[[g]][[k]][seq_len(L)], info = paste(g, k))
  }
})

# --- the independent streams ----------------------------------------------------------------------

test_that("portable: the stream draws are deterministic, keyed by (seed, gene, direction, permutation), and well formed", {
  # pinned values (checked against an independent implementation of src/stc_rng.h): a change of the key
  # derivation or of the generators changes every compareSpatial() result
  p <- .stc_stream_draws(0L, "Gpx1", 1L, 1L, 10L, 2L)
  expect_identical(p$key, "7ba0969ba4971c6d")
  expect_identical(p$perm, c(7L, 10L, 2L, 6L, 3L, 8L, 4L, 1L, 9L, 5L))
  expect_equal(p$noise[1:3, 1], c(-1.5893806325655002, -0.32467893449899676, -0.62614724774548847), tolerance = 1e-14)
  expect_equal(p$noise[1:2, 2], c(0.55560043724370078, -0.18545684579292143), tolerance = 1e-14)
  d <- .stc_stream_draws(0L, "Gpx1", 1L, 1L, 50L, 3L)
  expect_identical(d$key, p$key)  # the key does not depend on the permutation
  expect_identical(d, .stc_stream_draws(0L, "Gpx1", 1L, 1L, 50L, 3L))
  expect_identical(sort(d$perm), 1:50)
  expect_identical(dim(d$noise), c(50L, 3L))
  differ <- function(a, b) !identical(a$perm, b$perm) && !any(a$noise == b$noise)
  expect_true(differ(d, .stc_stream_draws(0L, "Gpx1", 1L, 2L, 50L, 3L)))
  expect_true(differ(d, .stc_stream_draws(0L, "Gpx1", 2L, 1L, 50L, 3L)))
  expect_true(differ(d, .stc_stream_draws(1L, "Gpx1", 1L, 1L, 50L, 3L)))
  expect_true(differ(d, .stc_stream_draws(0L, "Gpx2", 1L, 1L, 50L, 3L)))
  expect_false(identical(d$key, .stc_stream_draws(-1L, "Gpx1", 1L, 1L, 5L, 0L)$key))
  # a noise block does not depend on how many blocks are drawn; the names are taken as UTF-8
  expect_identical(.stc_stream_draws(0L, "Gpx1", 1L, 1L, 50L, 1L)$noise[, 1], d$noise[, 1])
  latin1 <- iconv("G\u00e9ne", "UTF-8", "latin1")
  expect_identical(Encoding(latin1), "latin1")
  expect_identical(.stc_stream_draws(0L, latin1, 1L, 1L, 4L, 0L)$key, .stc_stream_draws(0L, "G\u00e9ne", 1L, 1L, 4L, 0L)$key)
  # standard normal noise and uniform permutations (Lemire's unbiased bounded integers)
  z <- as.vector(.stc_stream_draws(5L, "g", 2L, 3L, 20000L, 2L)$noise)
  expect_lt(abs(mean(z)), 0.03)
  expect_lt(abs(stats::sd(z) - 1), 0.02)
  expect_gt(stats::ks.test(z, "pnorm")$p.value, 1e-3)
  first <- vapply(1:3000, function(b) .stc_stream_draws(0L, "g", 1L, b, 6L, 0L)$perm[1], 0L)
  expect_gt(stats::chisq.test(table(factor(first, 1:6)))$p.value, 1e-3)
})

test_that("portable: with the same permutations and noise, the independent-stream engine computes what the legacy engine computes", {
  P <- rf$pairs$aki
  g <- "Upk2"
  X <- P$X[g, ]
  Y <- P$Y[g, ]
  N <- length(X)
  grid <- c(0.05, 0.2, 0.5, 0.8)
  B <- 12L
  X1 <- matrix(X, dimnames = list(NULL, g))
  Y1 <- matrix(Y, dimnames = list(NULL, g))
  st <- .stc_engine_correlate(X1, Y1, P$pos, deltaX = grid, deltaY = rev(grid), nPermutations = B, seed = 5,
                              streams = "independent", returnPermutations = TRUE)
  expect_identical(st$status, "ok")
  s <- attr(st, "state")
  # a legacy session fed with the draws of each direction's stream reproduces that direction bit for bit
  for (dir in 1:2) {
    draws <- lapply(seq_len(B), function(b) .stc_stream_draws(5L, g, dir, b, N, length(grid)))
    perm <- vapply(draws, `[[`, integer(N), "perm")
    noise <- array(vapply(draws, `[[`, numeric(N * length(grid)), "noise"), c(N, length(grid), B))
    legacy <- .stc_engine_new(P$pos[, 2], P$pos[, 1], seq_len(N) - 1L, s$plan, s$deltas, as.integer(.stc_cor_mode()),
                              0, TRUE, 2L, 5L, TRUE)
    .stc_engine_define(legacy, cbind(X, Y), c(0L, 1L),
                       list(match(grid, s$deltas) - 1L, match(rev(grid), s$deltas) - 1L), list(1L, 0L),
                       list(abs(s$r$X), abs(s$r$Y)), c(0L, 0L), TRUE)
    .stc_engine_run(legacy, 0L, 1L, B, perm, noise, Inf, B)
    lt <- .stc_engine_task_results(legacy, dir - 1L, TRUE)[[1]]
    d <- c("X", "Y")[dir]
    expect_identical(lt$nulls[, 1], st[[paste0("null", d)]][[1]], info = d)
    expect_identical((if (dir == 1) grid else rev(grid))[lt$dstar], st[[paste0("deltaStar", d)]][[1]], info = d)
    expect_identical(lt$surrogates, st[[paste0("permutations", d)]][[1]], info = d)
  }
  # the directions of a gene draw from different streams
  expect_false(identical(st$nullX[[1]], st$nullY[[1]]))
})

test_that("portable: a progress callback runs on the main thread while the workers run; an error in it stops them and becomes an R error", {
  # 20 items of 50 ms on 2 threads take about 0.5 s: the main thread reports at least once (every 200 ms)
  seen <- numeric(0)
  w <- .stc_parallel_selftest(20L, 2L, -1L, 50L, progress = function(done) seen <<- c(seen, done))
  expect_length(w, 20L)
  expect_true(length(seen) >= 1L && all(diff(seen) >= 0) && all(seen <= 20))
  # 40 items of 50 ms on 1 thread would take 2 s; the error at the first report (after about 0.2 s) stops them
  t0 <- Sys.time()
  expect_error(.stc_parallel_selftest(40L, 1L, -1L, 50L, progress = function(done) stop("progress failed")),
               "progress failed")
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1.5)
  expect_length(.stc_parallel_selftest(10L, 2L), 10L)  # the pool still works afterwards
})

# --- calibration and controls -----------------------------------------------------------------------

test_that("portable, slow: calibration: on independent simulated fields the share of p < 0.05 is not above a binomial bound; correlated fields are detected", {
  skip_if_not_slow()
  cf <- fx_read("calibration_fixture.rds")
  jobs <- cf$test_jobs
  spe <- function(v, px, coords) {
    SpatialExperiment::SpatialExperiment(assays = list(pixelval = matrix(v, 1L, dimnames = list("g", px))),
                                         spatialCoords = coords)
  }
  run1 <- function(k) {
    i <- jobs$i[k]
    j <- jobs$j[k]
    sh <- which(!is.na(cf$fields[i, ]) & !is.na(cf$fields[j, ]))
    X <- unname(cf$fields[i, sh])
    Y <- unname(cf$fields[j, sh])
    if (jobs$rho[k] != 0) Y <- stc_mix(X, Y, jobs$rho[k], cf$mix_mu)
    px <- colnames(cf$fields)[sh]
    coords <- cf$coords[sh, , drop = FALSE]
    cs_run(list(spe(X, px, coords), spe(Y, px, coords)), tests = "correlation", nPermutations = 999,
           nThreads = fx_threads, seed = k)
  }
  t0 <- Sys.time()
  res <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(k) cbind(cs_table(run1(k)), rho = jobs$rho[k])))
  null <- res[res$rho == 0, ]
  message(sprintf("calibration: %d compareSpatial() calls in %.1f s; null: p < 0.05 in %d of %d (median %d permutations); rho = 0.6: %d of %d",
                  nrow(res), as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(null$p < 0.05), nrow(null),
                  as.integer(stats::median(null$nPermutations)), sum(res$p[res$rho != 0] < 0.05), sum(res$rho != 0)))
  expect_lte(sum(null$p < 0.05), stats::qbinom(0.999, nrow(null), 0.05))
  expect_true(all(null$status == "ok"))
  expect_lt(stats::median(null$nPermutations), 200)  # null genes stop early
  # the naive test is anti-conservative on these fields; correlated fields (rho = 0.6) are mostly detected
  expect_gt(mean(null$pNaive < 0.05), 0.05)
  expect_gte(mean(res$p[res$rho != 0] < 0.05), 0.7)
})

test_that("portable: speKidney: A vs B (negative) and A vs C (positive) are significant with 99 permutations", {
  rk <- fx_speKidney_raster()
  ab <- cs_run(list(A = rk$A, B = rk$B), nPermutations = 99)
  ac <- cs_run(list(A = rk$A, C = rk$C), nPermutations = 99)
  expect_identical(attr(ab, "params")$samples, c("A", "B"))
  expect_lt(ab$r, -0.9)
  expect_gt(ac$r, 0.9)
  for (res in list(ab, ac)) {
    expect_identical(res$p, 1 / 100)
    expect_identical(res$stop, "limit")
    expect_lt(res$padj, 0.05)
  }
  # opposite patterns are dissimilar, similar patterns at another level are dissimilar too: C is higher
  expect_gt(ac$dissimilarityY, 0.9)
})

# --- errors, dropped deltas and failed genes ---------------------------------------------------------------

test_that("portable: pixels whose coordinates differ between the samples are an error", {
  input <- cs_input(c("Gpx1", "Upk2"))
  y <- input[[2]]
  xy <- SpatialExperiment::spatialCoords(y)
  xy[17, 1] <- xy[17, 1] + 0.5
  SpatialExperiment::spatialCoords(y) <- xy
  expect_error(compareSpatial(input[[1]], y, verbose = FALSE),
               paste0("1 of the 311 pixels shared by x and y have different coordinates in the two samples \\(for example ",
                      colnames(y)[17], ", which is 0.5 apart\\).*rasterize both samples in one call"))
  # a shift below 1e-6 of the coordinate range is tolerated
  xy[17, 1] <- SpatialExperiment::spatialCoords(input[[2]])[17, 1] + 1e-9
  SpatialExperiment::spatialCoords(y) <- xy
  expect_s3_class(cs_run(input[[1]], y, tests = "similarity"), "STcompareResult")
  # samples rasterized separately share pixel names at other locations
  B <- STcompare::speKidney$B
  B <- B[, SpatialExperiment::spatialCoords(B)[, 2] > 1.5]
  ras <- function(s) SEraster::rasterizeGeneExpression(s, assay_name = "counts", resolution = 0.2, fun = "mean",
                                                        square = FALSE, BPPARAM = BiocParallel::SerialParam())
  expect_error(compareSpatial(ras(STcompare::speKidney$A), ras(B), verbose = FALSE), "different coordinates")
})

test_that("portable: deltas too small for the number of pixels are dropped with a message; none left is an error", {
  P <- rf$pairs$aki
  keep <- seq_len(150)
  small <- P
  small$X <- P$X[, keep]
  small$Y <- P$Y[, keep]
  small$pos <- P$pos[keep, ]
  small$pixel <- P$pixel[keep]
  input <- fx_spe_pair(small, c("Gpx1", "Upk2"))
  m <- cs_messages(compareSpatial(input, tests = "correlation", nPermutations = 9, progress = FALSE))
  expect_match(m$messages[1], "^compareSpatial: delta 0.01 dropped: with 150 shared pixels the smoothing window must hold at least 2 pixels \\(delta >= 0.0133\\)\n$")
  expect_identical(attr(m$value, "params")$deltaGrid, c(0.05, seq(0.1, 0.9, 0.1)))
  expect_error(cs_run(input, delta = c(0.001, 0.01), nPermutations = 9), "no delta is usable with 150 shared pixels")
})

test_that("portable: constant and rarely detected genes are skipped with a message, failed genes warn; the others are unaffected", {
  P <- rf$pairs$aki
  genes <- c("Gpx1", "Ech1", "Upk2", "Tspan8", "Agt")
  P$X["Ech1", ] <- 3              # constant in x: not tested, the similarity is still computed
  P$Y["Tspan8", 12] <- NA         # a missing value: nothing is computed
  P$X["Agt", ] <- 0               # detected in 5 pixels of x: fewer than sqrt(311), so not tested with gaussian surrogates
  P$X["Agt", 1:5] <- c(4, 1, 7, 2, 9)
  input <- fx_spe_pair(P, genes)
  # the sqrt(N) filter is the default of gaussian surrogates (with remapped ones, the default, see below)
  run <- function(...) cs_run(input, nPermutations = 30, exceedances = 3, delta = cs_delta, surrogate = "gaussian", ...)
  w <- fx_collect_warnings(cs_messages(compareSpatial(input, nPermutations = 30, exceedances = 3, delta = cs_delta,
                                                      surrogate = "gaussian", progress = FALSE, nThreads = 1L)))
  expect_length(w$warnings, 1L)
  expect_match(w$warnings, "^compareSpatial: 1 of 5 genes failed .*: Tspan8 \\(missing or infinite values in y\\)$")
  expect_match(w$value$messages[1], paste0("^compareSpatial: 2 of 5 genes are not tested for correlation \\(status \"skipped\"\\): ",
                                           "1 constant in a sample and 1 detected in fewer than 18 of the 311 shared pixels"))
  res <- w$value$value
  expect_identical(res$status, c("ok", "skipped", "ok", "failed", "skipped"))
  expect_identical(res$stop[c(2, 4, 5)], c("skipped", "failed", "skipped"))
  expect_identical(res$message[2], "constant in x (zero variance): not tested")
  kY <- sum(P$Y["Agt", ] > min(P$Y["Agt", ]))
  expect_identical(res$message[5], sprintf("detected in 5 (x) and %d (y) of the 311 shared pixels; minDetected asks for 18 in each: not tested", kY))
  for (col in c("p", "padj", "pX", "pY", "nPermutations", "deltaStarMedianX", "deltaGridEdge")) {
    expect_true(all(is.na(res[c("Ech1", "Tspan8", "Agt"), col])), info = col)
  }
  expect_true(is.na(res["Tspan8", "r"]) && is.na(res["Tspan8", "similarity"]))
  expect_false(is.na(res["Ech1", "similarity"]))
  expect_false(is.na(res["Agt", "r"]))  # the observed correlation is still reported
  # the adjustment counts only the tested genes; the tested genes' results do not depend on the others
  expect_identical(res$padj[c(1, 3)], stats::p.adjust(res$p[c(1, 3)], "BH"))
  ok <- cs_table(cs_run(input, genes = c("Gpx1", "Upk2"), nPermutations = 30, exceedances = 3, delta = cs_delta,
                        surrogate = "gaussian"))
  expect_identical(cs_table(res)[c("Gpx1", "Upk2"), ], ok)
  # minDetected = 0 tests every gene that is not constant; a larger share skips more genes
  expect_identical(suppressWarnings(run(minDetected = 0))$status, c("ok", "skipped", "ok", "failed", "ok"))
  expect_identical(attr(suppressWarnings(run(minDetected = 0.5)), "params")$minDetectedPixels, 156L)
  # with remapped surrogates (the default) there is no detection filter: Agt is tested, the message counts only the
  # constant gene, and the Agt row is the minDetected = 0 one
  w <- fx_collect_warnings(cs_messages(compareSpatial(input, nPermutations = 30, exceedances = 3, delta = cs_delta,
                                                      progress = FALSE, nThreads = 1L)))
  expect_length(w$warnings, 1L)
  expect_match(w$value$messages[1], "^compareSpatial: 1 of 5 genes are not tested for correlation \\(status \"skipped\"\\): 1 constant in a sample\n$")
  dflt <- w$value$value
  expect_identical(dflt$status, c("ok", "skipped", "ok", "failed", "ok"))
  expect_identical(attr(dflt, "params")$minDetectedPixels, 0L)
  expect_false(is.na(dflt["Agt", "p"]))
  expect_identical(cs_table(dflt), cs_table(suppressWarnings(cs_run(input, nPermutations = 30, exceedances = 3,
                                                                    delta = cs_delta, minDetected = 0))))
})

test_that("portable: invalid inputs and arguments give errors that say what to do", {
  input <- cs_input(c("Gpx1", "Upk2"))
  x <- input[[1]]
  y <- input[[2]]
  expect_error(compareSpatial(x), "two SpatialExperiment objects")
  expect_error(compareSpatial(list(x, y, y)), "two SpatialExperiment objects")
  expect_error(compareSpatial(x, SummarizedExperiment::assay(y)), "y must be a SpatialExperiment")
  expect_error(compareSpatial(x, y, assay = "logcounts"), "x has no assay \"logcounts\"; its assays are: counts")
  expect_error(compareSpatial(x, y, assay = 2), "x has no assay 2")
  expect_error(compareSpatial(x, y, genes = c("Gpx1", "Nope")), "1 of the genes are not in x, for example Nope")
  expect_error(compareSpatial(x, y, genes = c("Gpx1", "Gpx1")), "duplicates")
  dup <- x
  rownames(dup) <- c("Gpx1", "Gpx1")
  expect_error(compareSpatial(dup, y), "gene names \\(row names\\) of x must be unique")
  expect_error(compareSpatial(x[, 1:2], y), "share 2 pixel")
  expect_error(compareSpatial(x, y, tests = "power"), "tests must be")
  expect_error(compareSpatial(x, y, nPermutations = 0), "nPermutations must be a positive integer")
  expect_error(compareSpatial(x, y, exceedances = 0), "exceedances must be a positive whole number or Inf")
  expect_error(compareSpatial(x, y, exceedances = 2.5), "exceedances")
  expect_error(compareSpatial(x, y, delta = c(0.1, 1.5)), "delta must be a numeric vector of values in \\(0, 1\\]")
  expect_error(compareSpatial(x, y, maxDistPrctile = 0), "maxDistPrctile")
  expect_error(compareSpatial(x, y, minQuantile = 2), "minQuantile")
  expect_error(compareSpatial(x, y, adjustMethod = "bh"), "adjustMethod must be one of")
  expect_error(compareSpatial(x, y, seed = 1.5), "seed must be a single integer")
  expect_error(compareSpatial(x, y, nThreads = 0), "positive integer")
  expect_error(compareSpatial(x, y, nThreads = 1.5), "positive integer")
  expect_error(compareSpatial(x, y, nThreads = "2"), "positive integer")
  expect_error(compareSpatial(x, y, minDetected = 2), "minDetected must be a single number in \\[0, 1\\]")
  expect_error(compareSpatial(x, y, keepNulls = NA), "keepNulls must be TRUE or FALSE")
  expect_error(compareSpatial(x, y, maxDistPrctile = 1e-6, verbose = FALSE),
               "variograms cannot be computed .* Try a larger maxDistPrctile")
  # genes present in one sample only are dropped with a message
  m <- cs_messages(compareSpatial(x, y[c("Upk2"), ], tests = "similarity", progress = FALSE))
  expect_match(m$messages[1], "comparing the 1 genes in both samples \\(1 genes only in x and 0 only in y are dropped\\)")
  expect_identical(m$value$gene, "Upk2")
})

# --- similarity ------------------------------------------------------------------------------------------

test_that("portable: the similarity columns equal spatialSimilarity()'s on the fixtures and on random inputs", {
  cmp <- function(res, ss, info) {
    ss <- ss$similarityTable
    expect_identical(res$similarity, as.double(ss$percentSimilarity), info = info)
    expect_identical(res$dissimilarityX, as.double(ss$percentDissimilarityX), info = info)
    expect_identical(res$dissimilarityY, as.double(ss$percentDissimilarityY), info = info)
    expect_identical(res$thresholdX, as.double(ss$t1), info = info)
    expect_identical(res$thresholdY, as.double(ss$t2), info = info)
    expect_identical(res$nPixelsSimilarity, ss$numPixelInThresh, info = info)
  }
  P <- rf$pairs$aki
  input <- fx_spe_pair(P)
  cmp(cs_run(input, tests = "similarity"), spatialSimilarity(input), "AKI")
  cmp(cs_run(input, tests = "similarity", foldChange = 0.5, minQuantile = 0.5, minPixels = 0.6),
      spatialSimilarity(input, foldChange = 0.5, minQuantile = 0.5, minPixels = 0.6), "AKI, other settings")
  rk <- fx_speKidney_raster()
  for (s in c("B", "C")) cmp(cs_run(rk$A, rk[[s]], tests = "similarity"), spatialSimilarity(list(rk$A, rk[[s]])), s)
  # random inputs: zeros, sparse genes below minPixels, ties
  local_default_rng()
  set.seed(17)
  N <- 120
  G <- 12
  m <- function() {
    v <- matrix(stats::rexp(N * G, 0.2) * stats::rbinom(N * G, 1, rep(rep(c(0.9, 0.5, 0.02), length.out = G), each = N)), N, G)
    v[, 4] <- round(v[, 4])
    v
  }
  X <- m()
  Y <- m()
  dimnames(X) <- dimnames(Y) <- list(sprintf("p%03d", seq_len(N)), sprintf("g%02d", seq_len(G)))
  pos <- cbind(stats::runif(N), stats::runif(N))
  rownames(pos) <- rownames(X)
  spe <- function(M) SpatialExperiment::SpatialExperiment(assays = list(v = t(M)), spatialCoords = pos)
  ss <- spatialSimilarity(list(spe(X), spe(Y)))
  expect_true(any(is.na(ss$similarityTable$percentSimilarity)))
  cmp(cs_run(spe(X), spe(Y), tests = "similarity"), ss, "random")
  # the helper also takes fixed thresholds, as spatialSimilarity(t1, t2)
  h <- .stc_similarity(X, Y, t1 = 2, t2 = 3)
  ss2 <- spatialSimilarity(list(spe(X), spe(Y)), t1 = 2, t2 = 3)$similarityTable
  expect_identical(h$similarity, as.double(ss2$percentSimilarity))
  expect_identical(h$dissimilarityY, as.double(ss2$percentDissimilarityY))
  # negative values: spatialSimilarity() stops; compareSpatial() gives NA similarity with one warning, and the
  # other genes and the correlation test are not affected
  X[, 5] <- X[, 5] - 2
  expect_error(spatialSimilarity(list(spe(X), spe(Y))), "negative values in the first object, for example g05")
  w <- fx_collect_warnings(cs_run(spe(X), spe(Y), genes = c("g01", "g05"), nPermutations = 9, exceedances = 3))
  expect_length(w$warnings, 1L)
  expect_match(w$warnings, "the similarity compares fold changes, .* it is NA for 1 gene with negative values: g05$")
  res <- w$value
  expect_true(all(is.na(unlist(res["g05", c("similarity", "dissimilarityX", "dissimilarityY", "nPixelsSimilarity")]))))
  expect_identical(res$message, c("", "negative values in x: no similarity"))
  expect_identical(res$status, c("ok", "ok"))
  expect_false(is.na(res["g05", "p"]))
  expect_identical(res["g01", "similarity"], cs_run(spe(X), spe(Y), genes = "g01", tests = "similarity")$similarity)
})

# --- methods and progress --------------------------------------------------------------------------------

test_that("portable: print(), summary(), [ and as.data.frame() work; progress = TRUE gives the same result", {
  input <- cs_input(c("Gpx1", "Agt", "Upk2"))
  res <- cs_run(input, nPermutations = 30, exceedances = 3, delta = cs_delta)
  out <- utils::capture.output(print(res))
  expect_match(out[1], "^STcompareResult: 3 genes on 311 shared pixels; x vs y \\(assay counts\\)")
  expect_match(out[2], "adaptive permutation p-values \\(stop at 3 exceedances, at most 30 permutations\\); 5 deltas from 0.05 to 0.9")
  expect_true(any(grepl("padj < 0.05: [0-9]+ genes? \\([0-9]+ with r > 0, [0-9]+ with r < 0\\)", out)))
  expect_true(any(grepl("stopped early: [0-9]+; reached nPermutations: [0-9]+; skipped: 0; failed: 0", out)))
  expect_true(any(grepl("Gpx1", out)))
  s <- summary(res)
  expect_s3_class(s, "summary.STcompareResult")
  expect_identical(s$genes, 3L)
  expect_identical(s$permutations[["total"]], as.double(sum(res$nPermutations)))
  expect_output(print(s), "3 ok, 0 skipped, 0 failed")
  expect_output(print(s), "permutations per gene: median")
  df <- as.data.frame(res)
  expect_identical(class(df), "data.frame")
  expect_null(attr(df, "params"))
  expect_identical(names(df), names(res))
  expect_output(print(cs_run(input, tests = "similarity")), "Similarity \\(\\|log2\\(y / x\\)\\| <= 1\\)")
  # rows selected: still a result with its settings; columns selected: a plain data frame
  rows <- res[c("Upk2", "Gpx1"), ]
  expect_s3_class(rows, "STcompareResult")
  expect_identical(attr(rows, "params"), attr(res, "params"))
  expect_output(print(rows), "STcompareResult: 2 genes on 311 shared pixels")
  expect_identical(class(res[, c("gene", "p")]), "data.frame")
  expect_identical(class(res["p"]), "data.frame")
  expect_identical(res[, "p"], res$p)
  expect_identical(summary(res[, c("gene", "p")]), summary(df[, c("gene", "p")]))
  # a result whose essential columns were removed is printed as a data frame
  res2 <- res
  res2$status <- NULL
  expect_output(print(res2), "^ +gene")
  # a huge exceedances is printed as it is
  big <- cs_run(input, genes = "Upk2", nPermutations = 5, exceedances = 1e10, delta = cs_delta)
  expect_output(print(big), "stop at 1e\\+10 exceedances")
  # nThreads defaults to the option STcompare.nThreads
  old <- options(STcompare.nThreads = 2L)
  on.exit(options(old))
  expect_identical(attr(compareSpatial(input, tests = "similarity", progress = FALSE, verbose = FALSE), "params")$nThreads, 2L)
})

test_that("portable: the progress line fits the console, ends with the permutations of the result, and works for fixed B", {
  input <- cs_input(c("Gpx1", "Agt", "Upk2"))
  res <- cs_run(input, nPermutations = 30, exceedances = 3, delta = cs_delta)
  old <- options(width = 120L)
  on.exit(options(old))
  shown <- cs_messages(compareSpatial(input, nPermutations = 30, exceedances = 3, delta = cs_delta, progress = TRUE,
                                      verbose = TRUE, nThreads = 2L))
  m <- shown$messages
  expect_identical(cs_table(shown$value), cs_table(res))
  total <- sum(res$nPermutations)
  expect_true(any(grepl(sprintf("^\rcompareSpatial: 100%% \\| genes done 3/3 \\| %d permutations \\| elapsed 0:0[0-9] \\| ETA 0:00 *\n$",
                                total), m)))
  expect_true(any(grepl(sprintf("^compareSpatial: 3 genes on 311 shared pixels in [0-9.]+ s \\(%d permutations in total", total), m)))
  expect_silent(compareSpatial(input, nPermutations = 30, exceedances = 3, delta = cs_delta, progress = FALSE,
                               verbose = FALSE))
  # a narrow console: the elapsed time is left out, then the line is cut
  options(width = 70L)
  progress_lines <- function() {
    m <- cs_messages(compareSpatial(input, nPermutations = 30, exceedances = 3, delta = cs_delta, progress = TRUE,
                                    verbose = FALSE))$messages
    sub("\n$", "", sub("^\r", "", m))
  }
  narrow <- progress_lines()
  expect_true(all(nchar(narrow) <= 69L))
  expect_true(any(grepl(sprintf("^compareSpatial: 100%% \\| genes done 3/3 \\| %d permutations \\| ETA 0:00 *$", total),
                        narrow)))
  options(width = 40L)
  expect_true(all(nchar(progress_lines()) <= 39L))
  options(old)
  # fixed B (exceedances = Inf) with fewer permutations than a batch
  fixed <- cs_messages(compareSpatial(input, nPermutations = 19, exceedances = Inf, delta = cs_delta, progress = TRUE,
                                      verbose = FALSE))
  expect_true(any(grepl("^\rcompareSpatial: 100% \\| 3 genes x 19 permutations \\|", fixed$messages)))
  expect_identical(fixed$value$nPermutations, rep(19L, 3))
})

test_that("portable: an R error raised by the interrupt check (a time limit) stops the workers and stays an R error", {
  res <- tryCatch({
    setTimeLimit(elapsed = 0.5, transient = TRUE)
    .stc_parallel_selftest(60L, 1L, -1L, 50L)
  }, error = function(e) e, interrupt = function(i) i, finally = setTimeLimit(elapsed = Inf))
  expect_s3_class(res, "error")
  expect_match(conditionMessage(res), "time limit")
  expect_length(.stc_parallel_selftest(10L, 2L), 10L)  # the pool still works afterwards
})
