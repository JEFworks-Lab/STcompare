# Tests of the compiled engine (src/stc_engine.*, src/stc_engine_rcpp.cpp, R/engine.R; dev/engine-spec.md
# sections 2 to 6) and of the exported functions built on it: NA rows and warnings, inputs that crashed the
# original R implementation, determinism (threads, sub-chunks, batches), adaptive stopping, exceedance ties,
# the global RNG state, session continuation (iterative rounds), the within-sample mode, threads, interrupts
# and invalid inputs. The agreement with the original R implementation is tested against its stored outputs
# (test-reference-*.R); nothing here runs it. Labels are explained in helper-fixtures.R.

fx <- fx_read("kernel_fixture.rds")
rf <- fx_read("realistic_fixture.rds")

# --- helpers ------------------------------------------------------------------------------------

# The engine with the arguments of spatialCorrelation(); X and Y may hold several genes (columns).
eng <- function(X, Y, pos, B = 10, deltaX = NULL, deltaY = NULL, p = 0.25, seed = 0, ...) {
  .stc_engine_correlate(X, Y, pos, deltaX = deltaX, deltaY = deltaY, nPermutations = B, seed = seed,
                        maxDistPrctile = p, ...)
}

# A result without its session (for identical() comparisons).
strip <- function(res) {
  attr(res, "state") <- NULL
  res
}

# The AKI genes used by the determinism and adaptive tests (311 pixels, 11 deltas including 0.01 and 0.05).
aki_genes <- c("Gpx1", "Ech1", "Upk2", "Tspan8", "Mrps6", "Slc5a3")
aki_input <- function(genes = aki_genes) {
  P <- rf$pairs$aki
  list(X = t(P$X[genes, , drop = FALSE]), Y = t(P$Y[genes, , drop = FALSE]), pos = P$pos, delta = P$params$delta)
}

# Messages of an expression (muffled) and its value.
collect_messages <- function(expr) {
  m <- character(0)
  value <- withCallingHandlers(expr, message = function(cnd) {
    m <<- c(m, conditionMessage(cnd))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = m)
}

# --- NA rows and crashes -------------------------------------------------------------------------

test_that("portable: NA rows exactly where the R implementation gave them, and results elsewhere", {
  inp <- fx$cases$kidney_AB_jitter$input
  N <- length(inp$X)
  base <- list(X = inp$X, Y = inp$Y, pos = inp$pos, dX = 0.3, dY = 0.3, p = 0.25)
  # outcome of the original R implementation for each input (established by comparing with it)
  cases <- list(
    constant_zero = list(X = rep(0, N), na = TRUE),
    constant = list(X = rep(0.1, N), na = TRUE),
    na_X = list(X = replace(inp$X, 5, NA), na = TRUE),
    inf_X = list(X = replace(inp$X, 9, Inf), na = TRUE),
    na_Y = list(Y = replace(inp$Y, 7, NA), na = TRUE),
    delta_N_below_1 = list(dX = c(0.3, 0.5 / N), na = TRUE),       # k = 0: locfit failed
    delta_zero = list(dY = c(0, 0.3), na = TRUE),
    delta_nan = list(dX = c(0.3, NaN), na = TRUE),
    no_pair_below_max_dist = list(p = 0, na = TRUE),
    invalid_prctile = list(p = 2, na = TRUE),                       # quantile() failed
    one_bin = list(p = 1e-4, na = TRUE),                            # lm() with 1 bin: NA slope
    no_bins = list(p = 5e-5, na = TRUE),                            # lm() with 0 bins failed
    two_bins = list(p = 2e-4, na = FALSE),
    delta_above_1 = list(dX = c(0.3, 1.5), na = FALSE),
    duplicated_coordinates = list(pos = rbind(inp$pos, inp$pos[1:3, ]), X = c(inp$X, inp$X[1:3] + 1),
                                  Y = c(inp$Y, inp$Y[4:6]), na = FALSE))  # nugget bin
  for (nm in names(cases)) {
    a <- utils::modifyList(base, cases[[nm]])
    e <- eng(a$X, a$Y, a$pos, B = 2, deltaX = a$dX, deltaY = a$dY, p = a$p)
    if (a$na) {
      expect_identical(e$status, "failed", info = nm)
      expect_true(nzchar(e$message), info = nm)
      expect_true(all(is.na(c(e$pX, e$pY, e$p, e$L, e$bX, e$bY))), info = nm)
      for (col in c("nullX", "nullY", "deltaStarX", "deltaStarY")) expect_identical(e[[col]][[1]], NA, info = paste(nm, col))
    } else {
      expect_identical(e$status, "ok", info = paste(nm, e$message))
      for (d in c("X", "Y")) {
        null <- e[[paste0("null", d)]][[1]]
        expect_true(all(is.finite(null)), info = nm)
        expect_identical(e[[paste0("p", d)]], stc_empirical_p(null, e$r), info = nm)
        expect_true(all(e[[paste0("deltaStar", d)]][[1]] %in% a[[paste0("d", d)]]), info = nm)
      }
    }
  }
})

test_that("portable: inputs that crashed the R implementation give an NA row and a reason, never a crash", {
  inp <- fx$cases$kidney_AB_jitter$input
  N <- length(inp$X)
  # 1 <= delta * N < 2: k = 1, where locfit killed R
  e <- eng(inp$X, inp$Y, inp$pos, B = 2, deltaX = c(0.3, 1.5 / N), deltaY = 0.3)
  expect_identical(e$status, "failed")
  expect_match(e$message, "permuting X: delta = .*N \\* delta < 2")
  # delta < 0: locfit segfaulted
  e <- eng(inp$X, inp$Y, inp$pos, B = 2, deltaX = 0.3, deltaY = c(0.2, -0.1))
  expect_identical(e$status, "failed")
  expect_match(e$message, "permuting Y: delta = -0.1: invalid delta")
  # every location twice with k = 2: locfit's tree refined without end and overflowed the C stack; the
  # engine's tree runs out of vertex space instead (status 4)
  local_default_rng()
  set.seed(9)
  u <- cbind(stats::runif(150), stats::runif(150))
  xy <- u[rep(seq_len(150), each = 2), ]
  e <- eng(stats::rnorm(300), stats::rnorm(300), xy, B = 2, deltaX = c(0.3, 2 / 300), deltaY = 0.3)
  expect_identical(e$status, "failed")
  expect_match(e$message, "out of vertex space")
  info <- .stc_engine_deltas(attr(e, "state")$session)
  expect_identical(info$status[info$delta == 2 / 300], 4L)
})

test_that("portable: a failing gene leaves the other genes' results unchanged; spatialCorrelationGeneExp() warns once for it", {
  a <- aki_input(c("Gpx1", "Ech1", "Upk2"))
  X <- a$X
  X[, 2] <- 0  # constant: an NA row
  e <- strip(eng(X, a$Y, a$pos, B = 9, deltaX = a$delta, deltaY = a$delta))
  expect_identical(e$status, c("ok", "failed", "ok"))
  expect_match(e$message[2], "permuting X: the permuted values are constant")
  for (j in c(1, 3)) {
    e1 <- strip(eng(X[, j, drop = FALSE], a$Y[, j, drop = FALSE], a$pos, B = 9, deltaX = a$delta, deltaY = a$delta))
    rownames(e1) <- rownames(e)[j]
    expect_identical(e[j, ], e1, info = j)
  }
  # the exported function: one warning naming the gene, NA in every permutation column
  P <- rf$pairs$aki
  P$X["Ech1", ] <- 0
  delta <- rep(list(P$params$delta), 3)
  w <- fx_collect_warnings(spatialCorrelationGeneExp(fx_spe_pair(P, c("Gpx1", "Ech1", "Upk2")), nPermutations = 9,
                                                     deltaX = delta, deltaY = delta, verbose = FALSE))
  expect_length(w$warnings, 1L)
  expect_match(w$warnings, "^spatialCorrelationGeneExp: gene Ech1: no permutation p-values \\(NA row\\): permuting X: the permuted values are constant")
  o <- w$value
  expect_true(all(is.na(c(o["Ech1", "correlationCoef"], o["Ech1", "pValuePermuteX"], o["Ech1", "deltaStarMedianY"]))))
  expect_identical(o$nullCorrelationsY[2], I(list(NA)))
  expect_identical(unname(as.numeric(o["Gpx1", "nullCorrelationsX"][[1]])), e$nullX[[1]])
})

test_that("portable: sparse (dgCMatrix) assays and single-gene inputs give the same results", {
  P <- rf$pairs$aki
  genes <- c("Upk2", "Gpx1")
  delta <- rep(list(P$params$delta[3:6]), 2)
  run <- function(input, d = delta) {
    spatialCorrelationGeneExp(input, nPermutations = 6, deltaX = d, deltaY = d, verbose = FALSE, adjustMethod = "none")
  }
  dense <- fx_spe_pair(P, genes)
  sparse <- lapply(dense, function(s) {
    SummarizedExperiment::assay(s) <- methods::as(SummarizedExperiment::assay(s), "CsparseMatrix")
    s
  })
  expect_s4_class(SummarizedExperiment::assay(sparse[[1]]), "dgCMatrix")
  od <- run(dense)
  expect_identical(run(sparse), od)
  o1 <- run(lapply(dense, function(s) s["Upk2", ]), delta[1])
  expect_identical(o1, od["Upk2", ])
})

# --- determinism -----------------------------------------------------------------------------------

test_that("portable: results are identical for any number of threads, sub-chunk size and batch size", {
  a <- aki_input(c("Gpx1", "Ech1", "Upk2", "Tspan8"))
  run <- function(...) strip(eng(a$X, a$Y, a$pos, B = 19, deltaX = a$delta, deltaY = rev(a$delta), seed = 3,
                                 returnPermutations = TRUE, ...))
  ref <- run(nThreads = 1)
  expect_true(all(ref$status == "ok"))
  for (args in list(list(nThreads = 2), list(nThreads = 2, chunk = 1), list(nThreads = 2, chunk = 5),
                    list(nThreads = 2, chunk = 64), list(nThreads = 1, batch = 1), list(nThreads = 2, batch = 7, chunk = 3))) {
    expect_identical(do.call(run, args), ref, info = paste(names(args), unlist(args), collapse = " "))
  }
})

test_that("portable, slow: results are identical on 4 and 8 threads (more threads than the default suite uses)", {
  skip_if_not_slow()
  a <- aki_input()
  run <- function(...) strip(eng(a$X, a$Y, a$pos, B = 37, deltaX = a$delta, deltaY = rev(a$delta), seed = 5,
                                 returnPermutations = TRUE, ...))
  ref <- run(nThreads = 1)
  expect_identical(run(nThreads = 4), ref)
  expect_identical(run(nThreads = 8, chunk = 3), ref)
  expect_identical(run(nThreads = 4, batch = 10, chunk = 7), ref)
})

# --- adaptive (Besag-Clifford) stopping ---------------------------------------------------------------

# Besag-Clifford computed offline from complete sequences of nulls (dev/engine-spec.md 6.3).
bc_offline <- function(nullX, nullY, rX, rY, h, n_max) {
  cx <- cumsum(abs(nullX) >= abs(rX))
  cy <- cumsum(abs(nullY) >= abs(rY))
  hit <- which(cx >= h | cy >= h)
  stopped <- length(hit) > 0
  L <- if (stopped) hit[1] else n_max
  alone <- function(cs) if (any(cs >= h)) h / which(cs >= h)[1] else (cs[n_max] + 1) / (n_max + 1)
  list(L = L, stop = if (stopped) "exceedances" else "cap", bX = cx[L], bY = cy[L],
       p = if (stopped) max(cx[L], cy[L]) / L else (max(cx[L], cy[L]) + 1) / (n_max + 1),
       p_max_alone = max(alone(cx), alone(cy)))
}

test_that("portable: adaptive stopping with h = Inf equals the fixed-B run", {
  a <- aki_input(c("Gpx1", "Upk2", "Mrps6"))
  fixed <- strip(eng(a$X, a$Y, a$pos, B = 30, deltaX = a$delta, deltaY = a$delta))
  for (fb in c(1L, 7L, 64L)) {
    ad <- strip(.stc_engine_correlate(a$X, a$Y, a$pos, deltaX = a$delta, deltaY = a$delta,
                                      adaptive = list(h = Inf, n_max = 30, first_batch = fb, growth = 2)))
    expect_identical(ad, fixed, info = fb)
  }
})

test_that("portable: adaptive stopping equals a brute-force run with offline Besag-Clifford, for any schedule and threads", {
  # with h = 2 these genes stop at the cap (Gpx1), at once (Atp5a1), or late in X (Ndufb4, Agt, Pin1) or in Y
  # (Slc6a20b)
  a <- aki_input(c("Gpx1", "Atp5a1", "Ndufb4", "Agt", "Slc6a20b", "Pin1"))
  h <- 2
  n_max <- 40L
  brute <- strip(eng(a$X, a$Y, a$pos, B = n_max, deltaX = a$delta, deltaY = a$delta, returnPermutations = TRUE))
  schedules <- list(list(first_batch = 8, growth = 2, nThreads = 1L), list(first_batch = 1, growth = 1, nThreads = 2L),
                    list(first_batch = 64, growth = 2, nThreads = 2L))
  runs <- lapply(schedules, function(s) {
    strip(.stc_engine_correlate(a$X, a$Y, a$pos, deltaX = a$delta, deltaY = a$delta, returnPermutations = TRUE,
                                adaptive = list(h = h, n_max = n_max, first_batch = s$first_batch, growth = s$growth),
                                nThreads = s$nThreads, chunk = 5L))
  })
  for (k in seq_along(runs)[-1]) expect_identical(runs[[k]], runs[[1]], info = k)
  ad <- runs[[1]]
  expect_identical(ad$stop, c("cap", rep("exceedances", 5)))  # the gene set covers both outcomes
  expect_true(any(ad$bX >= h & ad$bY < h) && any(ad$bY >= h & ad$bX < h) && any(ad$L > 30 & ad$stop == "exceedances"))
  for (g in rownames(ad)) {
    rX <- stats::cor(a$X[, g], a$Y[, g])  # the observed r of each direction, as the method computes it
    rY <- stats::cor(a$Y[, g], a$X[, g])
    o <- bc_offline(brute[g, "nullX"][[1]], brute[g, "nullY"][[1]], rX, rY, h, n_max)
    L <- o$L
    expect_identical(ad[g, "L"], as.integer(L), info = g)
    expect_identical(ad[g, "stop"], o$stop, info = g)
    expect_identical(c(ad[g, "bX"], ad[g, "bY"]), as.integer(c(o$bX, o$bY)), info = g)
    expect_identical(ad[g, "p"], o$p, info = g)
    expect_identical(ad[g, "p"], o$p_max_alone, info = g)  # combined p = max of the standalone Besag-Clifford p-values
    for (col in c("nullX", "nullY", "deltaStarX", "deltaStarY")) {
      expect_identical(ad[g, col][[1]], brute[g, col][[1]][seq_len(L)], info = paste(g, col))
    }
    expect_identical(ad[g, "permutationsX"][[1]], brute[g, "permutationsX"][[1]][, seq_len(L), drop = FALSE], info = g)
    # the counts are the exceedances of the truncated nulls
    expect_identical(ad[g, "bX"], sum(abs(ad[g, "nullX"][[1]]) >= abs(rX)), info = g)
    expect_identical(ad[g, "bY"], sum(abs(ad[g, "nullY"][[1]]) >= abs(rY)), info = g)
  }
})

# --- exceedances and p-values --------------------------------------------------------------------------

test_that("portable: a null equal to |r| counts as an exceedance (|null| >= |r|), and p = (b + 1) / (B + 1)", {
  cs <- fx$cases$quakes_irregular$input
  e <- eng(cs$X, cs$Y, cs$pos, B = 6, deltaX = c(0.2, 0.6), deltaY = 0.4)
  st <- attr(e, "state")
  nullX <- e$nullX[[1]]
  nullY <- e$nullY[[1]]
  expect_identical(e$pX, (sum(abs(nullX) >= abs(e$r)) + 1) / (6 + 1))
  # a second session with the same plan, smoothers and streams, whose observed |r| are the absolute values
  # of the third null of X and the fifth null of Y: those nulls tie with |r| exactly
  s2 <- .stc_engine_new(cs$pos[, 2], cs$pos[, 1], as.integer(st$stream$ids) - 1L, st$plan, st$deltas,
                        as.integer(.stc_cor_mode()), 0, FALSE, 1L, 16L, FALSE)
  .stc_engine_define(s2, cbind(cs$X, cs$Y), c(0L, 1L), list(match(c(0.2, 0.6), st$deltas) - 1L, match(0.4, st$deltas) - 1L),
                     list(1L, 0L), list(abs(nullX[3]), abs(nullY[5])), c(0L, 0L), TRUE)
  .stc_engine_run(s2, 0L, 1L, 6L, .stc_stream_perms(st$stream, 1, 6), NULL, Inf, 6L)
  expect_identical(.stc_engine_task_results(s2, 0:1, FALSE)[[1]]$nulls[, 1], nullX)
  counts <- .stc_engine_units(s2)$counts[[1]]
  expect_identical(counts, c(sum(abs(nullX) >= abs(nullX[3])), sum(abs(nullY) >= abs(nullY[5]))))
  expect_identical(counts - 1L, c(sum(abs(nullX) > abs(nullX[3])), sum(abs(nullY) > abs(nullY[5]))))
})

# --- the global RNG state --------------------------------------------------------------------------

test_that("portable: the engine leaves the global RNG state as it was (no seed, a seed, non-default kinds)", {
  cs <- fx$cases$quakes_irregular$input
  run <- function() invisible(eng(cs$X, cs$Y, cs$pos, B = 3, deltaX = c(0.2, 0.6)))
  local_default_rng()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  run()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(123)
  before <- .Random.seed
  run()
  expect_identical(.Random.seed, before)
  suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rounding"))
  set.seed(77)
  before <- .Random.seed
  kinds <- RNGkind()
  run()
  expect_identical(RNGkind(), kinds)
  expect_identical(.Random.seed, before)
  # after an error too
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(5)
  before <- .Random.seed
  expect_error(eng(cs$X, cs$Y[-1], cs$pos, B = 3), "same dimensions")
  expect_error(.stc_engine_correlate(state = attr(eng(cs$X, cs$Y, cs$pos, B = 2), "state"), nPermutations = 3,
                                     adaptive = list(h = 0, n_max = 3)), "at least 1")
  expect_identical(.Random.seed, before)
})

test_that("portable: every exported function leaves the global RNG state unchanged", {
  q <- fx$cases$quakes_irregular$input
  m <- rbind(depth = q$X, mag = q$Y, mix = q$X / 100 + q$Y)
  colnames(m) <- sprintf("p%03d", seq_along(q$X))
  pos <- q$pos
  rownames(pos) <- colnames(m)
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = m), spatialCoords = pos)
  d2 <- rep(list(c(0.2, 0.6)), 3)
  calls <- list(
    viladomatCorrelation = function() viladomatCorrelation(cbind(q$X, q$Y, q$pos), c(0.2, 0.6), 0.25, 2),
    spatialCorrelation = function() spatialCorrelation(q$X, q$Y, q$pos, nPermutations = 2, deltaX = 0.3, deltaY = 0.3),
    spatialCorrelationGeneExp = function() {
      spatialCorrelationGeneExp(list(spe, spe[3:1, ]), nPermutations = 2, deltaX = d2, deltaY = d2, verbose = FALSE)
    },
    spatialCorrelationGeneExpIterPermutations = function() {
      spatialCorrelationGeneExpIterPermutations(list(spe, spe[3:1, ]), alpha = 1, nPermutations = c(2, 3),
                                                deltaX = d2, deltaY = d2, verbose = FALSE)
    },
    spatialCorrelationGeneExpWithinSample = function() {
      spatialCorrelationGeneExpWithinSample(spe, nPermutations = 2, delta = d2, verbose = FALSE)
    })
  local_default_rng()
  for (nm in names(calls)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    invisible(calls[[nm]]())
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE), info = nm)
    set.seed(42)
    before <- .Random.seed
    invisible(calls[[nm]]())
    expect_identical(.Random.seed, before, info = nm)
  }
})

test_that("portable: the noise does not depend on the session's normal.kind; the permutations follow its RNG kinds", {
  # The noise is always the L'Ecuyer-CMRG / Inversion stream (as inside the BiocParallel tasks of the original
  # code); the permutations are sample() draws in the session, under its kinds.
  cs <- fx$cases$quakes_irregular$input
  N <- length(cs$X)
  run <- function() strip(eng(cs$X, cs$Y, cs$pos, B = 3, deltaX = c(0.2, 0.6)))
  local_default_rng()
  ref <- run()
  for (nk in c("Box-Muller", "Kinderman-Ramage", "Ahrens-Dieter")) {
    suppressWarnings(RNGkind("Mersenne-Twister", nk, "Rejection"))
    expect_identical(run(), ref, info = nk)
    expect_identical(RNGkind()[2], nk)
  }
  for (k in list(c("Mersenne-Twister", "Inversion", "Rounding"), c("Knuth-TAOCP-2002", "Inversion", "Rejection"))) {
    suppressWarnings(RNGkind(k[1], k[2], k[3]))
    other <- run()
    expect_false(identical(other$nullX, ref$nullX), info = k[1])
    s <- .stc_legacy_stream(N, 0)
    perm <- .stc_stream_perms(s, 1, 3)
    set.seed(0)
    expect_identical(perm, vapply(1:3, function(b) sample.int(N, N), integer(N)), info = paste(k, collapse = "/"))
    expect_identical(RNGkind(), k)
  }
})

test_that("portable: noise drawn in R gives the same results as the C++ noise", {
  a <- aki_input(c("Gpx1", "Upk2"))
  e1 <- strip(eng(a$X, a$Y, a$pos, B = 12, deltaX = a$delta, deltaY = a$delta[3:6], noise = "cpp"))
  e2 <- strip(eng(a$X, a$Y, a$pos, B = 12, deltaX = a$delta, deltaY = a$delta[3:6], noise = "R", batch = 5))
  expect_identical(e2, e1)
  # the stream object: batches drawn in any order equal one draw, and the caller's RNG state is untouched
  local_default_rng()
  set.seed(1)
  before <- .Random.seed
  s <- .stc_legacy_stream(1200L, 3)
  p1 <- .stc_stream_perms(s, 1, 10)
  p2 <- cbind(.stc_stream_perms(s, 1, 4), .stc_stream_perms(s, 5, 10))
  p3 <- .stc_stream_perms(s, 6, 10)
  expect_identical(p2, p1)
  expect_identical(p3, p1[, 6:10])
  expect_identical(.Random.seed, before)
  set.seed(3)
  expect_identical(s$ids, sample(1200, 1000))
  expect_identical(p1[, 1], sample.int(1200, 1200))
})

# --- continuation (iterative rounds) ----------------------------------------------------------------

test_that("portable: continuing a session for some genes equals a single longer run (prefix property)", {
  a <- aki_input(c("Gpx1", "Ech1", "Upk2", "Tspan8"))
  run <- function(B, ...) eng(a$X, a$Y, a$pos, B = B, deltaX = a$delta, deltaY = a$delta, returnPermutations = TRUE, ...)
  r7 <- run(7)
  st <- attr(r7, "state")
  ext <- strip(.stc_engine_correlate(state = st, units = c("Gpx1", "Upk2"), nPermutations = 19))
  expect_identical(ext, strip(run(19))[c("Gpx1", "Upk2"), ])
  # a unit left behind is extended later, from an earlier permutation than the stream's position
  ext2 <- strip(.stc_engine_correlate(state = st, units = 2L, nPermutations = 12, nThreads = 2))
  expect_identical(ext2, strip(run(12))["Ech1", ])
  # the untouched unit still has its 7 permutations
  r <- strip(.stc_engine_correlate(state = st, units = 4L, nPermutations = 7))
  expect_identical(r, strip(r7)["Tspan8", ])
  # an adaptive run extended with a fixed number of permutations equals the fixed run
  ad <- eng(a$X, a$Y, a$pos, deltaX = a$delta, deltaY = a$delta, returnPermutations = TRUE,
            adaptive = list(h = 2, n_max = 15, first_batch = 4))
  ext3 <- strip(.stc_engine_correlate(state = attr(ad, "state"), nPermutations = 21))
  expect_identical(ext3, strip(run(21)))
})

test_that("portable: spatialCorrelationGeneExpIterPermutations() extends the promoted genes (= fresh runs at the larger B), never promotes NA rows, and adjusts at the end", {
  P <- rf$pairs$aki
  genes <- c("Gpx1", "Ech1", "Upk2", "Tspan8", "Rbp4", "Igf1", "Dcn")
  P$Y["Upk2", ] <- 7  # constant: an NA row (the R implementation crashed on it in round 2)
  input <- fx_spe_pair(P, genes)
  delta <- rep(list(P$params$delta[2:8]), length(genes))
  run <- function(f, ...) {
    w <- fx_collect_warnings(f(input, deltaX = delta, deltaY = delta, verbose = FALSE, nThreads = fx_threads, ...))
    expect_length(w$warnings, 1L)
    expect_match(w$warnings, "gene Upk2: no permutation p-values")
    w$value
  }
  it <- run(spatialCorrelationGeneExpIterPermutations, nPermutations = c(25, 10), alpha = 0.04, adjustMethod = "none",
            returnPermutations = TRUE)
  f10 <- run(spatialCorrelationGeneExp, nPermutations = 10, adjustMethod = "none", returnPermutations = TRUE)
  f25 <- run(spatialCorrelationGeneExp, nPermutations = 25, adjustMethod = "none", returnPermutations = TRUE)
  # screening after the 10-permutation round (nPermutations is sorted): both p < 100 * 0.04 / 10
  promoted <- !is.na(f10$pValuePermuteX) & f10$pValuePermuteX < 0.4 & f10$pValuePermuteY < 0.4
  expect_true(any(promoted))
  expect_true(any(!promoted & !is.na(f10$pValuePermuteX)))
  expect_identical(it[promoted, ], f25[promoted, ])
  expect_identical(it[!promoted, ], f10[!promoted, ])
  expect_identical(lengths(it$nullCorrelationsX), ifelse(promoted, 25L, ifelse(is.na(f10$pValuePermuteX), 1L, 10L)))
  # the final adjustment is applied across all genes, to each column separately
  bh <- run(spatialCorrelationGeneExpIterPermutations, nPermutations = c(10, 25), alpha = 0.04)
  expect_identical(bh$pValuePermuteX, stats::p.adjust(it$pValuePermuteX, method = "BH"))
  expect_identical(bh$pValuePermuteY, stats::p.adjust(it$pValuePermuteY, method = "BH"))
  expect_identical(bh$nullCorrelationsY, it$nullCorrelationsY)
  # a single round equals spatialCorrelationGeneExp()
  one <- run(spatialCorrelationGeneExpIterPermutations, nPermutations = 10)
  expect_identical(one, run(spatialCorrelationGeneExp, nPermutations = 10))
})

# --- within-sample mode ------------------------------------------------------------------------------

test_that("portable: spatialCorrelationGeneExpWithinSample() equals spatialCorrelation() on every pair; pairs with a constant gene are NA rows", {
  q <- fx$cases$quakes_irregular$input
  local_default_rng()
  set.seed(31)
  m <- rbind(depth = q$X, mag = q$Y, mix = q$X / 100 + q$Y + stats::rnorm(length(q$X), sd = 0.3), flat = 2)
  colnames(m) <- sprintf("p%03d", seq_along(q$X))
  pos <- q$pos
  rownames(pos) <- colnames(m)
  delta <- list(c(0.1, 0.5), c(0.2, 0.6, 0.9), c(0.15, 0.45), 0.3)
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = m), spatialCoords = pos)
  w <- fx_collect_warnings(spatialCorrelationGeneExpWithinSample(spe, nPermutations = 4, delta = delta, verbose = FALSE,
                                                                 nThreads = 2, returnPermutations = TRUE))
  o <- w$value
  pairs <- utils::combn(4, 2)
  expect_identical(o$first, rownames(m)[pairs[1, ]])
  expect_identical(o$second, rownames(m)[pairs[2, ]])
  expect_identical(rownames(o), c("cor", paste0("cor", 1:5)))
  expect_identical(names(o)[13:14], c("first", "second"))
  expect_length(w$warnings, 3L)  # every pair with the constant gene
  expect_true(all(grepl("^spatialCorrelationGeneExpWithinSample: genes [a-z]+ and flat: no permutation p-values \\(NA row\\): flat: the permuted values are constant",
                        w$warnings)))
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]
    j <- pairs[2, k]
    s <- suppressWarnings(spatialCorrelation(m[i, ], m[j, ], q$pos, nPermutations = 4, deltaX = delta[[i]],
                                             deltaY = delta[[j]], returnPermutations = TRUE))
    rownames(s) <- rownames(o)[k]
    expect_identical(o[k, 1:12], s, info = paste(i, j))
  }
  expect_error(.stc_engine_correlate(t(m), pos = q$pos, mode = "within", adaptive = list(h = 3, n_max = 9)),
               "only available for gene-wise")
})

test_that("portable: the naive r and p-values are cor.test()'s for every pair (finite, constant, missing and infinite values)", {
  # They are computed for all pairs at once (not by cor.test() itself) when both genes are finite; this
  # pins them to cor.test(), first for .stc_cor_tests() on many pairs, then through the within-sample function.
  local_default_rng()
  set.seed(32)
  f <- stats::rnorm(50)
  A <- sapply(seq_len(24), function(j) f * (j - 12) / 6 + stats::rnorm(50, sd = j / 8))
  A <- cbind(A, 3, A[, 1] * 2, -A[, 2] / 7, replace(A[, 3], 4, NA), replace(A[, 4], 9, -Inf))
  pairs <- utils::combn(ncol(A), 2L)
  for (r_given in c(FALSE, TRUE)) {
    rp <- if (r_given) suppressWarnings(stats::cor(A, A))[t(pairs)] else NULL
    nt <- .stc_cor_tests(A, A, pairs[1L, ], pairs[2L, ], r_pairs = rp)
    ref <- lapply(seq_len(ncol(pairs)), function(k) suppressWarnings(stats::cor.test(A[, pairs[1L, k]], A[, pairs[2L, k]])))
    expect_identical(nt$r, vapply(ref, function(ct) unname(ct$estimate), 0))
    expect_identical(nt$p, vapply(ref, function(ct) ct$p.value, 0))
    expect_identical(nt$msg, rep("", ncol(pairs)))
  }
  few <- .stc_cor_tests(A[1:2, 1:3], A[1:2, 1:3], c(1L, 1L), c(2L, 3L))
  expect_identical(few$msg, rep("cor.test() failed: not enough finite observations", 2L))
  expect_identical(few$p, c(NA_real_, NA_real_))

  q <- fx$cases$quakes_irregular$input
  m <- rbind(depth = q$X, mag = q$Y, mix = q$X / 100 + q$Y + stats::rnorm(length(q$X), sd = 0.3), flat = 2,
             twice = 2 * q$X, gap = replace(q$Y, 3, NA), inf = replace(q$X, 5, Inf))
  colnames(m) <- sprintf("p%03d", seq_along(q$X))
  pos <- q$pos
  rownames(pos) <- colnames(m)
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = m), spatialCoords = pos)
  o <- suppressWarnings(spatialCorrelationGeneExpWithinSample(spe, nPermutations = 2, delta = rep(list(0.3), nrow(m)),
                                                              verbose = FALSE))
  pairs <- utils::combn(nrow(m), 2L)
  for (k in seq_len(ncol(pairs))) {
    ct <- suppressWarnings(stats::cor.test(m[pairs[1L, k], ], m[pairs[2L, k], ]))
    expect_identical(o$correlationCoef[k], unname(ct$estimate), info = k)
    expect_identical(o$pValueNaive[k], ct$p.value, info = k)
  }
})

# --- verbose -------------------------------------------------------------------------------------

test_that("portable: verbose = TRUE gives one start and one end message", {
  P <- rf$pairs$aki
  genes <- c("Gpx1", "Upk2")
  input <- fx_spe_pair(P, genes)
  d <- rep(list(c(0.2, 0.5)), 2)
  m1 <- collect_messages(spatialCorrelationGeneExp(input, nPermutations = 3, deltaX = d, deltaY = d))$messages
  m2 <- collect_messages(spatialCorrelationGeneExpIterPermutations(input, nPermutations = c(3, 6), alpha = 1,
                                                                   deltaX = d, deltaY = d))$messages
  m3 <- collect_messages(spatialCorrelationGeneExpWithinSample(input[[1]], nPermutations = 3, delta = d))$messages
  for (m in list(m1, m2, m3)) {
    expect_length(m, 2L)
    expect_match(m[1], "thread")
    expect_match(m[2], "done in")
  }
  expect_match(m1[1], "^spatialCorrelationGeneExp: 2 gene\\(s\\) x 2 directions on 311 shared pixels, 3 permutations, 1 thread")
  expect_match(m2[2], "genes tested per round: 2 at B = 3, 2 at B = 6")
  expect_length(collect_messages(spatialCorrelationGeneExp(input, nPermutations = 3, deltaX = d, deltaY = d,
                                                           verbose = FALSE))$messages, 0L)
})

# --- threads, interrupts and errors -------------------------------------------------------------------

test_that("portable: .stc_threads() takes the workers of BPPARAM, otherwise nThreads", {
  expect_identical(.stc_threads(3), 3L)
  expect_identical(.stc_threads(2, BiocParallel::SerialParam()), 1L)
  # (at most 2 workers: BiocParallel refuses more under R CMD check's _R_CHECK_LIMIT_CORES_)
  expect_identical(.stc_threads(1, BiocParallel::SnowParam(workers = 2)), 2L)
  expect_error(.stc_threads(0), "positive integer")
  expect_error(.stc_threads(NA), "positive integer")
})

test_that("portable: the thread pool runs every item, and worker exceptions become R errors after the join", {
  w <- .stc_parallel_selftest(1000L, 2L)
  expect_length(w, 1000L)
  expect_true(all(w %in% 0:1))
  expect_identical(.stc_parallel_selftest(0L, 2L), integer(0))
  expect_identical(.stc_parallel_selftest(5L, 1L), rep(0L, 5))
  expect_error(.stc_parallel_selftest(100L, 2L, fail_item = 37L), "item 37 failed on a worker")
  expect_error(.stc_parallel_selftest(100L, 1L, fail_item = 0L), "item 0 failed on a worker")
  # the pool still works afterwards
  expect_length(.stc_parallel_selftest(10L, 2L), 10L)
})

test_that("portable: a user interrupt stops the workers and returns to R, which stays usable", {
  skip_on_cran()
  skip_on_os("windows")
  # A forked child inherits R's SIGINT handler; it waits in the pool (400 items of 50 ms on 2 threads, 10 s)
  # until the signal, which the main thread sees within about 100 ms.
  job <- parallel::mcparallel({
    t0 <- Sys.time()
    r <- tryCatch({
      .stc_parallel_selftest(400L, 2L, -1L, 50L)
      "finished"
    }, interrupt = function(e) "interrupted")
    list(r = r, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
         after = length(.stc_parallel_selftest(20L, 2L)))
  })
  Sys.sleep(1.5)
  tools::pskill(job$pid, tools::SIGINT)
  res <- parallel::mccollect(job, wait = TRUE, timeout = 60)[[1]]
  expect_identical(res$r, "interrupted")
  expect_lt(res$seconds, 8)
  expect_identical(res$after, 20L)
})

test_that("portable: invalid inputs and sessions give R errors", {
  cs <- fx$cases$kidney_AB$input
  expect_error(eng(cs$X, cs$Y, cs$pos[-1, ], B = 2), "one row per pixel")
  expect_error(eng(cs$X, cs$Y, replace(cs$pos, 3, NA), B = 2), "finite")
  expect_error(eng(cs$X, cs$Y, cs$pos, B = 0), "positive integer")
  expect_error(eng(cs$X, cs$Y, cs$pos, B = 2, seed = NA), "seed")
  expect_error(eng(cs$X, cs$Y, cs$pos, B = 2, seed = .Machine$integer.max), "seed")
  expect_error(eng(cs$X, cs$Y, cs$pos, B = 2, deltaX = list(0.1, 0.2)), "one element per gene")
  expect_error(eng(cs$X, cs$Y, cs$pos, B = 2, chunk = 0), "positive integer")
  e <- eng(cs$X, cs$Y, cs$pos, B = 2)
  st <- attr(e, "state")
  expect_error(.stc_engine_correlate(state = st, units = "nope", nPermutations = 3), "units")
  expect_error(.stc_engine_run(st$session, 0L, 5L, 5L, matrix(seq_along(cs$X), ncol = 1), NULL, Inf, 5L),
               "must start at permutation 3")
  expect_error(.stc_engine_run(st$session, 0L, 3L, 3L, matrix(0L, length(cs$X), 1), NULL, Inf, 3L), "out of range")
  expect_error(.stc_engine_units(e), "not an engine session")
  expect_error(.stc_engine_define(st$session, matrix(0, length(cs$X), 1), 0L, list(0L), list(0L), list(0), 0L, TRUE),
               "already defined")
  # results are unchanged by the failed calls
  expect_identical(strip(.stc_engine_correlate(state = st, nPermutations = 2)), strip(e))
})

test_that("portable: invalid arguments of the exported functions give R errors", {
  cs <- fx$cases$kidney_AB$input
  expect_error(spatialCorrelation(cs$X, cs$Y[-1], cs$pos), "same length")
  expect_error(spatialCorrelation(cs$X[1:2], cs$Y[1:2], cs$pos[1:2, ]), "at least 3")
  expect_error(spatialCorrelation(cs$X, cs$Y, cs$pos, nPermutations = 0), "positive integer")
  expect_error(spatialCorrelation(cs$X, cs$Y, cs$pos, nPermutations = 2.5), "positive integer")
  expect_error(spatialCorrelation(cs$X, cs$Y, cs$pos, nThreads = 0), "positive integer")
  expect_error(spatialCorrelation(as.character(cs$X), cs$Y, cs$pos), "numeric")
  expect_error(viladomatCorrelation(cbind(cs$X, cs$Y), 0.3, 0.25, 2), "4 columns")
  rk <- fx_speKidney_raster()
  expect_error(spatialCorrelationGeneExp(list(rk$A, rk$B), deltaX = list(0.1, 0.2), verbose = FALSE), "one element")
  expect_error(spatialCorrelationGeneExp(list(rk$A), verbose = FALSE), "list of two")
  expect_error(spatialCorrelationGeneExpIterPermutations(list(rk$A, rk$B), nPermutations = c(10, -1), verbose = FALSE),
               "must be positive")
  expect_error(spatialCorrelationGeneExpWithinSample(rk$A, verbose = FALSE), "at least 2 rows")
})
