# Tests of the rank-remapped surrogates, compareSpatial(surrogate = "remap") (the default) and
# .stc_engine_correlate(surrogate = "remap") (dev/surrogate-remap-spec.md, section 5; EngineTask::remap in
# src/stc_engine.h): every remapped surrogate is a rearrangement of the permuted gene's values in the rank order
# of the gaussian surrogate; for genes without spatial structure the test is the exact permutation test, ties
# included; the results do not depend on threads, batches or chunks; minDetected = NULL means 0 pixels with
# remapped surrogates and sqrt(N) with gaussian ones, and minDetected = 0 skips only constant genes; and
# surrogate = "gaussian" is the behaviour of the version before the surrogate argument existed. Labels are
# explained in helper-fixtures.R.

rf <- fx_read("realistic_fixture.rds")

# compareSpatial() without progress output or messages, and its table without the attributes.
rm_run <- function(input, ..., nThreads = 1L) {
  compareSpatial(input, ..., nThreads = nThreads, progress = FALSE, verbose = FALSE)
}
rm_table <- function(res) as.data.frame(res)

# The grid of the compareSpatial() tests: with 60 permutations and exceedances = 3, Upk2, Agt and Calm1 stop
# after 3, 40 and 47 permutations with gaussian surrogates (test-compareSpatial.R).
rm_grid <- c(0.05, 0.1, 0.3, 0.6, 0.9)

# A pixels x 1 matrix with a gene name (the independent streams key the gene name).
rm_col <- function(v, g = "g") matrix(unname(v), dimnames = list(NULL, g))

# The exceedance rule of remapped tasks: |null| >= |r| within 1e-9 relative (src/stc_engine.h, kRemapTieRel).
rm_exceedances <- function(null, r) sum(abs(null) >= abs(r) * (1 - 1e-9))

test_that("portable: every remapped surrogate has the values of the permuted gene, in the rank order of the gaussian surrogate", {
  P <- rf$pairs$aki
  genes <- c("Upk2", "Gpx1")  # detected in 19 and in all 311 pixels of x
  X <- t(P$X[genes, ])
  Y <- t(P$Y[genes, ])
  grid <- c(0.05, 0.2, 0.5, 0.8)
  B <- 8L
  run <- function(surrogate) {
    .stc_engine_correlate(X, Y, P$pos, deltaX = grid, deltaY = rev(grid), nPermutations = B, seed = 5,
                          streams = "independent", returnPermutations = TRUE, nThreads = 2L, surrogate = surrogate)
  }
  g <- run("gaussian")
  r <- run("remap")
  expect_identical(r$status, c("ok", "ok"))
  for (k in 1:2) {
    for (d in c("X", "Y")) {
      info <- paste(genes[k], d)
      src <- unname((if (d == "X") X else Y)[, k])
      tgt <- unname((if (d == "X") Y else X)[, k])
      Sg <- g[[paste0("permutations", d)]][[k]]
      Sr <- r[[paste0("permutations", d)]][[k]]
      expect_identical(dim(Sr), c(length(src), B), info = info)
      for (b in seq_len(B)) {
        expect_identical(sort(Sr[, b]), sort(src), info = info)          # the gene's values, zeros included
        expect_identical(Sr[order(Sg[, b]), b], sort(src), info = info)  # in the rank order of the gaussian surrogate
        expect_false(identical(sort(Sg[, b]), sort(src)), info = info)
      }
      # delta* is chosen before the remapping; the nulls are the correlations of the remapped surrogates
      expect_identical(r[[paste0("deltaStar", d)]][[k]], g[[paste0("deltaStar", d)]][[k]], info = info)
      expect_equal(r[[paste0("null", d)]][[k]], as.vector(stats::cor(Sr, tgt)), tolerance = 1e-12, info = info)
      expect_false(identical(r[[paste0("null", d)]][[k]], g[[paste0("null", d)]][[k]]), info = info)
    }
  }
  # the same with the legacy streams: the remapping does not depend on where the draws come from
  e <- .stc_engine_correlate(X, Y, P$pos, deltaX = grid, deltaY = grid, nPermutations = 3L, surrogate = "remap",
                             returnPermutations = TRUE)
  expect_identical(e$status, c("ok", "ok"))
  expect_identical(apply(e$permutationsY[[1]], 2L, sort), matrix(sort(unname(Y[, 1])), nrow(Y), 3L))
})

test_that("portable: for genes without spatial structure the remapped null is the plain permutation null, ties included", {
  P <- rf$pairs$aki
  N <- nrow(P$pos)
  local_default_rng()
  set.seed(11)
  y <- stats::rnorm(N)   # genes without spatial structure
  x <- stats::rnorm(N)
  forward <- function(x, y, B, seed) {
    .stc_engine_correlate(rm_col(x), rm_col(y), P$pos, deltaX = rm_grid, nPermutations = B, seed = seed,
                          mode = "forward", surrogate = "remap", nThreads = 2L)
  }
  # the remapped nulls of 300 permutations against the plain permutation null (Kolmogorov-Smirnov)
  plain <- replicate(20000, stats::cor(sample(x), y))
  expect_gt(stats::ks.test(forward(x, y, 300L, 1)$nullX[[1]], plain)$p.value, 0.05)
  # a gene detected in one pixel: every remapped surrogate puts its value at one pixel, so the test is the exact
  # permutation test, whose p-value is the share of pixels at least as far from mean(y) as the detected pixel
  B <- 1500L
  for (i0 in c(50L, 250L)) {
    spike <- numeric(N)
    spike[i0] <- 7
    p_exact <- mean(abs(y - mean(y)) >= abs(y[i0] - mean(y)))
    e <- forward(spike, y, B, i0)
    expect_true(abs(e$pX - p_exact) < 3 * sqrt(p_exact * (1 - p_exact) / B) + 1 / B,
                info = sprintf("pixel %d: exact p %.4f, pX %.4f", i0, p_exact, e$pX))
  }
  # against a sparse y most rearrangements give exactly the observed correlation (the one-pixel gene lands on a
  # zero of y): the exact p-value is close to 1, and the ties count although cor() rounds each of them differently
  ys <- numeric(N)
  ys[sample(N, 40)] <- stats::rexp(40) * 10
  spike <- numeric(N)
  spike[50] <- 7
  p_exact <- mean(abs(ys - mean(ys)) >= abs(ys[50] - mean(ys)))
  expect_gt(p_exact, 0.9)
  e <- forward(spike, ys, 500L, 3)
  null <- e$nullX[[1]]
  expect_gt(mean(abs(abs(null) - abs(e$r)) <= 1e-9 * abs(e$r)), 0.5)  # most nulls tie with r within rounding
  expect_true(any(abs(null) < abs(e$r)))                                   # and some of them round below it
  expect_identical(e$bX, rm_exceedances(null, e$r))
  expect_lt(abs(e$pX - p_exact), 3 * sqrt(p_exact * (1 - p_exact) / 500) + 1 / 500)
})

test_that("portable: in remap mode the results do not depend on threads, batches or chunks, and the global RNG state is unchanged", {
  genes <- c("Gpx1", "Upk2", "Agt", "Calm1")
  input <- fx_spe_pair(rf$pairs$aki, genes)
  run <- function(exceedances = 3, ...) {
    rm_table(rm_run(input, tests = "correlation", nPermutations = 60, exceedances = exceedances, delta = rm_grid,
                    surrogate = "remap", ...))
  }
  local_default_rng()
  set.seed(42)
  before <- .Random.seed
  ref <- run()
  expect_identical(.Random.seed, before)
  expect_true(all(ref$status == "ok"))
  expect_setequal(unique(ref$stop), c("exceedances", "limit"))
  expect_identical(run(nThreads = 2L), ref)
  for (s in list(list(first_batch = 1, growth = 1, max_batch = 7, chunk = 3),
                 list(first_batch = 5, growth = 3, max_batch = 1000, chunk = 1))) {
    old <- options(STcompare.compare_schedule = s)
    expect_identical(run(nThreads = 2L), ref, info = paste(names(s), unlist(s), collapse = " "))
    options(old)
  }
  fixed <- run(exceedances = Inf)
  old <- options(STcompare.compare_schedule = list(max_batch = 13, chunk = 4))
  expect_identical(run(exceedances = Inf, nThreads = 2L), fixed)
  options(old)
  # the stored nulls are the remapped ones (they differ from the gaussian nulls), and the table follows from them
  # with the tie rule; the result records and prints the mode
  kept <- rm_run(input, tests = "correlation", nPermutations = 60, exceedances = 3, delta = rm_grid, surrogate = "remap",
                 keepNulls = TRUE)
  gauss <- rm_run(input, tests = "correlation", nPermutations = 60, exceedances = 3, delta = rm_grid, keepNulls = TRUE,
                  surrogate = "gaussian")
  expect_identical(rm_table(kept), ref)
  expect_identical(attr(kept, "params")$surrogate, "remap")
  expect_identical(attr(gauss, "params")$surrogate, "gaussian")
  det <- attr(kept, "details")
  for (g in genes) {
    L <- kept[g, "nPermutations"]
    expect_false(identical(det[[g]]$nullX, attr(gauss, "details")[[g]]$nullX[seq_len(L)]), info = g)
    bX <- rm_exceedances(det[[g]]$nullX, kept[g, "r"])
    bY <- rm_exceedances(det[[g]]$nullY, kept[g, "r"])
    expect_identical(kept[g, "pX"], (bX + 1) / (L + 1), info = g)
    expect_identical(kept[g, "p"], if (kept[g, "stop"] == "exceedances") max(bX, bY) / L else (max(bX, bY) + 1) / (L + 1),
                     info = g)
  }
  expect_output(print(kept), "rank-remapped surrogates")
  expect_output(print(summary(kept)), "rank-remapped surrogates")
  expect_error(rm_run(input, surrogate = "aaft"), "should be one of")
})

test_that("portable: minDetected = NULL means 0 pixels with remapped surrogates and sqrt(N) with gaussian ones; 0 tests every gene that is not constant", {
  P <- rf$pairs$aki
  genes <- c("Gpx1", "Ech1", "Upk2")
  P$X["Ech1", ] <- 3                 # constant in x: skipped whatever minDetected is
  P$X["Upk2", ] <- 0                 # detected in one pixel of x: tested unless minDetected asks for more
  P$X["Upk2", 100] <- 5
  input <- fx_spe_pair(P, genes)
  run <- function(...) rm_run(input, tests = "correlation", nPermutations = 20, exceedances = Inf, delta = rm_grid, ...)
  for (mode in c("gaussian", "remap")) {
    res <- run(minDetected = 0, surrogate = mode)
    expect_identical(res$status, c("ok", "skipped", "ok"), info = mode)
    expect_identical(res$message[2], "constant in x (zero variance): not tested", info = mode)
    expect_identical(attr(res, "params")$minDetectedPixels, 0L, info = mode)
    expect_false(anyNA(res$p[c(1, 3)]), info = mode)
    expect_identical(res$nPermutations[c(1, 3)], c(20L, 20L), info = mode)
  }
  # the default of compareSpatial() is remap with no filter: the same as minDetected = 0 above
  default <- run()
  prm <- attr(default, "params")
  expect_identical(prm$surrogate, "remap")
  expect_null(prm$minDetected)
  expect_identical(prm$minDetectedPixels, 0L)
  expect_identical(rm_table(default), rm_table(run(minDetected = 0, surrogate = "remap")))
  expect_identical(default$status, c("ok", "skipped", "ok"))
  expect_output(print(default), "rank-remapped surrogates")
  # with gaussian surrogates the default filter is ceiling(sqrt(311)) = 18 pixels
  gauss <- run(surrogate = "gaussian")
  prm <- attr(gauss, "params")
  expect_identical(prm$surrogate, "gaussian")
  expect_null(prm$minDetected)
  expect_identical(prm$minDetectedPixels, 18L)
  expect_identical(gauss$status, c("ok", "skipped", "skipped"))
  expect_match(gauss$message[3], "^detected in 1 \\(x\\) and [0-9]+ \\(y\\) of the 311 shared pixels; minDetected asks for 18 in each: not tested$")
  cols <- setdiff(names(gauss), "padj")  # (padj is adjusted over the tested genes: 1 here, 2 with minDetected = 0)
  expect_identical(rm_table(gauss)[1, cols], rm_table(run(minDetected = 0, surrogate = "gaussian"))[1, cols])
  # an explicit minDetected applies in both modes
  for (mode in c("gaussian", "remap")) {
    expect_identical(attr(run(minDetected = 0.5, surrogate = mode), "params")$minDetectedPixels, 156L, info = mode)
  }
  # the brain grid: sqrt(2170) rounds up to 47 pixels with gaussian surrogates, 0 with remapped ones
  brain <- fx_spe_pair(rf$pairs$brain, c("Slc17a7", "Efemp1"))
  expect_identical(attr(rm_run(brain, tests = "similarity"), "params")$minDetectedPixels, NA_integer_)  # no correlation test
  brain_run <- function(...) {
    attr(rm_run(brain, tests = "correlation", nPermutations = 2, exceedances = Inf, delta = c(0.3, 0.6), ...), "params")
  }
  expect_identical(brain_run(surrogate = "gaussian")$minDetectedPixels, 47L)
  expect_identical(brain_run()$minDetectedPixels, 0L)
})

test_that("portable: surrogate = \"gaussian\" gives the results of the version before the surrogate argument existed; the default is remap", {
  genes <- c("Gpx1", "Upk2", "Agt")
  input <- fx_spe_pair(rf$pairs$aki, genes)
  ref <- rm_run(input, nPermutations = 60, exceedances = 3, delta = rm_grid, keepNulls = TRUE, surrogate = "gaussian")
  expect_identical(attr(ref, "params")$surrogate, "gaussian")
  expect_identical(ref$nPermutations, c(60L, 3L, 40L))  # as before the surrogate argument existed
  expect_identical(ref$p[1], 1 / 61)
  expect_output(print(ref), "5 deltas from 0.05 to 0.9; gaussian surrogates\n")
  expect_output(print(summary(ref)), "gaussian surrogates")
  # the exceedances of gaussian surrogates follow the exact rule |null| >= |r|, without the tie tolerance
  det <- attr(ref, "details")
  for (g in genes) {
    L <- ref[g, "nPermutations"]
    expect_identical(ref[g, "pX"], (sum(abs(det[[g]]$nullX) >= abs(ref[g, "r"])) + 1) / (L + 1), info = g)
  }
  # the default is remap, which gives other nulls and p-values
  dflt <- rm_run(input, nPermutations = 60, exceedances = 3, delta = rm_grid, keepNulls = TRUE)
  expect_identical(attr(dflt, "params")$surrogate, "remap")
  expect_false(identical(attr(dflt, "details")$Gpx1$nullX, det$Gpx1$nullX))
  expect_identical(rm_table(dflt), rm_table(rm_run(input, nPermutations = 60, exceedances = 3, delta = rm_grid,
                                                   surrogate = "remap")))
  # the internal interface keeps gaussian as its default (the legacy functions): it equals the explicit mode, and
  # the exceedances follow the exact rule
  P <- rf$pairs$aki
  X <- t(P$X[genes, ])
  Y <- t(P$Y[genes, ])
  strip <- function(res) {
    attr(res, "state") <- NULL
    res
  }
  e1 <- strip(.stc_engine_correlate(X, Y, P$pos, deltaX = rm_grid, deltaY = rm_grid, nPermutations = 7L))
  e2 <- strip(.stc_engine_correlate(X, Y, P$pos, deltaX = rm_grid, deltaY = rm_grid, nPermutations = 7L, surrogate = "gaussian"))
  expect_identical(e2, e1)
  expect_identical(e1$bX, vapply(seq_along(genes), function(k) sum(abs(e1$nullX[[k]]) >= abs(e1$r[k])), 0L))
})
