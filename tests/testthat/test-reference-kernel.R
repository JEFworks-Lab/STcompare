# Tier 0: reference tests of the legacy R implementation (fixtures/kernel_fixture.rds, built by
# data-raw/build_test_fixtures.R tier0; see data-raw/README.md). Labels ("exact", "portable", "canary",
# "fixture integrity") are explained in helper-fixtures.R.
#
# Cases (B = 10 unless noted): kidney_AB (hexagonal lattice, many pixel pairs tied at max.dist), kidney_AC,
# kidney_AB_jitter (jittered coordinates), quakes_irregular (datasets::quakes, coordinates on a 0.01-degree
# grid), brain_Oprk1_subsample (N = 2170 > 1000: variogram subsample; B = 5), kidney_AB_jitter_seed17
# (seed = 17, maxDistPrctile = 0.3, deltaX != deltaY; B = 3), brain_Oprk1_N1000 (N = 1000 exactly: no
# subsample; B = 2). Lines marked "C++:" (helper-fixtures.R) are where a compiled backend's kernels go.

fx <- fx_read("kernel_fixture.rds")
KF <- "kernel_fixture.rds"
tol <- 1e-12
dirs <- c(X = "forward", Y = "reverse")
dir_delta <- function(cs, dir) if (dir %in% c("forward", "X")) cs$input$deltaX else cs$input$deltaY

test_that("fixture integrity: kernel fixture (not regression coverage)", {
  expect_true(all(c("kidney_AB", "kidney_AC", "kidney_AB_jitter", "quakes_irregular", "brain_Oprk1_subsample",
                    "kidney_AB_jitter_seed17", "brain_Oprk1_N1000") %in% names(fx$cases)))
  for (cs in fx$cases) {
    expect_identical(cs$intermediate$forward$perm_index, cs$intermediate$reverse$perm_index, info = cs$name)
    for (dir in dirs) {
      ed <- cs$expected[[dir]]
      d <- cs$intermediate[[dir]]$detail
      expect_identical(ed$pValue, ed$nExtreme / cs$input$nPermutations, info = cs$name)
      expect_true(all(ed$deltaStar %in% dir_delta(cs, dir)), info = cs$name)
      if (!is.null(d$hat_star)) expect_identical(d$per_delta[[d$delta_star_id]]$hat, d$hat_star, info = cs$name)
    }
    expect_false(is.null(fx$meta$platform_signature$sets[[cs$name]]), info = cs$name)
  }
})

test_that("canary: tier-0 inputs are reproducible from data(speKidney), SEraster and datasets::quakes", {
  rk <- fx_speKidney_raster()
  coords <- function(s) SpatialExperiment::spatialCoords(s)
  for (nm in c("kidney_AB", "kidney_AC")) {
    inp <- fx$cases[[nm]]$input
    other <- if (nm == "kidney_AB") rk$B else rk$C
    sh <- intersect(rownames(coords(rk$A)), rownames(coords(other)))
    expect_identical(sh, inp$pixel, info = nm)
    expect_equal(as.numeric(SummarizedExperiment::assay(rk$A)[1, sh]), inp$X, tolerance = tol, info = nm)
    expect_equal(as.numeric(SummarizedExperiment::assay(other)[1, sh]), inp$Y, tolerance = tol, info = nm)
    expect_equal(unname(coords(rk$A)[sh, ]), inp$pos, tolerance = tol, info = nm)
  }
  ab <- fx$cases$kidney_AB$input
  jit <- fx$cases$kidney_AB_jitter$input
  local_default_rng()
  set.seed(7)
  expect_equal(jit$pos, ab$pos + matrix(runif(2 * length(ab$X), -0.02, 0.02), ncol = 2), tolerance = tol)
  expect_identical(jit$X, ab$X)
  s17 <- fx$cases$kidney_AB_jitter_seed17$input
  expect_identical(s17[c("X", "Y", "pos")], jit[c("X", "Y", "pos")])
  q <- datasets::quakes
  q <- q[!duplicated(cbind(q$lat, q$long)), ][1:300, ]
  qi <- fx$cases$quakes_irregular$input
  expect_identical(qi$X, as.numeric(q$depth))
  expect_identical(qi$Y, as.numeric(q$mag))
  expect_identical(qi$pos, unname(cbind(q$lat, q$long)))
  br <- fx$cases$brain_Oprk1_subsample$input
  b1k <- fx$cases$brain_Oprk1_N1000$input
  expect_identical(b1k$X, br$X[1:1000])
  expect_identical(b1k$Y, br$Y[1:1000])
  expect_identical(b1k$pos, br$pos[1:1000, ])
})

test_that("canary: RNG streams on every platform (subsample, permutation order, noise of permutation 1)", {
  local_default_rng()
  for (cs in fx$cases) {
    N <- length(cs$input$X)
    B <- cs$input$nPermutations
    set.seed(cs$input$seed)
    ids <- if (N > 1000) sample(N, 1000) else seq_len(N)                          # no subsample at N = 1000
    idx <- vapply(seq_len(B), function(i) sample.int(N, N), integer(N))
    for (dir in dirs) {
      expect_identical(ids, cs$intermediate[[dir]]$ids, info = paste(cs$name, dir))
      expect_identical(idx, cs$intermediate[[dir]]$perm_index, info = paste(cs$name, dir))
      expect_identical(cs$intermediate[[dir]]$detail$noise_seed, cs$input$seed + 1, info = paste(cs$name, dir))
    }
    # the package draws sample(X, length(X)), which consumes the same stream
    set.seed(cs$input$seed)
    if (N > 1000) invisible(sample(N, 1000))
    expect_identical(sample(cs$input$X, N), cs$input$X[cs$intermediate$forward$perm_index[, 1]], info = cs$name)
    # noise: rnorm(N) per delta after set.seed(seed + 1) under L'Ecuyer-CMRG; 1e-14 absorbs the last-bit
    # differences of R's qnorm() between platforms
    pd <- cs$intermediate$forward$detail$per_delta
    if (!is.null(pd[[1]]$noise)) {
      noise <- ref_noise(N, cs$input$seed + 1, length(pd))                       # C++: RNG (or injected)
      for (k in seq_along(pd)) expect_equal(noise[, k], pd[[k]]$noise, tolerance = 1e-14, info = paste(cs$name, k))
    }
  }
})

for (nm in names(fx$cases)) {
  test_that(sprintf("canary (exact): %s: permutation 1 replayed with the reference kernels matches the fixture", nm), {
    skip_if_not_exact(KF, nm)
    cs <- fx$cases[[nm]]
    lat <- cs$input$pos[, 1]
    long <- cs$input$pos[, 2]
    for (dir in dirs) {
      cap <- cs$intermediate[[dir]]
      d <- cap$detail
      rp <- fx_replay(cs, dir)
      info <- paste(nm, dir)
      expect_equal(rp$prctile, cap$prctile, tolerance = tol, info = info)
      expect_identical(as.numeric(rp$target$n), as.numeric(cap$target_variog$n), info = info)
      expect_equal(rp$target$u, cap$target_variog$u, tolerance = tol, info = info)
      expect_equal(rp$target$v, cap$target_variog$v, tolerance = tol, info = info)
      expect_equal(rp$target$bins.lim, cap$target_variog$bins.lim, tolerance = tol, info = info)
      for (k in seq_along(rp$delta)) {
        p <- d$per_delta[[k]]
        r <- rp$per[[k]]
        info_k <- sprintf("%s %s delta=%g", nm, dir, rp$delta[k])
        expect_identical(as.numeric(r$n1), as.numeric(cap$target_variog$n), info = info_k)
        expect_equal(r$v1, p$variog_fitted_v, tolerance = tol, info = info_k)
        expect_equal(r$bet, p$lm_coef, tolerance = tol, info = info_k)
        expect_equal(r$v2, p$variog_hat_v, tolerance = tol, info = info_k)
        expect_equal(r$rss, p$rss, tolerance = tol, info = info_k)
        if (!is.null(p$fitted)) {   # complete vectors: kidney_AB, kidney_AB_jitter, quakes_irregular (forward)
          expect_equal(r$xd, p$fitted, tolerance = tol, info = info_k)
          expect_equal(r$hat, p$hat, tolerance = tol, info = info_k)
          expect_equal(ref_smooth(rp$xr, long, lat, rp$delta[k], exact = TRUE), p$fitted_exact_evdat,
                       tolerance = tol, info = info_k)
        }
      }
      expect_equal(rp$rss, d$residus, tolerance = tol, info = info)
      expect_equal(rp$rss, cap$residus[, d$i], tolerance = tol, info = info)
      expect_identical(rp$argmin, d$delta_star_id, info = info)
    }
  })
}

for (nm in names(fx$cases)) {
  test_that(sprintf("exact: %s: matchingVariograms() and spatialCorrelation() reproduce the stored RSS, nulls, deltaStar and raw p-values", nm), {
    skip_if_not_exact(KF, nm)
    cs <- fx$cases[[nm]]
    for (dir in dirs) {
      mv <- fx_mv(cs, dir)
      d <- cs$intermediate[[dir]]$detail
      expect_equal(mv$residus, d$residus, tolerance = tol, info = paste(nm, dir))
      expect_identical(mv$delta.star.id, d$delta_star_id, info = paste(nm, dir))
      if (!is.null(d$hat_star)) expect_equal(as.numeric(mv$hat.X.delta.star), d$hat_star, tolerance = tol, info = paste(nm, dir))
    }
    o <- fx_case_run(cs)
    e <- cs$expected
    for (dir in names(dirs)) {
      ed <- e[[dirs[[dir]]]]
      info <- paste(nm, dir)
      null <- as.numeric(o[[paste0("nullCorrelations", dir)]][[1]])
      expect_equal(null, ed$nullCor, tolerance = tol, info = info)
      expect_identical(as.numeric(o[[paste0("deltaStar", dir)]][[1]]), ed$deltaStar, info = info)
      expect_identical(o[[paste0("deltaStarMedian", dir)]], ed$deltaStarMedian, info = info)
      expect_identical(sum(abs(null) > abs(ed$r_obs)), ed$nExtreme, info = info)
      expect_equal(o[[paste0("pValuePermute", dir)]], stc_empirical_p(ed$nullCor, ed$r_obs), info = info)
      P <- o[[paste0("permutations", dir)]][[1]]
      expect_equal(perm_fingerprint(P), ed$perm_fingerprint, tolerance = tol, info = info)
      if (!is.null(ed$perm1)) expect_equal(as.numeric(P[, 1]), ed$perm1, tolerance = tol, info = info)
    }
  })
}

for (nm in names(fx$cases)) {
  test_that(sprintf("portable: %s: matchingVariograms() equals the reference kernels on this machine; spatialCorrelation() follows the p-value definition and agrees loosely with the fixture", nm), {
    cs <- fx$cases[[nm]]
    for (dir in dirs) {
      mv <- fx_mv(cs, dir)
      rp <- fx_replay(cs, dir)
      expect_equal(mv$residus, rp$rss, tolerance = tol, info = paste(nm, dir))
      expect_identical(mv$delta.star.id, rp$argmin, info = paste(nm, dir))
      expect_equal(as.numeric(mv$hat.X.delta.star), rp$hat_star, tolerance = tol, info = paste(nm, dir))
    }
    o <- fx_case_run(cs)
    e <- cs$expected
    r <- unname(o$correlationCoef)
    expect_equal(r, e$spatialCorrelation$correlationCoef, tolerance = tol)
    expect_equal(o$pValueNaive, e$spatialCorrelation$pValueNaive, tolerance = 1e-10)
    for (dir in names(dirs)) {
      info <- paste(nm, dir)
      null <- as.numeric(o[[paste0("nullCorrelations", dir)]][[1]])
      ds <- as.numeric(o[[paste0("deltaStar", dir)]][[1]])
      expect_length(null, cs$input$nPermutations)
      expect_equal(o[[paste0("pValuePermute", dir)]], stc_empirical_p(null, r), info = info)
      expect_true(all(ds %in% dir_delta(cs, dir)), info = info)
      expect_identical(o[[paste0("deltaStarMedian", dir)]], stats::median(ds), info = info)
      expect_identical(dim(o[[paste0("permutations", dir)]][[1]]), as.integer(c(length(cs$input$X), cs$input$nPermutations)), info = info)
      expect_nulls_close(null, e[[dirs[[dir]]]]$nullCor, info = info)
    }
  })
}

test_that("portable: swapping X and Y swaps the two directions, and 3 permutations are a prefix of 10 (quakes)", {
  cs <- fx$cases$quakes_irregular
  o10 <- fx_case_run(cs)
  local_default_rng()
  o3 <- quiet_locfit(spatialCorrelation(cs$input$Y, cs$input$X, cs$input$pos, nPermutations = 3))
  expect_equal(unname(o3$correlationCoef), unname(o10$correlationCoef), tolerance = 1e-15)
  expect_equal(as.numeric(o3$nullCorrelationsX[[1]]), as.numeric(o10$nullCorrelationsY[[1]])[1:3], tolerance = tol)
  expect_equal(as.numeric(o3$nullCorrelationsY[[1]]), as.numeric(o10$nullCorrelationsX[[1]])[1:3], tolerance = tol)
  expect_identical(as.numeric(o3$deltaStarX[[1]]), as.numeric(o10$deltaStarY[[1]])[1:3])
  expect_identical(as.numeric(o3$deltaStarY[[1]]), as.numeric(o10$deltaStarX[[1]])[1:3])
})

test_that("portable: with returnPermutations = FALSE, pValuePermuteX/Y come from their own directions (quakes)", {
  cs <- fx$cases$quakes_irregular
  ot <- fx_case_run(cs)
  of <- fx_memo("sc_false:quakes_irregular", fx_sc_case(cs))                      # default returnPermutations
  r <- unname(of$correlationCoef)
  nx <- as.numeric(of$nullCorrelationsX[[1]])
  ny <- as.numeric(of$nullCorrelationsY[[1]])
  expect_false(any(c("permutationsX", "permutationsY") %in% names(of)))
  expect_identical(nx, as.numeric(ot$nullCorrelationsX[[1]]))
  expect_identical(ny, as.numeric(ot$nullCorrelationsY[[1]]))
  # this case has pX != pY (0 and 0.4 at B = 10 on every platform tested), so a swap is detected
  expect_equal(of$pValuePermuteX, stc_empirical_p(nx, r))
  expect_equal(of$pValuePermuteY, stc_empirical_p(ny, r))
})

test_that("exact: quakes with returnPermutations = FALSE reproduces the stored raw p-values", {
  skip_if_not_exact(KF, "quakes_irregular")
  cs <- fx$cases$quakes_irregular
  of <- fx_memo("sc_false:quakes_irregular", fx_sc_case(cs))
  expect_equal(of$pValuePermuteX, stc_empirical_p(cs$expected$forward$nullCor, cs$expected$forward$r_obs))
  expect_equal(of$pValuePermuteY, stc_empirical_p(cs$expected$reverse$nullCor, cs$expected$reverse$r_obs))
})

test_that("portable: the empirical p-value counts only strictly larger |null| (ties do not count)", {
  q <- fx$cases$quakes_irregular$input
  X <- sort(q$X)
  Y <- q$Y
  local_default_rng()
  # every null field is X itself, so every |null correlation| equals |r| exactly
  o <- testthat::with_mocked_bindings(
    viladomatCorrelation(cbind(X, Y, q$pos[, 1], q$pos[, 2]), 0.5, 0.25, 4, BPPARAM = BiocParallel::SerialParam()),
    matchingVariograms = function(X.randomized, ...) {
      list(residus = 0, delta.star.id = 1L, hat.X.delta.star = sort(X.randomized))
    },
    .package = "STcompare")
  expect_identical(abs(as.numeric(o$nullCorGlobal)), rep(abs(as.vector(stats::cor(X, Y))), 4))
  expect_identical(o$pValueGlobal, 0)
})

test_that("portable: results do not depend on the number of workers (nThreads = 2)", {
  skip_on_cran()
  skip_on_os("windows")
  cs <- fx$cases$kidney_AB_jitter
  o1 <- fx_case_run(cs)
  o2 <- fx_sc_case(cs, returnPermutations = TRUE, nThreads = 2)
  for (k in c("nullCorrelationsX", "nullCorrelationsY", "deltaStarX", "deltaStarY", "permutationsX", "permutationsY")) {
    expect_identical(o2[[k]][[1]], o1[[k]][[1]], info = k)
  }
  expect_identical(c(o2$pValuePermuteX, o2$pValuePermuteY), c(o1$pValuePermuteX, o1$pValuePermuteY))
})

test_that("portable: spatialCorrelationGeneExp() forwards seed, deltas and maxDistPrctile (seed = 17, B = 3)", {
  rk <- fx_speKidney_raster()
  sh <- intersect(rownames(SpatialExperiment::spatialCoords(rk$A)), rownames(SpatialExperiment::spatialCoords(rk$B)))
  local_default_rng()
  og <- quiet_locfit(spatialCorrelationGeneExp(list(rk$A, rk$B), nPermutations = 3, verbose = FALSE, seed = 17))
  os <- quiet_locfit(spatialCorrelation(as.numeric(SummarizedExperiment::assay(rk$A)[1, sh]),
                                        as.numeric(SummarizedExperiment::assay(rk$B)[1, sh]),
                                        SpatialExperiment::spatialCoords(rk$A)[sh, ], nPermutations = 3, seed = 17))
  expect_identical(nrow(og), 1L)
  expect_equal(unname(og$correlationCoef), unname(os$correlationCoef), tolerance = tol)
  for (k in c("nullCorrelationsX", "nullCorrelationsY", "deltaStarX", "deltaStarY")) {
    expect_identical(as.numeric(og[[k]][[1]]), as.numeric(os[[k]][[1]]), info = k)
  }
  # a single gene: a multiple-testing adjustment across genes leaves its p-values unchanged
  expect_identical(c(og$pValuePermuteX, og$pValuePermuteY), c(os$pValuePermuteX, os$pValuePermuteY))
})

test_that("exact: spatialCorrelationGeneExp() with default arguments reproduces kidney_AB", {
  # nulls, deltaStar and r only: the wrapper's p-value adjustment is applied per gene (a known legacy issue)
  skip_if_not_exact(KF, "kidney_AB")
  rk <- fx_speKidney_raster()
  cs <- fx$cases$kidney_AB
  local_default_rng()
  o <- quiet_locfit(spatialCorrelationGeneExp(list(rk$A, rk$B), nPermutations = cs$input$nPermutations, verbose = FALSE))
  expect_equal(unname(o$correlationCoef), cs$expected$spatialCorrelation$correlationCoef, tolerance = tol)
  expect_equal(as.numeric(o$nullCorrelationsX[[1]]), cs$expected$forward$nullCor, tolerance = tol)
  expect_equal(as.numeric(o$nullCorrelationsY[[1]]), cs$expected$reverse$nullCor, tolerance = tol)
  expect_identical(as.numeric(o$deltaStarX[[1]]), cs$expected$forward$deltaStar)
  expect_identical(as.numeric(o$deltaStarY[[1]]), cs$expected$reverse$deltaStar)
})

test_that("portable: legacy error paths return NA p-values instead of failing (B = 2)", {
  # (never the 1 <= delta * N < 2 inputs, where locfit kills the R session)
  inp <- fx$cases$kidney_AB_jitter$input
  # runs spatialCorrelation(), silencing the printed error and muffling (and counting) only the warning
  # expected for that input; any other warning reaches testthat
  run <- function(X, delta, expected_warning = NULL) {
    n_expected <- 0L
    out <- NULL
    withCallingHandlers(
      utils::capture.output(out <- quiet_locfit(spatialCorrelation(X, inp$Y, inp$pos, nPermutations = 2,
                                                                   deltaX = delta, deltaY = delta))),
      warning = function(w) {
        if (!is.null(expected_warning) && grepl(expected_warning, conditionMessage(w), fixed = TRUE)) {
          n_expected <<- n_expected + 1L
          invokeRestart("muffleWarning")
        }
      })
    list(out = out, n_expected = n_expected)
  }
  local_default_rng()
  # constant X: cor.test() warns that the standard deviation is zero; r, naive p and both empirical p are NA
  r1 <- run(rep(1, length(inp$X)), 0.3, "standard deviation is zero")
  o1 <- r1$out
  expect_gt(r1$n_expected, 0)
  expect_true(all(is.na(c(o1$correlationCoef, o1$pValueNaive, o1$pValuePermuteX, o1$pValuePermuteY))))
  # delta * N < 1: locfit warns "procv: no points with non-zero weight" and fails; the error is caught and
  # the empirical p-values are NA
  r2 <- run(inp$X, 0.001, "procv: no points with non-zero weight")
  o2 <- r2$out
  expect_gt(r2$n_expected, 0)
  expect_equal(unname(o2$correlationCoef), stats::cor(inp$X, inp$Y), tolerance = tol)
  expect_true(all(is.na(c(o2$pValuePermuteX, o2$pValuePermuteY))))
  # one NA in X: r and the naive p-value use the complete pairs; the empirical p-values are NA
  Xn <- inp$X
  Xn[5] <- NA
  o3 <- run(Xn, 0.3)$out
  expect_equal(unname(o3$correlationCoef), stats::cor(inp$X[-5], inp$Y[-5]), tolerance = tol)
  expect_equal(o3$pValueNaive, stats::cor.test(inp$X[-5], inp$Y[-5])$p.value, tolerance = 1e-10)
  expect_true(all(is.na(c(o3$pValuePermuteX, o3$pValuePermuteY))))
})

test_that("portable: spatialSimilarity() thresholds, pseudo-counts, fold-change boundaries and minPixels (hand-computed)", {
  # 20 pixels, 3 genes. Every gene has at least two zeros in x and in y, so the 5% quantile thresholds are
  # t1 = t2 = 0, and a pixel is kept only if x > 0 or y > 0 (strictly). Zeros become 1e-4 before y / x.
  px <- paste0("px", 1:20)
  x <- rbind(g1 = c(0, 0, 0, 1, 2, 4, 4, 2, 8, 3, 6, 1, 1, 5, 5, 2, 2, 7, 7, 9),
             g2 = c(rep(0, 19), 3),
             g3 = c(rep(0, 18), 3, 5))
  y <- rbind(g1 = c(0, 1, 3, 2, 1, 2, 8, 4, 4, 6, 3, 0, 1, 5, 10, 1, 4, 7, 14, 9),
             g2 = rep(0, 20),
             g3 = c(rep(0, 18), 6, 5))
  colnames(x) <- colnames(y) <- px
  mk <- function(m) SpatialExperiment::SpatialExperiment(assays = list(counts = m), spatialCoords = cbind(x = seq_len(20), y = 0))
  s <- spatialSimilarity(list(mk(x), mk(y)))
  st <- s$similarityTable
  expect_identical(st$gene, c("g1", "g2", "g3"))
  expect_identical(st$t1, c(0, 0, 0))
  expect_identical(st$t2, c(0, 0, 0))
  # g1: px1 = (0, 0) is out; 19 in. Ratios y/x of px2..px20; log2 = +-1 exactly (ratio 2 or 1/2) is similar
  # (inclusive), px2 and px3 (x = 0 -> 1e-4) are dissimilar towards Y, px12 (y = 0 -> 1e-4) towards X.
  g1 <- st[st$gene == "g1", ]
  expect_equal(g1$percentSimilarity, 16 / 19, tolerance = tol)
  expect_equal(g1$percentDissimilarityX, 1 / 19, tolerance = tol)
  expect_equal(g1$percentDissimilarityY, 2 / 19, tolerance = tol)
  expect_identical(as.numeric(g1$numPixelInThresh), 19)
  expect_identical(as.numeric(g1$numPixelOutThresh), 1)
  expect_identical(g1$pixelIDOutThresh[[1]], "px1")
  expect_identical(g1$pixelIDInThresh[[1]], px[2:20])
  expect_setequal(g1$similarPixelID[[1]], px[c(4:11, 13:20)])
  expect_identical(g1$dissimilarPixelIDX[[1]], "px12")
  expect_identical(g1$dissimilarPixelIDY[[1]], c("px2", "px3"))
  l2 <- s$pixelLogTransformation
  expect_identical(l2$gene, c("g1", "g3"))
  expect_equal(l2$log[[1]], log2(c(1 / 1e-4, 3 / 1e-4, 2 / 1, 1 / 2, 2 / 4, 8 / 4, 4 / 2, 4 / 8, 6 / 3, 3 / 6,
                                   1e-4 / 1, 1 / 1, 5 / 5, 10 / 5, 1 / 2, 4 / 2, 7 / 7, 14 / 7, 9 / 9)), tolerance = tol)
  # g2: one pixel passes, fewer than minPixels * 20 = 2: no similarity values (numPixelInThresh is not
  # checked: in this branch the legacy code reports the row count of its one-row summary)
  g2 <- st[st$gene == "g2", ]
  expect_true(all(is.na(c(g2$percentSimilarity, g2$percentDissimilarityX, g2$percentDissimilarityY))))
  expect_identical(as.numeric(g2$numPixelOutThresh), 19)
  # g3: exactly minPixels * 20 = 2 pixels pass (not fewer), so it is scored: ratios 2 and 1, both similar
  g3 <- st[st$gene == "g3", ]
  expect_identical(c(g3$percentSimilarity, g3$percentDissimilarityX, g3$percentDissimilarityY), c(1, 0, 0))
  expect_identical(as.numeric(c(g3$numPixelInThresh, g3$numPixelOutThresh)), c(2, 18))
  expect_equal(l2$log[[2]], c(1, 0), tolerance = tol)
  # given thresholds t1 = t2 = 2: a pixel with x == 2 or y == 2 (and the other not above 2) is out
  s2 <- spatialSimilarity(list(mk(x), mk(y)), t1 = 2, t2 = 2)
  st2 <- s2$similarityTable
  g1b <- st2[st2$gene == "g1", ]
  expect_identical(g1b$pixelIDInThresh[[1]], px[c(3, 6:11, 14:15, 17:20)])
  expect_equal(c(g1b$percentSimilarity, g1b$percentDissimilarityX, g1b$percentDissimilarityY), c(12, 0, 1) / 13, tolerance = tol)
  expect_identical(c(st2$t1, st2$t2), rep(2, 6))
  expect_true(is.na(st2$percentSimilarity[st2$gene == "g2"]))
  expect_identical(st2$percentSimilarity[st2$gene == "g3"], 1)
})

test_that("portable: spatialSimilarity() recomputed from data(speKidney) matches the fixture", {
  rk <- fx_speKidney_raster()
  for (nm in c("kidney_AB", "kidney_AC", "kidney_AC_foldChange2")) {
    ref <- fx$similarity[[nm]]
    second <- if (nm == "kidney_AB") rk$B else rk$C
    s <- do.call(spatialSimilarity, c(list(list(rk$A, second)), ref$args))
    st <- s$similarityTable
    expect_identical(st$gene, ref$table$gene, info = nm)
    for (col in c("percentSimilarity", "percentDissimilarityX", "percentDissimilarityY", "t1", "t2")) {
      expect_equal(st[[col]], ref$table[[col]], tolerance = tol, info = paste(nm, col))
    }
    expect_identical(as.numeric(st$numPixelInThresh), as.numeric(ref$table$numPixelInThresh), info = nm)
    expect_identical(as.numeric(st$numPixelOutThresh), as.numeric(ref$table$numPixelOutThresh), info = nm)
    expect_setequal(st$similarPixelID[[1]], ref$similarPixelID)
    expect_setequal(st$dissimilarPixelIDX[[1]], ref$dissimilarPixelIDX)
    expect_setequal(st$dissimilarPixelIDY[[1]], ref$dissimilarPixelIDY)
    expect_setequal(st$pixelIDInThresh[[1]], ref$pixelIDInThresh)
    expect_setequal(st$pixelIDOutThresh[[1]], ref$pixelIDOutThresh)
    l2 <- setNames(s$pixelLogTransformation$log[[1]], st$pixelIDInThresh[[1]])
    l2_ref <- setNames(ref$log2ratio, ref$pixelIDInThresh)
    expect_equal(l2[sort(names(l2))], l2_ref[sort(names(l2_ref))], tolerance = tol, info = nm)
  }
})
