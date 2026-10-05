# Tier 0: the exported functions (computed by the compiled engine) against the stored outputs of the
# original R implementation (fixtures/kernel_fixture.rds, built by data-raw/build_test_fixtures.R tier0; see
# data-raw/README.md). Labels ("exact", "portable", "canary", "fixture integrity") are explained in
# helper-fixtures.R.
#
# Cases (B = 10 unless noted): kidney_AB (hexagonal lattice, many pixel pairs tied at max.dist), kidney_AC,
# kidney_AB_jitter (jittered coordinates), quakes_irregular (datasets::quakes, coordinates on a 0.01-degree
# grid), brain_Oprk1_subsample (N = 2170 > 1000: variogram subsample; B = 5), kidney_AB_jitter_seed17
# (seed = 17, maxDistPrctile = 0.3, deltaX != deltaY; B = 3), brain_Oprk1_N1000 (N = 1000 exactly: no
# subsample; B = 2).

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

test_that("canary: R's RNG streams match the fixture on every platform (subsample, permutation order, noise of permutation 1)", {
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
    # the original code drew sample(X, length(X)), which consumes the same stream
    set.seed(cs$input$seed)
    if (N > 1000) invisible(sample(N, 1000))
    expect_identical(sample(cs$input$X, N), cs$input$X[cs$intermediate$forward$perm_index[, 1]], info = cs$name)
    # noise: rnorm(N) per delta after set.seed(seed + 1) under L'Ecuyer-CMRG; 1e-14 absorbs the last-bit
    # differences of R's qnorm() between platforms
    pd <- cs$intermediate$forward$detail$per_delta
    if (!is.null(pd[[1]]$noise)) {
      noise <- ref_noise(N, cs$input$seed + 1, length(pd))
      for (k in seq_along(pd)) expect_equal(noise[, k], pd[[k]]$noise, tolerance = 1e-14, info = paste(cs$name, k))
      # the C++ generator of the engine draws the same normals
      expect_equal(.stc_legacy_noise(as.integer(cs$input$seed + 1), N, length(pd)), noise, tolerance = 1e-14, info = cs$name)
    }
  }
})

for (nm in names(fx$cases)) {
  test_that(sprintf("exact: %s: spatialCorrelation() reproduces the stored nulls (1e-9), deltaStar, raw counts and permuted fields", nm), {
    skip_if_not_exact(KF, nm)
    cs <- fx$cases[[nm]]
    o <- fx_case_run(cs)
    for (dir in names(dirs)) {
      ed <- cs$expected[[dirs[[dir]]]]
      info <- paste(nm, dir)
      null <- as.numeric(o[[paste0("nullCorrelations", dir)]][[1]])
      expect_close(null, ed$nullCor, info = info)
      expect_identical(as.numeric(o[[paste0("deltaStar", dir)]][[1]]), ed$deltaStar, info = info)
      expect_identical(o[[paste0("deltaStarMedian", dir)]], ed$deltaStarMedian, info = info)
      expect_identical(sum(abs(null) > abs(ed$r_obs)), ed$nExtreme, info = info)
      expect_identical(o[[paste0("pValuePermute", dir)]], stc_empirical_p(ed$nullCor, ed$r_obs), info = info)
      P <- o[[paste0("permutations", dir)]][[1]]
      expect_close(perm_fingerprint(P), ed$perm_fingerprint, info = paste(info, "fingerprint"))
      if (!is.null(ed$perm1)) expect_close(P[, 1], ed$perm1, info = paste(info, "permutation 1"))
    }
  })
}

for (nm in names(fx$cases)) {
  test_that(sprintf("portable: %s: spatialCorrelation() output structure, the p-value definition and loose agreement with the fixture", nm), {
    cs <- fx$cases[[nm]]
    o <- fx_case_run(cs)
    e <- cs$expected
    N <- length(cs$input$X)
    B <- cs$input$nPermutations
    expect_s3_class(o, "data.frame")
    expect_identical(dim(o), c(1L, 12L))
    expect_identical(rownames(o), "cor")
    expect_identical(names(o), c("correlationCoef", "pValueNaive", "pValuePermuteX", "pValuePermuteY",
                                 "deltaStarMedianX", "deltaStarMedianY", "deltaStarX", "deltaStarY",
                                 "nullCorrelationsX", "nullCorrelationsY", "permutationsX", "permutationsY"))
    r <- o$correlationCoef
    expect_equal(r, e$spatialCorrelation$correlationCoef, tolerance = tol)
    expect_equal(o$pValueNaive, e$spatialCorrelation$pValueNaive, tolerance = 1e-10)
    for (dir in names(dirs)) {
      info <- paste(nm, dir)
      nc <- o[[paste0("nullCorrelations", dir)]]
      expect_s3_class(nc, "AsIs")
      expect_identical(dim(nc[[1]]), c(as.integer(B), 1L), info = info)
      expect_null(dimnames(nc[[1]]))
      null <- as.numeric(nc[[1]])
      ds <- o[[paste0("deltaStar", dir)]][[1]]
      expect_type(ds, "double")
      expect_length(ds, B)
      expect_identical(o[[paste0("pValuePermute", dir)]], stc_empirical_p(null, r), info = info)
      expect_true(all(ds %in% dir_delta(cs, dir)), info = info)
      expect_identical(o[[paste0("deltaStarMedian", dir)]], stats::median(ds), info = info)
      P <- o[[paste0("permutations", dir)]][[1]]
      expect_identical(dim(P), as.integer(c(N, B)), info = info)
      expect_null(dimnames(P))
      expect_nulls_close(null, e[[dirs[[dir]]]]$nullCor, info = info)
    }
  })
}

test_that("portable: swapping X and Y swaps the two directions, and 3 permutations are a prefix of 10 (quakes)", {
  cs <- fx$cases$quakes_irregular
  o10 <- fx_case_run(cs)
  local_default_rng()
  o3 <- spatialCorrelation(cs$input$Y, cs$input$X, cs$input$pos, nPermutations = 3)
  expect_equal(o3$correlationCoef, o10$correlationCoef, tolerance = 1e-15)
  expect_identical(as.numeric(o3$nullCorrelationsX[[1]]), as.numeric(o10$nullCorrelationsY[[1]])[1:3])
  expect_identical(as.numeric(o3$nullCorrelationsY[[1]]), as.numeric(o10$nullCorrelationsX[[1]])[1:3])
  expect_identical(o3$deltaStarX[[1]], o10$deltaStarY[[1]][1:3])
  expect_identical(o3$deltaStarY[[1]], o10$deltaStarX[[1]][1:3])
  expect_identical(o3$permutationsX, NULL)
})

test_that("portable: with returnPermutations = FALSE, pValuePermuteX/Y come from their own directions (quakes)", {
  cs <- fx$cases$quakes_irregular
  ot <- fx_case_run(cs)
  of <- fx_memo("sc_false:quakes_irregular", fx_sc_case(cs))                      # default returnPermutations
  r <- of$correlationCoef
  nx <- as.numeric(of$nullCorrelationsX[[1]])
  ny <- as.numeric(of$nullCorrelationsY[[1]])
  expect_identical(names(of), names(ot)[1:10])
  expect_identical(nx, as.numeric(ot$nullCorrelationsX[[1]]))
  expect_identical(ny, as.numeric(ot$nullCorrelationsY[[1]]))
  # this case has pX != pY at B = 10 on every platform tested, so a swap is detected
  expect_false(identical(of$pValuePermuteX, of$pValuePermuteY))
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

test_that("portable: results do not depend on the number of threads or on BPPARAM (kidney_AB_jitter)", {
  cs <- fx$cases$kidney_AB_jitter
  o1 <- fx_case_run(cs)
  o2 <- fx_sc_case(cs, returnPermutations = TRUE, nThreads = 2)
  expect_identical(o2, o1)
  # BPPARAM sets the number of threads (bpnworkers()); nothing is forked
  o3 <- fx_sc_case(cs, returnPermutations = TRUE, nThreads = 7, BPPARAM = BiocParallel::SerialParam())
  expect_identical(o3, o1)
})

test_that("portable: viladomatCorrelation() is the X direction of spatialCorrelation(), with its own output structure", {
  cs <- fx$cases$kidney_AB_jitter_seed17
  inp <- cs$input
  o <- fx_case_run(cs)
  local_default_rng()
  v <- viladomatCorrelation(cbind(inp$X, inp$Y, inp$pos), inp$deltaX, inp$maxDistPrctile, inp$nPermutations,
                            seed = inp$seed)
  expect_identical(names(v), c("deltaStarMedian", "deltaStar", "pValueGlobal", "nullCorGlobal", "permutations"))
  expect_identical(v$deltaStar, o$deltaStarX[[1]])
  expect_identical(v$deltaStarMedian, o$deltaStarMedianX)
  expect_identical(v$nullCorGlobal, o$nullCorrelationsX[[1]])
  expect_identical(v$permutations, o$permutationsX[[1]])
  expect_identical(v$pValueGlobal, o$pValuePermuteX)
  # the reverse direction is viladomatCorrelation() with X and Y swapped
  local_default_rng()
  vr <- viladomatCorrelation(data.frame(inp$Y, inp$X, inp$pos), inp$deltaY, inp$maxDistPrctile, inp$nPermutations,
                             seed = inp$seed)
  expect_identical(vr$nullCorGlobal, o$nullCorrelationsY[[1]])
  expect_identical(vr$deltaStar, o$deltaStarY[[1]])
})

test_that("portable: viladomatCorrelation() with an NA or a constant Y keeps the permutations; a constant X gives NA and a warning", {
  inp <- fx$cases$kidney_AB_jitter$input
  local_default_rng()
  ref <- viladomatCorrelation(cbind(inp$X, inp$Y, inp$pos), c(0.2, 0.6), 0.25, 3)
  for (Yb in list(replace(inp$Y, 4, NA), rep(2, length(inp$Y)))) {
    w <- fx_collect_warnings(viladomatCorrelation(cbind(inp$X, Yb, inp$pos), c(0.2, 0.6), 0.25, 3))
    expect_length(w$warnings, 1L)
    expect_match(w$warnings, "^viladomatCorrelation: nullCorGlobal or pValueGlobal is NA")
    expect_identical(w$value$deltaStar, ref$deltaStar)
    expect_identical(w$value$permutations, ref$permutations)
    expect_identical(w$value$nullCorGlobal, matrix(NA_real_, 3, 1))
    expect_identical(w$value$pValueGlobal, NA_real_)
  }
  w <- fx_collect_warnings(viladomatCorrelation(cbind(1, inp$Y, inp$pos), c(0.2, 0.6), 0.25, 3))
  expect_length(w$warnings, 1L)
  expect_match(w$warnings, "constant")
  expect_true(all(is.na(unlist(w$value))))
})

test_that("portable: spatialCorrelationGeneExp() forwards seed, deltas and maxDistPrctile (seed = 17, B = 3)", {
  rk <- fx_speKidney_raster()
  sh <- intersect(rownames(SpatialExperiment::spatialCoords(rk$A)), rownames(SpatialExperiment::spatialCoords(rk$B)))
  local_default_rng()
  og <- spatialCorrelationGeneExp(list(rk$A, rk$B), nPermutations = 3, verbose = FALSE, seed = 17,
                                  deltaX = list(c(0.1, 0.4)), deltaY = list(0.3), maxDistPrctile = 0.3)
  os <- spatialCorrelation(as.numeric(SummarizedExperiment::assay(rk$A)[1, sh]),
                           as.numeric(SummarizedExperiment::assay(rk$B)[1, sh]),
                           SpatialExperiment::spatialCoords(rk$A)[sh, ], nPermutations = 3, seed = 17,
                           deltaX = c(0.1, 0.4), deltaY = 0.3, maxDistPrctile = 0.3)
  expect_identical(nrow(og), 1L)
  expect_identical(rownames(og), "Gene")
  # a single gene: a multiple-testing adjustment across genes leaves its p-values unchanged
  rownames(og) <- "cor"
  expect_identical(og, os)
})

test_that("portable: spatialCorrelationGeneExp() adjusts p-values across genes with adjustMethod", {
  rk <- fx_speKidney_raster()
  sh <- intersect(rownames(SpatialExperiment::spatialCoords(rk$A)), rownames(SpatialExperiment::spatialCoords(rk$C)))
  a <- as.numeric(SummarizedExperiment::assay(rk$A)[1, sh])
  cc <- as.numeric(SummarizedExperiment::assay(rk$C)[1, sh])
  # three "genes" on the A-C pixels: X is always A; Y is C (positive), max(C) - C (negative) and C in
  # reversed pixel order (unrelated)
  mk <- function(rows) {
    m <- do.call(rbind, rows)
    dimnames(m) <- list(c("pos", "neg", "rev"), sh)
    SpatialExperiment::SpatialExperiment(assays = list(counts = m),
                                         spatialCoords = SpatialExperiment::spatialCoords(rk$A)[sh, ])
  }
  x <- mk(list(a, a, a))
  y <- mk(list(cc, max(cc) - cc, rev(cc)))
  local_default_rng()
  run <- function(method) spatialCorrelationGeneExp(list(x, y), nPermutations = 3, verbose = FALSE,
                                                    adjustMethod = method)
  raw <- run("none")
  expect_identical(rownames(raw), c("pos", "neg", "rev"))
  for (dir in c("X", "Y")) {
    col <- paste0("pValuePermute", dir)
    expect_equal(raw[[col]], mapply(stc_empirical_p, raw[[paste0("nullCorrelations", dir)]], raw$correlationCoef),
                 info = dir)
    for (method in c("BH", "bonferroni")) {
      expect_equal(run(method)[[col]], stats::p.adjust(raw[[col]], method = method), info = paste(dir, method))
    }
  }
  # an invalid method fails before any permutation is computed
  expect_error(spatialCorrelationGeneExp(list(x, y), nPermutations = 3, verbose = FALSE, adjustMethod = "nope"),
               "adjustMethod must be one of")
})

test_that("exact: spatialCorrelationGeneExp() with default arguments reproduces kidney_AB", {
  # nulls, deltaStar and r only; p-values are checked in the tests above
  skip_if_not_exact(KF, "kidney_AB")
  rk <- fx_speKidney_raster()
  cs <- fx$cases$kidney_AB
  local_default_rng()
  o <- spatialCorrelationGeneExp(list(rk$A, rk$B), nPermutations = cs$input$nPermutations, verbose = FALSE)
  expect_equal(o$correlationCoef, cs$expected$spatialCorrelation$correlationCoef, tolerance = tol)
  expect_close(as.numeric(o$nullCorrelationsX[[1]]), cs$expected$forward$nullCor)
  expect_close(as.numeric(o$nullCorrelationsY[[1]]), cs$expected$reverse$nullCor)
  expect_identical(o$deltaStarX[[1]], cs$expected$forward$deltaStar)
  expect_identical(o$deltaStarY[[1]], cs$expected$reverse$deltaStar)
})

test_that("portable: inputs where the R implementation gave NA p-values give an NA row and one warning (B = 2)", {
  inp <- fx$cases$kidney_AB_jitter$input
  N <- length(inp$X)
  run <- function(X, delta, pos = inp$pos) {
    fx_collect_warnings(spatialCorrelation(X, inp$Y, pos, nPermutations = 2, deltaX = delta, deltaY = delta))
  }
  expect_na_row <- function(w, pattern) {
    o <- w$value
    expect_length(w$warnings, 1L)
    expect_match(w$warnings, paste0("^spatialCorrelation: no permutation p-values \\(NA row\\): .*", pattern))
    expect_identical(names(o), c("correlationCoef", "pValueNaive", "pValuePermuteX", "pValuePermuteY",
                                 "deltaStarMedianX", "deltaStarMedianY", "deltaStarX", "deltaStarY",
                                 "nullCorrelationsX", "nullCorrelationsY"))
    expect_identical(c(o$pValuePermuteX, o$pValuePermuteY, o$deltaStarMedianX, o$deltaStarMedianY), rep(NA_real_, 4))
    for (col in c("deltaStarX", "deltaStarY", "nullCorrelationsX", "nullCorrelationsY")) {
      expect_identical(o[[col]], I(list(NA)), info = col)
    }
    o
  }
  local_default_rng()
  # constant X: r and naive p are NA too (cor.test()'s own warning is not repeated)
  o1 <- expect_na_row(run(rep(1, N), 0.3), "permuting X: the permuted values are constant")
  expect_true(all(is.na(c(o1$correlationCoef, o1$pValueNaive))))
  # delta * N < 1 (the R implementation failed in locfit)
  o2 <- expect_na_row(run(inp$X, 0.001), "N \\* delta < 2")
  expect_equal(o2$correlationCoef, stats::cor(inp$X, inp$Y), tolerance = tol)
  # one NA in X: r and the naive p-value use the complete pairs
  Xn <- replace(inp$X, 5, NA)
  o3 <- expect_na_row(run(Xn, 0.3), "contain NA")
  expect_equal(o3$correlationCoef, stats::cor(inp$X[-5], inp$Y[-5]), tolerance = tol)
  expect_equal(o3$pValueNaive, stats::cor.test(inp$X[-5], inp$Y[-5])$p.value, tolerance = 1e-10)
  # 1 <= delta * N < 2 (locfit killed the R session), and too few pairs for cor.test() (it failed with
  # "object 'corDF' not found")
  expect_na_row(run(inp$X, 1.5 / N), "N \\* delta < 2")
  o5 <- expect_na_row(run(c(inp$X[1:2], rep(NA, N - 2)), 0.3), "cor.test\\(\\) failed")
  expect_true(all(is.na(c(o5$correlationCoef, o5$pValueNaive))))
  # every location twice with k = 2 (locfit overflowed the C stack): out of vertex space
  local_default_rng()
  set.seed(9)
  u <- cbind(stats::runif(150), stats::runif(150))
  w <- fx_collect_warnings(spatialCorrelation(stats::rnorm(300), stats::rnorm(300), u[rep(seq_len(150), each = 2), ],
                                              nPermutations = 2, deltaX = c(0.3, 2 / 300), deltaY = 0.3))
  expect_length(w$warnings, 1L)
  expect_match(w$warnings, "out of vertex space")
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
