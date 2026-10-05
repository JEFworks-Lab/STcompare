# Tier 2: realistic regression tests against the published results (fixtures/realistic_fixture.rds, built
# by data-raw/build_test_fixtures.R tier2; see data-raw/README.md). Labels are explained in helper-fixtures.R.
#
# Inputs are real gene subsets: AKI kidney Visium (35 genes on 311 shared pixels, 11 deltas including 0.01 and
# 0.05) and MERFISH vs Visium brain (30 genes on 2170 pixels > 1000, so the variogram uses a 1000-pixel
# subsample; 9 deltas). The golden nulls and deltaStar are the first 100 of the published results in
# bench/published (computed by the authors with the original R implementation and 1000 permutations); by the
# prefix property a B-permutation run reproduces their first B values. They reproduce exactly where the AKI and
# brain pixel pairs are binned as on the build machine (macOS arm64), not on Linux, so the exact tests are
# gated (helper-platform.R).
#
# Default suite: all 35 AKI genes at B = 100 in one spatialCorrelationGeneExp() call, and 4 brain genes at
# B = 25. Slow tier (STCOMPARE_SLOW_TESTS=true): all 30 brain genes at B = 100, the 5 engineered negatives, and
# the published iterative protocol (100 then 1000 permutations) on the 35 AKI genes.

rf <- fx_read("realistic_fixture.rds")
RF <- "realistic_fixture.rds"
aki <- rf$pairs$aki
brain <- rf$pairs$brain
test_genes <- list(aki = c(positive = "Gpx1", negative = "Ech1", null_zero_inflated = "Upk2"),
                   brain = c(positive = "Slc17a7", null = "Efemp1", sparse_Y = "Top2a"))
brain_default <- c("Slc17a7", "Efemp1", "Top2a", "Gpr34")
B_brain <- 25L

# The published values of genes of one pair against a spatialCorrelationGeneExp() result `o` (unadjusted
# p-values) with B permutations: nulls 1e-9 relative, deltaStar identical, p-values from the published nulls.
expect_published <- function(o, P, genes, B, sign = 1) {
  for (g in genes) {
    G <- P$golden[[g]]
    for (dir in c("X", "Y")) {
      info <- paste(g, dir)
      expect_close(as.numeric(o[g, paste0("nullCorrelations", dir)][[1]]), sign * G[[paste0("null", dir)]][seq_len(B)],
                   info = info)
      expect_identical(o[g, paste0("deltaStar", dir)][[1]], G[[paste0("deltaStar", dir)]][seq_len(B)], info = info)
      expect_identical(o[g, paste0("pValuePermute", dir)],
                       stc_empirical_p(G[[paste0("null", dir)]][seq_len(B)], G$correlationCoef), info = info)
    }
  }
}

test_that("fixture integrity: realistic fixture (not regression coverage)", {
  for (pair in names(rf$pairs)) {
    P <- rf$pairs[[pair]]
    expect_identical(rownames(P$X), P$genes)
    expect_identical(dim(P$X), dim(P$Y))
    expect_identical(nrow(P$pos), ncol(P$X))
    expect_identical(length(P$pixel), ncol(P$X))
    expect_false(is.null(rf$meta$platform_signature$sets[[pair]]))
    for (g in P$genes) {
      G <- P$golden[[g]]
      expect_length(G$nullX, 100)
      expect_length(G$deltaStarY, 100)
      expect_true(all(c(G$deltaStarX, G$deltaStarY) %in% P$params$delta), info = g)
      expect_identical(G$pRawX_first100, stc_legacy_p(G$nullX, G$correlationCoef), info = g)
    }
  }
  expect_true(all(unlist(lapply(names(test_genes), function(p) test_genes[[p]] %in% rf$pairs[[p]]$genes))))
})

test_that("exact: AKI: all 35 genes at B = 100 reproduce the published first 100 nulls and deltaStar (one spatialCorrelationGeneExp() call)", {
  skip_if_not_exact(RF, "aki")
  o <- fx_genes_run("aki", aki$genes, 100L)
  expect_identical(rownames(o), aki$genes)
  expect_published(o, aki, aki$genes, 100L)
})

test_that(sprintf("exact: brain: %d genes at B = %d reproduce the published first nulls and deltaStar (N = 2170, variogram subsample)",
                  length(brain_default), B_brain), {
  skip_if_not_exact(RF, "brain")
  o <- fx_genes_run("brain", brain_default, B_brain)
  expect_published(o, brain, brain_default, B_brain)
})

test_that("portable: AKI and brain test genes: r, naive p, the p-value definition and loose agreement with the published nulls", {
  runs <- list(aki = fx_genes_run("aki", aki$genes, 100L), brain = fx_genes_run("brain", brain_default, B_brain))
  for (pair in names(test_genes)) {
    P <- rf$pairs[[pair]]
    o <- runs[[pair]]
    B <- if (pair == "aki") 100L else B_brain
    for (g in test_genes[[pair]]) {
      G <- P$golden[[g]]
      r <- o[g, "correlationCoef"]
      expect_equal(r, G$correlationCoef, tolerance = 1e-12, info = g)
      expect_equal(o[g, "pValueNaive"], G$pValueNaive, tolerance = 1e-10, info = g)
      for (dir in c("X", "Y")) {
        info <- paste(pair, g, dir)
        null <- as.numeric(o[g, paste0("nullCorrelations", dir)][[1]])
        expect_length(null, B)
        expect_identical(o[g, paste0("pValuePermute", dir)], stc_empirical_p(null, r), info = info)
        expect_true(all(o[g, paste0("deltaStar", dir)][[1]] %in% P$params$delta), info = info)
        expect_nulls_close(null, G[[paste0("null", dir)]][seq_len(B)], info = info)
      }
    }
  }
})

test_that("portable: engineered negative Y -> max(Y) - Y gives r -> -r, nullX -> -nullX and the same deltaStarX (brain Slc17a7)", {
  g <- brain$engineered_negatives$genes[1]
  o <- fx_genes_run("brain", brain_default, B_brain)
  of <- fx_genes_run("brain", g, 5L, flip = TRUE)   # B = 5: a prefix of the B = 25 run
  nx <- as.numeric(o[g, "nullCorrelationsX"][[1]])[1:5]
  expect_equal(of[g, "correlationCoef"], -o[g, "correlationCoef"], tolerance = 1e-12)
  expect_equal(as.numeric(of[g, "nullCorrelationsX"][[1]]), -nx, tolerance = 1e-12)
  expect_identical(of[g, "deltaStarX"][[1]], o[g, "deltaStarX"][[1]][1:5])
  expect_equal(of[g, "pValuePermuteX"], stc_empirical_p(nx, o[g, "correlationCoef"]))
  # permuting Yflip instead of Y changes the null fields, so pValuePermuteY is only statistically equal
})

test_that("exact: engineered negative reproduces -1 x the published nullX (brain Slc17a7, B = 5)", {
  skip_if_not_exact(RF, "brain")
  g <- brain$engineered_negatives$genes[1]
  G <- brain$golden[[g]]
  of <- fx_genes_run("brain", g, 5L, flip = TRUE)
  expect_equal(of[g, "correlationCoef"], -G$correlationCoef, tolerance = 1e-12)
  expect_close(as.numeric(of[g, "nullCorrelationsX"][[1]]), -G$nullX[1:5])
  expect_identical(of[g, "deltaStarX"][[1]], G$deltaStarX[1:5])
  expect_equal(of[g, "pValuePermuteX"], stc_empirical_p(G$nullX[1:5], G$correlationCoef))
})

# Slow (STCOMPARE_SLOW_TESTS=true): every brain gene at B = 100 against the 100 stored published nulls, and the
# 5 engineered negatives at B = 20 (the comparison of data-raw/build_test_fixtures.R verify).
test_that("exact, slow: brain: all 30 genes at B = 100 and the 5 engineered negatives at B = 20 reproduce the published results", {
  skip_if_not_slow()
  skip_if_not_exact(RF, "brain")
  t0 <- proc.time()
  o <- fx_genes_run("brain", brain$genes, 100L)
  expect_identical(rownames(o), brain$genes)
  expect_published(o, brain, brain$genes, 100L)
  neg <- brain$engineered_negatives$genes
  of <- fx_genes_run("brain", neg, 20L, flip = TRUE)
  for (g in neg) {
    G <- brain$golden[[g]]
    expect_equal(of[g, "correlationCoef"], -G$correlationCoef, tolerance = 1e-12, info = g)
    expect_close(as.numeric(of[g, "nullCorrelationsX"][[1]]), -G$nullX[1:20], info = g)
    expect_identical(of[g, "deltaStarX"][[1]], G$deltaStarX[1:20], info = g)
  }
  message(sprintf("brain: 30 genes at B = 100 and 5 engineered negatives at B = 20 on %d threads in %.1f s",
                  fx_threads, (proc.time() - t0)[["elapsed"]]))
})

# Slow: the published protocol of the kidney analysis (spatialCorrelationGeneExpIterPermutations(), 100 then
# 1000 permutations, alpha = 0.05) on the 35 AKI genes: the same screening decisions as the published table
# (permutations per gene), and the published first 100 nulls and deltaStar of every gene.
test_that("exact, slow: AKI: spatialCorrelationGeneExpIterPermutations() makes the published screening decisions (100, then 1000 permutations)", {
  skip_if_not_slow()
  skip_if_not_exact(RF, "aki")
  delta <- rep(list(aki$params$delta), length(aki$genes))
  local_default_rng()
  o <- spatialCorrelationGeneExpIterPermutations(fx_spe_pair(aki), nPermutations = aki$params$nPermutations_published,
                                                 alpha = aki$params$alpha, deltaX = delta, deltaY = delta,
                                                 seed = aki$params$seed, nThreads = fx_threads, verbose = FALSE)
  expect_identical(unname(lengths(o$nullCorrelationsX)),
                   as.integer(aki$selection$nperm[match(aki$genes, aki$selection$gene)]))
  expect_true(any(lengths(o$nullCorrelationsX) == 1000L) && any(lengths(o$nullCorrelationsX) == 100L))
  for (g in aki$genes) {
    G <- aki$golden[[g]]
    expect_close(as.numeric(o[g, "nullCorrelationsX"][[1]])[1:100], G$nullX, info = g)
    expect_close(as.numeric(o[g, "nullCorrelationsY"][[1]])[1:100], G$nullY, info = g)
    expect_identical(o[g, "deltaStarX"][[1]][1:100], G$deltaStarX, info = g)
    expect_identical(o[g, "deltaStarY"][[1]][1:100], G$deltaStarY, info = g)
  }
})
