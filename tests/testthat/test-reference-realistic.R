# Tier 2: realistic regression tests against the published results (fixtures/realistic_fixture.rds, built
# by data-raw/build_test_fixtures.R tier2; see data-raw/README.md). Labels are explained in helper-fixtures.R.
#
# Inputs are real gene subsets: AKI kidney Visium (311 shared pixels, 11 deltas including 0.01 and 0.05) and
# MERFISH vs Visium brain (2170 pixels > 1000, so the variogram uses a 1000-pixel subsample; 9 deltas). The
# golden nulls/deltaStar come from inst/extdata (computed by the authors with 1000 permutations); by the
# prefix property a 10-permutation run reproduces their first 10 values. They reproduce exactly where geoR
# bins the AKI and brain pixel pairs as on the build machine (macOS arm64), not on Linux, so the exact tests
# are gated (helper-platform.R). All tests are skipped on CRAN except those of AKI Upk2 (zero-inflated, small
# deltas, about 2 s), so that tier 2 is not absent from a plain R CMD check.

rf <- fx_read("realistic_fixture.rds")
RF <- "realistic_fixture.rds"
B <- 10
test_genes <- list(aki = c(positive = "Gpx1", negative = "Ech1", null_zero_inflated = "Upk2"),
                   brain = c(positive = "Slc17a7", null = "Efemp1", sparse_Y = "Top2a"))

# spatialCorrelation() of one gene at B = 10, computed once per session.
fx_gene_run <- function(pair, g, Bg = B, flip = FALSE) {
  fx_memo(sprintf("gene:%s:%s:%d:%s", pair, g, Bg, flip), {
    P <- rf$pairs[[pair]]
    Y <- if (flip) max(P$Y[g, ]) - P$Y[g, ] else P$Y[g, ]
    local_default_rng()
    quiet_locfit(spatialCorrelation(P$X[g, ], Y, P$pos, nPermutations = Bg, deltaX = P$params$delta,
                                    deltaY = P$params$delta, maxDistPrctile = P$params$maxDistPrctile,
                                    seed = P$params$seed))
  })
}

test_that("fixture integrity: realistic fixture (not regression coverage)", {
  skip_on_cran()
  for (pair in names(rf$pairs)) {
    P <- rf$pairs[[pair]]
    expect_identical(rownames(P$X), P$genes)
    expect_identical(dim(P$X), dim(P$Y))
    expect_identical(nrow(P$pos), ncol(P$X))
    expect_false(is.null(rf$meta$platform_signature$sets[[pair]]))
    for (g in P$genes) {
      G <- P$golden[[g]]
      expect_length(G$nullX, 100)
      expect_length(G$deltaStarY, 100)
      expect_true(all(c(G$deltaStarX, G$deltaStarY) %in% P$params$delta), info = g)
      expect_identical(G$pRawX_first100, stc_empirical_p(G$nullX, G$correlationCoef), info = g)
    }
  }
  expect_true(all(unlist(lapply(names(test_genes), function(p) test_genes[[p]] %in% rf$pairs[[p]]$genes))))
})

for (pair in names(test_genes)) {
  for (cls in names(test_genes[[pair]])) {
    g <- test_genes[[pair]][[cls]]
    on_cran_too <- identical(g, "Upk2")
    test_that(sprintf("exact: %s %s (%s): first %d nulls and deltaStar reproduce the published results", pair, g, cls, B), {
      if (!on_cran_too) skip_on_cran()
      skip_if_not_exact(RF, pair)
      G <- rf$pairs[[pair]]$golden[[g]]
      o <- fx_gene_run(pair, g)
      expect_equal(as.numeric(o$nullCorrelationsX[[1]]), G$nullX[1:B], tolerance = 1e-12)
      expect_equal(as.numeric(o$nullCorrelationsY[[1]]), G$nullY[1:B], tolerance = 1e-12)
      expect_identical(as.numeric(o$deltaStarX[[1]]), G$deltaStarX[1:B])
      expect_identical(as.numeric(o$deltaStarY[[1]]), G$deltaStarY[1:B])
      expect_equal(o$pValuePermuteX, stc_empirical_p(G$nullX[1:B], G$correlationCoef))
      expect_equal(o$pValuePermuteY, stc_empirical_p(G$nullY[1:B], G$correlationCoef))
    })
    test_that(sprintf("portable: %s %s (%s): r, naive p, p-value definition and loose agreement with the published nulls", pair, g, cls), {
      if (!on_cran_too) skip_on_cran()
      P <- rf$pairs[[pair]]
      G <- P$golden[[g]]
      o <- fx_gene_run(pair, g)
      r <- unname(o$correlationCoef)
      expect_equal(r, G$correlationCoef, tolerance = 1e-12)
      expect_equal(o$pValueNaive, G$pValueNaive, tolerance = 1e-10)
      for (dir in c("X", "Y")) {
        null <- as.numeric(o[[paste0("nullCorrelations", dir)]][[1]])
        expect_length(null, B)
        expect_equal(o[[paste0("pValuePermute", dir)]], stc_empirical_p(null, r), info = dir)
        expect_true(all(as.numeric(o[[paste0("deltaStar", dir)]][[1]]) %in% P$params$delta), info = dir)
        expect_nulls_close(null, G[[paste0("null", dir)]][1:B], info = paste(pair, g, dir))
      }
    })
  }
}

test_that("portable: engineered negative Y -> max(Y) - Y gives r -> -r, nullX -> -nullX and the same deltaStarX (brain Slc17a7)", {
  skip_on_cran()
  g <- rf$pairs$brain$engineered_negatives$genes[1]
  o <- fx_gene_run("brain", g)                       # B = 10
  of <- fx_gene_run("brain", g, Bg = 5, flip = TRUE) # B = 5: also a prefix of the B = 10 run
  nx <- as.numeric(o$nullCorrelationsX[[1]])[1:5]
  expect_equal(unname(of$correlationCoef), -unname(o$correlationCoef), tolerance = 1e-12)
  expect_equal(as.numeric(of$nullCorrelationsX[[1]]), -nx, tolerance = 1e-12)
  expect_identical(as.numeric(of$deltaStarX[[1]]), as.numeric(o$deltaStarX[[1]])[1:5])
  expect_equal(of$pValuePermuteX, stc_empirical_p(nx, unname(o$correlationCoef)))
  # permuting Yflip instead of Y changes the null fields, so pValuePermuteY is only statistically equal
})

test_that("exact: engineered negative reproduces -1 x the published nullX (brain Slc17a7, B = 5)", {
  skip_on_cran()
  skip_if_not_exact(RF, "brain")
  g <- rf$pairs$brain$engineered_negatives$genes[1]
  G <- rf$pairs$brain$golden[[g]]
  of <- fx_gene_run("brain", g, Bg = 5, flip = TRUE)
  expect_equal(unname(of$correlationCoef), -G$correlationCoef, tolerance = 1e-12)
  expect_equal(as.numeric(of$nullCorrelationsX[[1]]), -G$nullX[1:5], tolerance = 1e-12)
  expect_identical(as.numeric(of$deltaStarX[[1]]), G$deltaStarX[1:5])
  expect_equal(of$pValuePermuteX, stc_empirical_p(G$nullX[1:5], G$correlationCoef))
})

# Slow (STCOMPARE_SLOW_TESTS=true): every gene of the fixture at B = 100 against the 100 stored published
# nulls, and the 5 engineered negatives at B = 20 (the same comparison as data-raw/build_test_fixtures.R
# verify). About 2500 CPU-seconds; parallel over genes on STCOMPARE_TEST_WORKERS processes.
for (pair in names(rf$pairs)) {
  test_that(sprintf("exact, slow: %s: all %d genes at B = 100 reproduce the published first 100 nulls and deltaStar", pair, length(rf$pairs[[pair]]$genes)), {
    skip_if_not_slow()
    skip_if_not_exact(RF, pair)
    P <- rf$pairs[[pair]]
    jobs <- data.frame(gene = P$genes, flip = FALSE, Bj = 100L)
    if (!is.null(P$engineered_negatives)) jobs <- rbind(jobs, data.frame(gene = P$engineered_negatives$genes, flip = TRUE, Bj = 20L))
    run1 <- function(k) {
      g <- jobs$gene[k]
      G <- P$golden[[g]]
      flip <- jobs$flip[k]
      Bj <- jobs$Bj[k]
      Y <- if (flip) max(P$Y[g, ]) - P$Y[g, ] else P$Y[g, ]
      local_default_rng()
      o <- quiet_locfit(spatialCorrelation(P$X[g, ], Y, P$pos, nPermutations = Bj, deltaX = P$params$delta,
                                           deltaY = P$params$delta, maxDistPrctile = P$params$maxDistPrctile,
                                           BPPARAM = BiocParallel::SerialParam(), seed = P$params$seed))
      s <- if (flip) -1 else 1
      data.frame(gene = g, flip = flip, B = Bj,
                 dr = abs(unname(o$correlationCoef) - s * G$correlationCoef),
                 dnullX = max(abs(as.numeric(o$nullCorrelationsX[[1]]) - s * G$nullX[1:Bj])),
                 dnullY = if (flip) 0 else max(abs(as.numeric(o$nullCorrelationsY[[1]]) - G$nullY[1:Bj])),
                 deltaStarX = identical(as.numeric(o$deltaStarX[[1]]), G$deltaStarX[1:Bj]),
                 deltaStarY = flip || identical(as.numeric(o$deltaStarY[[1]]), G$deltaStarY[1:Bj]))
    }
    t0 <- Sys.time()
    res <- do.call(rbind, BiocParallel::bplapply(seq_len(nrow(jobs)), run1, BPPARAM = fx_bpparam()))
    message(sprintf("%s: %d jobs on %d worker(s) in %.0f s; max |dnull| = %.2g", pair, nrow(res), fx_workers(),
                    as.numeric(difftime(Sys.time(), t0, units = "secs")), max(res$dnullX, res$dnullY)))
    expect_true(all(res$dr <= 1e-12), info = paste(res$gene[res$dr > 1e-12], collapse = " "))
    bad <- res$gene[res$dnullX > 1e-12 | res$dnullY > 1e-12 | !res$deltaStarX | !res$deltaStarY]
    expect_identical(bad, character(0))
  })
}
