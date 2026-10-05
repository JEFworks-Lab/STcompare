# Tier 1: statistical calibration and power (fixtures/calibration_fixture.rds, built by
# data-raw/build_test_fixtures.R tier1; see data-raw/README.md). Labels are explained in helper-fixtures.R.
# Skipped unless STCOMPARE_SLOW_TESTS=true. 60 calls of spatialCorrelation() with B = 100 (~7 s each on one
# core); they run in parallel over pairs on STCOMPARE_TEST_WORKERS processes (default
# BiocParallel::multicoreWorkers()), and the results do not depend on the number of workers.
#
# Null pairs: 40 disjoint pairs of the 100 independent simulated fields of data(simRanPatternRasts), so the
# 40 p-values are independent. Power pairs: the first 20 of them with Y = stc_mix(f_i, f_j, rho = 0.6), a field
# with the same covariance model and population correlation 0.6 with X = f_i.

CF <- "calibration_fixture.rds"

calibration_results <- function() {
  fx_memo("calibration_results", {
    cf <- fx_read(CF)
    jobs <- cf$test_jobs
    B <- cf$reference_test_jobs$B
    run1 <- function(k) {
      i <- jobs$i[k]
      j <- jobs$j[k]
      rho <- jobs$rho[k]
      sh <- which(!is.na(cf$fields[i, ]) & !is.na(cf$fields[j, ]))
      X <- unname(cf$fields[i, sh])
      Y <- if (rho == 0) unname(cf$fields[j, sh]) else stc_mix(X, unname(cf$fields[j, sh]), rho, cf$mix_mu)
      # BiocParallel workers run L'Ecuyer-CMRG; the permutation stream is defined for R's default generators
      local_default_rng()
      o <- quiet_locfit(spatialCorrelation(X, Y, unname(cf$coords[sh, , drop = FALSE]), nPermutations = B,
                                           BPPARAM = BiocParallel::SerialParam(), seed = 0))
      list(row = data.frame(i = i, j = j, rho = rho, N = length(sh), r = unname(o$correlationCoef), pNaive = o$pValueNaive,
                            pX = o$pValuePermuteX, pY = o$pValuePermuteY,
                            nExtremeX = sum(abs(o$nullCorrelationsX[[1]]) > abs(o$correlationCoef)),
                            nExtremeY = sum(abs(o$nullCorrelationsY[[1]]) > abs(o$correlationCoef)), B = B,
                            deltaStarMedianX = o$deltaStarMedianX, deltaStarMedianY = o$deltaStarMedianY),
           deltaStarX = as.numeric(o$deltaStarX[[1]]), deltaStarY = as.numeric(o$deltaStarY[[1]]))
    }
    t0 <- Sys.time()
    out <- BiocParallel::bplapply(seq_len(nrow(jobs)), run1, BPPARAM = fx_bpparam())
    res <- do.call(rbind, lapply(out, `[[`, "row"))
    message(sprintf("calibration: %d spatialCorrelation() calls (B = %d) on %d worker(s) in %.1f s",
                    nrow(res), B, fx_workers(), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    for (rho in unique(res$rho)) {
      d <- res[res$rho == rho, ]
      message(sprintf("  rho = %.1f (n = %d): mean r = %.3f; naive p < 0.05: %.3f; pX < 0.05: %.3f; pY < 0.05: %.3f; max(pX, pY) < 0.05: %.3f",
                      rho, nrow(d), mean(d$r), mean(d$pNaive < 0.05), mean(d$pX < 0.05), mean(d$pY < 0.05),
                      mean(pmax(d$pX, d$pY) < 0.05)))
    }
    list(results = res, deltaStarX = sapply(out, `[[`, "deltaStarX"), deltaStarY = sapply(out, `[[`, "deltaStarY"))
  })
}

test_that("fixture integrity: the helper's stc_mix() is the recipe the calibration fixture was built with", {
  cf <- fx_read(CF)
  chk <- cf$mix_check
  expect_false(any(vapply(cf, is.function, NA)))
  expect_equal(stc_mix(chk$fi, chk$fj, chk$rho, cf$mix_mu), chk$value, tolerance = 1e-15)
})

test_that("portable: null pairs: rejection rate at alpha = 0.05 is not above a one-sided 99.9% binomial bound", {
  skip_if_not_slow()
  res <- calibration_results()$results
  d <- res[res$rho == 0, ]
  bound <- stats::qbinom(0.999, nrow(d), 0.05)
  expect_lte(sum(d$pX < 0.05), bound)
  expect_lte(sum(d$pY < 0.05), bound)
  # the naive test is badly anti-conservative on these autocorrelated fields (positive control for the test)
  expect_gt(mean(d$pNaive < 0.05), 0.05)
})

test_that("portable: power: mix(rho = 0.6) pairs are mostly significant in both directions", {
  skip_if_not_slow()
  res <- calibration_results()$results
  d <- res[res$rho == 0.6, ]
  expect_gt(mean(d$r), 0.4)
  expect_gte(mean(pmax(d$pX, d$pY) < 0.05), 0.7)
})

test_that("portable: p-values agree in distribution with the stored reference and the bandwidth selection is not degenerate", {
  # Backend-agnostic acceptance check for a statistically (not bit-) equivalent backend. Paired over the
  # 60 jobs: no systematic shift of the p-values (paired Wilcoxon test, alpha = 0.001; independent Monte
  # Carlo noise alone gives mean |dp| of about 0.05 at B = 100), and the selected deltas use the grid like
  # the legacy code (share of permutations at the smallest and at the largest delta within 0.25 of the
  # reference's, at least 3 distinct values).
  skip_if_not_slow()
  cr <- calibration_results()
  ref <- fx_read(CF)$reference_test_jobs
  expect_identical(cr$results[, c("i", "j", "rho")], ref$results[, c("i", "j", "rho")])
  grid <- seq(0.1, 0.9, 0.1)
  for (dir in c("X", "Y")) {
    p_new <- cr$results[[paste0("p", dir)]]
    # the reference stores tail counts; convert them with the package's p-value definition
    p_ref <- stc_p_from_count(ref$results[[paste0("nExtreme", dir)]], ref$results$B)
    dp <- p_new - p_ref
    expect_lte(mean(abs(dp)), 0.1)
    if (any(dp != 0)) {
      expect_gte(stats::wilcox.test(p_new, p_ref, paired = TRUE, exact = FALSE)$p.value, 0.001)
    }
    ds_new <- cr[[paste0("deltaStar", dir)]]
    ds_ref <- ref[[paste0("deltaStar", dir)]]
    expect_identical(dim(ds_new), dim(ds_ref))
    expect_lte(abs(mean(ds_new == min(grid)) - mean(ds_ref == min(grid))), 0.25)
    expect_lte(abs(mean(ds_new == max(grid)) - mean(ds_ref == max(grid))), 0.25)
    expect_gte(length(unique(as.vector(ds_new))), 3)
  }
})

test_that("exact (legacy backend): calibration results equal the stored reference", {
  # Only jobs whose pixel coordinates geoR bins as on the build machine are compared (all of them with
  # STCOMPARE_EXACT_TESTS=true); a statistically equivalent backend is judged by the test above instead.
  skip_if_not_slow()
  cr <- calibration_results()
  ref <- fx_read(CF)$reference_test_jobs
  k <- seq_len(nrow(ref$results))
  if (fx_exact_mode() == "false") skip("exact reference checks disabled (STCOMPARE_EXACT_TESTS=false)")
  if (fx_exact_mode() == "auto") k <- k[vapply(k, function(j) fx_exact_status(CF, sprintf("job%02d", j))$exact, NA)]
  if (!length(k)) skip("geoR bins every calibration coordinate set differently on this machine than on the build machine")
  message(sprintf("exact calibration comparison on %d of %d jobs", length(k), nrow(ref$results)))
  res <- cr$results
  expect_identical(res[k, c("i", "j", "rho", "N")], ref$results[k, c("i", "j", "rho", "N")])
  expect_equal(res$r[k], ref$results$r[k], tolerance = 1e-12)
  expect_identical(res$nExtremeX[k], ref$results$nExtremeX[k])
  expect_identical(res$nExtremeY[k], ref$results$nExtremeY[k])
  expect_identical(cr$deltaStarX[, k], ref$deltaStarX[, k])
  expect_identical(cr$deltaStarY[, k], ref$deltaStarY[, k])
})
