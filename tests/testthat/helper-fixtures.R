# Helpers for the reference tests. testthat (and devtools::load_all()) sources this file before the tests;
# it only defines functions. The fixtures in tests/testthat/fixtures/ are built from the legacy R
# implementation by data-raw/build_test_fixtures.R (see data-raw/README.md for what each tier pins down).
#
# Test labels used in the test files:
#   "exact: ..."    compare with the stored legacy values at the documented tolerances (1e-12; integers, bin
#                   counts and deltaStar identical; RNG draws 1e-14). They run only for coordinate sets that
#                   geoR bins exactly as on the build machine (see helper-platform.R) and are skipped
#                   elsewhere, naming the differing probes. STCOMPARE_EXACT_TESTS=true forces them on every
#                   machine, STCOMPARE_EXACT_TESTS=false skips them.
#   "portable: ..." hold on every platform: STcompare against the reference kernels run on the same machine,
#                   invariants (prefix property, X/Y symmetry, worker independence, the p-value definition,
#                   error paths), hand-computed spatialSimilarity() values, and loose agreement with the
#                   fixtures.
#   "canary: ..."   dependency canaries and C++ kernel slots: R's RNG, SEraster, geoR, locfit and lm (through
#                   the ref_*() helpers) against the fixture. They call no STcompare code.
#   "fixture integrity: ..." checks of the fixture files themselves; they are not regression coverage.

fx_cache <- new.env(parent = emptyenv())

# Read a fixture once per session.
fx_read <- function(name) {
  if (is.null(fx_cache[[name]])) {
    path <- testthat::test_path("fixtures", name)
    if (!file.exists(path)) stop("Missing test fixture ", path, "; rebuild it with data-raw/build_test_fixtures.R")
    fx_cache[[name]] <- readRDS(path)
  }
  fx_cache[[name]]
}

# Evaluate `expr` once per session and keep the value under `key`.
fx_memo <- function(key, expr) {
  if (!exists(key, envir = fx_cache, inherits = FALSE)) assign(key, expr, envir = fx_cache)
  get(key, envir = fx_cache, inherits = FALSE)
}

# Slow tests (statistical calibration, all realistic genes) run only when STCOMPARE_SLOW_TESTS=true.
skip_if_not_slow <- function() {
  testthat::skip_if_not(identical(tolower(Sys.getenv("STCOMPARE_SLOW_TESTS")), "true"),
                        "slow test; set STCOMPARE_SLOW_TESTS=true to run it")
}

# Worker processes for the slow tests (STCOMPARE_TEST_WORKERS, default BiocParallel::multicoreWorkers()).
fx_workers <- function() {
  w <- Sys.getenv("STCOMPARE_TEST_WORKERS", "")
  if (nzchar(w)) max(1L, as.integer(w)) else max(1L, BiocParallel::multicoreWorkers())
}

fx_bpparam <- function() {
  w <- fx_workers()
  if (w > 1 && .Platform$OS.type != "windows") BiocParallel::MulticoreParam(workers = w) else BiocParallel::SerialParam()
}

# --- RNG ---------------------------------------------------------------------------------------------------
# The caller's RNG kind and seed are put back when `env` (a test or a function frame) exits.
fx_defer_rng_restore <- function(env) {
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  restore <- function() {
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE), envir = env)
  invisible(NULL)
}

# R's default generators (the permutation stream is defined for them) until the calling test or function
# exits; then the caller's RNG kind and seed are restored.
local_default_rng <- function(env = parent.frame()) {
  fx_defer_rng_restore(env)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  invisible(NULL)
}

# Evaluate `code` under L'Ecuyer-CMRG (Inversion, Rejection) seeded with `seed`: the generator a BiocParallel
# task uses when matchingVariograms() calls set.seed(seed + i). The caller's RNG kind and seed are restored.
with_lecuyer <- function(seed, code) {
  fx_defer_rng_restore(environment())
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  code
}

# --- Warnings ----------------------------------------------------------------------------------------------
# Muffle only locfit's expected "Estimated rdf < 1.0; not estimating variance" (small deltas); any other
# warning reaches testthat and is reported.
quiet_locfit <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Estimated rdf < 1.0", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

# --- Exactness gate ----------------------------------------------------------------------------------------
fx_signature <- function(name) {
  fx_memo(paste0("signature:", name), stc_platform_signature(stc_signature_sets(name, fx_read(name))))
}

# Can coordinate set `set` of fixture `name` be compared exactly on this machine?
fx_exact_status <- function(name, set) {
  ref <- fx_read(name)$meta$platform_signature
  if (is.null(ref)) return(list(exact = FALSE, diff = "fixture has no platform signature", global = character(0), built_on = "?"))
  cur <- fx_signature(name)
  d <- stc_set_diff(ref, cur, set)
  list(exact = !length(d), diff = d, global = stc_global_diff(ref, cur),
       built_on = sprintf("%s, %s", ref$info$platform, ref$info$R))
}

fx_exact_mode <- function() {
  m <- tolower(Sys.getenv("STCOMPARE_EXACT_TESTS", "auto"))
  if (m %in% c("true", "false")) m else "auto"
}

skip_if_not_exact <- function(name, set) {
  mode <- fx_exact_mode()
  if (mode == "true") return(invisible(TRUE))
  if (mode == "false") testthat::skip("exact reference checks disabled (STCOMPARE_EXACT_TESTS=false)")
  st <- fx_exact_status(name, set)
  if (!st$exact) {
    testthat::skip(sprintf(paste0("%s [%s] was built on %s; geoR bins this coordinate set differently here ",
                                  "(differs: %s; global probes differing: %s). Portable checks still run; ",
                                  "STCOMPARE_EXACT_TESTS=true forces the exact checks"),
                           name, set, st$built_on, paste(st$diff, collapse = ", "),
                           if (length(st$global)) paste(st$global, collapse = ", ") else "none"))
  }
  invisible(TRUE)
}

# Loose agreement with stored null correlations, for machines where geoR bins differently: at least 60%
# (rounded down, at least one) of the nulls within 0.01 of the stored ones. A different bin count moves every
# null by about 1e-4, and the few permutations whose deltaStar changes by up to about 0.2; a wrong seed,
# permutation order, delta grid or smoother moves nearly all of them by more.
expect_nulls_close <- function(null, ref, info = NULL, tol = 0.01, frac = 0.6) {
  k <- sum(abs(null - ref) <= tol)
  need <- max(1L, floor(frac * length(ref)))
  testthat::expect(length(null) == length(ref) && k >= need,
                   sprintf("%d of %d null correlations within %g of the stored values (need %d)%s", k, length(ref),
                           tol, need, if (is.null(info)) "" else paste0(" [", info, "]")))
  invisible(null)
}

# --- Definitions shared by tests and fixtures --------------------------------------------------------------
# Empirical p-value as defined by the package: (b + 1) / (B + 1), where b counts the null correlations
# whose absolute value is at least |r|. Tests compare package p-values with this function applied to raw
# nulls; the fixtures keep raw nulls, so changing the definition needs no rebuild.
stc_empirical_p <- function(null, r) (sum(abs(null) >= abs(r)) + 1) / (length(null) + 1)
stc_p_from_count <- function(b, B) (b + 1) / (B + 1)

# The legacy definition (b / B with strict ">", so it can be exactly 0). The fixtures were built with it, so
# fixture-integrity checks of stored p-values use this; package outputs are checked with stc_empirical_p().
stc_legacy_p <- function(null, r) sum(abs(null) > abs(r)) / length(null)

# Positive-control recipe of the calibration tier (calibration_fixture.rds$mix_recipe): a field with the same
# covariance model as f_i and f_j whose population correlation with f_i is rho.
stc_mix <- function(fi, fj, rho, mu = 10) rho * (fi - mu) + sqrt(1 - rho^2) * (fj - mu) + mu

# Column fingerprint of an N x B matrix of permuted fields (the fixtures store this, not the fields).
perm_fingerprint <- function(P) {
  rbind(mean = colMeans(P), sd = apply(P, 2, stats::sd), min = apply(P, 2, min), max = apply(P, 2, max))
}

fx_default_delta <- seq(0.1, 0.9, 0.1)

# spatialCorrelation() on a kernel case. Arguments equal to the documented defaults (deltaX = deltaY =
# seq(0.1, 0.9, 0.1), maxDistPrctile = 0.25, seed = 0) are omitted so the defaults themselves are exercised;
# the others are passed explicitly.
fx_sc_case <- function(cs, B = cs$input$nPermutations, ...) {
  inp <- cs$input
  args <- list(inp$X, inp$Y, inp$pos, nPermutations = B, ...)
  if (!identical(inp$deltaX, fx_default_delta)) args$deltaX <- inp$deltaX
  if (!identical(inp$deltaY, fx_default_delta)) args$deltaY <- inp$deltaY
  if (!identical(inp$maxDistPrctile, 0.25)) args$maxDistPrctile <- inp$maxDistPrctile
  if (!identical(inp$seed, 0)) args$seed <- inp$seed
  local_default_rng()
  quiet_locfit(do.call(spatialCorrelation, args))
}

# The returnPermutations = TRUE run of a kernel case, computed once per session.
fx_case_run <- function(cs) fx_memo(paste0("sc:", cs$name), fx_sc_case(cs, returnPermutations = TRUE))

# data(speKidney) rasterized exactly as the tier-0 fixture (computed once per session).
fx_speKidney_raster <- function() {
  fx_memo("speKidney_raster", SEraster::rasterizeGeneExpression(
    STcompare::speKidney, assay_name = "counts", resolution = 0.2, fun = "mean", square = FALSE))
}

# --- Reference kernels -------------------------------------------------------------------------------------
# One function per step of viladomatCorrelation() / matchingVariograms(), written exactly as the package
# calls them (lat = pos[, 1], long = pos[, 2]; variogram coordinates are cbind(long, lat)). A compiled backend
# is validated step by step by swapping its kernel in here; the tests compare every step with the fixture
# (on the build machine) and with STcompare itself (on every machine).

# max.dist: quantile of the pairwise distances of the (subsampled) points; note dist(cbind(lat, long)).
ref_prctile <- function(lat, long, p) unname(stats::quantile(stats::dist(cbind(lat, long)), probs = p))

ref_variog <- function(z, long, lat, max_dist) {
  geoR::variog(data = z, coords = cbind(long, lat), max.dist = max_dist, option = "bin", messages = FALSE)
}

# locfit local-constant Gaussian-kernel smoother with nearest-neighbour bandwidth nn = delta. The package
# uses locfit's default adaptive kd-tree with interpolation; exact = TRUE evaluates at every point
# (ev = dat()), the reference for a smoother that is not bit-compatible with the tree.
ref_smooth <- function(x, long, lat, delta, exact = FALSE) {
  fit <- if (exact) {
    quiet_locfit(locfit::locfit(x ~ locfit::lp(long, lat, nn = delta, deg = 0), kern = "gauss", maxk = 300,
                                ev = locfit::dat()))
  } else {
    quiet_locfit(locfit::locfit(x ~ locfit::lp(long, lat, nn = delta, deg = 0), kern = "gauss", maxk = 300))
  }
  as.numeric(stats::fitted(fit))
}

# Least squares of the target variogram on the candidate variogram: c(intercept, slope).
ref_lm <- function(target_v, v) as.numeric(stats::lm(target_v ~ 1 + v)$coefficients)

# Noise of permutation i: n_delta draws of rnorm(N) after set.seed(seed + i) under L'Ecuyer-CMRG (N x n_delta).
ref_noise <- function(N, seed_i, n_delta) with_lecuyer(seed_i, vapply(seq_len(n_delta), function(k) stats::rnorm(N), numeric(N)))

ref_rescale <- function(xd, noise, bet) xd * sqrt(abs(bet[2])) + noise * sqrt(abs(bet[1]))

# Permutation 1 of a kernel case and direction ("forward" permutes X, "reverse" permutes Y), replayed step by
# step with the reference kernels on this machine (computed once per session). Uses the stored permutation
# indices and variogram subsample, which the RNG canary checks on every platform.
fx_replay <- function(cs, dir) fx_memo(paste0("replay:", cs$name, ":", dir), {
  cap <- cs$intermediate[[dir]]
  d <- cap$detail
  z <- if (dir == "forward") cs$input$X else cs$input$Y
  delta <- if (dir == "forward") cs$input$deltaX else cs$input$deltaY
  lat <- cs$input$pos[, 1]
  long <- cs$input$pos[, 2]
  ids <- cap$ids
  prctile <- ref_prctile(lat[ids], long[ids], cs$input$maxDistPrctile)          # C++: distance quantile
  tv <- ref_variog(z[ids], long[ids], lat[ids], prctile)                         # C++: variogram
  xr <- z[cap$perm_index[, d$i]]
  noise <- ref_noise(length(z), d$noise_seed, length(delta))                    # C++: RNG (or injected)
  per <- lapply(seq_along(delta), function(k) {
    xd <- ref_smooth(xr, long, lat, delta[k])                                    # C++: smoother
    v1 <- ref_variog(xd[ids], long[ids], lat[ids], prctile)
    bet <- ref_lm(tv$v, v1$v)                                                    # C++: least squares
    hat <- ref_rescale(xd, noise[, k], bet)
    v2 <- ref_variog(hat[ids], long[ids], lat[ids], prctile)
    list(xd = xd, n1 = v1$n, v1 = v1$v, bet = bet, hat = hat, n2 = v2$n, v2 = v2$v, rss = sum((v2$v - tv$v)^2))
  })
  rss <- vapply(per, function(p) p$rss, 0)
  list(prctile = prctile, target = tv, xr = xr, noise = noise, delta = delta, per = per, rss = rss,
       argmin = which.min(rss), hat_star = per[[which.min(rss)]]$hat)
})

# STcompare's matchingVariograms() on the same permutation, target variogram and noise seed as fx_replay(),
# under L'Ecuyer-CMRG as inside viladomatCorrelation() (computed once per session).
fx_mv <- function(cs, dir) fx_memo(paste0("mv:", cs$name, ":", dir), {
  rp <- fx_replay(cs, dir)
  cap <- cs$intermediate[[dir]]
  with_lecuyer(cap$detail$noise_seed, quiet_locfit(matchingVariograms(
    rp$xr, cs$input$pos[, 2], cs$input$pos[, 1], rp$delta, rp$target, rp$prctile, cap$ids, cap$detail$i,
    seed = cap$detail$noise_seed)))
})
