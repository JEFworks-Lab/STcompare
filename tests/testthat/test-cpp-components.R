# Component tests of the compiled engine's building blocks (src/stc_*.cpp, R/engine.R) against the R
# code they replace: geoR::variog(), fitted(locfit()), rnorm() under L'Ecuyer-CMRG, cor() and lm()
# (dev/engine-spec.md, section 5). No exported function uses these blocks yet. Labels ("exact",
# "portable") are explained in helper-fixtures.R.
#
# Coordinate conventions: lat = pos[, 1] and long = pos[, 2]. The variogram takes geoR's column order
# (long, lat) and the smoother locfit's order lp(long, lat), so x1 = long and x2 = lat in both.

fx <- fx_read("kernel_fixture.rds")
KF <- "kernel_fixture.rds"
cpp_sets <- c("kidney_AB", "kidney_AB_jitter", "quakes_irregular", "brain_Oprk1_subsample")

# --- helpers ------------------------------------------------------------------------------------

# geoR::variog() as STcompare calls it. capture.output() silences geoR's "co-locatted data" message,
# which geoR prints with cat() even with messages = FALSE.
geoR_variog <- function(z, long, lat, max_dist) {
  out <- NULL
  utils::capture.output(out <- ref_variog(z, long, lat, max_dist))
  out
}

# The variogram inputs of a kernel case as viladomatCorrelation() builds them.
case_variog_input <- function(cs) {
  ids <- cs$intermediate$forward$ids
  lat <- cs$input$pos[ids, 1]
  long <- cs$input$pos[ids, 2]
  list(ids = ids, lat = lat, long = long, max_dist = ref_prctile(lat, long, cs$input$maxDistPrctile))
}

# u, n, v and bins.lim of the C++ variogram identical() to geoR's, for every column of Z.
expect_variog_identical <- function(plan, Z, long, lat, max_dist, info) {
  V <- .stc_variog_eval(plan, Z)
  for (b in seq_len(ncol(Z))) {
    g <- geoR_variog(Z[, b], long, lat, max_dist)
    expect_identical(plan$u, g$u, info = info)
    expect_identical(plan$n, g$n, info = info)
    expect_identical(V[, b], g$v, info = info)
    expect_identical(plan$bins_lim, g$bins.lim, info = info)
  }
}

# A hexagonal lattice patch of N points with spacing sp (the layout SEraster's hexagonal pixels have).
hex_patch <- function(N, sp) {
  side <- ceiling(sqrt(N)) + 1
  g <- expand.grid(i = 0:(side - 1), j = 0:(side - 1))
  xy <- cbind((g$i + (g$j %% 2) / 2) * sp, g$j * sp * sqrt(3) / 2)
  xy[order(rowSums((xy - colMeans(xy))^2))[seq_len(N)], , drop = FALSE]
}

# Random coordinate sets: uniform, clustered, hexagonal and square lattices with offsets, and one
# with duplicated points (call with the default RNG set).
random_layouts <- function(n_sets) {
  lapply(seq_len(n_sets), function(r) {
    N <- sample(c(30L, 60L, 150L, 400L), 1)
    sp <- sample(c(0.055, 1, 50), 1)
    off <- sample(c(0, -1000, 1000), 1)
    type <- c("uniform", "clustered", "hex", "square", "duplicated")[(r - 1) %% 5 + 1]
    xy <- switch(type,
      uniform = cbind(stats::runif(N), stats::runif(N)) * sp * sqrt(N),
      clustered = {
        centres <- matrix(stats::runif(8), 4) * sp * sqrt(N)
        k <- sample(4, N, replace = TRUE)
        centres[k, ] + matrix(stats::rnorm(2 * N, sd = sp), N)
      },
      hex = hex_patch(N, sp),
      square = as.matrix(expand.grid(seq_len(ceiling(sqrt(N))), seq_len(ceiling(sqrt(N)))))[seq_len(N), ] * sp,
      duplicated = {
        xy <- cbind(stats::runif(N), stats::runif(N)) * sp * sqrt(N)
        rbind(xy, xy[c(2, 5, 9, 9), ])
      })
    list(type = type, xy = unname(xy) + off, p = sample(c(0.05, 0.1, 0.25, 0.5, 0.9), 1))
  })
}

# locfit's tree and fitted values for one (coordinate set, y, delta), computed once per session. xev holds
# locfit's vertex coordinates in its vertex order (one row per vertex).
lf_reference <- function(key, y, x1, x2, delta) fx_memo(paste0("cpp:lf:", key), {
  fit <- quiet_locfit(locfit::locfit(y ~ locfit::lp(x1, x2, nn = delta, deg = 0), kern = "gauss", maxk = 300))
  list(fitted = as.numeric(stats::fitted(fit)), nv = unname(fit$nvc["nv"]), nvm = unname(fit$nvc["nvm"]),
       xev = matrix(fit$eva$xev, ncol = 2, byrow = TRUE))
})

# The smoother on the fixture sets: permuted X (as the legacy code smooths), the default delta grid
# plus delta = 1.5 (fewer deltas for the N = 2170 set).
smoother_fixture_runs <- function() fx_memo("cpp:smoother_fixture_runs", {
  local_default_rng()
  set.seed(5)
  out <- list()
  for (nm in cpp_sets) {
    cs <- fx$cases[[nm]]
    x1 <- cs$input$pos[, 2]
    x2 <- cs$input$pos[, 1]
    y <- sample(cs$input$X)
    deltas <- if (nm == "brain_Oprk1_subsample") c(0.1, 0.5, 0.9, 1.5) else c(seq(0.1, 0.9, 0.1), 1.5)
    for (d in deltas) {
      key <- sprintf("%s delta = %g", nm, d)
      out[[key]] <- list(info = key, y = y, ref = lf_reference(key, y, x1, x2, d),
                         op = .stc_smoother(x1, x2, d), exact = .stc_smoother_fitted_exact(x1, x2, y, d))
    }
  }
  out
})

# The smoother on random layouts (uniform, clustered, lattices, duplicated points), N = 30 to 400.
smoother_random_runs <- function() fx_memo("cpp:smoother_random_runs", {
  local_default_rng()
  set.seed(77)
  out <- list()
  for (L in random_layouts(10)) {
    x1 <- L$xy[, 1]
    x2 <- L$xy[, 2]
    y <- stats::rnorm(length(x1), 5)
    for (d in sample(c(0.1, 0.3, 0.7, 1.2), 2)) {
      key <- sprintf("%s N = %d delta = %g", L$type, length(x1), d)
      out[[key]] <- list(info = key, y = y, x1 = x1, x2 = x2, ref = lf_reference(key, y, x1, x2, d),
                         op = .stc_smoother(x1, x2, d), exact = .stc_smoother_fitted_exact(x1, x2, y, d))
    }
  }
  out
})

# Exact comparisons with locfit need locfit to round like the replica. Builds of locfit whose
# compiler fuses multiply-adds (GCC on Linux arm64) do not, and neither does the replica match them.
# The exact tests are skipped where locfit rounds differently both from the replica (on a small fixed
# problem) and from the build machine (the global locfit probe of the platform signature). On the build
# machine they therefore always run, so a change in the replica cannot hide behind this gate.
skip_if_locfit_rounds_differently <- function() {
  mode <- fx_exact_mode()
  if (mode == "true") return(invisible(TRUE))
  if (mode == "false") testthat::skip("exact reference checks disabled (STCOMPARE_EXACT_TESTS=false)")
  st <- fx_memo("cpp:locfit_gate", {
    g <- expand.grid(a = seq(0, 1, length.out = 12), b = seq(0, 1, length.out = 12))
    px <- g$a + sin(17 * g$b) / 50
    py <- g$b + cos(13 * g$a) / 50
    z <- sin(3 * px) + cos(5 * py)
    like_replica <- identical(.stc_smoother_fitted_exact(px, py, z, 0.3)$fitted, ref_smooth(z, px, py, 0.3))
    like_build <- identical(stc_global_probe()$locfit, fx$meta$platform_signature$global$locfit)
    list(like_replica = like_replica, like_build = like_build)
  })
  if (!st$like_replica && !st$like_build) {
    testthat::skip(paste0("locfit on this machine (", R.version$platform, ") rounds differently from the replica and ",
                          "from the build machine (for example, a build with fused multiply-adds); the portable ",
                          "tests compare to 1e-13. STCOMPARE_EXACT_TESTS=true forces the exact checks"))
  }
  invisible(TRUE)
}

# --- variogram ----------------------------------------------------------------------------------

test_that("portable: the variogram plan and evaluator are identical() to geoR::variog on the fixture coordinate sets", {
  local_default_rng()
  set.seed(11)
  for (nm in cpp_sets) {
    cs <- fx$cases[[nm]]
    vi <- case_variog_input(cs)
    plan <- .stc_variog_plan(vi$long, vi$lat, vi$max_dist)
    expect_true(plan$ok, info = nm)
    M <- length(vi$ids)
    Z <- cbind(cs$input$X[vi$ids], cs$input$Y[vi$ids], stats::rnorm(M), stats::rpois(M, 2))
    expect_variog_identical(plan, Z, vi$long, vi$lat, vi$max_dist, nm)
    u <- as.vector(stats::dist(cbind(vi$long, vi$lat)))
    expect_identical(plan$umax, max(u[u < vi$max_dist]), info = nm)
    expect_identical(sum(plan$n), as.numeric(length(plan$i)), info = nm)
    expect_true(all(diff(plan$ptr) == plan$n), info = nm)
  }
})

for (nm in cpp_sets) {
  test_that(sprintf("exact: %s: the C++ variogram reproduces the stored target variogram", nm), {
    skip_if_not_exact(KF, nm)
    cs <- fx$cases[[nm]]
    vi <- case_variog_input(cs)
    plan <- .stc_variog_plan(vi$long, vi$lat, vi$max_dist)
    for (dir in c("forward", "reverse")) {
      tv <- cs$intermediate[[dir]]$target_variog
      z <- if (dir == "forward") cs$input$X else cs$input$Y
      expect_identical(plan$n, as.numeric(tv$n), info = dir)
      expect_identical(plan$u, tv$u, info = dir)
      expect_identical(as.vector(.stc_variog_eval(plan, matrix(z[vi$ids]))), tv$v, info = dir)
    }
  })
}

test_that("portable: the variogram is identical() to geoR on random layouts, duplicated points and exact bin-edge ties", {
  local_default_rng()
  set.seed(20261004)
  for (L in random_layouts(15)) {
    long <- L$xy[, 1]
    lat <- L$xy[, 2]
    info <- sprintf("%s N = %d p = %g", L$type, length(long), L$p)
    md <- ref_prctile(lat, long, L$p)
    plan <- .stc_variog_plan(long, lat, md)
    expect_true(plan$ok, info = info)
    expect_identical(plan$nugget, L$type == "duplicated", info = info)
    Z <- cbind(stats::rnorm(length(long)), stats::rpois(length(long), 3))
    expect_variog_identical(plan, Z, long, lat, md, info)
  }
  # integer grid: with max.dist = 13 + 1e-9, umax is exactly 13, the bin edges are exactly 1, ..., 12 and
  # thousands of pairs lie exactly on an edge (geoR counts them in the upper bin) or at umax (dropped);
  # with max.dist = 13, umax = sqrt(164) and the pairs at 13 are beyond it
  g <- as.matrix(expand.grid(1:20, 1:20)) + 0
  Zg <- cbind(stats::rnorm(400), (1:400) %% 7)
  for (md in c(13 + 1e-9, 13)) {
    plan <- .stc_variog_plan(g[, 1], g[, 2], md)
    expect_variog_identical(plan, Zg, g[, 1], g[, 2], md, sprintf("integer grid, max.dist = %.9f", md))
  }
  plan <- .stc_variog_plan(g[, 1], g[, 2], 13 + 1e-9)
  expect_identical(plan$umax, 13)
  expect_identical(plan$bins_lim, c(1e-12, 1:13))
  # tiles of any width give identical columns (the per-column summation order does not change)
  Zt <- matrix(stats::rnorm(400 * 37), 400)
  ref <- .stc_variog_eval(plan, Zt, tile = 1L)
  for (tile in c(3L, 4L, 7L, 16L, 33L, 1000L)) expect_identical(.stc_variog_eval(plan, Zt, tile = tile), ref, info = tile)
})

test_that("portable: variogram inputs on which geoR fails give ok = FALSE or an R error, never a crash", {
  x <- c(0, 1, 2, 3)
  y <- c(0, 0, 0, 0)
  # no pair closer than max.dist: geoR's max() is -Inf and seq() fails
  p <- .stc_variog_plan(x, y, 0.5)
  expect_false(p$ok)
  expect_identical(p$nbins, 0L)
  expect_error(suppressWarnings(geoR_variog(1:4 + 0, x, y, 0.5)))
  expect_error(.stc_variog_eval(p, matrix(c(1, 2, 3, 4))), "no pair table")
  # co-located points and only the nugget bin left: geoR fails at variogram.R:255
  p2 <- .stc_variog_plan(c(x, 0, 1), c(y, 0, 0), 0.5)
  expect_false(p2$ok)
  expect_error(geoR_variog(1:6 + 0, c(x, 0, 1), c(y, 0, 0), 0.5))
  # a plan without bins (its one pair has fewer than pairs.min = 2): geoR returns empty u, v and n, and so does
  # the evaluator, for any number of columns and tile width (and without touching an empty output buffer)
  p0 <- .stc_variog_plan(c(0, 1, 3), c(0, 0, 0), 2.5)
  expect_true(p0$ok)
  expect_identical(p0$nbins, 0L)
  for (tile in c(1L, 16L, 64L)) expect_identical(dim(.stc_variog_eval(p0, matrix(seq_len(51) / 7, 3), tile = tile)), c(0L, 17L))
  pp <- .stc_variog_plan(x, y, 2.5)
  expect_identical(dim(.stc_variog_eval(pp, matrix(numeric(0), 4, 0))), c(pp$nbins, 0L))
  # NA in the data: geoR fails with "NA/NaN/Inf in foreign function call", and so does the evaluator
  plan <- .stc_variog_plan(x, y, 2.5)
  expect_true(plan$ok)
  expect_error(.stc_variog_eval(plan, matrix(c(1, NA, 3, 4))), "NA/NaN/Inf")
  # corrupted tables and limits are rejected before any memory access
  bad <- plan
  bad$i[1] <- 99L
  expect_error(.stc_variog_eval(bad, matrix(c(1, 2, 3, 4))), "out of range")
  bad <- plan
  bad$ptr[2] <- 1e6L
  expect_error(.stc_variog_eval(bad, matrix(c(1, 2, 3, 4))), "invalid pair table")
  expect_error(.stc_variog_pairs(x, y, 2.5, c(0, 1e-12, 1, 2), FALSE), "invalid input")
  expect_error(.stc_variog_pairs(x, y, 2.5, c(-1, 1e-12, 2, 1), FALSE), "invalid input")
  expect_error(.stc_variog_plan(c(1, NA), c(1, 2), 1), "finite")
  expect_error(.stc_variog_plan(1, 1, 1), "at least 2 points")
})

# --- smoother -----------------------------------------------------------------------------------

test_that("portable: the smoother replica and the factored operator agree with fitted(locfit()) to 1e-13 (fixture sets, delta > 1)", {
  for (r in smoother_fixture_runs()) {
    expect_identical(r$op$status, 0L, info = r$info)
    expect_identical(r$exact$status, 0L, info = r$info)
    # same tree size and vertex capacity as locfit (atree_guessnv with maxk = 300)
    expect_identical(c(r$op$nv, r$op$nvm), as.integer(c(r$ref$nv, r$ref$nvm)), info = r$info)
    scale <- max(abs(r$ref$fitted))
    expect_lt(max(abs(r$exact$fitted - r$ref$fitted)) / scale, 1e-13)
    ap <- as.vector(.stc_smoother_apply(r$op, matrix(r$y)))
    expect_lt(max(abs(ap - r$ref$fitted)) / scale, 1e-13)
    # rows of Wn and of M are weights that sum to 1
    expect_lt(max(abs(rowSums(r$op$Wn) - 1)), 1e-13)
    expect_lt(max(abs(rowsum(r$op$M$val, rep(seq_len(r$op$n), diff(r$op$M$ptr)))[, 1] - 1)), 1e-13)
  }
})

test_that("exact: the smoother replica is identical() to fitted(locfit()) (fixture sets, default delta grid, delta > 1)", {
  skip_if_locfit_rounds_differently()
  for (r in smoother_fixture_runs()) {
    expect_identical(r$exact$fitted, r$ref$fitted, info = r$info)
    # the cells are visited in locfit's depth-first order: same vertices, created in the same order
    expect_identical(r$op$tree$xev, r$ref$xev, info = r$info)
  }
})

test_that("portable: the smoother agrees with locfit on random layouts and duplicated points; block and row-subset application", {
  for (r in smoother_random_runs()) {
    expect_identical(r$op$status, 0L, info = r$info)
    expect_identical(c(r$op$nv, r$op$nvm), as.integer(c(r$ref$nv, r$ref$nvm)), info = r$info)
    scale <- max(abs(r$ref$fitted))
    expect_lt(max(abs(r$exact$fitted - r$ref$fitted)) / scale, 1e-13)
    expect_lt(max(abs(as.vector(.stc_smoother_apply(r$op, matrix(r$y))) - r$ref$fitted)) / scale, 1e-13)
  }
  # a block of columns gives the columns' own results, and rows = ids gives those rows
  r <- smoother_random_runs()[[3]]
  local_default_rng()
  set.seed(3)
  Y <- matrix(stats::rnorm(length(r$y) * 33), length(r$y))
  full <- .stc_smoother_apply(r$op, Y)
  for (b in c(1, 17, 33)) expect_identical(full[, b], as.vector(.stc_smoother_apply(r$op, Y[, b, drop = FALSE])), info = b)
  rows <- sort(sample(length(r$y), 20))
  expect_identical(.stc_smoother_apply(r$op, Y, rows = rows), full[rows, , drop = FALSE])
})

test_that("exact: the smoother replica is identical() to fitted(locfit()) on random layouts and duplicated points", {
  skip_if_locfit_rounds_differently()
  for (r in smoother_random_runs()) {
    expect_identical(r$exact$fitted, r$ref$fitted, info = r$info)
    expect_identical(r$op$tree$xev, r$ref$xev, info = r$info)
  }
})

test_that("portable: the factored operator is offset-invariant (1e3 + N(0, 1) and 1e6 + N(0, 1) against fitted(locfit()))", {
  # Rows of Wn and M sum to 1 only to rounding, so M (Wn y) on raw values is off by about eps * |mean(y)|. The
  # operator centres y first (as locfit does); what is left is one rounding of the result. Bound: 1e-14 of
  # the spread of the data plus 2 eps max|fitted|. The uncentred operator exceeds it 2 to 4 times here.
  eps <- .Machine$double.eps
  for (nm in c("kidney_AB_jitter", "quakes_irregular")) {
    cs <- fx$cases[[nm]]
    x1 <- cs$input$pos[, 2]
    x2 <- cs$input$pos[, 1]
    local_default_rng()
    set.seed(12)
    z <- stats::rnorm(length(x1))
    for (d in c(0.1, 0.5, 1.5)) {
      op <- .stc_smoother(x1, x2, d)
      for (off in c(1e3, 1e6)) {
        y <- off + z
        ref <- as.numeric(stats::fitted(quiet_locfit(locfit::locfit(y ~ locfit::lp(x1, x2, nn = d, deg = 0), kern = "gauss", maxk = 300))))
        got <- as.vector(.stc_smoother_apply(op, matrix(y)))
        info <- sprintf("%s delta = %g offset = %g", nm, d, off)
        expect_lte(max(abs(got - ref)), 1e-14 * diff(range(y)) + 2 * eps * max(abs(ref)), label = info)
        # the same column without the offset is the same smooth, shifted
        expect_lt(max(abs(got - off - as.vector(.stc_smoother_apply(op, matrix(z))))), 4 * eps * off, label = info)
      }
    }
  }
})

# Every location twice, k = 2: status, k, nv = nvm and depth of the tree for N points.
duplicated_tree <- function(N) {
  u <- cbind(stats::runif(N / 2), stats::runif(N / 2))
  xy <- u[rep(seq_len(N / 2), each = 2), ]
  r <- .stc_smoother(xy[, 2], xy[, 1], 2 / N, build_operator = FALSE)
  list(xy = xy, r = r)
}

test_that("portable: exactly duplicated coordinates with a small k end with status 4 (out of vertex space), never a crash", {
  # With every location present twice and N * delta = 2 (k = 2), the cells around a duplicated point refine
  # without end until the vertex capacity is used up. locfit dies there with a C stack overflow; the tree
  # builder visits the cells from a heap stack (no recursion), so it returns the overflow status. A recursive
  # build overflowed R's 8 MB C stack from N = 2000 (N = 5000 is in the slow tier).
  local_default_rng()
  set.seed(9)
  d <- duplicated_tree(2000L)
  expect_identical(c(d$r$status, d$r$k), c(4L, 2L))
  expect_identical(d$r$nv, d$r$nvm)  # it failed when the capacity was exhausted
  expect_gt(d$r$depth, 1000L)        # far deeper than locfit's C stack allows (about 60 levels)
  d <- duplicated_tree(200L)
  expect_identical(.stc_smoother_fitted_exact(d$xy[, 2], d$xy[, 1], stats::rnorm(200), 2 / 200)$status, 4L)
  # a single duplicated pair among distinct points does the same for k = 2; with k = 3 the tree is ordinary
  u <- cbind(stats::runif(1000), stats::runif(1000))
  xy <- rbind(u, u[1, , drop = FALSE])
  expect_identical(.stc_smoother(xy[, 2], xy[, 1], 2.5 / 1001, build_operator = FALSE)$status, 4L)
  r3 <- .stc_smoother(xy[, 2], xy[, 1], 3.5 / 1001, build_operator = FALSE)
  expect_identical(r3$status, 0L)
  expect_lt(r3$depth, 40L)
  # the default grid on duplicated points: k is far above the number of copies
  r <- .stc_smoother(xy[, 2], xy[, 1], 0.1)
  expect_identical(r$status, 0L)
  expect_lt(r$depth, 20L)
})

test_that("portable, slow: duplicated coordinates at N = 5000 with k = 2 end with status 4", {
  skip_if_not_slow()
  local_default_rng()
  set.seed(10)
  d <- duplicated_tree(5000L)
  expect_identical(c(d$r$status, d$r$k), c(4L, 2L))
  expect_identical(d$r$nv, d$r$nvm)
  expect_gt(d$r$depth, 1000L)
})

test_that("portable: the smoother returns error statuses instead of failing or crashing", {
  cs <- fx$cases$kidney_AB_jitter
  x1 <- cs$input$pos[, 2]
  x2 <- cs$input$pos[, 1]
  N <- length(x1)
  # 1 <= delta * N < 2 gives k = 1, where locfit kills the R session (never call locfit with it)
  for (d in c(1, 1.5, 1.999) / N) {
    r <- .stc_smoother(x1, x2, d)
    expect_identical(c(r$status, r$k), c(3L, 1L), info = d * N)
    expect_null(r$Wn)
    expect_identical(.stc_smoother_fitted_exact(x1, x2, cs$input$X, d)$status, 3L)
  }
  # delta * N < 1 gives k = 0 (locfit: "procv: no points with non-zero weight")
  expect_identical(.stc_smoother(x1, x2, 0.5 / N)$status, 3L)
  # k = 2 works
  expect_identical(.stc_smoother(x1, x2, 2 / N)$status, 0L)
  # delta must be a positive number
  for (d in c(0, -0.1, NaN, Inf, NA)) expect_identical(.stc_smoother(x1, x2, d)$status, 2L, info = d)
  # out of vertex space: quakes_irregular needs 249 vertices at delta = 0.1; locfit fails for maxk = 78
  # (capacity 246) and succeeds for maxk = 79 (capacity 249), and so does the C++ tree
  q <- fx$cases$quakes_irregular$input
  qy <- q$X
  qx1 <- q$pos[, 2]
  qx2 <- q$pos[, 1]
  for (mk in c(78L, 79L)) {
    lf_ok <- tryCatch({
      quiet_locfit(locfit::locfit(qy ~ locfit::lp(qx1, qx2, nn = 0.1, deg = 0), kern = "gauss", maxk = mk))
      TRUE
    }, error = function(e) {
      expect_match(conditionMessage(e), "out of vertex space")
      FALSE
    })
    r <- .stc_smoother(qx1, qx2, 0.1, maxk = mk)
    expect_identical(r$status, if (lf_ok) 0L else 4L, info = mk)
    expect_identical(lf_ok, mk == 79L)
  }
  # invalid input
  expect_identical(.stc_smoother(1, 1, 0.5)$status, 1L)
  expect_identical(.stc_smoother(c(1, NA, 3), c(1, 2, 3), 0.9)$status, 1L)
  expect_identical(.stc_smoother(c(1, Inf, 3), c(1, 2, 3), 0.9)$status, 1L)
  expect_identical(.stc_smoother(x1, x2, 0.5, maxk = 0L)$status, 1L)
  expect_error(.stc_smoother(1:3, 1:2, 0.5), "same length")
  expect_error(.stc_smoother_apply(.stc_smoother(x1, x2, 1 / N), matrix(1, N, 1)), "no operator")
  op <- .stc_smoother(x1, x2, 0.5)
  expect_error(.stc_smoother_apply(op, matrix(1, 3, 1)), "nrow")
  expect_error(.stc_smoother_apply(op, matrix(1, N, 1), rows = N + 1L), "out of range")
})

test_that("portable: the smoother's results survive garbage collection while they are returned (gctorture)", {
  # .stc_smoother() and .stc_smoother_fitted_exact() used to keep the tree, Wn, M and the fitted values in raw
  # SEXP variables after the Rcpp objects that protected them had gone out of scope. A garbage collection
  # while the result list was built freed them; R then handed their memory out again (4 of 5 rounds here).
  local_default_rng()
  set.seed(4)
  N <- 40
  x1 <- stats::runif(N)
  x2 <- stats::runif(N)
  y <- stats::rnorm(N)
  ref_f <- .stc_smoother_fitted_exact(x1, x2, y, 0.5)$fitted
  ref_s <- .stc_smoother(x1, x2, 0.5)
  for (round in 1:3) {
    tryCatch({
      gctorture(TRUE)
      f <- .stc_smoother_fitted_exact(x1, x2, y, 0.5)
      s <- .stc_smoother(x1, x2, 0.5)
    }, finally = gctorture(FALSE))
    # allocations of the same sizes reuse any memory that was freed
    junk <- lapply(1:2000, function(i) list(rep(-1.5, N), matrix(-2.5, nrow(ref_s$Wn), ncol(ref_s$Wn)), 1:6))
    rm(junk)
    info <- paste("round", round)
    expect_identical(f$fitted, ref_f, info = info)
    expect_identical(s$tree, ref_s$tree, info = info)
    expect_identical(s$Wn, ref_s$Wn, info = info)
    expect_identical(s$M, ref_s$M, info = info)
  }
})

# --- RNG ----------------------------------------------------------------------------------------

test_that("portable: the C++ L'Ecuyer-CMRG stream is R's: set.seed() state, runif() and rnorm() for seeds 0..1000, negative and large seeds", {
  seeds <- c(0:1000, -1L, -2L, -1000L, -123456789L, .Machine$integer.max, -.Machine$integer.max, 123456789L, 2L^30)
  as_unsigned <- function(s) as.numeric(s) %% 2^32
  same <- vapply(seeds, function(s) {
    with_lecuyer(s, c(state = identical(.stc_lecuyer_seed(s), as_unsigned(.Random.seed[2:7])),
                      rnorm = identical(.stc_lecuyer_rnorm(s, 100L), stats::rnorm(100))))
  }, c(state = NA, rnorm = NA))
  expect_length(seeds[!same["state", ]], 0)
  expect_length(seeds[!same["rnorm", ]], 0)
  unif <- vapply(seeds[c(1:50, 1002:1009)], function(s) with_lecuyer(s, identical(.stc_lecuyer_runif(s, 50L), stats::runif(50))), NA)
  expect_true(all(unif))
})

test_that("portable: C++ normals for n = 1e5 and the legacy noise blocks are identical() to R's", {
  for (s in c(0L, 42L, -7L, .Machine$integer.max)) {
    expect_identical(with_lecuyer(s, stats::rnorm(1e5)), .stc_lecuyer_rnorm(s, 100000L), info = s)
  }
  # the noise of permutation b: one rnorm(N) per delta after set.seed(seed + b) (matchingVariograms())
  for (s in c(1L, 2L, 18L)) expect_identical(.stc_legacy_noise(s, 273L, 9L), ref_noise(273, s, 9), info = s)
  expect_identical(dim(.stc_legacy_noise(1L, 0L, 3L)), c(0L, 3L))
})

test_that("portable: the C++ generator leaves R's RNG state alone and rejects an NA seed", {
  local_default_rng()
  set.seed(1)
  before <- .Random.seed
  invisible(.stc_lecuyer_rnorm(5L, 1000L))
  invisible(.stc_legacy_noise(5L, 100L, 3L))
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  expect_error(.stc_lecuyer_rnorm(NA_integer_, 3L), "not a valid integer")
})

test_that("portable: no compiled entry point creates or changes .Random.seed", {
  # dev/engine-spec.md 1.5. Rcpp wraps an exported function in GetRNGstate()/PutRNGstate() unless its
  # attribute says rng = false, and GetRNGstate() creates a time-seeded .Random.seed when there is none.
  local_default_rng()
  x <- c(0, 1, 3, 4, 7)
  y <- c(0, 0, 1, 2, 1)
  z <- c(1, 2, 3, 5, 8)
  plan <- .stc_variog_plan(x, y, 10)
  op <- .stc_smoother(x, y, 0.9)
  calls <- list(
    variog_pairs = function() .stc_variog_pairs(x, y, 10, c(-1, 1e-12, 2, 4, 8), FALSE),
    variog_eval = function() .stc_variog_eval(plan, matrix(z)),
    smoother = function() .stc_smoother(x, y, 0.9),
    smoother_fitted_exact = function() .stc_smoother_fitted_exact(x, y, z, 0.9),
    smoother_apply = function() .stc_smoother_apply(op, matrix(z)),
    lecuyer_seed = function() .stc_lecuyer_seed(1L),
    lecuyer_runif = function() .stc_lecuyer_runif(1L, 3L),
    lecuyer_rnorm = function() .stc_lecuyer_rnorm(1L, 3L),
    legacy_noise = function() .stc_legacy_noise(1L, 4L, 2L),
    cor_cols = function() .stc_cor_cols(matrix(z), c(2, 1, 4, 3, 5), 0L),
    ols = function() .stc_ols(c(1, 2, 3), c(2, 4, 7)),
    variog_plan = function() .stc_variog_plan(x, y, 10))
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

# --- correlation --------------------------------------------------------------------------------

# Test inputs: random columns of many lengths, legacy surrogates (the stored hat vectors of kidney_AB
# against Y), and the NA / NaN / Inf / constant / clamping cases.
cor_inputs <- function() fx_memo("cpp:cor_inputs", {
  local_default_rng()
  set.seed(8)
  out <- lapply(c(2:12, 31, 64, 65, 273, 1000, 2170), function(n) {
    list(X = matrix(stats::rnorm(n * 8, 10, 3), n), y = stats::rnorm(n) * 2 + 1)
  })
  cs <- fx$cases$kidney_AB
  out[[length(out) + 1]] <- list(X = vapply(cs$intermediate$forward$detail$per_delta, function(p) p$hat, numeric(273)),
                                 y = cs$input$Y)
  y <- c(2, 1, 5, 3, 8, 4, 7, 6, 9, 11)
  out[[length(out) + 1]] <- list(X = cbind(1:10, c(1:9, NA), rep(3, 10), c(1:9, Inf), c(1:9, NaN), 3 * y + 1, -2 * y), y = y)
  out
})

test_that("portable: the C++ correlation is within 1e-15 of cor() (selected mode and the double mode), with R's NA cases", {
  mode <- .stc_cor_mode()
  expect_true(mode %in% 0:3)
  for (inp in cor_inputs()) {
    r <- suppressWarnings(as.vector(stats::cor(inp$X, inp$y)))
    for (m in unique(c(as.integer(mode), 3L))) {
      got <- .stc_cor_cols(inp$X, inp$y, m)
      expect_identical(is.na(got), is.na(r), info = m)
      expect_identical(is.nan(got), is.nan(r), info = m)
      ok <- !is.na(r)
      expect_lt(max(abs(got[ok] - r[ok]), 0), 1e-15)
    }
  }
  # NA in y, n = 1, n = 0 (no "subscript out of bounds" warning from Rcpp either)
  expect_identical(.stc_cor_cols(matrix(1:4 + 0, 4), c(1, NA, 2, 3), mode), NA_real_)
  expect_identical(.stc_cor_cols(matrix(1, 1, 2), 1, mode), c(NA_real_, NA_real_))
  expect_silent(r0 <- .stc_cor_cols(matrix(numeric(0), 0, 2), numeric(0), mode))
  expect_identical(r0, c(NA_real_, NA_real_))
  expect_error(.stc_cor_cols(matrix(1, 3, 1), 1:2 + 0), "nrow")
})

test_that("portable: .stc_cor_mode() selects an exact mode whenever one of the four modes reproduces cor()", {
  # mode 3 (double accumulators) is R's arithmetic when R is built with --disable-long-double
  r <- lapply(cor_inputs(), function(inp) suppressWarnings(as.vector(stats::cor(inp$X, inp$y))))
  exact <- vapply(0:3, function(m) {
    all(mapply(function(inp, ri) identical(.stc_cor_cols(inp$X, inp$y, m), ri), cor_inputs(), r))
  }, NA)
  mode <- .stc_cor_mode()
  if (any(exact)) {
    expect_true(isTRUE(attr(mode, "exact")))
    expect_true(exact[as.integer(mode) + 1L])
  } else {
    expect_false(isTRUE(attr(mode, "exact")))
  }
})

test_that("exact: the C++ correlation is identical() to cor() in the mode .stc_cor_mode() selects", {
  mode <- .stc_cor_mode()
  if (!isTRUE(attr(mode, "exact"))) {
    testthat::skip(sprintf("no mode reproduces this R build's cor() bit for bit (using mode %d, within 1e-15)", mode))
  }
  for (inp in cor_inputs()) {
    expect_identical(.stc_cor_cols(inp$X, inp$y, mode), suppressWarnings(as.vector(stats::cor(inp$X, inp$y))))
  }
})

# --- least squares ------------------------------------------------------------------------------

# (target, candidate) variogram pairs as matchingVariograms() fits them: the pairs stored in the
# fixture, and pairs from permutations smoothed with the C++ operator.
ols_pairs <- function() fx_memo("cpp:ols_pairs", {
  pairs <- list()
  for (cs in fx$cases) for (dir in c("forward", "reverse")) {
    tv <- cs$intermediate[[dir]]$target_variog$v
    for (p in cs$intermediate[[dir]]$detail$per_delta) {
      if (!is.null(p$variog_fitted_v)) pairs[[length(pairs) + 1]] <- list(y = tv, x = p$variog_fitted_v)
    }
  }
  for (nm in c("kidney_AB_jitter", "quakes_irregular")) {
    cs <- fx$cases[[nm]]
    vi <- case_variog_input(cs)
    plan <- .stc_variog_plan(vi$long, vi$lat, vi$max_dist)
    tv <- as.vector(.stc_variog_eval(plan, matrix(cs$input$X[vi$ids])))
    P <- cs$intermediate$forward$perm_index
    Xr <- matrix(cs$input$X[P], nrow(P))
    for (d in seq(0.1, 0.9, 0.2)) {
      S <- .stc_smoother_apply(.stc_smoother(cs$input$pos[, 2], cs$input$pos[, 1], d), Xr, rows = vi$ids)
      V <- .stc_variog_eval(plan, S)
      for (b in seq_len(ncol(V))) pairs[[length(pairs) + 1]] <- list(y = tv, x = V[, b])
    }
  }
  pairs
})

test_that("portable: closed-form least squares agrees with lm() on realistic variogram pairs (slope 1e-14, intercept 1e-12 relative)", {
  rel <- vapply(ols_pairs(), function(p) {
    l <- ref_lm(p$y, p$x)
    o <- .stc_ols(p$x, p$y)
    c(status = o$status, b0 = abs(o$coefficients[1] - l[1]) / abs(l[1]), b1 = abs(o$coefficients[2] - l[2]) / abs(l[2]))
  }, c(status = 0, b0 = 0, b1 = 0))
  expect_true(all(rel["status", ] == 0))
  expect_gt(ncol(rel), 150)
  expect_lt(max(rel["b1", ]), 1e-14)
  expect_lt(max(rel["b0", ]), 1e-12)
})

test_that("portable: least squares makes lm()'s NA-slope decisions on near-degenerate and degenerate inputs", {
  # x = c + e * z with ||x - mean(x)|| / ||x|| = ratio; lm() drops the slope below about 1e-7
  for (K in c(2L, 3L, 5L, 10L, 13L)) for (ratio in c(1.2e-7, 1.01e-7, 0.99e-7, 9e-8, 1e-8, 0)) for (cc in c(1, 1e-3, 1e4)) {
    z <- if (K == 2L) c(-1, 1) else sin(seq_len(K) * 1.7) - mean(sin(seq_len(K) * 1.7))
    z <- z / sqrt(sum(z^2))
    x <- cc + ratio * sqrt(K) * cc / sqrt(1 - ratio^2) * z
    y <- 1 + seq_len(K) / 3
    l <- stats::lm(y ~ 1 + x)$coefficients
    o <- .stc_ols(x, y)
    info <- sprintf("K = %d ratio = %g c = %g", K, ratio, cc)
    expect_identical(is.na(o$coefficients[2]), is.na(unname(l[2])), info = info)
    expect_identical(o$status, if (is.na(l[2])) 1L else 0L, info = info)
    if (is.na(l[2])) expect_equal(o$coefficients[1], unname(l[1]), tolerance = 1e-12, info = info)
  }
  # an all-zero or constant candidate variogram: NA slope, intercept mean(y), as lm()
  y <- c(1, 2, 4, 3)
  for (x in list(c(0, 0, 0, 0), c(5, 5, 5, 5))) {
    o <- .stc_ols(x, y)
    expect_identical(o$status, 1L)
    expect_equal(o$coefficients, unname(stats::lm(y ~ 1 + x)$coefficients), tolerance = 1e-14)
  }
  # one bin: lm() returns y and an NA slope; no bins: lm() fails
  o <- .stc_ols(3, 5)
  expect_identical(o$status, 2L)
  expect_identical(o$coefficients, c(5, NA_real_))
  expect_identical(.stc_ols(numeric(0), numeric(0))$status, 2L)
  expect_identical(.stc_ols(c(1, NA, 3), c(1, 2, 3))$status, 3L)
  expect_identical(.stc_ols(c(1, 2, 3), c(1, Inf, 3))$coefficients, c(NA_real_, NA_real_))
  expect_error(.stc_ols(1:3 + 0, 1:2 + 0), "same length")
})

# --- the blocks together ------------------------------------------------------------------------

test_that("portable: the blocks chained replay permutation 1 like the reference kernels (kidney_AB_jitter: RSS 1e-10, same deltaStar)", {
  cs <- fx$cases$kidney_AB_jitter
  lat <- cs$input$pos[, 1]
  long <- cs$input$pos[, 2]
  for (dir in c("forward", "reverse")) {
    rp <- fx_replay(cs, dir)  # locfit, geoR and lm on this machine
    cap <- cs$intermediate[[dir]]
    ids <- cap$ids
    z <- if (dir == "forward") cs$input$X else cs$input$Y
    plan <- .stc_variog_plan(long[ids], lat[ids], rp$prctile)
    tv <- as.vector(.stc_variog_eval(plan, matrix(z[ids])))
    expect_identical(tv, rp$target$v, info = dir)
    noise <- .stc_legacy_noise(cap$detail$noise_seed, length(z), length(rp$delta))
    expect_identical(noise, rp$noise, info = dir)
    rss <- vapply(seq_along(rp$delta), function(k) {
      xd <- as.vector(.stc_smoother_apply(.stc_smoother(long, lat, rp$delta[k]), matrix(rp$xr)))
      v1 <- as.vector(.stc_variog_eval(plan, matrix(xd[ids])))
      bet <- .stc_ols(v1, tv)$coefficients
      hat <- xd * sqrt(abs(bet[2])) + noise[, k] * sqrt(abs(bet[1]))
      v2 <- as.vector(.stc_variog_eval(plan, matrix(hat[ids])))
      sum((v2 - tv)^2)
    }, 0)
    expect_equal(rss, rp$rss, tolerance = 1e-10, info = dir)
    expect_identical(which.min(rss), rp$argmin, info = dir)
  }
})
