# Internal helpers of the compiled engine (dev/engine-spec.md).
#
# This file holds the R side of the engine's exact building blocks. Their C++ entry points are in
# R/RcppExports.R (generated from src/stc_rcpp.cpp). All names start with a dot, so NAMESPACE's
# exportPattern("^[[:alpha:]]+") does not export them, and no exported function uses them yet.
#
#   .stc_variog_plan()      R part of geoR::variog() (umax, nugget flag, bin limits) + C++ pair table
#   .stc_variog_eval()      variograms of many data vectors with a plan (C++)
#   .stc_smoother()         locfit's tree and the factored smoothing operator (C++)
#   .stc_smoother_apply()   apply the operator to a block of columns (C++)
#   .stc_smoother_fitted_exact()  exact replica of fitted(locfit(...)), for tests (C++)
#   .stc_lecuyer_rnorm(), .stc_legacy_noise()  R-compatible L'Ecuyer-CMRG normals (C++)
#   .stc_cor_cols()         R's cor(X, y), column by column (C++); .stc_cor_mode() picks its mode
#   .stc_ols()              lm(y ~ 1 + x) coefficients by the closed form (C++)

# Session cache of the engine (results of one-off probes).
.stc_engine_cache <- new.env(parent = emptyenv())

# Pair table of
#   geoR::variog(coords = cbind(long, lat), data = z, max.dist = max_dist, option = "bin")
# (every other argument at its default; geoR 1.9-6) for fixed coordinates, so that
# .stc_variog_eval(plan, Z) returns the same u, n and v as geoR for every column of Z.
#
# long, lat: coordinates of the (subsampled) points in geoR's column order, as STcompare passes
#   them (coords = cbind(long[ids], lat[ids]) with lat = pos[, 1] and long = pos[, 2]).
# max_dist: the max.dist STcompare passes, quantile(dist(cbind(lat, long)), maxDistPrctile).
#
# Everything geoR computes with R's own dist() is computed here with dist(), so it follows the
# running R build (dist() uses a fused multiply-add on some platforms, and on lattices many pairs
# tie with umax to the last bit; dev/investigation/02-geoR-variogram.md):
#   u <- dist(coords); nugget <- min(u) < 1e-12; umax <- max(u[u < max.dist])  (variogram.R:93-119)
#   bins.lim <- seq(0, umax, length.out = 14); bins.lim <- c(0, 1e-12, bins.lim[bins.lim > 1e-12])
#   u (centres) <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])          (.define.bins())
#   bins.lim[1] <- -1                                                           (variogram.R:127)
# The C++ part then runs geoR's binit() loop with hypot() (src/geoR.c:370-422).
#
# Returns a list:
#   ok, reason   ok = FALSE where geoR::variog() itself fails (no pair closer than max.dist; or
#                co-located points with fewer than two bins kept, where variogram.R:255 fails);
#   nbins, u, n  the bins geoR returns (n >= 2 pairs; the nugget bin only with co-located points);
#   bins_lim     geoR's $bins.lim;
#   ptr, i, j    the pairs of each bin in geoR's summation order (CSR, 0-based point indices);
#   bin          geoR's bin index of each returned bin (0 = nugget bin);
#   umax, max_dist, nugget, lims, npoints, npairs_le_maxdist.
# Unlike geoR, nothing is printed for co-located points.
.stc_variog_plan <- function(long, lat, max_dist) {
  long <- as.double(long)
  lat <- as.double(lat)
  max_dist <- as.double(max_dist)  # quantile() names its value; geoR passes as.double(max.dist)
  if (length(long) != length(lat)) stop("long and lat must have the same length")
  if (length(long) < 2L) stop("the variogram needs at least 2 points")
  if (length(max_dist) != 1L || is.na(max_dist)) stop("max_dist must be a single number")
  if (!all(is.finite(long)) || !all(is.finite(lat))) stop("coordinates must be finite")
  fail <- function(reason) {
    list(ok = FALSE, reason = reason, nbins = 0L, u = numeric(0), n = numeric(0),
         npoints = length(long), max_dist = max_dist)
  }
  u <- as.vector(stats::dist(cbind(long, lat)))
  nugget <- min(u) < 1e-12
  below <- u[u < max_dist]
  if (!length(below)) return(fail("no pair of points is closer than max.dist"))
  umax <- max(below)
  bins_lim <- seq(0, umax, length.out = 14)
  bins_lim <- c(0, 1e-12, bins_lim[bins_lim > 1e-12])
  centres <- 0.5 * (bins_lim[-1] + bins_lim[-length(bins_lim)])
  lims <- bins_lim
  if (lims[1] < 1e-16) lims[1] <- -1
  tab <- .stc_variog_pairs(long, lat, max_dist, lims, nugget, 2L)
  u_out <- centres[tab$bin + 1L]
  if (nugget) {
    # variogram.R:252-257: with co-located points geoR replaces the second centre when the first two
    # are below 1e-11, and fails (an NA condition) when fewer than two bins are left
    first2 <- all(u_out[1:2] < 1e-11)
    if (is.na(first2)) return(fail("co-located points and fewer than two variogram bins"))
    if (first2) u_out[2] <- sum(bins_lim[2:3]) / 2
  }
  list(ok = TRUE, reason = "", nbins = length(tab$n), u = u_out, n = tab$n,
       bins_lim = if (nugget) bins_lim else bins_lim[-1],
       ptr = tab$ptr, i = tab$i, j = tab$j, bin = tab$bin,
       umax = umax, max_dist = max_dist, nugget = nugget, lims = lims,
       npoints = length(long), npairs_le_maxdist = tab$npairs_le_maxdist)
}

# The arithmetic mode of .stc_cor_cols() for this machine (see src/stc_stats.h): 0 R's code as
# written (long double), 1 every accumulation fused, 2 the pattern of CRAN's R for macOS arm64, 3 double
# accumulators (R built with --disable-long-double). It compares cor() with the C++ code on fixed probes
# once per session and returns the first of modes 0, 2, 1 and 3 that reproduces cor() bit for bit
# (attribute exact = TRUE). If none does, it returns 0 when that agrees with cor() to 1e-14, and
# otherwise 3, the mode without x87 arithmetic (Rosetta 2 mis-executes the long double code; exact =
# FALSE in both cases; mode 3 then agrees with cor() to a few 1e-15).
.stc_cor_mode <- function() {
  if (is.null(.stc_engine_cache$cor_mode)) {
    probe <- function(n) {
      k <- seq_len(n)
      y <- sin(0.731 * k) * 3 + k / 17
      X <- vapply(seq_len(12), function(j) cos(k * (0.37 + j / 11)) + sin(j * k / 5) / (j + 1) + k * j / 1000,
                  numeric(n))
      list(X = X, y = y, r = as.vector(stats::cor(X, y)))
    }
    probes <- lapply(c(69L, 64L, 7L), probe)
    exact <- function(mode) all(vapply(probes, function(p) identical(.stc_cor_cols(p$X, p$y, mode), p$r), NA))
    maxdiff <- function(mode) max(vapply(probes, function(p) max(abs(.stc_cor_cols(p$X, p$y, mode) - p$r)), 0))
    mode <- NA_integer_
    for (m in c(0L, 2L, 1L, 3L)) {
      if (exact(m)) {
        mode <- m
        break
      }
    }
    is_exact <- !is.na(mode)
    if (!is_exact) mode <- if (isTRUE(maxdiff(0L) <= 1e-14)) 0L else 3L
    .stc_engine_cache$cor_mode <- structure(mode, exact = is_exact)
  }
  .stc_engine_cache$cor_mode
}
