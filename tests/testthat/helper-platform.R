# Platform arithmetic signature. Sourced by testthat before the tests and by data-raw/build_test_fixtures.R,
# which stores the signature of the build machine in every fixture (meta$platform_signature). It uses only
# base R, stats, geoR and locfit, so both sides compute it with the same code.
#
# Why it exists: the legacy results reproduce at the 1e-12 level only where geoR bins the pixel pairs exactly
# as on the build machine. geoR::variog() takes umax = max(u[u < max.dist]) from R's dist() and counts a pair
# only if libm hypot(dx, dy) < umax, so in EVERY data set the pair that defines umax sits exactly on the last
# bin edge, and lattice coordinates put many tied pairs there. Whether those pairs are counted depends on the
# last bit of dist() (compiled with or without fused multiply-add) and of the libm hypot(). One different bin
# count changes the target variogram and so every null correlation (by about 1e-4) and some deltaStar values
# (those nulls then move by up to about 0.2). The fixtures were built on macOS arm64. On Linux (glibc 2.39)
# the bin counts differ for kidney_AB, quakes_irregular and the brain sets on x86-64, and for the brain sets on
# arm64; the AKI and brain sets of the realistic tier differ on both.
#
# Other last-bit differences (R's qnorm() behind rnorm(), long double accumulation in sum()/mean()/cor(),
# locfit and lm() arithmetic) move results by about 1e-16 relative, which the tests absorb (tolerance 1e-12;
# RNG draws 1e-14). They are recorded as global probes and reported, but they do not gate anything.
#
# The signature has:
#   info:   build platform and R version (reported only);
#   global: qnorm() at fixed probabilities, the first rnorm() draws under L'Ecuyer-CMRG, a long double
#           accumulation probe, dist() and hypot() of fixed points, lm() coefficients and locfit fitted values of
#           fixed small problems, as exact hexadecimal strings (reported only);
#   sets:   for every coordinate set a fixture uses, the probes that decide exactness: max.dist (the dist()
#           quantile the package computes), geoR's umax, the number of pairs tied at umax and geoR's bin counts
#           (which depend only on the coordinates, so they are the bin counts of every variogram of that set).

stc_hex <- function(x) sprintf("%a", as.numeric(x))

# Run `code` with the given RNG kinds and restore the caller's kinds and seed afterwards.
stc_with_rng <- function(kind, code) {
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  RNGkind(kind[1], kind[2], kind[3])
  code
}

# Variogram subsample of viladomatCorrelation(): sample(N, 1000) right after set.seed(seed) when N > 1000
# (Mersenne-Twister, Rejection sampling), otherwise all points.
stc_subsample_ids <- function(N, seed) {
  if (N <= 1000) return(seq_len(N))
  stc_with_rng(c("Mersenne-Twister", "Inversion", "Rejection"), { set.seed(seed); sample(N, 1000) })
}

# Probes of one coordinate set as viladomatCorrelation() sees it: lat = pos[, 1], long = pos[, 2], `ids` the
# variogram subsample. max.dist is the quantile of dist(cbind(lat, long)); geoR computes u from
# coords = cbind(long, lat).
stc_coord_probe <- function(pos, ids, maxDistPrctile) {
  lat <- as.numeric(pos[ids, 1])
  long <- as.numeric(pos[ids, 2])
  prctile <- unname(stats::quantile(stats::dist(cbind(lat, long)), probs = maxDistPrctile))
  u <- as.vector(stats::dist(cbind(long, lat)))
  umax <- max(u[u < prctile])
  v <- geoR::variog(data = seq_along(lat), coords = cbind(long, lat), max.dist = prctile, option = "bin",
                    messages = FALSE)
  list(max_dist = stc_hex(prctile), umax = stc_hex(umax), n_pairs_at_umax = sum(u == umax),
       bin_n = as.integer(v$n))
}

# Probes that depend on no fixture (reported, not gating). The caller's RNG kind and seed are restored.
stc_global_probe <- function() {
  p <- c(1e-300, 1e-10, 1e-4, 0.02425, 0.1, 0.2, 0.3, 0.4, 0.45, 0.55, 0.6, 0.7, 0.8, 0.9, 0.97575, 1 - 1e-4)
  rn <- stc_with_rng(c("L'Ecuyer-CMRG", "Inversion", "Rejection"), { set.seed(1); stats::rnorm(64) })
  x <- seq(0.05, 1, length.out = 13)
  y <- 0.3 + 2 * x + sin(7 * x) / 10
  g <- expand.grid(a = seq(0, 1, length.out = 12), b = seq(0, 1, length.out = 12))
  px <- g$a + sin(17 * g$b) / 50
  py <- g$b + cos(13 * g$a) / 50
  z <- sin(3 * px) + cos(5 * py)
  fit <- suppressWarnings(locfit::locfit(z ~ locfit::lp(px, py, nn = 0.3, deg = 0), kern = "gauss", maxk = 300))
  list(qnorm = stc_hex(stats::qnorm(p)), rnorm_lecuyer_seed1 = stc_hex(rn),
       long_double_sum = stc_hex(sum(c(1, 2^-60, -1))),
       dist = stc_hex(stats::dist(cbind(px[1:40], py[1:40]))),
       hypot = stc_hex(Mod(complex(real = px[1:40] - px[41:80], imaginary = py[1:40] - py[41:80]))),
       lm = stc_hex(stats::lm(y ~ 1 + x)$coefficients), locfit = stc_hex(stats::fitted(fit)))
}

# sets: named list of list(pos, ids, maxDistPrctile).
stc_platform_signature <- function(sets) {
  list(info = list(platform = R.version$platform, R = R.version.string,
                   long_double = unname(capabilities("long.double"))),
       global = stc_global_probe(),
       sets = lapply(sets, function(s) stc_coord_probe(s$pos, s$ids, s$maxDistPrctile)))
}

# Probes of coordinate set `set` that differ between two signatures (these decide exactness).
stc_set_diff <- function(ref, cur, set) {
  a <- ref$sets[[set]]
  b <- cur$sets[[set]]
  if (is.null(a) || is.null(b)) return("no signature for this set")
  names(a)[!vapply(names(a), function(k) identical(a[[k]], b[[k]]), NA)]
}

# Global probes that differ (reported only).
stc_global_diff <- function(ref, cur) {
  names(ref$global)[!vapply(names(ref$global), function(k) identical(ref$global[[k]], cur$global[[k]]), NA)]
}

# Coordinate sets of each fixture, read from the fixture (the builder calls this on the fixture it writes).
stc_signature_sets <- function(name, fx) {
  if (name == "kernel_fixture.rds") {
    return(lapply(fx$cases, function(cs) list(pos = cs$input$pos, ids = cs$intermediate$forward$ids,
                                             maxDistPrctile = cs$input$maxDistPrctile)))
  }
  if (name == "realistic_fixture.rds") {
    return(lapply(fx$pairs, function(P) list(pos = P$pos, ids = stc_subsample_ids(nrow(P$pos), P$params$seed),
                                            maxDistPrctile = P$params$maxDistPrctile)))
  }
  if (name == "calibration_fixture.rds") {
    jobs <- fx$test_jobs
    sets <- lapply(seq_len(nrow(jobs)), function(k) {
      sh <- which(!is.na(fx$fields[jobs$i[k], ]) & !is.na(fx$fields[jobs$j[k], ]))
      list(pos = unname(fx$coords[sh, , drop = FALSE]), ids = stc_subsample_ids(length(sh), 0), maxDistPrctile = 0.25)
    })
    names(sets) <- sprintf("job%02d", seq_len(nrow(jobs)))
    return(sets)
  }
  stop("unknown fixture ", name)
}
