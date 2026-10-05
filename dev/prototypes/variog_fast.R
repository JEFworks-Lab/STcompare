## Pure-R, vectorised re-implementation of
##   geoR::variog(coords = coords, data = z, max.dist = max.dist, option = "bin", messages = FALSE)
## with every other argument at its default (geoR 1.9-6), for MANY data vectors at once.
##
## Split into
##   prep  <- variog_prep_r(coords, max.dist)   # depends on coordinates + max.dist only (once)
##   out   <- variog_compute_r(prep, Z)         # Z: N x B matrix, one variogram per column
## Reproduces geoR bit-for-bit (u, v, n) on this platform; see verify script.
##
## Source references (geoR 1.9-6 tarball):
##   R/variogram.R:103      u <- as.vector(dist(as.matrix(coords)))       (R's dist(), FMA on arm64 CRAN R)
##   R/variogram.R:104-113  nugget.tolerance 1e-12 ; nt.ind <- min(u) < 1e-12
##   R/variogram.R:132      umax <- max(u[u < max.dist])                   (strict <)
##   R/variogram.R:1023-28  .define.bins: seq(0, umax, l = 14); c(0, 1e-12, >1e-12); midpoints
##   R/variogram.R:137      bins.lim[1] <- -1
##   src/geoR.c:370-422     binit(): hypot(dx,dy) <= max.dist ; left-closed bins ; d >= umax dropped
##   R/variogram.R:151-168  pairs.min = 2 filter, drop nugget bin, keep.NA = FALSE -> omit bins

variog_prep_r <- function(coords, max.dist, uvec = 13, pairs.min = 2,
                          nugget.tolerance = 1e-12) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  n <- nrow(coords)
  max.dist <- as.double(max.dist)          # .C(as.double(max.dist)) drops the "25%" name
  ## --- R side of geoR::variog: distances from R's dist() (NOT hypot) ---------------------
  u_r <- as.vector(dist(coords))           # variogram.R:103
  nt.ind <- min(u_r) < nugget.tolerance    # variogram.R:112-113
  umax <- max(u_r[u_r < max.dist])         # variogram.R:132
  ## --- .define.bins (variogram.R:1023-1028) ------------------------------------------------
  bl <- seq(0, umax, length.out = uvec + 1)               # 0, k*(umax/13) (k=1..12), umax
  bins.lim <- c(0, nugget.tolerance, bl[bl > nugget.tolerance])
  uvec_mid <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
  nbins <- length(bins.lim) - 1L
  lims <- bins.lim
  if (lims[1] < 1e-16) lims[1] <- -1                      # variogram.R:137
  ## --- C side (binit): pairs in loop order j < i, j outer == order of as.vector(dist()) ----
  jj <- rep.int(seq_len(n - 1L), (n - 1L):1L)
  ii <- sequence((n - 1L):1L, from = 2:n)
  dx <- coords[ii, 1] - coords[jj, 1]
  dy <- coords[ii, 2] - coords[jj, 2]
  d_b <- Mod(complex(real = dx, imaginary = dy))          # == C hypot(dx, dy) bit-for-bit (checked)
  rm(dx, dy)
  bin <- findInterval(d_b, lims)                          # = binit's while loop (#edges <= d)
  keep <- (d_b <= max.dist) & (bin <= nbins)              # bin == nbins+1  <=>  d >= umax: dropped
  ii <- ii[keep]; jj <- jj[keep]; bin <- bin[keep]
  cnt <- tabulate(bin, nbins)
  indp <- cnt >= pairs.min                                # variogram.R:151
  out_bins <- seq_len(nbins)
  if (!nt.ind) { indp[1] <- FALSE }                       # variogram.R:154-159 (nugget bin removed)
  out_bins <- out_bins[indp]
  list(ii = ii, jj = jj, bin = bin, nbins = nbins, out_bins = out_bins,
       u = uvec_mid[out_bins], n = as.double(cnt[out_bins]),
       bins.lim = if (nt.ind) bins.lim else bins.lim[-1],
       umax = umax, max.dist = max.dist, nt.ind = nt.ind,
       npairs_total = length(d_b), npairs_maxdist = sum(d_b <= max.dist),
       npairs_used = length(ii))
}

## Z: N x B numeric matrix (or a vector). Returns list(u, n, v) with v an (nbins_out x B) matrix.
## Sums are accumulated with rowsum(), which adds rows sequentially in pair order into a double
## accumulator -- the same order and precision as binit, hence bit-identical.
variog_compute_r <- function(prep, Z, chunk = 64L) {
  Z <- as.matrix(Z)
  B <- ncol(Z)
  V <- matrix(NA_real_, length(prep$out_bins), B)
  for (s in seq(1L, B, by = chunk)) {
    cols <- s:min(B, s + chunk - 1L)
    D <- Z[prep$ii, cols, drop = FALSE] - Z[prep$jj, cols, drop = FALSE]
    D <- (D * D) / 2                                      # binit: v = (v*v)/2.0
    S <- rowsum(D, prep$bin, reorder = TRUE)              # sequential per-bin sums
    S <- S[match(prep$out_bins, as.integer(rownames(S))), , drop = FALSE]
    V[, cols] <- S / prep$n                               # binit: vbin[j] / cbin[j]
  }
  list(u = prep$u, n = prep$n, v = V)
}
