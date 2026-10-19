## Exact factorisation S_delta = W %*% V of locfit's deg=0 gauss smoother (fitted values)
##  V: Nadaraya-Watson weights at the adaptive-tree vertices (bandwidth = k-NN distance of the vertex,
##     untruncated Gaussian exp(-(2.5 d/h)^2/2), k = floor(N*nn + 1e-12))
##  W: (multi)linear interpolation weights from fitted vertices to data points, extracted from
##     locfit's own 'sfitted' routine by feeding unit coefficient vectors (pseudo-vertices are
##     resolved inside locfit). locfit stores vertex values relative to a parametric component
##     (global mean for deg=0); that term cancels because interpolation weights sum to one.
factorS <- function(long, lat, nn) {
  N <- length(long)
  y0 <- rnorm(N)
  fit <- suppressWarnings(locfit::locfit(y0 ~ locfit::lp(long, lat, nn = nn, deg = 0), kern = "gauss", maxk = 300))
  nv <- fit$nvc[4]
  xev <- matrix(fit$eva$xev, ncol = 2, byrow = TRUE)[seq_len(nv), , drop = FALSE]
  h <- fit$eva$coef[seq_len(nv), 13]
  pv <- fit$cell$s[seq_len(nv)] == 1
  X <- cbind(long, lat)
  D2 <- outer(xev[,1], X[,1], "-")^2 + outer(xev[,2], X[,2], "-")^2
  Wk <- exp(-(2.5 * sqrt(D2) / h)^2 / 2)
  V <- Wk / rowSums(Wk)
  dat <- data.frame(y0 = y0, long = long, lat = lat)
  coef0 <- fit$eva$coef
  fz <- fit; fz$eva$coef[, 1] <- 0
  base <- fitted(fz, data = dat)            # parametric-component contribution only
  W <- matrix(0, N, nv)
  for (v in which(!pv)) {
    fz$eva$coef[, 1] <- 0; fz$eva$coef[v, 1] <- 1
    W[, v] <- fitted(fz, data = dat) - base
  }
  list(W = W[, !pv, drop = FALSE], V = V[!pv, , drop = FALSE], nv = nv, npv = sum(pv), h = h,
       vertex_check = max(abs((V[!pv, ] %*% y0 - mean(y0)) - coef0[which(!pv), 1])))
}
