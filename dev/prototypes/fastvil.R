## Prototype: exact, vectorised re-implementation of STcompare's
## viladomatCorrelation()/matchingVariograms() using
##  (a) a precomputed linear smoother operator per delta (shared by all genes / both directions)
##  (b) a precomputed pair -> bin structure for the geoR "bin" variogram
##  (c) closed-form OLS
## It reproduces the package's RNG stream exactly (sample() for permutations,
## set.seed(seed + i) then rnorm(N) per delta for noise).
suppressPackageStartupMessages({ library(Rcpp) })

## --- exact replica of geoR::variog(option="bin") bin assignment (binit uses hypot) ---
cppFunction('
List pairBins(NumericVector xc, NumericVector yc, NumericVector lims, double maxdist) {
  int n = xc.size(), nbins = lims.size() - 1;
  std::vector<int> I, J, K;
  for (int j = 0; j < n; j++) for (int i = j + 1; i < n; i++) {
    double dist = hypot(xc[i] - xc[j], yc[i] - yc[j]);
    if (dist <= maxdist) {
      int ind = 0;
      while (ind < nbins && dist >= lims[ind]) ind++;
      if (dist < lims[ind]) { I.push_back(i + 1); J.push_back(j + 1); K.push_back(ind); }
    }
  }
  return List::create(_["i"] = I, _["j"] = J, _["bin"] = K);
}')

## Precompute the variogram structure for a coordinate subsample (mirrors geoR::variog defaults:
## uvec = 13 bins on [0, umax], nugget bin, pairs.min = 2)
precomputeVariog <- function(coords, max.dist) {
  u <- as.vector(dist(coords))
  umax <- max(u[u < max.dist])
  nt <- 1e-12
  bl <- seq(0, umax, l = 14)
  bins.lim <- c(0, nt, bl[bl > nt])
  nt.ind <- min(u) < nt
  lims <- bins.lim; if (lims[1] < 1e-16) lims[1] <- -1
  pb <- pairBins(coords[,1], coords[,2], lims, max.dist)
  nb <- length(lims) - 1
  cnt <- tabulate(pb$bin, nbins = nb)
  keep <- cnt >= 2
  if (!nt.ind) keep[1] <- FALSE
  list(i = pb$i, j = pb$j, bin = pb$bin, nbins = nb, cnt = cnt, keep = keep)
}

## Binned classical variogram for each column of Z (rows = subsample points)
variogCols <- function(Z, vs) {
  Z <- as.matrix(Z)
  D <- (Z[vs$i, , drop = FALSE] - Z[vs$j, , drop = FALSE])
  D <- D * D / 2
  S <- rowsum(D, vs$bin, reorder = TRUE)          # bins present, sorted
  present <- sort(unique(vs$bin))
  V <- matrix(NA_real_, vs$nbins, ncol(Z)); V[present, ] <- S
  V <- V / vs$cnt
  V[vs$keep, , drop = FALSE]
}

## Smoother operator for one delta, built from linearity of locfit's fitted values
buildS_unit <- function(long, lat, nn) {
  N <- length(long)
  S <- matrix(0, N, N)
  for (i in seq_len(N)) {
    e <- numeric(N); e[i] <- 1
    fit <- suppressWarnings(locfit::locfit(e ~ locfit::lp(long, lat, nn = nn, deg = 0), kern = "gauss", maxk = 300))
    S[, i] <- fitted(fit)
  }
  S
}

## Fast one-direction replica of viladomatCorrelation()
fastViladomat <- function(X, Y, coords, delta, Slist, vs, ids, prctile, B, seed = 0, rowsOnlyForSelection = TRUE) {
  N <- length(X); K <- length(delta)
  set.seed(seed)
  if (N > 1000) { ids2 <- sample(N, 1000) } else ids2 <- 1:N
  stopifnot(identical(ids2, ids))
  target <- variogCols(X[ids], vs)[, 1]
  Xr <- vapply(1:B, function(i) sample(X, size = N, replace = FALSE), numeric(N))  # N x B
  ## noise exactly as matchingVariograms: set.seed(seed+i); rnorm(N) for k = 1..K
  ## NB: BiocParallel::bplapply() evaluates matchingVariograms() with RNGkind("L'Ecuyer-CMRG"),
  ## so set.seed(seed + i) seeds L'Ecuyer-CMRG there (not Mersenne-Twister).
  E <- array(0, c(N, K, B))
  oldkind <- RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG")
  for (i in 1:B) { set.seed(seed + i); E[, , i] <- matrix(rnorm(N * K), N, K) }
  RNGkind(oldkind)
  tbar <- mean(target)
  RSS <- matrix(NA_real_, K, B); b0 <- b1 <- RSS
  for (k in 1:K) {
    Xd_s <- Slist[[k]][ids, , drop = FALSE] %*% Xr            # smoothed values at subsample only
    G <- variogCols(Xd_s, vs)                                  # nb x B
    gbar <- colMeans(G)
    Gc <- sweep(G, 2, gbar)
    b1[k, ] <- colSums(Gc * (target - tbar)) / colSums(Gc * Gc)
    b0[k, ] <- tbar - b1[k, ] * gbar
    H <- sweep(Xd_s, 2, sqrt(abs(b1[k, ])), `*`) + sweep(E[ids, k, ], 2, sqrt(abs(b0[k, ])), `*`)
    RSS[k, ] <- colSums((variogCols(H, vs) - target)^2)
  }
  dstar <- apply(RSS, 2, which.min)
  ## full-length smoothing only for the chosen delta of each permutation
  P <- matrix(0, N, B)
  for (k in unique(dstar)) {
    cols <- which(dstar == k)
    Xd <- Slist[[k]] %*% Xr[, cols, drop = FALSE]
    P[, cols] <- sweep(Xd, 2, sqrt(abs(b1[k, cols])), `*`) + sweep(matrix(E[, k, cols], N), 2, sqrt(abs(b0[k, cols])), `*`)
  }
  r0 <- as.vector(cor(X, Y)); rn <- cor(P, Y)
  list(deltaStarMedian = median(delta[dstar]), deltaStar = delta[dstar],
       pValueGlobal = sum(abs(rn) > abs(r0)) / B, nullCorGlobal = rn, permutations = P, RSS = RSS)
}

## Same as fastViladomat but with factorised operators Flist[[k]] = list(W = N x m, V = m x N)
fastViladomatF <- function(X, Y, delta, Flist, vs, ids, B, seed = 0) {
  N <- length(X); K <- length(delta)
  set.seed(seed)
  if (N > 1000) { ids2 <- sample(N, 1000) } else ids2 <- 1:N
  stopifnot(identical(ids2, ids))
  target <- variogCols(X[ids], vs)[, 1]
  Xr <- vapply(1:B, function(i) sample(X, size = N, replace = FALSE), numeric(N))
  E <- array(0, c(N, K, B))
  oldkind <- RNGkind()[1]; RNGkind("L'Ecuyer-CMRG")
  for (i in 1:B) { set.seed(seed + i); E[, , i] <- matrix(rnorm(N * K), N, K) }
  RNGkind(oldkind)
  tbar <- mean(target)
  RSS <- matrix(NA_real_, K, B); b0 <- b1 <- RSS
  VX <- vector("list", K)
  for (k in 1:K) {
    VX[[k]] <- Flist[[k]]$V %*% Xr                               # m x B (vertex fits)
    Xd_s <- Flist[[k]]$W[ids, , drop = FALSE] %*% VX[[k]]        # smoothed values at subsample
    G <- variogCols(Xd_s, vs)
    gbar <- colMeans(G); Gc <- sweep(G, 2, gbar)
    b1[k, ] <- colSums(Gc * (target - tbar)) / colSums(Gc * Gc)
    b0[k, ] <- tbar - b1[k, ] * gbar
    H <- sweep(Xd_s, 2, sqrt(abs(b1[k, ])), `*`) + sweep(matrix(E[ids, k, ], length(ids)), 2, sqrt(abs(b0[k, ])), `*`)
    RSS[k, ] <- colSums((variogCols(H, vs) - target)^2)
  }
  dstar <- apply(RSS, 2, which.min)
  P <- matrix(0, N, B)
  for (k in unique(dstar)) {
    cols <- which(dstar == k)
    Xd <- Flist[[k]]$W %*% VX[[k]][, cols, drop = FALSE]
    P[, cols] <- sweep(Xd, 2, sqrt(abs(b1[k, cols])), `*`) + sweep(matrix(E[, k, cols], N), 2, sqrt(abs(b0[k, cols])), `*`)
  }
  r0 <- as.vector(cor(X, Y)); rn <- cor(P, Y)
  list(deltaStarMedian = median(delta[dstar]), deltaStar = delta[dstar],
       pValueGlobal = sum(abs(rn) > abs(r0)) / B, nullCorGlobal = rn, permutations = P, RSS = RSS, b0 = b0, b1 = b1)
}

## one-time per coordinate set: subsample ids, variogram pair structure, factorised operators
prepareCoords <- function(coords, delta, maxDistPrctile = 0.25, seed = 0) {
  N <- nrow(coords)
  lat <- coords[,1]; long <- coords[,2]       # as in viladomatCorrelation (lat=data[,3], long=data[,4])
  set.seed(seed)
  ids <- if (N > 1000) sample(N, 1000) else 1:N
  dists <- dist(cbind(lat[ids], long[ids]))
  prctile <- quantile(dists, probs = maxDistPrctile)
  vs <- precomputeVariog(cbind(long[ids], lat[ids]), prctile)
  Flist <- lapply(delta, function(nn) factorS(long, lat, nn))
  list(ids = ids, prctile = prctile, vs = vs, Flist = Flist)
}
