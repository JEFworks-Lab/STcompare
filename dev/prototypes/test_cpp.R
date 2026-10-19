suppressPackageStartupMessages({ library(Rcpp); library(RcppArmadillo) })
source("fastvil.R"); source("factorS.R")
sourceCpp("fastgene.cpp")
delta <- seq(0.1, 0.9, 0.1)
d <- readRDS("sim_N5000.rds"); coords <- d$pos; N <- nrow(coords)
prep <- prepareCoords(coords, delta)
Ws <- lapply(prep$Flist, `[[`, "W"); Vs <- lapply(prep$Flist, `[[`, "V")
runCpp <- function(X, Y, B, seed = 0) {
  K <- length(delta); ids <- prep$ids; vs <- prep$vs
  t0 <- Sys.time()
  set.seed(seed); if (N > 1000) ids2 <- sample(N, 1000) else ids2 <- 1:N
  Xr <- vapply(1:B, function(i) sample(X, size = N, replace = FALSE), numeric(N))
  E <- array(0, c(N, K, B)); ok <- RNGkind()[1]; RNGkind("L'Ecuyer-CMRG")
  for (i in 1:B) { set.seed(seed + i); E[, , i] <- matrix(rnorm(N * K), N, K) }
  RNGkind(ok)
  t1 <- Sys.time()
  target <- variogCols(X[ids], vs)[, 1]
  out <- fastGeneCpp(Xr, E, Ws, Vs, ids, vs$i, vs$j, vs$bin, vs$cnt, vs$keep, target)
  t2 <- Sys.time()
  rn <- cor(out$permutations, Y); r0 <- cor(X, Y)
  list(p = sum(abs(rn) > abs(r0)) / B, deltaStar = delta[out$dstar], permutations = out$permutations,
       t_rng = as.numeric(t1 - t0, units = "secs"), t_cpp = as.numeric(t2 - t1, units = "secs"))
}
X <- d$X[2, ]; Y <- d$Y[2, ]
B <- 100
ref <- fastViladomatF(X, Y, delta, prep$Flist, prep$vs, prep$ids, B, seed = 0)   # exact R replica (validated vs package)
cc <- runCpp(X, Y, B)
cat(sprintf("N=%d B=%d one direction: C++ core %.3fs + R RNG %.3fs | deltaStar identical to exact R replica: %s | p %g vs %g | max|perm diff| %.3g\n",
    N, B, cc$t_cpp, cc$t_rng, identical(cc$deltaStar, ref$deltaStar), cc$p, ref$pValueGlobal, max(abs(cc$permutations - ref$permutations))))
## full gene: both directions
t0 <- Sys.time(); a <- runCpp(X, Y, 100); b <- runCpp(Y, X, 100); tg <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("full gene (2 directions x B=100): %.2fs total (C++ %.2fs, RNG in R %.2fs)\n", tg, a$t_cpp + b$t_cpp, a$t_rng + b$t_rng))
## also against the package directly at B=10 (independent check)
suppressPackageStartupMessages(devtools::load_all("/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare", quiet = TRUE))
pk <- viladomatCorrelation(data.frame(X = X, Y = Y, x = coords[,1], y = coords[,2]), delta = delta, maxDistPrctile = 0.25, nPermutations = 10, BPPARAM = BiocParallel::SerialParam(), seed = 0)
c10 <- runCpp(X, Y, 10)
cat(sprintf("vs package (B=10): deltaStar identical %s, p %g vs %g, max|perm diff| %.3g\n", identical(pk$deltaStar, c10$deltaStar), pk$pValueGlobal, c10$p, max(abs(pk$permutations - c10$permutations))))
