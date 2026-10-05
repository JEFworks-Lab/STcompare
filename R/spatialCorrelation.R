#' Spatially autocorrelated surrogates and permutation p-value for one direction
#'
#' @description Function to calculate Pearson's correlation between two spatial
#'   datasets, X and Y. To replace the analytical p-value which results in a
#'   high false positive rate for autocorrelated spatial patterns, it calculates
#'   empirical p-values from empirical null distributions generated from
#'   permuting dataset X by randomly shuffling the values and then smoothing to
#'   maintain the original degree of autocorrelation of X
#'
#' @details Each permutation shuffles the values of X across the locations,
#'   smooths the shuffled values with a Gaussian kernel whose nearest-neighbour
#'   bandwidth covers a proportion delta of the locations, and rescales the
#'   smoothed values and adds Gaussian noise so that their variogram matches the
#'   variogram of X (Viladomat et al. 2014). This is repeated for every delta,
#'   and the delta whose variogram matches best (delta star) gives the permuted
#'   field. If X has more than 1000 values, the variograms use a random
#'   subsample of 1000 locations.
#'
#'   The computation is done by compiled code (C++), which reimplements the
#'   local-constant kernel smoother of \code{locfit} (its adaptive k-d tree with
#'   vertex interpolation, \code{nn = delta}, \code{kern = "gauss"}) and the
#'   binned variogram of \code{geoR::variog()}. It gives the same results as the
#'   original R implementation: the same delta star for every permutation and
#'   the same null correlations up to floating-point rounding (about 1e-14 on
#'   the published analyses).
#'
#'   If the null cannot be computed (for example X is constant or has missing
#'   values, \eqn{N \times delta < 2}{N * delta < 2} for a delta, or the
#'   variogram has fewer than 2 bins), every element of the result is \code{NA}
#'   and a warning gives the reason. If Y has missing values or is constant, the
#'   permutations are returned, and \code{nullCorGlobal} and
#'   \code{pValueGlobal} are \code{NA} with a warning.
#'
#' @param data \code{matrix} An N x 4 matrix with the first column as the
#'   values of X, the second column as the values of Y, the third column as the
#'   x-coordinates, and the fourth column as the y-coordinates.
#'
#' @param delta \code{numeric vector} Given a point in dataset X, the percentage
#'   of neighbors which should be within the smoothing kernel. This can be a
#'   single numeric or a vector of candidate values; the best one is chosen for
#'   every permutation. Each delta must satisfy \eqn{N \times delta \ge 2}{N *
#'   delta >= 2}.
#'
#' @param maxDistPrctile \code{numeric}: Percentile of distances between pixels
#'   to use as max distance when calculating variograms. At greater distances
#'   the variogram is less precise because there are fewer pairs of points with
#'   that distance between them. Therefore, since the goal is to minimize the
#'   difference between the variogram of X and those of its permutations, the
#'   variogram should be subsetted to the percentile that is more robust.
#'
#' @param nPermutations \code{integer}: Number of permutations to generate to
#'   build the empirical null distribution. This number will determine the
#'   precision of the p-value. The smallest possible p-value is
#'   \eqn{1 / (nPermutations + 1)}, for example about 0.0099 when
#'   \code{nPermutations <- 100}
#'
#' @param nThreads \code{integer}: Number of threads of the compiled code.
#'   Default = 1. The permutations are distributed over the threads, and the
#'   results do not depend on the number of threads. We recommend the number of
#'   cores available (\code{parallel::detectCores(logical = FALSE)}).
#'
#' @param BPPARAM \code{BiocParallelParam}: Optional. If not \code{NULL}, its
#'   number of workers (\code{BiocParallel::bpnworkers(BPPARAM)}) is used as the
#'   number of threads instead of \code{nThreads}. No BiocParallel back-end is
#'   started (nothing is forked). Default is \code{NULL}.
#'
#' @param seed \code{integer}: Seed for the random number generator. Default
#'   \code{0}. The permutations are drawn after \code{set.seed(seed)} with the
#'   session's random number generator kinds (\code{RNGkind()}), and the noise
#'   of permutation \eqn{b} after \code{set.seed(seed + b)} with
#'   \code{"L'Ecuyer-CMRG"}. The results depend on the seed only, and the
#'   global random number generator state (\code{.Random.seed}) is left
#'   unchanged.
#'
#' @return The output is returned as a list.
#'  \describe{
#'   \item{\code{deltaStarMedian}}{numeric, the median of the deltas that minimize
#'   the residual sum of squares across each permutation}
#'   \item{\code{deltaStar}}{numeric vector of length \code{nPermutations}:
#'   the delta that minimizes the residual sum of squares for each permutation}
#'   \item{\code{pValueGlobal}}{numeric, empirical p-value for the Pearson's
#'   correlation of X and Y, computed as \eqn{(b + 1) / (B + 1)} where \eqn{b}
#'   is the number of null correlations whose absolute value is at least the
#'   absolute value of the observed correlation}
#'   \item{\code{nullCorGlobal}}{a B x 1 matrix, where B is \code{nPermutations}.
#'   This matrix is the correlation coefficients between the permutations and
#'   Y that compose that null distribution used to calculate the empirical p-value}
#'   \item{\code{permutations}}{an N x B matrix, where B is \code{nPermutations}.
#'   Each column is the resulting values of a permutation of X}
#'   }
#'
#' @references adapted from: Viladomat, Júlia et al. “Assessing the significance
#' of global and local correlations under spatial autocorrelation: a
#' nonparametric approach.” Biometrics vol. 70,2 (2014): 409-18.
#' doi:10.1111/biom.12139
#'
#' @export
#'
#' @examples
#'
#' data(quakes)
#'
#' #remove duplicated positions
#' quakes_data <- quakes[!duplicated(cbind(quakes$lat, quakes$long)),]
#'
#' data <- cbind(quakes_data$depth,
#'               quakes_data$mag,
#'               quakes_data$lat,
#'               quakes_data$long)
#'
#' #sequence of deltas
#' delta <- seq(0.1,0.9,0.1)
#'
#' # maximum distance for the variogram set at the 25% percentile of
#' # the distribution of pairs of distances:
#' maxDistPrctile <- 0.25
#'
#' #number of permutations
#' nPermutations <- 10
#'
#' resultsPermuteX <- viladomatCorrelation(data, delta, maxDistPrctile, nPermutations)
#' resultsPermuteX$pValueGlobal
#' resultsPermuteX$deltaStar
#'
viladomatCorrelation <- function(data, delta, maxDistPrctile, nPermutations,
                                 nThreads = 1, BPPARAM = NULL, seed = 0) {
  if (length(dim(data)) != 2L || ncol(data) < 4L) {
    stop("data must be a matrix or data frame with 4 columns: X, Y and the two coordinates")
  }
  X <- .stc_values(data[, 1], "data[, 1] (X)")
  Y <- .stc_values(data[, 2], "data[, 2] (Y)")
  pos <- cbind(data[, 3], data[, 4])
  nThreads <- .stc_threads(nThreads, BPPARAM)
  B <- .stc_check_nperm(nPermutations)

  e <- .stc_engine_correlate(X, Y, pos, deltaX = list(delta), nPermutations = B, seed = seed,
                             maxDistPrctile = maxDistPrctile, nThreads = nThreads,
                             returnPermutations = TRUE, mode = "forward")
  if (e$status[1] != "ok") {
    warning("viladomatCorrelation: no permutations (all results are NA): ", e$message[1], call. = FALSE)
    return(list(deltaStarMedian = NA_real_, deltaStar = NA_real_, pValueGlobal = NA_real_,
                nullCorGlobal = NA_real_, permutations = NA_real_))
  }
  if (nzchar(e$message[1])) {
    warning("viladomatCorrelation: nullCorGlobal or pValueGlobal is NA: ", e$message[1], call. = FALSE)
  }
  deltaStar <- e$deltaStarX[[1]]
  list(deltaStarMedian = stats::median(deltaStar),
       deltaStar = deltaStar,
       pValueGlobal = e$pX[1],
       nullCorGlobal = matrix(e$nullX[[1]], ncol = 1L),
       permutations = e$permutationsX[[1]])
}


#' Spatial correlation test of two vectors of values at the same locations
#'
#' @description Function to calculate Pearson's correlation between two spatial
#'   datasets. To replace the analytical p-value which results in a high false
#'   positive rate for autocorrelated spatial patterns, it calculates empirical
#'   p-values from empirical null distributions generated from permuting the
#'   datasets and then smoothing to maintain the original degree of
#'   autocorrelation
#'
#' @details The null distribution of the correlation is built twice: by
#'   permuting X while keeping Y fixed, and by permuting Y while keeping X fixed
#'   (see \code{\link{viladomatCorrelation}} for how a permutation keeps the
#'   autocorrelation). Both directions use the same seed, so they use the same
#'   shuffles of the locations. The empirical p-value of a direction is
#'   \eqn{(b + 1) / (B + 1)}, where \eqn{b} is the number of null correlations
#'   whose absolute value is at least the absolute value of the observed
#'   correlation and \eqn{B} is \code{nPermutations}; it is never 0.
#'
#'   The computation is done by compiled code (C++) that reimplements the
#'   \code{locfit} smoother and the \code{geoR} variogram of the original R
#'   implementation and gives the same results: the same delta star for every
#'   permutation, and the same null correlations up to floating-point rounding
#'   (about 1e-14 on the published analyses). Results do not depend on the
#'   number of threads, and the global random number generator state is left
#'   unchanged.
#'
#'   If the null cannot be computed in a direction, the columns computed from
#'   permutations are \code{NA} (an NA row) and a warning gives the reason.
#'   This happens when X or Y is constant or has missing values, when
#'   \eqn{N \times delta < 2}{N * delta < 2} for a delta (fewer than 2
#'   locations in the smoothing neighbourhood), when \code{maxDistPrctile} is
#'   so small that the variogram has fewer than 2 bins, or with exactly
#'   duplicated coordinates and a delta so small that the neighbourhood holds
#'   no more than the copies of a location. \code{correlationCoef} and
#'   \code{pValueNaive} are still computed by \code{cor.test()} where possible
#'   (on the complete pairs if there are missing values).
#'
#' @param X \code{numeric}: a numeric vector with N observations (a matrix
#'   with one row or one column is also accepted)
#'
#' @param Y \code{numeric}: a numeric vector with N observations, like
#'   \code{X}
#'
#' @param pos \code{matrix}: an N x 2 matrix of the spatial x,y coordinates of
#'   observations
#'
#' @param nPermutations \code{integer} or \code{double}: number of permutations
#'   to generate to build the empirical null distribution. This number will
#'   determine the precision of the p-value. Default is \code{100}, such that
#'   the smallest possible p-value is \eqn{1 / 101}, about 0.0099
#'
#' @param deltaX \code{numeric}: A single numeric or a numeric vector for
#'   controlling the degree of smoothing in permutations of X. Delta is a
#'   proportion calculated by dividing k neighbors by N total observations in X,
#'   where k is the number of neighbors in the permutation of X that should be
#'   within the radius smoothed by the Gaussian kernel to achieve the amount of
#'   autocorrelation present in the original X. If a single delta is not known,
#'   a sequence of deltas can be inputted and the best delta will be found such
#'   that it minimizes the sum of squares of the residuals between the variogram
#'   of the permutation generated from the delta and the variogram of the
#'   target. Default is \code{NULL}. If no value is supplied for \code{deltaX},
#'   \code{seq(0.1,0.9,0.1)}, the sequence of every 0.1 from 0.1 to 0.9, will be
#'   used to find the best delta for X.
#'
#' @param deltaY \code{numeric}: A single numeric or numeric vector for
#'   controlling the degree of smoothing in permutations of Y. \code{deltaY} is
#'   like \code{deltaX} but for observation in Y instead of X. Default is
#'   \code{NULL}. If no value is supplied for \code{deltaY},
#'   \code{seq(0.1,0.9,0.1)}, the sequence of every 0.1 from 0.1 to 0.9, will be
#'   used to find the best delta for permutations of Y.
#'
#' @param maxDistPrctile \code{numeric}: percentile of distances between pixels
#'   to use as max distance when calculating variograms. Default = 0.25. At
#'   greater distances the variogram is less precise because there are fewer
#'   pairs of points with that distance between them. Therefore, since the goal
#'   is to minimize the difference between the variogram of X and those of its
#'   permutations, the variogram should be subsetted to the percentile that is
#'   more robust.
#'
#' @param returnPermutations \code{logical}: \code{FALSE} (default) indicate
#'   whether the outputted dataframe will have a column with the values of the
#'   permutations used to calculate the null correlations and the empirical
#'   p-value.
#'
#' @param nThreads \code{integer}: Number of threads of the compiled code.
#'   Default = 1. The permutations are distributed over the threads, and the
#'   results do not depend on the number of threads. We recommend the number of
#'   cores available (\code{parallel::detectCores(logical = FALSE)}).
#'
#' @param BPPARAM \code{BiocParallelParam}: Optional. If not \code{NULL}, its
#'   number of workers (\code{BiocParallel::bpnworkers(BPPARAM)}) is used as the
#'   number of threads instead of \code{nThreads}. No BiocParallel back-end is
#'   started (nothing is forked). Default is \code{NULL}.
#'
#' @param seed \code{integer}: Seed for the random number generator. Default
#'   \code{0}. The permutations are drawn after \code{set.seed(seed)} with the
#'   session's random number generator kinds (\code{RNGkind()}), and the noise
#'   of permutation \eqn{b} after \code{set.seed(seed + b)} with
#'   \code{"L'Ecuyer-CMRG"}. The results depend on the seed only, and the
#'   global random number generator state (\code{.Random.seed}) is left
#'   unchanged.
#'
#' @return The output is returned as a \code{data.frame} with one row (named
#'   \code{"cor"}) containing the columns:
#' \describe{
#'   \item{\code{correlationCoef}}{Pearson's correlation coefficient.}
#'   \item{\code{pValueNaive}}{the analytical p-value naively assuming independent
#'   observations}
#'   \item{\code{pValuePermuteX}}{the empirical p-value when creating a null
#'   from permutations of observations in X, computed as \eqn{(b + 1) / (B + 1)}
#'   where \eqn{b} is the number of null correlations whose absolute value is
#'   at least the absolute value of the observed correlation}
#'   \item{\code{pValuePermuteY}}{the empirical p-value when creating a null from
#'   permutations of observations in Y, computed like \code{pValuePermuteX}}
#'   \item{\code{deltaStarMedianX}}{the median delta star (the delta which
#'   minimizes the difference between the variogram of the permutation and the
#'   variogram of observations) across permutations of X}
#'   \item{\code{deltaStarMedianY}}{the median delta star across permutations of Y}
#'   \item{\code{deltaStarX}}{list of delta star for all permutations of X}
#'   \item{\code{deltaStarY}}{list of delta star for all permutations of Y}
#'   \item{\code{nullCorrelationsX}}{list of a B x 1 matrix: the correlation
#'   coefficients for pairing Y and all permutations of X}
#'   \item{\code{nullCorrelationsY}}{list of a B x 1 matrix: the correlation
#'   coefficients for pairing X and all permutations of Y}
#'   \item{\code{permutationsX}}{(optional) an N x B matrix, where N is the
#'   length of X and B is \code{nPermutations}. Each column is the resulting values
#'   of a permutation of X}
#'   \item{\code{permutationsY}}{(optional) an N x B matrix, where N is the
#'   length of Y and B is \code{nPermutations}. Each column is the resulting values
#'   of a permutation of Y}
#'   }
#'
#' @export
#'
#' @examples
#'
#' data(quakes)
#'
#' #remove duplicated positions
#' quakes_data <- quakes[!duplicated(cbind(quakes$lat, quakes$long)),]
#'
#' cor <- spatialCorrelation(X = quakes_data$depth,
#'                           Y = quakes_data$mag,
#'                           pos = cbind(quakes_data$lat, quakes_data$long),
#'                           nThreads = 2)
#' cor
#'
#' # plot the delta star (the delta which minimizes the difference between the
#' # variogram of the permutation and the variogram of observations) for all
#' # permutations to see if clear peak found in the range inputted
#' hist(cor$deltaStarX[[1]])
#' hist(cor$deltaStarY[[1]])
#'
#' #plot null distribution of correlations to see if normally distributed
#' hist(cor$nullCorrelationsX[[1]])
#' hist(cor$nullCorrelationsY[[1]])
#'
#' #example of inputting specific range for deltaX and deltaY
#' cor2 <- spatialCorrelation(X = quakes_data$depth,
#'                           Y = quakes_data$mag,
#'                           pos = cbind(quakes_data$lat, quakes_data$long),
#'                           deltaX = seq(0.05, 0.9, 0.05),
#'                           deltaY = seq(0.02, 0.5, 0.02),
#'                           nThreads = 2)
#'
#' cor2
#'
#' hist(cor2$deltaStarX[[1]])
#' hist(cor2$deltaStarY[[1]])
#' hist(cor2$nullCorrelationsX[[1]])
#' hist(cor2$nullCorrelationsY[[1]])
#'
#' #visualizations of the spatial data to verify negative correlation
#' library(ggplot2)
#' p1 <- ggplot2::ggplot(quakes_data,
#'                       ggplot2::aes(x = long, y = lat, color = depth)) +
#' ggplot2::geom_point(size = 2) +
#'   ggplot2::scale_color_gradient(low = "lightblue", high = "blue") +
#'   ggplot2::labs(title = "Locations of Earthquakes off Fiji",
#'                 x = "Longitude", y = "Latitude",
#'                 color = "Depth (km)") +
#'   ggplot2::theme_minimal() +
#'   ggplot2::coord_quickmap()
#'
#' p2 <- ggplot2::ggplot(quakes_data,
#'                       ggplot2::aes(x = long, y = lat, color = mag)) +
#' ggplot2::geom_point(size = 2) +
#'   ggplot2::scale_color_gradient(low = "lightblue", high = "blue") +
#'   ggplot2::labs(title = "Locations of Earthquakes off Fiji",
#'                 x = "Longitude", y = "Latitude",
#'                 color = "Richter Magnitude") +
#'   ggplot2::theme_minimal() +
#'   ggplot2::coord_quickmap()
#'
#' p1
#' p2
spatialCorrelation <- function(X, Y, pos, nPermutations = 100,
                               deltaX = NULL, deltaY = NULL,
                               maxDistPrctile = 0.25,
                               returnPermutations = FALSE,
                               nThreads = 1, BPPARAM = NULL,
                               seed = 0){
  X <- .stc_values(X, "X")
  Y <- .stc_values(Y, "Y")
  if (length(X) != length(Y)) stop("X and Y must have the same length")
  if (length(X) < 3L) stop("X and Y must have at least 3 values")
  nThreads <- .stc_threads(nThreads, BPPARAM)
  B <- .stc_check_nperm(nPermutations)

  naive <- .stc_cor_tests(matrix(X), matrix(Y))
  e <- .stc_engine_correlate(X, Y, pos, deltaX = list(deltaX), deltaY = list(deltaY),
                             nPermutations = B, seed = seed, maxDistPrctile = maxDistPrctile,
                             nThreads = nThreads, returnPermutations = returnPermutations)
  .stc_warn_na_rows("spatialCorrelation", NULL, e, naive)
  .stc_result_table(naive, e, isTRUE(returnPermutations), "cor")
}

#' Spatial correlation test for every gene of two samples
#'
#' @description Function to calculate Pearson's correlation between assays from
#'   two SpatialExperiment datasets. To replace the analytical p-value which
#'   results in a high false positive rate for autocorrelated spatial patterns,
#'   it calculates empirical p-values from empirical null distributions
#'   generated from permuting the datasets and then smoothing to maintain the
#'   original degree of autocorrelation
#'
#' @details For every gene (row), the correlation of its values in the two
#'   datasets over their shared pixels is tested as in
#'   \code{\link{spatialCorrelation}}: with nulls from permutations of X and
#'   from permutations of Y. All genes are computed in one call of the compiled
#'   code, which shares the smoothing operators and the shuffles of the
#'   locations between genes and distributes the permutations of all genes
#'   over \code{nThreads} threads. The results are the same as for one gene at
#'   a time and do not depend on the number of threads.
#'
#'   Genes for which a null cannot be computed (see
#'   \code{\link{spatialCorrelation}}: for example a gene that is constant or
#'   has missing values on the shared pixels, or \eqn{N \times delta < 2}{N *
#'   delta < 2}) get \code{NA} in every column computed from permutations, and
#'   one warning per such gene gives the reason.
#'
#' @param input \code{list} List of two SpatialExperiment objects with matched
#'   spatial locations. The first element corresponds to the first
#'   SpatialExperiment (\code{X}), and the second to the second SpatialExperiment
#'   (\code{Y}). The SpatialCoords of the two SpatialExperiment objects should be on
#'   the same coordinate framework and observations at the same coordinate
#'   location in both datasets should have the same row names. If the
#'   SpatialExperiment objects do not have shared locations, use
#'   \code{SEraster::rasterizeGeneExpression()} to generate SpatialExperiment objects
#'   with shared pixel locations. See \code{assayName} parameter if the
#'   SpatialExperiment objects have more than one assay.
#'
#' @param nPermutations \code{integer} or \code{double}: number of permutations
#'   to generate to build the empirical null distribution. This number will
#'   determine the precision of the p-value. Default is \code{100}, such that
#'   the smallest possible unadjusted p-value is \eqn{1 / 101}, about 0.0099
#'
#' @param deltaX \code{list}: List of single numerics or list of numeric vectors
#'   to use for delta, the parameter controlling the degree of smoothing in
#'   permutations of X. The length of the list should be the same as the number of
#'   rows in the SpatialExperiment.  Delta is a proportion calculated by
#'   dividing k neighbors by N total observations (columns) in X, where k is the
#'   number of neighbors in the permutation of X that should be within the
#'   radius smoothed by the Gaussian kernel to achieve the amount of
#'   autocorrelation present in the original X. If a single delta is not known,
#'   a sequence of deltas can be inputted and the best delta will be found such
#'   that it minimizes the sum of squares of the residuals between the variogram
#'   of the permutation generated from the delta and the variogram of the
#'   target. Default is \code{NULL}. If no value is supplied for \code{deltaX},
#'   \code{seq(0.1,0.9,0.1)}, the sequence of every 0.1 from 0.1 to 0.9, will be
#'   used to find the best delta for each row (gene) in X.
#'
#' @param deltaY \code{list}: List of single numerics or list of numeric vectors
#'   to use for delta, the parameter controlling the degree of smoothing in
#'   permutations of Y. \code{deltaY} is like \code{deltaX} but for permuting
#'   data in Y instead of X. Default is \code{NULL}. If no value is supplied for
#'   \code{deltaY}, \code{seq(0.1,0.9,0.1)}, the sequence of every 0.1 from 0.1
#'   to 0.9, will be used to find the best delta for permutations for each row
#'   (gene) in Y.
#'
#' @param maxDistPrctile \code{numeric}: percentile of distances between pixels
#'   to use as max distance when calculating variograms. Default = 0.25. At
#'   greater distances the variogram is less precise because there are fewer
#'   pairs of points with that distance between them. Therefore, since the goal
#'   is to minimize the difference between the variogram of X and those of its
#'   permutations, the variogram should be subsetted to the percentile that is
#'   more robust.
#'
#' @param returnPermutations \code{logical}: indicate whether the dataframe
#'   returned as output will have a column with the values of the permutations
#'   used to calculate the null correlations and the empirical p-value. Default
#'   is \code{FALSE}
#'
#' @param assayName \code{character} or \code{integer} A character string or
#'   numeric specifying the assay in the SpatialExperiment to use. Default is
#'   \code{NULL}. If no value is supplied for \code{assayName}, then the first
#'   assay is used as a default
#'
#' @param nThreads \code{integer}: Number of threads of the compiled code.
#'   Default = 1. The permutations of all genes are distributed over the
#'   threads, and the results do not depend on the number of threads. We
#'   recommend the number of cores available
#'   (\code{parallel::detectCores(logical = FALSE)}).
#'
#' @param BPPARAM \code{BiocParallelParam}: Optional. If not \code{NULL}, its
#'   number of workers (\code{BiocParallel::bpnworkers(BPPARAM)}) is used as the
#'   number of threads instead of \code{nThreads}. No BiocParallel back-end is
#'   started (nothing is forked). Default is \code{NULL}.
#'
#' @param verbose \code{logical}: if \code{TRUE} (default), print a message
#'   when the computation starts (genes, permutations and threads) and when it
#'   ends (elapsed time).
#'
#' @param seed \code{integer}: Seed for the random number generator. Default
#'   \code{0}. The permutations are drawn after \code{set.seed(seed)} with the
#'   session's random number generator kinds (\code{RNGkind()}), and the noise
#'   of permutation \eqn{b} after \code{set.seed(seed + b)} with
#'   \code{"L'Ecuyer-CMRG"}; every gene uses the same seed. The results depend
#'   on the seed only, and the global random number generator state
#'   (\code{.Random.seed}) is left unchanged.
#'
#' @param adjustMethod \code{character}: multiple-testing correction method
#'   passed to \code{stats::p.adjust()}. It is applied across all genes,
#'   separately to the \code{pValuePermuteX} column and to the
#'   \code{pValuePermuteY} column. Must be one of \code{p.adjust.methods}; use
#'   \code{"none"} for unadjusted p-values. Default is \code{"BH"}.
#'
#' @return The output is returned as a \code{data.frame}. The rownames are the
#'   rownames of the SpatialExperiments. The names of the columns and their
#'   contents are as follows:
#' \describe{
#'   \item{\code{correlationCoef}}{Pearson's correlation coefficient.}
#'   \item{\code{pValueNaive}}{the analytical p-value naively assuming independent
#'   observations}
#'   \item{\code{pValuePermuteX}}{the empirical p-value when creating a null from
#'   permutations of observations in X, computed as \eqn{(b + 1) / (B + 1)} (see
#'   \code{\link{spatialCorrelation}}) and then adjusted across genes with
#'   \code{adjustMethod}}
#'   \item{\code{pValuePermuteY}}{the empirical p-value when creating a null from
#'   permutations of observations in Y, adjusted across genes like
#'   \code{pValuePermuteX}}
#'   \item{\code{deltaStarMedianX}}{the median delta star (the delta which
#'   minimizes the difference between the variogram of the permutation and the
#'   variogram of observations) across permutations of X}
#'   \item{\code{deltaStarMedianY}}{the median delta star across permutations of Y}
#'   \item{\code{deltaStarX}}{list of delta star for all permutations of X}
#'   \item{\code{deltaStarY}}{list of delta star for all permutations of Y}
#'   \item{\code{nullCorrelationsX}}{list of B x 1 matrices: the correlation
#'   coefficients for Y and all permutations of X}
#'   \item{\code{nullCorrelationsY}}{list of B x 1 matrices: the correlation
#'   coefficients for X and all permutations of Y}
#'   \item{\code{permutationsX}}{(optional) an N x B matrix, where N is the length of X and B is \code{nPermutations}.
#'   Each column is the resulting values of a permutation of X}
#'   \item{\code{permutationsY}}{(optional) an N x B matrix, where N is the length of Y and B is \code{nPermutations}.
#'   Each column is the resulting values of a permutation of Y}
#'   }
#'
#' @export
#'
#' @examples
#'
#' data(speKidney)
#'
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                square = FALSE)
#'
#' ##### Use STcompare to calculate Pearson's correlation coefficient #####
#' rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
#' rastGexpListAC <- list(A = rastKidney$A, C = rastKidney$C)
#'
#' negCorrelation <- spatialCorrelationGeneExp(rastGexpListAB, nThreads = 2)
#' posCorrelation <- spatialCorrelationGeneExp(rastGexpListAC, nThreads = 2)
#'
#' negCorrelation
#' posCorrelation

spatialCorrelationGeneExp <- function(input, nPermutations = 100,
                                      deltaX = NULL, deltaY = NULL,
                                      maxDistPrctile = 0.25,
                                      returnPermutations = FALSE,
                                      assayName = NULL,
                                      nThreads = 1, BPPARAM = NULL,
                                      verbose = TRUE,
                                      seed = 0,
                                      adjustMethod = "BH"){

  # correction method should be from p.adjust.methods
  if (!adjustMethod %in% stats::p.adjust.methods) {
    stop("adjustMethod must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "))
  }
  nThreads <- .stc_threads(nThreads, BPPARAM)
  B <- .stc_check_nperm(nPermutations)
  d <- .stc_pair_input(input, assayName)
  G <- length(d$genes)
  deltaX <- .stc_gene_deltas(deltaX, G, "deltaX")
  deltaY <- .stc_gene_deltas(deltaY, G, "deltaY")

  t0 <- proc.time()
  if (verbose) {
    message(sprintf("spatialCorrelationGeneExp: %d gene(s) x 2 directions on %d shared pixels, %d permutations, %d thread(s)",
                    G, nrow(d$pos), B, nThreads))
  }
  naive <- .stc_cor_tests(d$X, d$Y)
  e <- .stc_engine_correlate(d$X, d$Y, d$pos, deltaX = deltaX, deltaY = deltaY, nPermutations = B,
                             seed = seed, maxDistPrctile = maxDistPrctile, nThreads = nThreads,
                             returnPermutations = returnPermutations)
  n_na <- .stc_warn_na_rows("spatialCorrelationGeneExp", paste0("gene ", d$genes), e, naive)
  out <- .stc_result_table(naive, e, isTRUE(returnPermutations), make.unique(d$genes, sep = ""))

  # mht correct for pValuePermuteX and pValuePermuteY separately, across genes
  out$pValuePermuteX <- stats::p.adjust(out$pValuePermuteX, method = adjustMethod)
  out$pValuePermuteY <- stats::p.adjust(out$pValuePermuteY, method = adjustMethod)
  if (verbose) {
    message(sprintf("spatialCorrelationGeneExp: done in %s (%d of %d gene(s) with NA permutation p-values)",
                    .stc_elapsed(t0), n_na, G))
  }
  out
}

#' Spatial correlation test between pairs of genes of one sample
#'
#' @description Function to calculate Pearson's correlation between rows from
#'   one SpatialExperiment dataset. To replace the analytical p-value which
#'   results in a high false positive rate for autocorrelated spatial patterns,
#'   it calculates empirical p-values from empirical null distributions
#'   generated from permuting the data and then smoothing to maintain the
#'   original degree of autocorrelation
#'
#' @details Every pair of rows (genes) is tested as in
#'   \code{\link{spatialCorrelation}}, with the first gene of the pair as X and
#'   the second as Y. A gene's permutations do not depend on its partner, so
#'   the compiled code computes the permutations of every gene once and
#'   correlates them with every other gene, on \code{nThreads} threads. The
#'   results are the same as for one pair at a time and do not depend on the
#'   number of threads. Pairs with a gene for which a null cannot be computed
#'   (for example a constant gene) get \code{NA} in every column computed from
#'   permutations, and one warning per such pair gives the reason. The
#'   empirical p-values are not adjusted for multiple testing.
#'
#' @param input \code{SpatialExperiment} A SpatialExperiment object. See
#'   \code{assayName} parameter if the SpatialExperiment object has more than
#'   one assay.
#'
#' @param nPermutations \code{integer} or \code{double}: number of permutations
#'   to generate to build the empirical null distribution. This number will
#'   determine the precision of the p-value. Default is \code{100}, such that
#'   the smallest possible p-value is \eqn{1 / 101}, about 0.0099
#'
#' @param delta \code{list}: List of single numerics or list of numeric vectors
#'   to use for delta, the parameter controlling the degree of smoothing in
#'   permutations of each row (gene). The length of the list should be the same as
#'   the number of rows in the SpatialExperiment. Delta is a proportion
#'   calculated by dividing k neighbors by N total observations (columns), where
#'   k is the number of neighbors in the permutation that should be within the
#'   radius smoothed by the Gaussian kernel to achieve the amount of
#'   autocorrelation present in the original data. If a single delta is not
#'   known, a sequence of deltas can be inputted and the best delta will be
#'   found such that it minimizes the sum of squares of the residuals between
#'   the variogram of the permutation generated from the delta and the
#'   variogram of the target. Default is \code{NULL}. If no value is supplied
#'   for \code{delta}, \code{seq(0.1,0.9,0.1)}, the sequence of every 0.1 from
#'   0.1 to 0.9, will be used to find the best delta for each row (gene).
#'
#' @param maxDistPrctile \code{numeric}: percentile of distances between pixels
#'   to use as max distance when calculating variograms. Default = 0.25. At
#'   greater distances the variogram is less precise because there are fewer
#'   pairs of points with that distance between them. Therefore, since the goal
#'   is to minimize the difference between the variogram of X and those of its
#'   permutations, the variogram should be subsetted to the percentile that is
#'   more robust.
#'
#' @param returnPermutations \code{logical}: indicate whether the dataframe
#'   returned as output will have a column with the values of the permutations
#'   used to calculate the null correlations and the empirical p-value. Default
#'   is \code{FALSE}
#'
#' @param assayName \code{character} or \code{integer} A character string or
#'   numeric specifying the assay in the SpatialExperiment to use. Default is
#'   \code{NULL}. If no value is supplied for \code{assayName}, then the first
#'   assay is used as a default
#'
#' @param nThreads \code{integer}: Number of threads of the compiled code.
#'   Default = 1. The permutations of all genes are distributed over the
#'   threads, and the results do not depend on the number of threads. We
#'   recommend the number of cores available
#'   (\code{parallel::detectCores(logical = FALSE)}).
#'
#' @param BPPARAM \code{BiocParallelParam}: Optional. If not \code{NULL}, its
#'   number of workers (\code{BiocParallel::bpnworkers(BPPARAM)}) is used as the
#'   number of threads instead of \code{nThreads}. No BiocParallel back-end is
#'   started (nothing is forked). Default is \code{NULL}.
#'
#' @param verbose \code{logical}: if \code{TRUE} (default), print a message
#'   when the computation starts (genes, pairs, permutations and threads) and
#'   when it ends (elapsed time).
#'
#' @param seed \code{integer}: Seed for the random number generator. Default
#'   \code{0}. The permutations are drawn after \code{set.seed(seed)} with the
#'   session's random number generator kinds (\code{RNGkind()}), and the noise
#'   of permutation \eqn{b} after \code{set.seed(seed + b)} with
#'   \code{"L'Ecuyer-CMRG"}; every gene uses the same seed. The results depend
#'   on the seed only, and the global random number generator state
#'   (\code{.Random.seed}) is left unchanged.
#'
#' @return The output is returned as a \code{data.frame} with one row per pair
#'   of rows (genes), in the order of \code{combn(rownames(input), 2)}. The
#'   rownames are arbitrary (\code{"cor"}, \code{"cor1"}, ...). The names of the
#'   columns and their contents are as follows:
#' \describe{
#'   \item{\code{correlationCoef}}{Pearson's correlation coefficient.}
#'   \item{\code{pValueNaive}}{the analytical p-value naively assuming independent
#'   observations}
#'   \item{\code{pValuePermuteX}}{the p-value when creating an empirical null
#'   from permutations of the first gene of the pair, computed as
#'   \eqn{(b + 1) / (B + 1)} (see \code{\link{spatialCorrelation}})}
#'   \item{\code{pValuePermuteY}}{the p-value when creating an empirical null from
#'   permutations of the second gene of the pair}
#'   \item{\code{deltaStarMedianX}}{the median delta star (the delta which
#'   minimizes the difference between the variogram of the permutation and the
#'   variogram of observations) across permutations of the first gene}
#'   \item{\code{deltaStarMedianY}}{the median delta star across permutations
#'   of the second gene}
#'   \item{\code{deltaStarX}}{list of delta star for all permutations of the
#'   first gene}
#'   \item{\code{deltaStarY}}{list of delta star for all permutations of the
#'   second gene}
#'   \item{\code{nullCorrelationsX}}{list of B x 1 matrices: the correlation
#'   coefficients for the second gene and all permutations of the first}
#'   \item{\code{nullCorrelationsY}}{list of B x 1 matrices: the correlation
#'   coefficients for the first gene and all permutations of the second}
#'   \item{\code{permutationsX}}{(optional) an N x B matrix, where N is the
#'   number of pixels and B is \code{nPermutations}. Each column is the resulting
#'   values of a permutation of the first gene}
#'   \item{\code{permutationsY}}{(optional) the same for the second gene}
#'   \item{\code{first}}{the name of the first row in the pair}
#'   \item{\code{second}}{the name of the second row in the pair}
#'   }
#'
#' @export
#'
#' @examples
#'
#' data(speKidney)
#'
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                square = FALSE)
#'
#' # one SpatialExperiment whose three rows are the gene in samples A, C and
#' # B, on the pixels the three samples share
#' shared <- Reduce(intersect, lapply(rastKidney, colnames))
#' expr <- t(sapply(rastKidney, function(s) {
#'   SummarizedExperiment::assay(s)[1, shared]
#' }))
#' speACB <- SpatialExperiment::SpatialExperiment(
#'   assays = list(pixelval = expr),
#'   spatialCoords = SpatialExperiment::spatialCoords(rastKidney$A)[shared, ])
#'
#' sc_within_sample <- spatialCorrelationGeneExpWithinSample(speACB,
#'                                                           nThreads = 2)
#' sc_within_sample[, c("first", "second", "correlationCoef",
#'                      "pValuePermuteX", "pValuePermuteY")]
#'
spatialCorrelationGeneExpWithinSample <- function(input,
                                                  nPermutations = 100,
                                                  delta = NULL,
                                                  maxDistPrctile = 0.25,
                                                  returnPermutations = FALSE,
                                                  assayName = NULL,
                                                  nThreads = 1,
                                                  BPPARAM = NULL,
                                                  verbose = TRUE,
                                                  seed = 0){
  nThreads <- .stc_threads(nThreads, BPPARAM)
  B <- .stc_check_nperm(nPermutations)
  if (is.null(assayName)) {
    assayName <- 1
  }
  genes <- rownames(input)
  if (is.null(genes) || anyNA(genes)) stop("input must have row names (genes)")
  G <- length(genes)
  if (G < 2L) stop("input must have at least 2 rows (genes)")
  pos <- SpatialExperiment::spatialCoords(input)
  M <- t(as.matrix(SummarizedExperiment::assay(input, assayName)))
  storage.mode(M) <- "double"
  colnames(M) <- genes
  delta <- .stc_gene_deltas(delta, G, "delta")
  pairs <- utils::combn(G, 2L)

  t0 <- proc.time()
  if (verbose) {
    message(sprintf("spatialCorrelationGeneExpWithinSample: %d genes (%d pairs) on %d pixels, %d permutations, %d thread(s)",
                    G, ncol(pairs), nrow(M), B, nThreads))
  }
  e <- .stc_engine_correlate(M, pos = pos, deltaX = delta, nPermutations = B, seed = seed,
                             maxDistPrctile = maxDistPrctile, nThreads = nThreads,
                             returnPermutations = returnPermutations, mode = "within")
  # cor.test() of every pair, with the engine's cor(M[, i], M[, j]) as the estimate where cor.test()
  # would compute that same value
  naive <- .stc_cor_tests(M, M, pairs[1, ], pairs[2, ], r_pairs = attr(e, "state")$r[t(pairs)])
  labels <- sprintf("genes %s and %s", genes[pairs[1, ]], genes[pairs[2, ]])
  n_na <- .stc_warn_na_rows("spatialCorrelationGeneExpWithinSample", labels, e, naive)
  out <- .stc_result_table(naive, e, isTRUE(returnPermutations),
                           make.unique(rep("cor", ncol(pairs)), sep = ""))
  out$first <- genes[pairs[1, ]]
  out$second <- genes[pairs[2, ]]
  if (verbose) {
    message(sprintf("spatialCorrelationGeneExpWithinSample: done in %s (%d of %d pairs with NA permutation p-values)",
                    .stc_elapsed(t0), n_na, ncol(pairs)))
  }
  out
}

#' Plot a gene's values in two samples at the shared pixels
#'
#' @description Scatter plot of a gene's values in two SpatialExperiment
#'   objects on their shared pixels, with the correlation coefficient and the
#'   permutation p-value of the gene in the title.
#'
#' @details The pixels are matched by name, as in
#'   \code{\link{spatialCorrelationGeneExp}()}. Both axes start at 0, or at the
#'   smallest value if a value is negative, and end at the largest value, so
#'   every pixel with values in both objects is drawn.
#'
#' @param speList \code{list} List of two SpatialExperiment objects with matched
#'   spatial locations. The first element corresponds to the first
#'   SpatialExperiment (\code{X}), and the second to the second
#'   SpatialExperiment (\code{Y}). The names of the list label the axes.
#'
#' @param spatialCorrelation \code{dataframe}: the output of
#'   \code{\link{spatialCorrelationGeneExp}()} or
#'   \code{\link{spatialCorrelationGeneExpIterPermutations}()} for
#'   \code{speList}, or of \code{\link{compareSpatial}()} for the same objects.
#'
#' @param geneName \code{character}: The name of the gene (row) in both
#'   SpatialExperiment objects of \code{speList} and in
#'   \code{spatialCorrelation}.
#'
#' @param assayName \code{character} or \code{integer} A character string or
#'   numeric specifying the assay in the SpatialExperiment to use. Default is
#'   \code{NULL}. If no value is supplied for \code{assayName}, then the first
#'   assay is used as a default
#'
#' @return A ggplot object: a scatterplot with the values of the gene in the
#'   first SpatialExperiment, i.e
#'   \code{SummarizedExperiment::assay(speList[[1]], assayName)[geneName, ]},
#'   on the x-axis and values of the gene in the second SpatialExperiment on
#'   the y-axis. The title includes the gene name, the correlation coefficient
#'   (rounded to 3 decimals), and the empirical p-value (3 significant digits):
#'   \code{p_E}, the greater of \code{pValuePermuteX} and
#'   \code{pValuePermuteY} (\code{NA} if either is \code{NA}) for the output
#'   of the legacy functions, or \code{padj} for a \code{compareSpatial()}
#'   result.
#'
#'
#' @export
#'
#' @examples
#'
#' data(speKidney)
#'
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                square = FALSE)
#'
#' ##### Use STcompare to calculate Pearson's correlation coefficient #####
#' rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
#' rastGexpListAC <- list(A = rastKidney$A, C = rastKidney$C)
#'
#' negCorrelation <- spatialCorrelationGeneExp(rastGexpListAB, nThreads = 2)
#' posCorrelation <- spatialCorrelationGeneExp(rastGexpListAC, nThreads = 2)
#'
#' negCorrelation
#' posCorrelation
#'
#' expAB <- plotCorrelationGeneExp(rastGexpListAB, negCorrelation, "Gene")
#' expAC <- plotCorrelationGeneExp(rastGexpListAC, posCorrelation, "Gene")
#'
#' expAB
#' expAC
plotCorrelationGeneExp <- function(speList,
                                   spatialCorrelation,
                                   geneName,
                                   assayName = NULL){

  ## if name of assay to use in the SpatialExperiment object is not provided,
  # use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }
  if (!is.character(geneName) || length(geneName) != 1L || is.na(geneName)) {
    stop("geneName must be one gene name")
  }
  if (!geneName %in% rownames(spatialCorrelation)) {
    stop(sprintf("gene %s is not a row of spatialCorrelation", geneName))
  }

  #store names for labeling axes and title
  nameList <- names(speList)
  if (is.character(assayName)) {
    assayNameChar <- assayName
  }
  else{
    assayNameChar <- SummarizedExperiment::assayNames(speList[[1]])[assayName]
  }

  # the correlation and the p-value shown in the title: for the legacy functions, the greater of the two
  # empirical p-values (NA if either is NA); for compareSpatial(), the adjusted permutation p-value
  res <- spatialCorrelation[geneName, , drop = FALSE]
  if (all(c("r", "padj") %in% names(res))) {
    r <- res$r
    p <- res$padj
    pLabel <- "padj"
  } else {
    r <- res$correlationCoef
    p <- max(res$pValuePermuteX, res$pValuePermuteY)
    pLabel <- "p_E"
  }

  #Determine the positions of shared pixels between two rasterized spatial
  #experiments
  Y <- speList[[2]]
  X <- speList[[1]]
  for (s in list(X, Y)) {
    if (!geneName %in% rownames(s)) stop(sprintf("gene %s is not in both SpatialExperiment objects of speList", geneName))
  }
  sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(Y)),
                            rownames(SpatialExperiment::spatialCoords(X)))

  rastDf <- data.frame(YGexp = SummarizedExperiment::assay(Y, assayName)[geneName, sharedPixels],
                       XGexp = SummarizedExperiment::assay(X, assayName)[geneName, sharedPixels])
  # the same range on both axes, from 0 (or the smallest value, if negative) to the largest value
  lims <- range(0, rastDf$XGexp, rastDf$YGexp, na.rm = TRUE)

  pltGexp <- ggplot2::ggplot(data = rastDf, mapping = ggplot2::aes(x = .data$XGexp, y = .data$YGexp)) +
    ggplot2::geom_point(alpha = 0.5, size = 1, na.rm = TRUE) +
    ggplot2::ylim(lims) +
    ggplot2::xlim(lims) +
    ggplot2::theme_classic() + ggplot2::theme(legend.position="right") +
    ggplot2::labs(x = paste0(assayNameChar, " in ", nameList[1]),
                  y = paste0(assayNameChar, " in ", nameList[2]),
                  title = paste(geneName, " r = ",
                                round(r, 3),
                                paste0(" ", pLabel, " = "),
                                signif(p, 3)))
  pltGexp
}
