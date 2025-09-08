# matching function.
matchingVariograms <- function( X.randomized, long, lat, delta, target_variog, prctile, ids, i ) {
  
  #initialize vectors for fitting variograms for each delta
  variog.X.delta <- vector(mode = "list", length = length(delta))
  linear.fit <- variog.X.delta
  hat.X.delta <- variog.X.delta
  resid.sum.squares <- rep(0, length(delta))
  
  for (k in 1:length(delta)) { 
    # smooth X.randomized using locfit:
    fit <- locfit(X.randomized ~ lp(long, lat, nn = delta[k], deg = 0), kern = "gauss", maxk = 300)
    X.delta <- fitted(fit) 
    # variogram of X.delta:
    variog.X.delta[[k]] <- variog(data = X.delta[ids], coords = cbind(long[ids],lat[ids]), option = "bin", max.dist = prctile, messages = FALSE)
    
    # linear regression between the target and X.delta variograms:
    linear.fit[[k]] <- lm(target_variog$v ~ 1 + variog.X.delta[[k]]$v)  
    # least square estimates:
    bet.hat <- as.numeric(linear.fit[[k]]$coefficients)
    # transformed X.delta:
    hat.X.delta[[k]] <- X.delta * sqrt(abs(bet.hat[2])) + rnorm(length(X.delta)) * sqrt(abs(bet.hat[1]))
    variog.hat.X.delta <- variog(data = hat.X.delta[[k]][ids], coords = cbind(long[ids],lat[ids]), option = "bin", max.dist = prctile, messages = FALSE)
    
    # sum of squares of the residuals:
    resid.sum.squares[k] <- sum((variog.hat.X.delta$v-target_variog$v) ^ 2)
  }
  # delta that minimizes the residual sum of squares:
  delta.star.id <- which.min(resid.sum.squares)
  # permutation for which its variogram matches the variogram of X:
  hat.X.delta.star <- hat.X.delta[[delta.star.id]]
  return(list(residus = resid.sum.squares, delta.star.id = delta.star.id, hat.X.delta.star = hat.X.delta.star ))
}

viladomatCorrelation <- function(data, delta, maxDistPrctile, BPPARAM = BPPARAM, nPermutations) {
  
  suppressMessages(library(geoR))
  suppressMessages(library(locfit))
  library(foreach)
  #suppressMessages(library(doParallel)) #do i need this?
  #registerDoParallel(cores=workers) # do i need this?
  
  # load the data X, Y and the coordinates of the N data locations:
  X <- data[,1]
  Y <- data[,2]
  lat <- data[,3]
  long <- data[,4]
  N <- length(X)
  
  # If N is too big, we take a subsample of size N_s every time we calculate a variogram:
  N_s <- 1000
  if (length(X) > N_s) {
    ids <- sample(N, N_s)
    X_s <- X[ids]
    long_s <- long[ids]
    lat_s <- lat[ids]
  } else {
    X_s <- X
    long_s <- long
    lat_s <- lat
    ids <- 1:N
  }
  
  # maximum distance for the variogram set at the 25% percentile of 
  # the distribution of pairs of distances:
  dists <- dist(cbind(lat_s, long_s))
  prctile <- quantile(dists, probs = maxDistPrctile)
  
  # variogram of variable X that will be used as target when doing the matching:
  target_variog <- variog(data = X_s, coords = cbind(long_s,lat_s), max.dist = prctile, option = "bin", messages = FALSE)
  
  # ALGORITHM:
  # It returns B random fields with the same autocorrelation 
  # as X but independent of Y, stored in permutations. The basis
  # to calculate B realizations of the null we are interested in.
  
  B <- nPermutations
  permutations <- vector(mode = "list", length = B)
  
  
  # random permutation of the values of X across locations:
  X.randomized <- lapply(1:B, function(i) {
    sample(X, size = length(X), replace = FALSE)
  })
  
  output <- BiocParallel::bplapply(1:B, function(i) {
    # smoothing and scaling step to match the target variogram:
    output <- matchingVariograms(X.randomized[[i]], long, lat, delta, target_variog, prctile, ids, i) 
  }, BPPARAM=BPPARAM)
  
  
  permutations <- do.call(cbind, BiocParallel::bplapply(1:B, function(i) {
    # store permutations after smoothing and scaling to match the target variogram                        
    hat.X.delta.star <- output[[i]]$hat.X.delta.star
    c(hat.X.delta.star)
  }, BPPARAM=BPPARAM)
  )
  
  delta.star <- do.call(cbind, BiocParallel::bplapply(1:B, function(i) {
    # index of delta that minimizes the residual sum of squares between between 
    # the target and X.delta variograms for each permutation
    delta.star.id <- output[[i]]$delta.star.id
    # evaluate the delta at that index
    c(delta[delta.star.id])
  }, BPPARAM=BPPARAM)
  )
  
  # store median of the delta stars to return
  delta.star.median<- median(delta.star)
  
  # ASSESSING THE SINGLE PEARSON'S CORRELATION COEFFICIENT (GLOBAL CORRELATION):
  # observed global correlation:
  cor.global.obs <- as.vector(cor(X,Y))
  
  # null distribution for the global correlation:
  cor.global <- cor(permutations,Y)
  
  # p-value:
  extreme <- sum(abs(cor.global) > abs(cor.global.obs))
  p.value.global <- extreme / B
  
  return(list(deltaStarMedian = delta.star.median, 
              deltaStar = c(delta.star),
              pValueGlobal = p.value.global,
              nullCorGlobal = cor.global,
              permutations = permutations))
}


#' spatialCorrelation
#'
#' @description Function to calculate Pearson's correlation between two spatial
#'   datasets. To replace the analytical p-value which results in a high false
#'   positive rate for autocorrelated spatial patterns, it calculates empirical
#'   p-values from empirical null distributions generated from permuting the
#'   datasets and then smoothing to maintain the original degree of
#'   autocorrelation
#'
#' @param X \code{numeric} or \code{matrix}: a 1 x N numeric vector or matrix
#'   with N observations
#'
#' @param Y \code{numeric} or \code{matrix}: a 1 x N numeric vector or matrix
#'   with N observations
#'
#' @param pos \code{matrix}: a N x 2 matrix array of spatial x,y coordinates of
#'   observations
#'
#' @param nPermutations \code{integer} or \code{double}: number of permutations
#'   to generate to build the empirical null distribution. This number will
#'   determine the precision of the p-value. Default is \code{100}, such that
#'   the smallest p-value is 0.01
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
#'   to use as max distance in when calculating variograms. Default = 0.25. At
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
#' @param nThreads \code{integer}: Number of threads for parallelization.
#'   Default = 1. Inputting this argument when the \code{BPPARAM} argument is
#'   \code{NULL} would set parallel execution back-end to be
#'   \code{BiocParallel::MulticoreParam(workers = nThreads)}. We recommend
#'   setting this argument to be the number of cores available
#'   (\code{parallel::detectCores(logical = FALSE)}). If \code{BPPARAM} argument
#'   is not \code{NULL}, the \code{BPPARAM} argument would override
#'   \code{nThreads} argument.
#'
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for
#'   parallelization. This argument is provided for advanced users of
#'   \code{BiocParallel} for further flexibility for setting up
#'   parallel-execution back-end. Default is NULL. If provided, this is assumed
#'   to be an instance of \code{BiocParallelParam}.
#'
#' @return The output is returned as a \code{data.frame} containing the columns:
#' \itemize{
#'   \item{\code{correlationCoef}}{Pearson's correlation coefficient.}
#'   \item{\code{pValueNaive}}{the analytical p-value naively assuming independent
#'   observations}
#'   \item{\code{pValuePermuteX}}{the p-value when creating an empirical null from permutations
#'   of observations in X}
#'   \item{\code{pValuePermuteY}}{the p-value when creating an empirical null from
#'   permutations of observations in Y}
#'   \item{\code{deltaStarMedianX}}{the median delta star (the delta which
#'   minimizes the difference between the variogram of the permutation and the
#'   variogram of observations) across permutations of X}
#'   \item{\code{deltaStarMedianY}}{the median delta star across permutations of Y}
#'   \item{\code{deltaStarX}}{list of delta star for all permutations of X}
#'   \item{\code{deltaStarY}}{list of delta star for all permutations of Y}
#'   \item{\code{nullCorrelationsX}}{correlation coefficients for Y and all permuations of X}
#'   \item{\code{nullCorrelationsY}}{correlation coefficients for X and all permuations of Y}
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
#'                           pos = cbind(quakes_data$lat, quakes_data$long))
#' cor
#'
#' #plot the delta star (the delta which minimizes the difference between the
#' variogram of the permutation and the variogram of observations) for all
#' permutations to see if clear peak found in the range inputted
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
#'                           deltaY = seq(0.02, 0.5, 0.02))
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
#' p1 <- ggplot2::ggplot(quakes_data, ggplot2::aes(x = long, y = lat, color = depth)) +
#' ggplot2::geom_point(size = 2) + # Add points for each earthquake
#'   ggplot2::scale_color_gradient(low = "lightblue", high = "blue") + # Color gradient for depth
#'   ggplot2::labs(title = "Locations of Earthquakes off Fiji", x = "Longitude", y = "Latitude", color = "Depth (km)") +
#'   ggplot2::theme_minimal() + # Apply a minimal theme
#'   ggplot2::coord_map() # Use a map projection for better representation of the globe
#'
#' p2 <- ggplot2::ggplot(quakes_data, ggplot2::aes(x = long, y = lat, color = mag)) +
#' ggplot2::geom_point(size = 2) + # Add points for each earthquake
#'   ggplot2::scale_color_gradient(low = "lightblue", high = "blue") + # Color gradient for mag
#'   ggplot2::labs(title = "Locations of Earthquakes off Fiji", x = "Longitude", y = "Latitude", color = "Richter Magnitude") +
#'   ggplot2::theme_minimal() + # Apply a minimal theme
#'   ggplot2::coord_map() # Use a map projection for better representation of the globe
#'
#' p1
#' p2  
spatialCorrelation <- function(X, Y, pos, nPermutations = 100, deltaX = NULL, deltaY = NULL, maxDistPrctile = 0.25, returnPermutations = FALSE, nThreads = 1, BPPARAM = NULL){
  
  ## Set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nThreads)
  }
  
  ## Organize data into dataframes
  dataForward <- data.frame(X = X, 
                            Y = Y,
                            x = pos[,1],
                            y = pos[,2])
  
  dataReverse <- data.frame(X = dataForward$Y, 
                            Y = dataForward$X, 
                            x = dataForward$x,
                            y = dataForward$y)
  
  ## If deltas to test are not supplied, try 0.1 to 0.9 for each dataset 
  if (is.null(deltaX)){
    deltaX = seq(0.1,0.9,0.1)
  }
  
  if (is.null(deltaY)){
    deltaY = seq(0.1,0.9,0.1)
  }
  
  tryCatch({
    ## Calculate Pearson's correlation and return data frame with correlation
    #estimate and naive p-value assuming independence
    corDF <- cor.test(dataForward$X, dataForward$Y)
    
    ## Calculate corrected p-value for Pearson's correlation
    resultsPermuteX <- viladomatCorrelation(dataForward, delta = deltaX, maxDistPrctile = maxDistPrctile, BPPARAM = BPPARAM, nPermutations = nPermutations)
    resultsPermuteY <- viladomatCorrelation(dataReverse, delta = deltaY, maxDistPrctile = maxDistPrctile, BPPARAM = BPPARAM, nPermutations = nPermutations)
    
    ## Store correlation value, naive p-value, and corrected p-value for
    #permuting either source and target as dataframe with 1 row
    if(returnPermutations == TRUE){
      output <- data.frame(correlationCoef = corDF$estimate,
                           pValueNaive = corDF$p.value,
                           pValuePermuteX = resultsPermuteX[["pValueGlobal"]],
                           pValuePermuteY = resultsPermuteY[["pValueGlobal"]],
                           deltaStarMedianX = resultsPermuteX[["deltaStarMedian"]],
                           deltaStarMedianY = resultsPermuteY[["deltaStarMedian"]],
                           deltaStarX = I(list(resultsPermuteX[["deltaStar"]])),
                           deltaStarY = I(list(resultsPermuteY[["deltaStar"]])),
                           nullCorrelationsX = I(list(resultsPermuteX[["nullCorGlobal"]])),
                           nullCorrelationsY = I(list(resultsPermuteY[["nullCorGlobal"]])),
                           permutationsX = I(list(resultsPermuteX[["permutations"]])),
                           permutationsY = I(list(resultsPermuteY[["permutations"]]))
      )
    } else {
      output <- data.frame(correlationCoef = corDF$estimate,
                           pValueNaive = corDF$p.value,
                           pValuePermuteX = resultsPermuteX[["pValueGlobal"]],
                           pValuePermuteY = resultsPermuteY[["pValueGlobal"]],
                           deltaStarMedianX = resultsPermuteX[["deltaStarMedian"]],
                           deltaStarMedianY = resultsPermuteY[["deltaStarMedian"]],
                           deltaStarX = I(list(resultsPermuteX[["deltaStar"]])),
                           deltaStarY = I(list(resultsPermuteY[["deltaStar"]])),
                           nullCorrelationsX = I(list(resultsPermuteX[["nullCorGlobal"]])),
                           nullCorrelationsY = I(list(resultsPermuteY[["nullCorGlobal"]]))
      )
    } 
    
    return(output)
    
  }
  ,
  
  # warning = function(cond) {
  #   print(cond)
  #   #if get warning that cor cannot be calculated because standard deviation is zero and can't divide by zero
  #   #message(paste0("Expression of has 0 variance in X and/or target"))
  #   #return NA for correlation value, naive p-value, and corrected p-values (permuting source or permuting target) as dataframe with 1 row
  #   output <- data.frame(correlationCoef = NA,
  #                          pValueNaive = NA,
  #                          pValuePermuteX = NA,
  #                          pValuePermuteY = NA,
  #                          deltaStarMedianX = NA,
  #                          deltaStarMedianY = NA,
  #                          deltaStarX = NA,
  #                          deltaStarY = NA,
  #                          nullCorrelationsX = NA,
  #                          nullCorrelationsY = NA
  #     )
  # 
  #     return(output)
  #  }
  # 
  # ,
  error = function(cond) {
    print(cond)
    #if get error in the main correction function
    #return NA for corrected p-values (permuting source or permuting target) as dataframe with 1 row
    
    if(returnPermutations == TRUE){
      output <- data.frame(correlationCoef = corDF$estimate,
                           pValueNaive = corDF$p.value,
                           pValuePermuteX = NA,
                           pValuePermuteY = NA,
                           deltaStarMedianX = NA,
                           deltaStarMedianY = NA,
                           deltaStarX = NA,
                           deltaStarY = NA,
                           nullCorrelationsX = NA,
                           nullCorrelationsY = NA,
                           permutationsX = NA,
                           permutationsY = NA)
    } 
    else {
      output <- data.frame(correlationCoef = corDF$estimate,
                           pValueNaive = corDF$p.value,
                           pValuePermuteX = NA,
                           pValuePermuteY = NA,
                           deltaStarMedianX = NA,
                           deltaStarMedianY = NA,
                           deltaStarX = NA,
                           deltaStarY = NA,
                           nullCorrelationsX = NA,
                           nullCorrelationsY = NA)
    }
    
    return(output)
    
  }
  )
  
}

#' spatialCorrelationGeneExp
#'
#' @description Function to calculate Pearson's correlation between assays from
#'   two SpatialExperiment datasets. To replace the analytical p-value which
#'   results in a high false positive rate for autocorrelated spatial patterns,
#'   it calculates empirical p-values from empirical null distributions
#'   generated from permuting the datasets and then smoothing to maintain the
#'   original degree of autocorrelation
#'
#' @param input \code{list} List of two SpatialExperiment objects with matched
#'   spatial locations. The first element corresponds to the first
#'   SpatialExperiment (`X`), and the second to the second SpatialExperiment
#'   (`Y`). The SpatialCoords of the two SpatialExperiment objects should be on
#'   the same coordinate framework and observations at the same coordinate
#'   location in both datasets should have the same row names. If the
#'   SpatialExperiment objects do not have shared locations, use
#'   `SEraster::rasterizeGeneExpression()` to generate SpatialExperiment objects
#'   with shared pixel locations. See \code{assayName} parameter if the
#'   SpatialExperiment objects have more than one assay.
#'
#' @param nPermutations \code{integer} or \code{double}: number of permutations
#'   to generate to build the empirical null distribution. This number will
#'   determine the precision of the p-value. Default is \code{100}, such that
#'   the smallest p-value is 0.01
#'
#' @param deltaX \code{list}: List of single numerics or list of numeric vectors
#'   to use for delta, the parameter controlling the degree of smoothing in
#'   permutations of X. The length of the list should the same as the number of
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
#'   to use as max distance in when calculating variograms. Default = 0.25. At
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
#' @param nThreads \code{integer}: Number of threads for parallelization.
#'   Default = 1. Inputting this argument when the \code{BPPARAM} argument is
#'   \code{NULL} would set parallel execution back-end to be
#'   \code{BiocParallel::MulticoreParam(workers = nThreads)}. We recommend
#'   setting this argument to be the number of cores available
#'   (\code{parallel::detectCores(logical = FALSE)}). If \code{BPPARAM} argument
#'   is not \code{NULL}, the \code{BPPARAM} argument would override
#'   \code{nThreads} argument.
#'
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for
#'   parallelization. This argument is provided for advanced users of
#'   \code{BiocParallel} for further flexibility for setting up
#'   parallel-execution back-end. Default is NULL. If provided, this is assumed
#'   to be an instance of \code{BiocParallelParam}.
#'
#' @param verbose \code{logical}: indicate whether to print row number and name
#'   to show progress as the function iterates through the rows of the
#'   SpatialExperiments to calculate a correlation coefficient and empirical
#'   p-value for each row
#'
#' @return The output is returned as a \code{data.frame}. The rownames are the
#'   rownames of the SpatialExperiments. The names of the columns and their
#'   contents are as follows:
#' \itemize{
#'   \item{\code{correlationCoef}}{Pearson's correlation coefficient.}
#'   \item{\code{pValueNaive}}{the analytical p-value naively assuming independent
#'   observations}
#'   \item{\code{pValuePermuteX}}{the p-value when creating an empirical null from permutations
#'   of observations in X}
#'   \item{\code{pValuePermuteY}}{the p-value when creating an empirical null from
#'   permutations of observations in Y}
#'   \item{\code{deltaStarMedianX}}{the median delta star (the delta which
#'   minimizes the difference between the variogram of the permutation and the
#'   variogram of observations) across permutations of X}
#'   \item{\code{deltaStarMedianY}}{the median delta star across permutations of Y}
#'   \item{\code{deltaStarX}}{list of delta star for all permutations of X}
#'   \item{\code{deltaStarY}}{list of delta star for all permutations of Y}
#'   \item{\code{nullCorrelationsX}}{correlation coefficients for Y and all permuations of X}
#'   \item{\code{nullCorrelationsY}}{correlation coefficients for X and all permuations of Y}
#'   }
#'
#' @export
#'
#' @examples
#'
#' data(quakes)
#'
#' ##### Rasterize to get pixels at matched spatial locations #####
#'rastKidney <- SEraster::rasterizeGeneExpression(speKidney, assay_name = 'counts', resolution = 0.2, fun = "mean",BPPARAM = BiocParallel::MulticoreParam(), square = FALSE)
#'
#' cor <- spatialCorrelation(X = quakes_data$depth,
#'                           Y = quakes_data$mag,
#'                           pos = cbind(quakes_data$lat, quakes_data$long))
#' cor
#'
#' #plot the delta star (the delta which minimizes the difference between the
#' variogram of the permutation and the variogram of observations) for all
#' permutations to see if clear peak found in the range inputted
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
#'                           deltaY = seq(0.02, 0.5, 0.02))
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
#' p1 <- ggplot2::ggplot(quakes_data, ggplot2::aes(x = long, y = lat, color = depth)) +
#' ggplot2::geom_point(size = 2) + # Add points for each earthquake
#'   ggplot2::scale_color_gradient(low = "lightblue", high = "blue") + # Color gradient for depth
#'   ggplot2::labs(title = "Locations of Earthquakes off Fiji", x = "Longitude", y = "Latitude", color = "Depth (km)") +
#'   ggplot2::theme_minimal() + # Apply a minimal theme
#'   ggplot2::coord_map() # Use a map projection for better representation of the globe
#'
#' p2 <- ggplot2::ggplot(quakes_data, ggplot2::aes(x = long, y = lat, color = mag)) +
#' ggplot2::geom_point(size = 2) + # Add points for each earthquake
#'   ggplot2::scale_color_gradient(low = "lightblue", high = "blue") + # Color gradient for mag
#'   ggplot2::labs(title = "Locations of Earthquakes off Fiji", x = "Longitude", y = "Latitude", color = "Richter Magnitude") +
#'   ggplot2::theme_minimal() + # Apply a minimal theme
#'   ggplot2::coord_map() # Use a map projection for better representation of the globe
#'
#' p1
#' p2  
spatialCorrelationGeneExp <- function(input, nPermutations = 100, 
                                      deltaX = NULL, deltaY = NULL,  
                                      maxDistPrctile = 0.25, 
                                      returnPermutations = FALSE, 
                                      assayName = NULL,
                                      nThreads = 1, BPPARAM = NULL,
                                      verbose = TRUE){
  
  ## set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nThreads)
  }
  
  #Determine the positions of shared pixels between two rasterized spatial experiments
  source <- input[[1]]
  target <- input[[2]]
  shared_pixels <- intersect(rownames(SpatialExperiment::spatialCoords(source)),
                             rownames(SpatialExperiment::spatialCoords(target)))
  pos <- SpatialExperiment::spatialCoords(source)[shared_pixels,]
  
  #If lists of deltas to test are not supplied, try 0.1 to 0.9 for both datasets for each gene
  if (is.null(deltaX)){
    deltaX = replicate(length(rownames(source)), list(seq(0.1,0.9,0.1)), simplify = FALSE)
  }
  
  if (is.null(deltaY)){
    deltaY = replicate(length(rownames(source)), list(seq(0.1,0.9,0.1)), simplify = FALSE)
  }
  
  ## if name of assay to use in the SpatialExperiment object is not provided, use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }
  
  #calculate Pearson's correlation of expression between shared pixels in datasets for each gene,
  #naive p-value assuming independence and corrected p-value using empirical null from permutations
  correctedCorrelation <- do.call(rbind, lapply(1:length(rownames(source)), function(i) {
    
    #store name of gene
    g <- rownames(source)[i]
    
    #print number of iteration and name of gene
    if (verbose) {
      message(paste0(i, ': ', g))
    }
    
    #store gene expression matrices from SpatialExperiments for gene "g" 
    X <- SummarizedExperiment::assays(source)[[assayName]][g, shared_pixels]
    Y <- SummarizedExperiment::assays(target)[[assayName]][g, shared_pixels]
    
    #calculate correlation, naive p-value, corrected p-value using empirical null from permutations
    output <- spatialCorrelation(X, Y, pos, nPermutations = nPermutations, deltaX = deltaX[[i]], deltaY = deltaY[[i]],  maxDistPrctile = maxDistPrctile, returnPermutations = returnPermutations, nThreads = nThreads, BPPARAM = NULL)
    
    #name row of dataframe with gene name
    row.names(output) <- g
    
    return(output)
    
  }))
  
  return(correctedCorrelation)
}

spatialCorrelationGeneExpWithinSample <- function(rastGexp, nPermutations = 100, delta = NULL, returnPermutations = FALSE, assayName = NULL, nThreads = 1, BPPARAM = NULL, verbose = TRUE){
  
  ## set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nThreads)
  }
  
  #Store positions of pixels in the rasterized spatial experiment
  pos <- SpatialExperiment::spatialCoords(rastGexp)
  
  #If lists of deltas to test are not supplied, try 0.1 to 0.9 for each gene
  if (is.null(delta)){
    delta = replicate(length(rownames(rastGexp)), list(seq(0.1,0.9,0.1)))
  }
  names(delta) <- rownames(rastGexp)
  
  ## if name of assay to use in the SpatialExperiment object is not provided, use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }
  
  #identify all unique combinations of pairs of genes
  genePairs <- combn(rownames(rastGexp), 2)
  #store number of pairs of genes
  n <- dim(genePairs)[2]
  
  correctedCorrelation <- do.call(rbind, lapply(1:n, function(i) {
    
    #store names of genes in the ith pair
    g <- genePairs[1,i]
    g2 <- genePairs[2,i]
      
    #print number of iteration and names of genes in the pair
    if (verbose) {
      message(paste0(i, ': ', g, ' and ', g2))
    }
    
    #Assuming gene expression is set as the first assay element in the 
    #SpatialExperiments, store gene expression matrices for genes in the pair
    X <- assay(rastGexp)[[assayName]][g, ]
    Y <- assay(rastGexp)[[assayName]][g2, ]

    #calculate correlation, naive p-value, corrected p-value using empirical null from permutations
    output <- spatialCorrelation(X, Y, pos, nPermutations = nPermutations, deltaX = delta[[g]], deltaY = delta[[g2]], returnPermutations = returnPermutations, nThreads = nThreads, BPPARAM = BPPARAM)
      
    #add columns with gene names from the pair
    output$first <- g
    output$second <- g2
      
      
    return(output)
    }))
  
  return(correctedCorrelation)
}
  

# @examples

set.seed(0)
t <- Sys.time()
selectGenes <- sig.genes.naive.tp[1:5]
selectGenes <- good.genes[1:5]
input <- list( AKI_irl = rast_norm$AKI_irl[selectGenes ,],
               AKI_ctrl = rast_norm$AKI_ctrl[selectGenes ,])

test_output <-spatialCorrelationGeneExp(input, nThreads = 5
                                        #, delta = rep(0.1, 5)
)
tf <- Sys.time() - t
print(tf) #Time difference of 33.35451 secs


set.seed(0)
t <- Sys.time()
input <- list('AKI_ctrl'= AKI_ctrl_SE[, ],
              'AKI_irl'= AKI_irl_SE[, ])
rast_bad <- rasterizeGeneExpression(input,
                                      assay_name = 'libnorm',
                                      resolution = 3, fun = "mean",
                                      BPPARAM = BiocParallel::MulticoreParam(),
                                      square = FALSE)
input <- list( AKI_irl = rast_bad$AKI_irl[1 ,],
               AKI_ctrl = rast_bad$AKI_ctrl[1 ,])

test_output <-spatialCorrelationGeneExp(input, nThreads = 5
                                        #, delta = rep(0.1, 5)
                                        )
tf <- Sys.time() - t
print(tf) 

t <- Sys.time()
test_output2 <-spatialCorrelationGeneExpWithinSample(rastGexp = rast_norm$AKI_irl[sig.genes.naive.tp[1:5] ,], nThreads = 5)
tf <- Sys.time() - t
print(tf)

load("~/ST_compare/packageFunctions/ischemiaKidney.Rdata")
set.seed(0)
t <- Sys.time()
sc_geneExp <- spatialCorrelationGeneExp(ischemiaKidRastGexp, nThreads = 5)
save(sc_geneExp, file = "~/ST_compare/data/kidney_data/sc_geneExp.Rdata")
tf <- Sys.time() - t
print(tf) #Time difference of 23.98663 hours
