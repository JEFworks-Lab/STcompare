# Helper to run `spatialCorrelation()` for a subset of genes at a fixed
# permutation count while keeping the outer loop focused on the iterative
# rerun logic.
.run_spatial_correlation_iteration <- function(
    genes,
    source,
    target,
    shared_pixels,
    pos,
    assayName,
    nPermutations,
    deltaX,
    deltaY,
    maxDistPrctile,
    returnPermutations,
    nThreads,
    BPPARAM,
    verbose,
    seed) {

  if (length(genes) == 0) {
    return(NULL)
  }

  res_list <- lapply(seq_along(genes), function(j) {

    # store the gene name
    gene <- genes[j]
    gene_idx <- match(gene, rownames(source))

    if (verbose) {
      message(sprintf(
        "nPermutations=%s | gene %d/%d: %s",
        nPermutations, j, length(genes), gene
      ))
    }

    X <- SummarizedExperiment::assays(source)[[assayName]][gene, shared_pixels]
    Y <- SummarizedExperiment::assays(target)[[assayName]][gene, shared_pixels]

    out <- spatialCorrelation(
      X, Y, pos,
      nPermutations = nPermutations,
      deltaX = deltaX[[gene_idx]],
      deltaY = deltaY[[gene_idx]],
      maxDistPrctile = maxDistPrctile,
      returnPermutations = returnPermutations,
      nThreads = nThreads,
      BPPARAM = BPPARAM,
      seed = seed
    )

    rownames(out) <- gene
    out
  })

  do.call(rbind, res_list)
}

# Helper to decide which genes should be rerun with a larger number of
# permutations. Only genes with both empirical permutation p-values below the
# current screening threshold are carried forward.
.get_genes_to_repermute <- function(results_df, alpha, nPermutes) {
  t <- alpha / nPermutes
  keep <- results_df$pValuePermuteX < t & results_df$pValuePermuteY < t
  rownames(results_df)[keep]
}


#' spatialCorrelationGeneExpIterPermutations
#'
#' @description Function to calculate Pearson's correlation between assays from
#'   two SpatialExperiment datasets using an iterative permutation strategy. It
#'   first evaluates all genes with a smaller number of permutations, then only
#'   reruns genes that remain potentially significant with progressively larger
#'   permutation counts to refine their empirical p-values while accounting for
#'   original degree of autocorrelation.
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
#' @param alpha \code{numeric}: significance threshold used to decide which genes
#'   should be rerun at the next permutation level. After each iteration, a gene
#'   is carried forward only when both \code{pValuePermuteX} and
#'   \code{pValuePermuteY} are less than \code{alpha / nPermutations[k]}.
#'   Default is \code{0.05}.
#'
#' @param nPermutations \code{numeric vector}: numbers of permutations to use
#'   across iterative rounds. The vector is applied from smallest to largest.
#'   All genes are first tested with \code{nPermutations[1]}; genes passing the
#'   screening rule defined by \code{alpha} are rerun with
#'   \code{nPermutations[2]}, and so on. Default is \code{c(100, 1000)}.
#'
#' @param deltaX \code{list}: List of single numerics or list of numeric vectors
#'   to use for delta, the parameter controlling the degree of smoothing in
#'   permutations of X. The length of the list should be the same as the number
#'   of rows in the SpatialExperiment. Delta is a proportion calculated by
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
#'   returned as output will have columns with the values of the permutations
#'   used to calculate the null correlations and empirical p-values. Default is
#'   \code{FALSE}.
#'
#' @param assayName \code{character} or \code{integer} A character string or
#'   numeric specifying the assay in the SpatialExperiment to use. Default is
#'   \code{NULL}. If no value is supplied for \code{assayName}, then the first
#'   assay is used as a default.
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
#' @param verbose \code{logical}: indicate whether to print the current
#'   permutation level together with the row number and name to show progress as
#'   the function iterates through genes.
#'
#' @param seed \code{integer}: Seed for the random number generator used to
#'   generate noise in the variogram matching step. Ensures reproducibility of
#'   empirical p-values regardless of parallelization back-end. Default is
#'   \code{0}.
#'
#' @param adjustMethod \code{character}: multiple-testing correction method
#'   passed to \code{stats::p.adjust()} for the final \code{pValuePermuteX} and
#'   \code{pValuePermuteY} columns separately. Must be one of
#'   \code{p.adjust.methods}. Default is \code{"BH"}.
#'
#' @return The output is returned as a \code{data.frame}. The rownames are the
#'   rownames of the SpatialExperiments, and each row reflects the last
#'   permutation round in which that gene was evaluated. The columns and their
#'   contents are as follows:
#' \itemize{
#'   \item{\code{correlationCoef}}{Pearson's correlation coefficient.}
#'   \item{\code{pValueNaive}}{the analytical p-value naively assuming independent
#'   observations}
#'   \item{\code{pValuePermuteX}}{multiple-testing-adjusted p-value from an
#'   empirical null generated by permutations of observations in X}
#'   \item{\code{pValuePermuteY}}{multiple-testing-adjusted p-value from an
#'   empirical null generated by permutations of observations in Y}
#'   \item{\code{deltaStarMedianX}}{the median delta star (the delta which
#'   minimizes the difference between the variogram of the permutation and the
#'   variogram of observations) across permutations of X}
#'   \item{\code{deltaStarMedianY}}{the median delta star across permutations of Y}
#'   \item{\code{deltaStarX}}{list of delta star for all permutations of X}
#'   \item{\code{deltaStarY}}{list of delta star for all permutations of Y}
#'   \item{\code{nullCorrelationsX}}{correlation coefficients for Y and all
#'   permutations of X}
#'   \item{\code{nullCorrelationsY}}{correlation coefficients for X and all
#'   permutations of Y}
#'   \item{\code{permutationsX}}{(optional) a N x B matrix, where N is the
#'   length of X and B is the final \code{nPermutations} used for that gene.
#'   Each column is the resulting values of a permutation of X}
#'   \item{\code{permutationsY}}{(optional) a N x B matrix, where N is the
#'   length of Y and B is the final \code{nPermutations} used for that gene.
#'   Each column is the resulting values of a permutation of Y}
#'   }
#'
#' @export
#'
#' @examples
#' data(speKidney)
#' \dontrun{
#' rastKidney <- SEraster::rasterizeGeneExpression(
#'   speKidney,
#'   assay_name = "counts",
#'   resolution = 0.2,
#'   fun = "mean",
#'   BPPARAM = BiocParallel::MulticoreParam(),
#'   square = FALSE
#' )
#'
#' rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
#'
#' corr <- spatialCorrelationGeneExpIterPermutations(
#'   rastGexpListAB,
#'   nPermutations = c(100, 1000),
#'   nThreads = 5
#' )
#'
#' negCorrelation
#' }
spatialCorrelationGeneExpIterPermutations <- function(
    input,
    alpha = 0.05,
    nPermutations = c(1e2, 1e3),
    deltaX = NULL,
    deltaY = NULL,
    maxDistPrctile = 0.25,
    returnPermutations = FALSE,
    assayName = NULL,
    nThreads = 1,
    BPPARAM = NULL,
    verbose = TRUE,
    seed = 0,
    adjustMethod = "BH") {

  # correction method should be from p.adjust.methods
  if (!adjustMethod %in% p.adjust.methods) {
    stop("adjustMethod must be one of: ", paste(p.adjust.methods, collapse = ", "))
  }

  if (!is.numeric(nPermutations) || length(nPermutations) < 1) {
    stop("nPermutations must be a numeric vector of length >= 1.")
  }

  if (any(is.na(nPermutations)) || any(nPermutations <= 0)) {
    stop("All nPermutations values must be positive.")
  }

  if (is.unsorted(nPermutations, strictly = FALSE)) {
    nPermutations <- sort(nPermutations)
  }

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nThreads)
  }

  # Determine the positions of shared pixels between two rasterized spatial
  # experiments
  source <- input[[1]]
  target <- input[[2]]
  shared_pixels <- intersect(rownames(SpatialExperiment::spatialCoords(source)),
                             rownames(SpatialExperiment::spatialCoords(target)))
  pos <- SpatialExperiment::spatialCoords(source)[shared_pixels,]

  # If lists of deltas to test are not supplied, try 0.1 to 0.9 for both datasets
  # for each gene
  if (is.null(deltaX)){
    deltaX <- rep(list(seq(0.1,0.9,0.1)), length(rownames(source)))
  }

  if (is.null(deltaY)){
    deltaY <- rep(list(seq(0.1,0.9,0.1)), length(rownames(source)))
  }

  # if name of assay to use in the SpatialExperiment object is not provided,
  # use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }

  genes_all <- rownames(source)

  # creates an index for each gene name
  # storage for latest result per gene
  final_results <- NULL

  # genes to test in current iteration
  genes_current <- genes_all

  # go through all nPermutations
  for (k in seq_along(nPermutations)) {

    # skipping if no more genes are significant
    if (length(genes_current) == 0) {
      if (verbose) {
        message(sprintf(
          "Iteration %d/%d (nPermutations=%d): no genes left to test — skipping",
          k, length(nPermutations), nPermutations[k]
        ))
      }
      next
    }

    iter_res <- .run_spatial_correlation_iteration(
      genes = genes_current,
      source = source,
      target = target,
      shared_pixels = shared_pixels,
      pos = pos,
      assayName = assayName,
      nPermutations = nPermutations[k],
      deltaX = deltaX,
      deltaY = deltaY,
      maxDistPrctile = maxDistPrctile,
      returnPermutations = returnPermutations,
      nThreads = nThreads,
      BPPARAM = BPPARAM,
      verbose = verbose,
      seed = seed
    )

    # overwrite previous results for genes rerun at this iteration
    if (is.null(final_results)) {
      final_results <- iter_res
    } else {
      final_results[rownames(iter_res), colnames(iter_res)] <- iter_res
    }

    # if there are more iterations left
    # decide which genes needs more permutations
    if (k < length(nPermutations)) {
      genes_next <- .get_genes_to_repermute(iter_res, alpha = alpha, nPermutes = nPermutations[k])
      genes_current <- genes_next
    } 
  }

  # mht correct for pValuePermuteX and pValuePermuteY seperately
  final_results$pValuePermuteX <- stats::p.adjust(final_results$pValuePermuteX, method = adjustMethod)
  final_results$pValuePermuteY <- stats::p.adjust(final_results$pValuePermuteY, method = adjustMethod)

  # order rows back to original gene order
  final_results <- final_results[genes_all, , drop = FALSE]
  final_results
}










