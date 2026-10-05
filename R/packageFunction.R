#' Extract the values of one gene on the pixels shared by two SpatialExperiments
#'
#' @param x A SpatialExperiment object representing the first dataset.
#' @param y A SpatialExperiment object representing the second dataset.
#' @param gene A character string specifying the gene of interest.
#' @param assayName A character string or numeric specifying the assay to use. Default is \code{1}, the first
#'   assay.
#'
#' @return A data frame with one row per shared pixel (the column names in both objects, in the order of
#'   \code{x}) and three columns: \code{pixel}, the pixel name; \code{x}, the gene's value in \code{x}; and
#'   \code{y}, its value in \code{y}. Only the gene's row is extracted from each assay, so a sparse assay is
#'   not densified.
#'
#' @noRd
getGenePixelDF <- function(x, y, gene, assayName = 1) {
  pixels <- intersect(colnames(x), colnames(y))
  data.frame(
    pixel = pixels,
    x = as.vector(as.matrix(SummarizedExperiment::assay(x, assayName)[gene, pixels, drop = FALSE])),
    y = as.vector(as.matrix(SummarizedExperiment::assay(y, assayName)[gene, pixels, drop = FALSE]))
  )
}

#' Computes spatial similarity of gene expression between two
#' SpatialExperiments.
#'
#' This function calculates the spatial similarity of gene expression patterns
#' between two SpatialExperiment objects based on thresholding and fold-change
#' criteria. Similarity is defined as the proportion of pixels where the fold
#' change in gene expression falls within a specified range.
#'
#' @details The genes compared are those in both objects (row names), and the
#'   pixels compared are those whose names (column names) are in both objects:
#'   rasterize the two samples in one call, such as
#'   \code{SEraster::rasterizeGeneExpression(list(x, y), ...)}, so that a pixel
#'   name means the same location in both.
#'
#'   For each gene, a pixel is kept when its value in the first object is
#'   above \code{t1} or its value in the second object is above \code{t2}
#'   (strictly), so the pixels removed are those at or below the threshold in
#'   both objects. Zeros among the kept pixels are replaced by 1e-4. With
#'   \eqn{l = \log_2(y / x)}{l = log2(y / x)} for a kept pixel, where \eqn{x}
#'   and \eqn{y} are its values in the first and the second object, the pixel
#'   is similar when \eqn{-b \le l \le b}{-b <= l <= b} for
#'   \eqn{b} = \code{foldChange} (the default 1 means within two-fold, both
#'   ends included), higher in the first object when \eqn{l < -b}, and higher
#'   in the second object when \eqn{l > b}. If fewer than a share
#'   \code{minPixels} of the shared pixels are kept, the gene gets no
#'   similarity (\code{NA}).
#'
#'   The values must be finite and not negative, such as counts or normalized
#'   expression on a linear scale: a fold change is not defined for negative
#'   values, and genes with negative, missing or infinite values on the shared
#'   pixels stop the function with an error that names them. All genes are
#'   computed at once by the helper that also computes the similarity columns
#'   of \code{\link{compareSpatial}()}, which gives the same numbers.
#'
#' @param input List of two SpatialExperiment objects (only the first two
#'   elements are used). The first element corresponds to the first spatial
#'   experiment (\code{x}), and the second to the second spatial experiment
#'   (\code{y}). Similarity is computed as the proportion of pixels satisfying
#'   \code{-b <= log2(y/x) <= b} for each gene, where \code{b = foldChange}
#'   (default 1).
#'
#' @param t1 \code{numeric}: Gene expression threshold for the first spatial
#'   experiment (\code{x}), a single number used for every gene. Only pixels
#'   with values greater than \code{t1} in \code{x} or greater than \code{t2}
#'   in \code{y} are used to calculate the similarity score. Default is
#'   \code{NULL}: each gene's threshold is then its \code{minQuantile}
#'   quantile in \code{x}.
#' @param t2 \code{numeric}: Gene expression threshold for the second spatial
#'   experiment (\code{y}), like \code{t1}. Default is \code{NULL}: each gene's
#'   threshold is then its \code{minQuantile} quantile in \code{y}.
#' @param minQuantile \code{numeric}: The quantile (between 0 and 1) of each
#'   gene's values in each experiment used as its threshold when \code{t1} or
#'   \code{t2} is not supplied. Default is \code{0.05}.
#' @param minPixels \code{numeric}: The smallest proportion (between 0 and 1)
#'   of the shared pixels that must pass the thresholds. If fewer pixels pass,
#'   a spatial similarity score is not calculated for the gene (\code{NA}).
#'   Default is \code{0.1}.
#' @param foldChange \code{numeric}: The similarity band on the log2 scale:
#'   pixels with \eqn{|\log_2(y / x)| \le}{|log2(y / x)| <=} \code{foldChange}
#'   are similar. Default is \code{1}, that is, within two-fold.
#' @param assayName A character string or numeric specifying the assay in the
#'   Spatial Experiment to use. Default is \code{NULL}. If no value is supplied
#'   for \code{assayName}, then the first assay is used as a default. The assay
#'   used is recorded in the result, and \code{\link{linearRegression}()},
#'   \code{\link{pixelClass}()} and \code{\link{savePlots}()} use it by
#'   default.
#' @param verbose \code{logical}: if \code{TRUE}, print a message with the
#'   number of genes and shared pixels compared. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{similarityTable}}{A data frame with one row per gene and the
#'   following columns:
#'   \describe{
#'     \item{\code{gene}}{Gene name.}
#'     \item{\code{percentSimilarity}}{Proportion (between 0 and 1) of the
#'     kept pixels that are similar: the similarity score S.}
#'     \item{\code{percentDissimilarityX}}{Proportion of the kept pixels for
#'     which \code{log2(y/x) < -foldChange} (higher in the first experiment).}
#'     \item{\code{percentDissimilarityY}}{Proportion of the kept pixels for
#'     which \code{log2(y/x) > foldChange} (higher in the second experiment).}
#'     \item{\code{similarPixelID}}{List of pixel IDs classified as similar.}
#'     \item{\code{dissimilarPixelIDX}}{List of pixel IDs for
#'     \code{log2(y/x) < -foldChange}.}
#'     \item{\code{dissimilarPixelIDY}}{List of pixel IDs for
#'     \code{log2(y/x) > foldChange}.}
#'     \item{\code{numPixelInThresh}}{Number of pixels kept: above the
#'     threshold in either experiment.}
#'     \item{\code{pixelIDInThresh}}{List of the kept pixel IDs.}
#'     \item{\code{numPixelOutThresh}}{Number of pixels removed: at or below
#'     the threshold in both experiments.}
#'     \item{\code{pixelIDOutThresh}}{List of the removed pixel IDs.}
#'     \item{\code{t1}}{Threshold value for this gene in first spatial experiment.}
#'     \item{\code{t2}}{Threshold value for this gene in second spatial experiment.}
#'   }
#'   For a gene with fewer than \code{minPixels} kept pixels, the three
#'   proportions and the three lists of classified pixels are \code{NA}.}
#'   \item{\code{pixelLogTransformation}}{A data frame with one row per gene
#'   that has a similarity score: \code{gene}, and \code{log}, a list with the
#'   \code{log2(y/x)} values of its kept pixels (in the order of
#'   \code{pixelIDInThresh}).}
#'   \item{\code{parameters}}{A list with \code{foldChange}, \code{minPixels},
#'   the \code{input} objects, and the \code{assayName} (the assay used),
#'   \code{minQuantile}, \code{t1} and \code{t2} arguments.}
#' }
#'
#' @seealso \code{\link{compareSpatial}()}, which computes the similarity
#'   together with the spatial correlation test; \code{\link{linearRegression}()},
#'   \code{\link{pixelClass}()} and \code{\link{savePlots}()} to plot the
#'   result.
#'
#' @export
#'
#' @examples
#' data(speKidney)
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                square = FALSE)
#'
#' s <- spatialSimilarity(list(rastKidney$A, rastKidney$C))
#' s$similarityTable[, c("gene", "percentSimilarity", "percentDissimilarityX",
#'                       "percentDissimilarityY", "numPixelInThresh")]
#'
spatialSimilarity <- function(
    input,
    t1 = NULL,
    t2 = NULL,
    minQuantile = 0.05,
    minPixels = 0.1,
    foldChange = 1,
    assayName = NULL,
    verbose = FALSE
    ) {

  if (is.null(assayName)) {
    assayName <- 1
  }
  if (length(input) < 2L || !inherits(input[[1]], "SummarizedExperiment") ||
      !inherits(input[[2]], "SummarizedExperiment")) {
    stop("input must be a list of two SpatialExperiment objects")
  }
  for (a in c("t1", "t2")) {
    v <- get(a)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v))) {
      stop(sprintf("%s must be NULL or a single number", a))
    }
  }
  .stc_check_unit_interval(minQuantile, "minQuantile")
  .stc_check_unit_interval(minPixels, "minPixels")
  if (!is.numeric(foldChange) || length(foldChange) != 1L || is.na(foldChange) || foldChange < 0) {
    stop("foldChange must be a single non-negative number")
  }

  x <- input[[1]]
  y <- input[[2]]
  genes <- intersect(rownames(x), rownames(y))
  pixels <- intersect(colnames(x), colnames(y))
  if (!length(pixels)) {
    stop("the two SpatialExperiment objects share no pixel (column) names. Rasterize both samples in one call, ",
         "SEraster::rasterizeGeneExpression(list(x, y), ...), so that they have the same pixels")
  }
  if (verbose) {
    message(sprintf("spatialSimilarity: %d gene(s) on %d shared pixels", length(genes), length(pixels)))
  }

  # pixels x genes matrices of the shared genes on the shared pixels, each assay subset (and densified) once
  values <- function(s, label) {
    M <- t(as.matrix(SummarizedExperiment::assay(s, assayName)[genes, pixels, drop = FALSE]))
    storage.mode(M) <- "double"
    bad <- function(v) genes[v]
    nonfinite <- bad(colSums(!is.finite(M)) > 0L)
    if (length(nonfinite)) {
      stop(sprintf("spatialSimilarity() needs finite values: %d gene(s) have missing or infinite values in %s on the shared pixels, for example %s",
                   length(nonfinite), label, paste(utils::head(nonfinite, 5L), collapse = ", ")))
    }
    negative <- bad(colSums(M < 0) > 0L)
    if (length(negative)) {
      stop(sprintf(paste0("spatialSimilarity() compares fold changes, which need values that are not negative ",
                          "(counts or normalized expression on a linear scale): %d gene(s) have negative values in %s, for example %s"),
                   length(negative), label, paste(utils::head(negative, 5L), collapse = ", ")))
    }
    M
  }
  X <- values(x, "the first object")
  Y <- values(y, "the second object")

  # the similarity of every gene, with the kept pixels and their log2(y / x)
  s <- .stc_similarity(X, Y, minQuantile = minQuantile, minPixels = minPixels, foldChange = foldChange,
                       t1 = t1, t2 = t2, details = TRUE)
  d <- attr(s, "details")
  G <- length(genes)
  N <- length(pixels)
  # the kept pixels of each class, NA for genes without a similarity score
  classify <- function(in_class) {
    I(lapply(seq_len(G), function(j) {
      if (d$scored[j]) pixels[d$kept[[j]][in_class(d$log[[j]])]] else NA
    }))
  }
  output <- data.frame(gene = genes, percentSimilarity = s$similarity, percentDissimilarityX = s$dissimilarityX,
                       percentDissimilarityY = s$dissimilarityY, stringsAsFactors = FALSE)
  output$similarPixelID <- classify(function(l) l >= -foldChange & l <= foldChange)
  output$dissimilarPixelIDX <- classify(function(l) l < -foldChange)
  output$dissimilarPixelIDY <- classify(function(l) l > foldChange)
  output$numPixelInThresh <- s$nPixelsSimilarity
  output$pixelIDInThresh <- I(lapply(d$kept, function(k) pixels[k]))
  output$numPixelOutThresh <- N - s$nPixelsSimilarity
  output$pixelIDOutThresh <- I(lapply(d$kept, function(k) pixels[!(seq_len(N) %in% k)]))
  output$t1 <- s$thresholdX
  output$t2 <- s$thresholdY

  logTransGenes <- data.frame(gene = genes[d$scored], stringsAsFactors = FALSE)
  logTransGenes$log <- I(d$log[d$scored])

  return(list(
    similarityTable = output,
    pixelLogTransformation = logTransGenes,
    parameters = list(
      foldChange = foldChange,
      minPixels = minPixels,
      input = input,
      assayName = assayName,
      minQuantile = minQuantile,
      t1 = t1,
      t2 = t2
      )
    ))

}
