
library(SpatialExperiment)



#' Extract Gene Expression Values for Shared Pixels Between Two SpatialExperiments
#'
#' This function creates a data frame containing gene expression values for a
#' specified gene across shared spatial pixels in two SpatialExperiment objects.
#'
#' @param x A SpatialExperiment object representing the first dataset.
#' @param y A SpatialExperiment object representing the second dataset.
#' @param gene A character string specifying the gene of interest.
#'
#' @return A data frame with three columns:
#' \item{pixel}{Character vector of shared pixel IDs between both experiments.}
#' \item{y}{Numeric vector of gene expression values from `y` (second SpatialExperiment).}
#' \item{x}{Numeric vector of gene expression values from `x` (first SpatialExperiment).}
#'
#' @export
#'
#' @examples
#' getGenePixelDF(spatialExp1, spatialExp2, "Aqp4")
getGenePixelDF <- function (x, y, gene) {

  # pixels that are in both SpatialExperiments
  pixels <- intersect(colnames(x), colnames(y))

  # data frame of pixel values for only given gene
  geneDF <- data.frame(
    pixel = pixels,
    # y = as.vector(assay(y)[gene, pixels]),
    # x = as.vector(assay(x)[gene, pixels])

    # x <- as.vector(assay(x, "pixelval")[gene, pixels, drop = TRUE]),
    # y <- as.vector(assay(y, "pixelval")[gene, pixels, drop = TRUE])

    x = as.vector(as.matrix(assay(x, "pixelval"))[gene, pixels, drop = TRUE]),
    y = as.vector(as.matrix(assay(y, "pixelval"))[gene, pixels, drop = TRUE])

  )

  return(geneDF)

}

#' Threshold and preprocess gene expression data
#'
#' This function processes a gene-to-pixel data frame by removing missing values,
#' applying a threshold to filter out low-expression pixels, and approximating zero values.
#'
#' @param genePixelDF A data frame containing gene expression for one gene
#'        for each pixel.
#'        It should have columns "x" and "y" representing expression levels in two SpatialExperiments.
#' @param t1 Threshold value for gene expression in the second SpatialExperiment (y).
#' @param t2 Threshold value for gene expression in the first SpatialExperiment (x).
#'
#' @return A filtered and processed data frame where:
#'         - NA values are removed,
#'         - Pixels with both x and y values below their respective thresholds are excluded,
#'         - Zero values are approximated to 0.001 to avoid numerical instability.
#'
#' @export
#'
#' @examples
threshold <- function (genePixelDF, t1, t2) {


  # remove NA value
  df <- na.omit(genePixelDF)

  # threshold
  df <- df[!(df$y < t1 & df$x <t2), ]

  # list of pixels inside of
  pixInThresh <- df$pixel

  # list of pixels outside of threshold
  # pixels that are in genePixelDF, but not in df
  pixOutThresh <- setdiff(genePixelDF$pixel, df$pixel)

  # approximate 0s to 0.0001
  df$x[df$x == 0] <- 0.0001
  df$y[df$y == 0] <- 0.0001

  output <- data.frame(
    pixel = I(list(df$pixel)),
    x = I(list(df$x)),
    y = I(list(df$y)),
    numPixelOutThresh = length(pixOutThresh),
    pixelIDOutThresh = I(list(pixOutThresh))
  )

  return(output)

}

#' Computes spatial similarity of gene expression between two SpatialExperiments.
#'
#' This function calculates the spatial similarity of gene expression patterns
#' between two SpatialExperiment objects based on thresholding and fold-change criteria.
#' Similarity is defined as the proportion of pixels where the fold change in gene
#' expression falls within a specified range.
#'
#' @param input List of two SpatialExperiment objects. The first element corresponds
#' to the first spatial experiment (`y`), and the second to the second spatial experiment (`x`).
#' Similarity is computed as the proportion of pixels satisfying:
#' `-1 <= log2(y/x) <= 1` for each gene.
#'
#' @param t1 Numeric. Gene expression threshold for the first spatial experiment (`y`).
#' @param t2 Numeric. Gene expression threshold for the second spatial experiment (`x`).
#' @param foldChange Numeric. Fold-change threshold defining similarity in gene
#' expression across spatial locations.
#'
#' @return A list containing:
#' \item{similarityTable}{A data frame with the following columns:}
#'   \itemize{
#'     \item{\strong{gene}}{Gene name.}
#'     \item{\strong{percentSimilarity}}{Percentage of similar pixels for the gene.}
#'     \item{\strong{similarPixelID}}{List of pixel IDs classified as similar.}
#'     \item{\strong{numPixelInThresh}}{Number of pixels above the threshold in both experiments.}
#'     \item{\strong{pixelIDInThresh}}{List of pixel IDs above the threshold.}
#'     \item{\strong{numPixelOutThresh}}{Number of pixels below the threshold.}
#'     \item{\strong{pixelIDOutThresh}}{List of pixel IDs below the threshold.}
#'   }
#' \item{pixelLogTransformation}{A data frame containing log-transformed fold-change values for each gene.}
#' \item{parameters}{A list of input parameters used in the computation.}
#'
#' @export
#'
#' @examples

spatialSimilarity <- function (input, t1 = 1, t2 = 1, foldChange = 1) {

  output <- data.frame(
    gene = character(),
    percentSimilarity = numeric(),
    similarPixelID = I(list()),
    numPixelInThresh = numeric(),
    pixelIDInThresh = I(list()),
    numPixelOutThresh = numeric(),
    pixelIDOutThresh = I(list())
  )

  logTransGenes <- data.frame(
    gene = character(),
    log = I(list())
  )

  y = input[[1]]
  x = input[[2]]

  shared_genes <- intersect(rownames(x), rownames(y))

  for (gene in shared_genes) {

    # gene to pixel matrix
    genePixel <- getGenePixelDF(x = x, y = y, gene = gene);

    # threshold the gene to pixel matrix and
    # and remove NA values and make approximations for 0 numbers
    thresh <- threshold(genePixel, t1, t2)
    threshDF <- data.frame(
      pixel = thresh$pixel[[1]],
      x = thresh$x[[1]],
      y = thresh$y[[1]]
    )

    # skip gene if less than 10% of all gene pass the threshold
    if(dim(threshDF)[1] < (0.1 * dim(genePixel)[1])) {

      output <- rbind(output, data.frame(
        gene = gene,
        percentSimilarity = NA,
        similarPixelID = NA,
        numPixelInThresh = dim(thresh)[1],
        pixelIDInThresh = NA,
        numPixelOutThresh = thresh$numPixelOutThresh[[1]],
        pixelIDOutThresh = thresh$pixelIDOutThresh
      ))
      next
    }

    # calculate the slope for each pixel
    threshDF$slope <- threshDF$y / threshDF$x

    # log transformation
    logTrans <- data.frame(
      pixel = threshDF$pixel,
      log = log2(threshDF$slope)
    )

    # calculate similarity for gene based on fold change
    similarPixels <- logTrans[logTrans$log >= -foldChange & logTrans$log <= foldChange, ]
    geneSimilarity <- dim(similarPixels)[1] / dim(logTrans)[1]

    output <- rbind(output, data.frame(
      gene = gene,
      percentSimilarity = geneSimilarity,
      similarPixelID = I(list(similarPixels$pixel)),
      numPixelInThresh = dim(logTrans)[1],
      pixelIDInThresh = thresh$pixel,
      numPixelOutThresh = thresh$numPixelOutThresh[[1]],
      pixelIDOutThresh = thresh$pixelIDOutThresh
    ))

    logTransGenes <- rbind(logTransGenes, data.frame(
      gene = gene,
      log = I(list(logTrans$log))
    ))

  }

  return(list(
    similarityTable = output,
    pixelLogTransformation = logTransGenes,
    parameters = list(
      t1 = t1,
      t2 = t2,
      foldChange = foldChange,
      input = input
      )
    ))

}

