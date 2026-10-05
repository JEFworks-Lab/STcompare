
# 1. linear regression plot ###########################################

#' Plot a gene's values in two samples, coloured by similarity class
#'
#' This function creates a scatter plot comparing gene expression levels between
#' two spatial experiments for a specified gene. It colors data points based on
#' similarity classification and overlays fold-change threshold lines. No
#' regression line is fitted: the solid line is \eqn{y = x}. Both axes end at
#' the larger of the two 95\% quantiles of the gene's values, so the pixels
#' with the highest values are not drawn (with a warning from ggplot2).
#'
#' @param input A list. Results from \code{spatialSimilarity()}. This includes
#'   the similarity table, log-transformed pixel data, and analysis parameters.
#' @param gene Character. The name of the gene to visualize.
#' @param assayName A character string or numeric specifying the assay in the
#'   Spatial Experiment to use. Default is \code{NULL}: the assay that
#'   \code{spatialSimilarity()} used (the first assay for results of
#'   STcompare 0.1.0 and earlier, which did not record it).
#'
#' @return A ggplot2 scatter plot displaying gene expression values from two
#'   spatial experiments. Data points are colored as follows:
#' \describe{
#'   \item{\strong{blue}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{yellow}}{Pixels with greater expression in dataset X than Y.}
#'   \item{\strong{red}}{Pixels with greater expression in dataset Y than X.}
#'   \item{\strong{grey}}{Pixels with gene expression below the threshold in both experiments.}
#' }
#' The plot includes:
#' \describe{
#'   \item{\strong{Solid line:}}{y = x (equal expression).}
#'   \item{\strong{Dashed lines:}}{Fold-change similarity thresholds (upper and lower bounds).}
#' }
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
#' s <- spatialSimilarity(list(A = rastKidney$A, C = rastKidney$C))
#' linearRegression(s, "Gene")
#'
linearRegression <- function (input, gene, assayName = NULL) {

  ## the assay that spatialSimilarity() used, unless another one is given
  assayName <- .stc_similarity_assay(input, assayName)

  df <- getGenePixelDF(
    y = input$parameters$input[[2]],
    x = input$parameters$input[[1]],
    gene = gene,
    assayName = assayName
  )

  # similarity score from the similarity table
  s <- round(
    input$similarityTable[input$similarityTable$gene == gene, ]$percentSimilarity,
    digits = 3
  )

  # names for the axis from the input
  y_name <- names(input$parameters$input)[[2]]
  x_name <- names(input$parameters$input)[[1]]

  # plots the 95 percent quantile to avoid outliers
  y_max <- stats::quantile(df$y, 0.95)
  x_max <- stats::quantile(df$x, 0.95)

  max <- max(y_max, x_max)

  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh

  dissimilarPixelsX <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDX
  dissimilarPixelsY <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDY

  # Assign color:
  #   blue = similar,
  #   yellow = greater in dataset X than Y
  #   red = greater in dataset Y than X
  #   grey = below threshold
  df$color <- dplyr::case_when(df$pixel %in% similarPixels[[1]] ~ "blue",
                               df$pixel %in% dissimilarPixelsX[[1]] ~ viridis::plasma(3)[3],
                               df$pixel %in% dissimilarPixelsY[[1]] ~ "red",
                               df$pixel %in% belowThreshPixels[[1]] ~ "grey"
                               )


  plt <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$x, y = .data$y, color = .data$color
  )) +
    ggplot2::geom_point(alpha = 0.5, size=1) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = paste0("Expression of pixel in ", x_name),
      y = paste0("Expression of pixel in ", y_name),
      title = paste(gene, " S = ", s)
    ) +
    ggplot2::xlim(0, max) +
    ggplot2::ylim(0, max) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 18),
      legend.position = "none"
    ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 2^(input$parameters$foldChange), linetype = "dashed", linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1/(2^(input$parameters$foldChange)), linetype = "dashed", linewidth = 1)

  return(plt)
}

#' Map the similarity class of every pixel
#'
#' This function visualizes the spatial distribution of gene expression
#' similarity by classifying pixels into four categories: below threshold,
#' similar, higher in the first dataset and higher in the second dataset. The
#' plot is drawn with the pixel geometries of rasterized SpatialExperiment
#' objects (\code{colData()} column \code{geometry}, which needs the sf
#' package), otherwise with points at the spatial coordinates.
#'
#' @param input A list. Results from \code{spatialSimilarity()}. This includes
#'   the similarity table, log-transformed pixel data, and analysis parameters.
#' @param gene Character. The name of the gene to visualize.
#' @param assayName A character string or numeric specifying the assay in the
#'   Spatial Experiment to use. Default is \code{NULL}: the assay that
#'   \code{spatialSimilarity()} used. The classification comes from
#'   \code{input}, so the plot does not depend on it.
#'
#' @return A ggplot2 spatial plot displaying classified pixels, where:
#' \describe{
#'   \item{\strong{blue}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{yellow}}{Pixels with greater expression in dataset X than Y.}
#'   \item{\strong{red}}{Pixels with greater expression in dataset Y than X.}
#'   \item{\strong{grey}}{Pixels with gene expression below the threshold in both experiments.}
#' }
#'   The plot title includes the gene name and its similarity score.
#'
#' @export
#'
#' @examples
#' data(speKidney)
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                 assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                 square = FALSE)
#' s <- spatialSimilarity(list(A = rastKidney$A, B = rastKidney$B))
#' pixelClass(s, "Gene")
#'
pixelClass <- function (input, gene, assayName = NULL) {
  ## the assay that spatialSimilarity() used, unless another one is given
  assayName <- .stc_similarity_assay(input, assayName)

  df <- getGenePixelDF(
    y = input$parameters$input[[2]],
    x = input$parameters$input[[1]],
    gene = gene,
    assayName = assayName
  )

  # names for the axis from the input
  y_name <- names(input$parameters$input)[[2]]
  x_name <- names(input$parameters$input)[[1]]

  # assignFill creates a new column fill, fill, in df
  df <- assignFill(input = input, gene = gene, df = df)
  df$fill <- as.factor(df$fill)

  if ("geometry" %in% names(SummarizedExperiment::colData(input$parameters$input[[1]]))) {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("pixelClass() needs the sf package to draw the pixels of rasterized objects: install.packages(\"sf\")")
  }
  df_sf <- sf::st_sf(geometry = SummarizedExperiment::colData(input$parameters$input[[1]])[df$pixel, ]$geometry,
                     row.names = df$pixel)

  df_sf <- cbind(df_sf, fill = df$fill)
  df_sf$fill <- factor(df_sf$fill, levels = c(1, 2, 3, 4))

  # Assign color:
  #   blue = similar,
  #   yellow = greater in dataset X than Y
  #   red = greater in dataset Y than X
  #   grey = below threshold

  # geom_sf() draws with coord_sf(), which keeps the aspect ratio
  plt <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = df_sf, ggplot2::aes(fill = .data$fill)) +
    ggplot2::scale_fill_manual(
      values = c("blue", viridis::plasma(3)[3], "red", "grey"),
      labels = c("similar pixels", paste0(x_name, "+ pixels"), paste0(y_name, "+ pixels"), "below threshold"),
      drop = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  s <- round(
      input$similarityTable[input$similarityTable$gene == gene, ]$percentSimilarity,
      digits = 3
  )

  plt <- plt + ggplot2::ggtitle(paste0(gene, " S = ", s))

  }

  else {

    sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(input$parameters$input[[1]])),
                              rownames(SpatialExperiment::spatialCoords(input$parameters$input[[2]])))

    dfclass <- data.frame(SpatialExperiment::spatialCoords(input$parameters$input[[1]])[sharedPixels,])
    colnames(dfclass) <- c("X", "Y")
    dfclass <- cbind(df, dfclass)
    dfclass$fill <- factor(dfclass$fill, levels = c(1, 2, 3, 4))

    plt <- ggplot2::ggplot(dfclass, ggplot2::aes(x = .data$X, y = .data$Y, color = .data$fill)) +
      ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18)+
      ggplot2::scale_color_manual(
        values = c("blue", viridis::plasma(3)[3], "red", "grey"),
        labels = c("similar pixels", paste0(x_name, "+ pixels"), paste0(y_name, "+ pixels"), "below threshold"),
        drop = FALSE
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())

    s <- round(
      input$similarityTable[input$similarityTable$gene == gene, ]$percentSimilarity,
      digits = 3
    )

    plt <- plt + ggplot2::ggtitle(paste0(gene, " S = ", s))


  }

  return(plt)
}


#' Assigns classification labels to pixels based on gene expression similarity.
#'
#' Adds a column \code{fill} to \code{df}, a data frame with a column
#' \code{pixel}: 1 for similar pixels, 2 for pixels with greater expression in
#' the first dataset, 3 for greater expression in the second dataset and 4 for
#' pixels below the threshold in both, as classified by
#' \code{spatialSimilarity()} (\code{input}) for \code{gene}.
#'
#' @noRd
assignFill <- function (input, gene, df) {

  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh

  dissimilarPixelsX <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDX
  dissimilarPixelsY <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDY

  # Assign factors:
  #   1 = similar,
  #   2 = greater in dataset X than Y
  #   3 = greater in dataset Y than X
  #   4 = below threshold
  df$fill <- dplyr::case_when(df$pixel %in% similarPixels[[1]] ~ 1,
                               df$pixel %in% dissimilarPixelsX[[1]] ~ 2,
                               df$pixel %in% dissimilarPixelsY[[1]] ~ 3,
                               df$pixel %in% belowThreshPixels[[1]] ~ 4
  )

  return(df)
}

# The assay of a spatialSimilarity() result: assayName if given, otherwise the one the result records (results
# of STcompare 0.1.0 and earlier record none: the first assay).
.stc_similarity_assay <- function(input, assayName = NULL) {
  if (!is.null(assayName)) return(assayName)
  if (!is.null(input$parameters$assayName)) return(input$parameters$assayName)
  1
}

#' Maps, pixel classes and scatter plot of genes, optionally saved as PDF
#'
#' This function creates a multi-panel plot for each gene showing spatial expression
#' patterns, pixel classifications, and the expression scatter plot. Each gene gives
#' a four-panel figure, which can also be saved as a PDF file. It needs the
#' patchwork package.
#'
#' @param geneNames Character vector. Names of genes to visualize and save.
#' @param spatialSimilarity A list. Results from \code{spatialSimilarity()}
#'   containing similarity tables and analysis parameters.
#' @param rastGexp A list of two SpatialExperiment objects: the rasterized gene
#'   expression data that \code{spatialSimilarity()} compared, in the same
#'   order. The names of the list are the titles of the expression panels.
#' @param assayName A character string or numeric specifying the assay in the
#'   Spatial Experiment to use. Default is \code{NULL}: the assay that
#'   \code{spatialSimilarity()} used. It is used for all four panels.
#' @param filePath Character. Directory where a PDF file is saved for each
#'   gene. Default is \code{FALSE}: no files are written.
#'
#' @return A list containing the arranged plot objects for each gene, with gene names
#' as list element names. Each plot object is a four-panel arrangement showing:
#' \describe{
#'   \item{\strong{Panel 1:}}{Spatial expression plot for the first experiment.}
#'   \item{\strong{Panel 2:}}{Spatial expression plot for the second experiment.}
#'   \item{\strong{Panel 3:}}{Pixel classification plot showing similarity categories (\code{\link{pixelClass}()}).}
#'   \item{\strong{Panel 4:}}{Expression scatter plot comparing the experiments (\code{\link{linearRegression}()}).}
#' }
#'
#' @details If \code{filePath} is a directory, each gene's figure is saved there
#' as a PDF file named after the gene (\code{<gene>.pdf}), 17 inches wide by 5
#' inches tall.
#'
#' @export
#'
#' @examples
#' data(speKidney)
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                 assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                 square = FALSE)
#' rastAB <- list(A = rastKidney$A, B = rastKidney$B)
#' s <- spatialSimilarity(rastAB)
#' plts <- savePlots("Gene", s, rastAB)
#' plts
#'
savePlots <- function (geneNames, spatialSimilarity, rastGexp, assayName = NULL, filePath = FALSE) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("savePlots() needs the patchwork package to arrange the panels: install.packages(\"patchwork\")")
  }
  ## the assay that spatialSimilarity() used, unless another one is given
  assayName <- .stc_similarity_assay(spatialSimilarity, assayName)
  titles <- names(rastGexp)
  if (is.null(titles)) titles <- c("", "")

  sharedPixels <- intersect(colnames(rastGexp[[1]]), colnames(rastGexp[[2]]))

  output <- list()

  for (gene in geneNames) {

    a <- .stc_expression_panel(rastGexp[[1]], gene, assayName, sharedPixels, titles[[1]])
    b <- .stc_expression_panel(rastGexp[[2]], gene, assayName, sharedPixels, titles[[2]])
    pc <- pixelClass(spatialSimilarity, gene, assayName = assayName)
    c <- linearRegression(input = spatialSimilarity, gene = gene, assayName = assayName)

    plts <- patchwork::wrap_plots(a, b, pc, c, ncol = 4)

    output[[gene]] <- plts

    # Only save to file if filePath is given
    if (!isFALSE(filePath) && !is.null(filePath)) {
      gene_plot_file <- file.path(filePath, paste0(gene, ".pdf"))
      ggplot2::ggsave(gene_plot_file,
             plts,
             width = 17, height = 5, units = "in"
      )
    }

  }

  return (output)
}

# One sample's expression of gene for savePlots(), with the assay assayName: SEraster::plotRaster() for
# rasterized objects (colData column "geometry"), otherwise points at the spatial coordinates of the pixels
# shared by the two samples. Both have a 3-inch colour bar.
.stc_expression_panel <- function(spe, gene, assayName, sharedPixels, title, textSize = 18) {
  if (!nzchar(title)) title <- NULL
  if ("geometry" %in% names(SummarizedExperiment::colData(spe))) {
    # plotRaster() takes the assay by name (NULL: the first assay)
    assay_name <- if (is.character(assayName)) assayName else SummarizedExperiment::assayNames(spe)[assayName]
    if (length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
      if (!identical(as.numeric(assayName), 1)) stop("give assayName as a name: the assays have no names")
      assay_name <- NULL
    }
    plt <- SEraster::plotRaster(spe[gene, ], assay_name = assay_name, plotTitle = title, showAxis = TRUE)
  } else {
    coords <- SpatialExperiment::spatialCoords(spe)[sharedPixels, 1:2, drop = FALSE]
    df <- data.frame(x = coords[, 1], y = coords[, 2],
                     value = as.vector(as.matrix(SummarizedExperiment::assay(spe, assayName)[gene, sharedPixels, drop = FALSE])))
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$value)) +
      ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18) +
      viridis::scale_color_viridis() +
      ggplot2::coord_fixed() +
      ggplot2::theme_void()
    if (!is.null(title)) plt <- plt + ggplot2::ggtitle(title)
  }
  bar <- ggplot2::guide_colorbar(theme = ggplot2::theme(legend.key.height = ggplot2::unit(3, "in"),
                                                        legend.key.width = ggplot2::unit(0.25, "in")))
  plt +
    ggplot2::theme(legend.text = ggplot2::element_text(size = textSize)) +
    ggplot2::guides(fill = bar, colour = bar)
}
