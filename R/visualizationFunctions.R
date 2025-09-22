
# 1. linear regression plot ###########################################

#' Generates linear regression plot for a given gene.
#'
#' This function creates a scatter plot comparing gene expression levels between
#' two spatial experiments for a specified gene. It colors data points based on
#' similarity classification and overlays fold-change threshold lines.
#'
#' @param input A list. Results from `spatialSimilarity()`. This includes the
#'   similarity table, log-transformed pixel data, and analysis parameters.
#' @param gene Character. The name of the gene to visualize.
#' @param assayName A character string or numeric specifying the assay in the
#' Spatial Experiment to use. Default is \code{NULL}. If no value is supplied
#' for \code{assayName}, then the first assay is used as a default
#'
#' @return A ggplot2 scatter plot displaying gene expression values from two
#'   spatial experiments. Data points are colored as follows:
#' \itemize{
#'   \item{\strong{blue}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{yellow}}{Pixels with greater expression in dataset X than Y.}
#'   \item{\strong{red}}{Pixels with greater expression in dataset Y than X.}
#'   \item{\strong{grey}}{Pixels with gene expression below the threshold in both experiments.}
#' }
#' The plot includes:
#' \itemize{
#'   \item{\strong{Solid line:}}{y = x (perfect correlation).}
#'   \item{\strong{Dashed lines:}}{Fold-change similarity thresholds (upper and lower bounds).}
#' }
#'
#' @export
#'
#' @examples
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                BPPARAM = BiocParallel::MulticoreParam(), square = FALSE)
#'
#' s <- spatialSimilarity(list(rastKidney$A, rastKidney$C))
#' linearRegression(s, "Gene")
#'
linearRegression <- function (input, gene, assayName = NULL) {

  ## if name of assay to use in the SpatialExperiment object is not provided, use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }

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
  y_max <- quantile(df$y, 0.95)
  x_max <- quantile(df$x, 0.95)

  max <- max(y_max, x_max)

  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh

  dissimilarPixelsX <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDX
  dissimilarPixelsY <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDY

  # Assign color: black = below = threshold, blue = similar, red = not similar
  # df$color <- ifelse(df$pixel %in% similarPixels[[1]], "blue",
  #                    ifelse(df$pixel %in% belowThreshPixels[[1]], "grey", "red"))

  # Assign color:
  #   blue = similar,
  #   yellow = greater in dataset X than Y
  #   red = greater in dataset Y than X
  #   grey = below threshold
  df$color <- dplyr::case_when(df$pixel %in% similarPixels[[1]] ~ "blue",
                               df$pixel %in% dissimilarPixelsX[[1]] ~ viridis::plasma(3)[3],
                               df$pixel %in% dissimilarPixelsY[[1]] ~ "red",
                               df$pixel %in% belowThreshPixels[[1]] ~ "grey",
                               )


  plt <- ggplot2::ggplot(df, ggplot2::aes(
    x = x, y = y, color = color,
  )) +
    # geom_point(size=2) +
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

#' Generates a spatial plot of pixel classifications for a given gene.
#'
#' This function visualizes the spatial distribution of gene expression
#' similarity by classifying pixels into three categories: below threshold,
#' similar, and not similar. The plot is generated using spatial geometry from
#' the SpatialExperiment object.
#'
#' @param input A list. Results from `spatialSimilarity()`. This includes the
#'   similarity table, log-transformed pixel data, and analysis parameters.
#' @param gene Character. The name of the gene to visualize.
#' @param assayName A character string or numeric specifying the assay in the
#'   Spatial Experiment to use. Default is \code{NULL}. If no value is supplied
#'   for \code{assayName}, then the first assay is used as a default
#'
#' @return A ggplot2 spatial plot displaying classified pixels, where:
#' \itemize{
#'   \item{\strong{blue}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{yellow}}{Pixels with greater expression in dataset X than Y.}
#'   \item{\strong{red}}{Pixels with greater expression in dataset Y than X.}
#'   \item{\strong{grey}}{Pixels with gene expression below the threshold in both experiments.}
#' }
#'   The plot is generated using `sf` (simple features) for spatial
#'   representation and is overlaid with pixel classifications. The plot title
#'   includes the gene name and its similarity score.
#'
#' @export
#'
#' @examples
#' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                 assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                 BPPARAM = BiocParallel::MulticoreParam(), square = FALSE)
#' s <- spatialSimilarity(list(rastKidney$A, rastKidney$B))
#' pixelClass(s, "Gene")
#'
pixelClass <- function (input, gene, assayName = NULL) {
  ## if name of assay to use in the SpatialExperiment object is not provided, use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }

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


  df_sf <- sf::st_sf(geometry = SummarizedExperiment::colData(input$parameters$input[[1]])[df$pixel, ]$geometry,
                     row.names = df$pixel)

  df_sf <- cbind(df_sf, fill = df$fill)
  #df_sf$fill <- factor(df_sf$fill, levels = c(1, 2, 3))
  df_sf$fill <- factor(df_sf$fill, levels = c(1, 2, 3, 4))

  # Assign color:
  #   blue = similar,
  #   yellow = greater in dataset X than Y
  #   red = greater in dataset Y than X
  #   grey = below threshold

  plt <- ggplot2::ggplot() + ggplot2::coord_fixed() +
    ggplot2::geom_sf(data = df_sf, ggplot2::aes(fill = fill)) +
    # ggplot2::scale_fill_manual(
    #   values = c("grey", "blue", "red"),
    #   labels = c("below threshold", "similar pixels", " not similar pixels"),
    #   drop = FALSE
    # ) +
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

    plt <- ggplot2::ggplot(dfclass, ggplot2::aes(x = X, y = Y, color = fill)) +
      ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18)+
      # ggplot2::scale_color_manual(
      #   values = c("grey", "blue", "red"),
      #   labels = c("below threshold", "similar pixels", " not similar pixels"),
      #   drop = FALSE
      # ) +
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
#' This function categorizes pixels into four groups based on their similarity
#' in gene expression between two spatial experiments. The assigned classification
#' is stored in a new column, `fill`, within the input data frame.
#'
#' @param input A list. Results from `spatialSimilarity()`.
#' This includes the similarity table and classified pixel IDs.
#' @param gene Character. The name of the gene for which pixel classification is performed.
#' @param df Data frame containing pixel-level gene expression data.
#'
#' @return The input data frame with an additional column, `fill`, which assigns
#' classification labels to pixels:
#' \itemize{
#'   \item{\strong{1}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{2}}{Pixels with greater expression in dataset X than Y}
#'   \item{\strong{3}}{Pixels with greater expression in dataset Y than X}
#'   \item{\strong{4}}{Pixels with gene expression below the threshold in both experiments.}
#' }
#'
#'
assignFill <- function (input, gene, df) {

  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh

  dissimilarPixelsX <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDX
  dissimilarPixelsY <- input$similarityTable[input$similarityTable$gene == gene, ]$dissimilarPixelIDY

  # # Assign color: 1 = grey (below thresh), 2 = blue (similar), 3 = red (not similar)
  # df$fill <- ifelse(df$pixel %in% similarPixels[[1]], 2,
  #                    ifelse(df$pixel %in% belowThreshPixels[[1]], 1, 3))

  # Assign factors:
  #   1 = similar,
  #   2 = greater in dataset X than Y
  #   3 = greater in dataset Y than X
  #   4 = below threshold
  df$fill <- dplyr::case_when(df$pixel %in% similarPixels[[1]] ~ 1,
                               df$pixel %in% dissimilarPixelsX[[1]] ~ 2,
                               df$pixel %in% dissimilarPixelsY[[1]] ~ 3,
                               df$pixel %in% belowThreshPixels[[1]] ~ 4,
  )

  return(df)
}


#' Generates and saves comprehensive visualization plots for specified genes.
#'
#' This function creates a multi-panel plot for each gene showing spatial expression
#' patterns, pixel classifications, and correlation analysis. Each gene generates
#' a four-panel figure that is saved as a PDF file.
#'
#' @param geneNames Character vector. Names of genes to visualize and save.
#' @param spatialSimilarity A list. Results from `spatialSimilarity()` containing
#' similarity tables and analysis parameters.
#' @param rastGexp A list of two SpatialExperiment objects. Rasterized gene expression
#' data from the spatial experiments being compared.
#' @param assayName A character string or numeric specifying the assay in the
#'   Spatial Experiment to use. Default is \code{NULL}. If no value is supplied
#'   for \code{assayName}, then the first assay is used as a default
#' @param filePath Character. Directory path where PDF files will be saved.
#'
#' @return A list containing the arranged plot objects for each gene, with gene names
#' as list element names. Each plot object is a four-panel arrangement showing:
#' \itemize{
#'   \item{\strong{Panel 1:}}{Spatial expression plot for the first experiment.}
#'   \item{\strong{Panel 2:}}{Spatial expression plot for the second experiment.}
#'   \item{\strong{Panel 3:}}{Pixel classification plot showing similarity categories.}
#'   \item{\strong{Panel 4:}}{Linear regression scatter plot comparing expression between experiments.}
#' }
#'
#' @details The function saves each gene's visualization as a PDF file named
#' "{gene_name}.pdf" in the specified directory. Plot dimensions are set to
#' 17 inches wide by 5 inches tall at 300 DPI resolution.
#'
#' @export
#'
#' @examples
#' #' ##### Rasterize to get pixels at matched spatial locations #####
#' rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
#'                 assay_name = 'counts', resolution = 0.2, fun = "mean",
#'                 BPPARAM = BiocParallel::MulticoreParam(), square = FALSE)
#' s <- spatialSimilarity(list(rastKidney$A, rastKidney$B))
#' plts <- savePlots("Gene", s, rastKidney)
#' plts


savePlots <- function (geneNames, spatialSimilarity, rastGexp, assayName = NULL, filePath = FALSE) {

  library(ggplot2)
  library(gridExtra)
  library(patchwork)


  output <- list()

  for (gene in geneNames) {

    X <- 18 # size of fill

    pc <- pixelClass(spatialSimilarity, gene)

    if ("geometry" %in% names(SummarizedExperiment::colData(rastGexp[[1]]))) {
        a <- SEraster::plotRaster(rastGexp[[1]][gene, ],
                                  plotTitle = names(rastGexp)[[1]],
                                  showAxis = TRUE) +
          ggplot2::theme(legend.text=ggplot2::element_text(size=X)) +
          ggplot2::guides(fill = ggplot2::guide_colorbar(
            barheight = unit(3, "in"),
            barwidth = unit(0.25, "in")
          ))

        b <- SEraster::plotRaster(rastGexp[[2]][gene, ],
                                  plotTitle = names(rastGexp)[[2]],
                                  showAxis = TRUE) +
          ggplot2::theme(legend.text=ggplot2::element_text(size=X)) +
          ggplot2::guides(fill = ggplot2::guide_colorbar(
            barheight = unit(3, "in"),
            barwidth = unit(0.25, "in")
          ))
    }
    else{
      ## if name of assay to use in the SpatialExperiment object is not provided, use the first assay as a default
      if (is.null(assayName)) {
        assayName <- 1
      }

      sharedPixels <- intersect(rownames(SpatialExperiment::spatialCoords(rastGexp[[1]])),
                                rownames(SpatialExperiment::spatialCoords(rastGexp[[2]])))

      dfa <- data.frame(SpatialExperiment::spatialCoords(rastGexp[[1]])[sharedPixels,],
                        color = SummarizedExperiment::assay(rastGexp[[1]], assayName)[gene, sharedPixels])
      dfb <- data.frame(SpatialExperiment::spatialCoords(rastGexp[[2]])[sharedPixels,],
                        color = SummarizedExperiment::assay(rastGexp[[1]], assayName)[gene, sharedPixels])

      a <- ggplot2::ggplot(dfa, ggplot2::aes(x = x, y = y, color = color)) +
        ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18) +
        viridis::scale_color_viridis() +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.text=ggplot2::element_text(size=X)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(
          barheight = unit(3, "in"),
          barwidth = unit(0.25, "in")
        ))

      b <- ggplot2::ggplot(dfb, ggplot2::aes(x = x, y = y, color = color)) +
        ggplot2::geom_point(size = 0.7, alpha = 1, shape = 18) +
        viridis::scale_color_viridis() +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.text=ggplot2::element_text(size=X)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(
          barheight = unit(3, "in"),
          barwidth = unit(0.25, "in")
        ))


    }

    c <- linearRegression(input = spatialSimilarity, gene = gene)

    plts <- a + b + pc + c + plot_layout(ncol = 4)

    output[[gene]] <- plts

    # Only save to file if filePath is not FALSE
    if (filePath != FALSE) {
      gene_plot_file <- file.path(filePath, paste0(gene, ".pdf"))
      ggplot2::ggsave(gene_plot_file,
             plts,
             width = 17, height = 5, dpi = 300, units = "in"
      )
    }

  }

  return (output)
}






