
# 1. linear regression plot ###########################################

#' Generates linear regression plot for a given gene.
#'
#' This function creates a scatter plot comparing gene expression levels between 
#' two spatial experiments for a specified gene. It colors data points based on 
#' similarity classification and overlays fold-change threshold lines.
#'
#' @param input A list. Results from `spatialSimilarity()`. 
#' This includes the similarity table, log-transformed pixel data, and analysis parameters.
#' @param gene Character. The name of the gene to visualize.
#'
#' @return A ggplot2 scatter plot displaying gene expression values from two spatial 
#' experiments. Data points are colored as follows:
#' \itemize{
#'   \item{\strong{Blue}}{Pixels with similar gene expression (within the fold-change threshold).}
#'   \item{\strong{Red}}{Pixels outside the similarity threshold.}
#'   \item{\strong{Black}}{Pixels below the gene expression threshold in both experiments.}
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
linearRegression <- function (input, gene) {
  
  df <- getGenePixelDF(
    y = input$parameters$input[[1]],
    x = input$parameters$input[[2]],
    gene = gene
  )  
  
  # similarity score from the similarity table 
  s <- round(
    input$similarityTable[input$similarityTable$gene == gene, ]$percentSimilarity, 
    digits = 3
  )
  
  # names for the axis from the input 
  y_name <- names(input$parameters$input)[[1]]
  x_name <- names(input$parameters$input)[[2]]
  
  # plots the 95 percent quantile to avoid outliers  
  y_max <- quantile(df[[3]], 0.95)
  x_max <- quantile(df[[2]], 0.95)
  
  max <- max(y_max, x_max)
  
  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh
  
  # Assign color: black = below = threshold, blue = similar, red = not similar 
  df$color <- ifelse(df$pixel %in% similarPixels[[1]], "blue",
                     ifelse(df$pixel %in% belowThreshPixels[[1]], "black", "red"))
  
  
  plt <- ggplot2::ggplot(df, aes(
    x = x, y = y, color = color, 
  )) + 
    # geom_point(size=2) +
    geom_point(alpha = 0.5, size=1) +
    scale_color_identity() +  
    labs(
      y = y_name,
      x = x_name
      # title = paste(gene, " S = ", s)
    ) + 
    xlim(0, max) +
    ylim(0, max) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),     
      axis.text = element_text(size = 14),     
      plot.title = element_text(size = 18),     
      legend.position = "none"              
    ) + 
    geom_abline(intercept = 0, slope = 1, linetype = "solid", size = 1) +
    geom_abline(intercept = 0, slope = 2^(input$parameters$foldChange), linetype = "dashed", size = 1) +
    geom_abline(intercept = 0, slope = 1/(2^(input$parameters$foldChange)), linetype = "dashed", size = 1)
  
  return(plt)
}

linearRegression <- function (input, gene) {
  
  df <- getGenePixelDF(
    y = input$parameters$input[[1]],
    x = input$parameters$input[[2]],
    gene = gene
  )  
  
  # similarity score from the similarity table 
  s <- round(
    input$similarityTable[input$similarityTable$gene == gene, ]$percentSimilarity, 
    digits = 3
  )
  
  # names for the axis from the input 
  y_name <- names(input$parameters$input)[[1]]
  x_name <- names(input$parameters$input)[[2]]

  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh
  
  # Assign color: black = below = threshold, blue = similar, red = not similar 
  df$color <- ifelse(df$pixel %in% similarPixels[[1]], "blue",
                     ifelse(df$pixel %in% belowThreshPixels[[1]], "black", "red"))
  
  
  plt <- ggplot2::ggplot(df, aes(
    x = x, y = y, color = color, 
  )) + 
    # geom_point(size=2) +
    geom_point(alpha = 0.5, size=1) +
    scale_color_identity() +  
    labs(
      y = y_name,
      x = x_name,
      title = paste(gene, " S = ", s)
    ) + 
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),     
      axis.text = element_text(size = 14),     
      plot.title = element_text(size = 18),     
      legend.position = "none"              
    ) + 
    geom_abline(intercept = 0, slope = 1, linetype = "solid", size = 1) +
    geom_abline(intercept = 0, slope = 2^(input$parameters$foldChange), linetype = "dashed", size = 1) +
    geom_abline(intercept = 0, slope = 1/(2^(input$parameters$foldChange)), linetype = "dashed", size = 1)
  
  return(plt)
}

#' Generates a spatial plot of pixel classifications for a given gene.
#'
#' This function visualizes the spatial distribution of gene expression similarity 
#' by classifying pixels into three categories: below threshold, similar, and not similar. 
#' The plot is generated using spatial geometry from the SpatialExperiment object.
#'
#' @param input A list. Results from `spatialSimilarity()`. 
#' This includes the similarity table, log-transformed pixel data, and analysis parameters.
#' @param gene Character. The name of the gene to visualize.
#'
#' @return A ggplot2 spatial plot displaying classified pixels, where:
#' \itemize{
#'   \item{\strong{Grey}}{Pixels with gene expression below the threshold in both experiments.}
#'   \item{\strong{Blue}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{Red}}{Pixels classified as not similar.}
#' }
#' The plot is generated using `sf` (simple features) for spatial representation 
#' and is overlaid with pixel classifications. The plot title includes the gene 
#' name and its similarity score.
#'
#' @export
#'
#' @examples
pixelClass <- function (input, gene) {
  df <- getGenePixelDF(
    y = input$parameters$input[[1]],
    x = input$parameters$input[[2]],
    gene = gene
  ) 
  
  # assignFill creates a new column fill, fill, in df 
  df <- assignFill(input = input, gene = gene, df = df)
  df$fill <- as.factor(df$fill)
  
  
  df_sf <- sf::st_sf(geometry = colData(input$parameters$input[[1]])[df$pixel, ]$geometry, 
                     row.names = df$pixel)
  
  df_sf <- cbind(df_sf, fill = df$fill)
  df_sf$fill <- factor(df_sf$fill, levels = c(1, 2, 3))
  
  plt <- ggplot2::ggplot() + ggplot2::coord_fixed() +
    ggplot2::geom_sf(data = df_sf, ggplot2::aes(fill = fill)) +
    scale_fill_manual(
      values = c("grey", "blue", "red"),
      labels = c("below threshold", "similar pixels", " not similar pixels"),
      drop = FALSE
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) 

  s <- round(
      input$similarityTable[input$similarityTable$gene == gene, ]$percentSimilarity, 
      digits = 3
  )
  
  plt <- plt + ggplot2::ggtitle(paste0(gene, " S = ", s))
  
  return(plt)
}


#' Assigns classification labels to pixels based on gene expression similarity.
#'
#' This function categorizes pixels into three groups based on their similarity 
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
#'   \item{\strong{1 (Grey)}}{Pixels with gene expression below the threshold in both experiments.}
#'   \item{\strong{2 (Blue)}}{Pixels classified as similar (within the fold-change threshold).}
#'   \item{\strong{3 (Red)}}{Pixels classified as not similar.}
#' }
#'
#' @export
#'
#' @examples
assignFill <- function (input, gene, df) {
  
  similarPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$similarPixelID
  belowThreshPixels <- input$similarityTable[input$similarityTable$gene == gene, ]$pixelIDOutThresh
  
  # Assign color: 1 = grey (below thresh), 2 = blue (similar), 3 = red (not similar)
  df$fill <- ifelse(df$pixel %in% similarPixels[[1]], 2,
                     ifelse(df$pixel %in% belowThreshPixels[[1]], 1, 3))
  
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

savePlots <- function (geneNames, spatialSimilarity, rastGexp, filePath = FALSE) {
  
  library(ggplot2)
  library(gridExtra)
  library(patchwork)
  
  
  output <- list()
  
  for (gene in geneNames) {
    
    X <- 18 # size of fill 
    
    pc <- pixelClass(spatialSimilarity, gene) 
    
    a <- SEraster::plotRaster(rastGexp[[1]][gene, ], 
                              plotTitle = names(rastGexp)[[1]], 
                              showAxis = TRUE) + 
      theme(legend.text=element_text(size=X)) +  
      guides(fill = guide_colorbar(
        barheight = unit(3, "in"),  
        barwidth = unit(0.25, "in")   
      ))
    
    b <- SEraster::plotRaster(rastGexp[[2]][gene, ], 
                              plotTitle = names(rastGexp)[[2]], 
                              showAxis = TRUE) + 
      theme(legend.text=element_text(size=X)) + 
      guides(fill = guide_colorbar(
        barheight = unit(3, "in"), 
        barwidth = unit(0.25, "in")    
      ))
    
    c <- linearRegression(input = spatialSimilarity, gene = gene)
    
    plts <- a + b + pc + c + plot_layout(ncol = 4)
    
    output[[gene]] <- plts
    
    # Only save to file if filePath is not FALSE
    if (filePath != FALSE) {
      gene_plot_file <- file.path(filePath, paste0(gene, ".pdf"))
      ggsave(gene_plot_file,
             plts,
             width = 17, height = 5, dpi = 300, units = "in"
      )
    }
    
  }
  
  return (output)
}






