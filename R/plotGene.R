# ---------------------------------------------------------------------------------------------------
# Standard figures of compareSpatial() results: plotGene() for one gene and plotOverview() for all genes.
# Every panel is a ggplot2 plot; the figure is a patchwork of them, so that it can be printed, saved with
# ggplot2::ggsave() or changed with patchwork's operators. Every text size is relative to baseSize, so
# that a figure can be drawn at a journal's width (baseSize = 7 for 180 mm) and rescaled with
# `& theme(text = element_text(size = ...))`.
# ---------------------------------------------------------------------------------------------------

# Colours of the four similarity classes (Okabe-Ito: blue, orange, reddish purple; the pixels below the
# thresholds are grey) and of the significance classes of plotOverview() (vermillion, blue, grey). Every
# pair of coloured classes is separable under the common colour-vision deficiencies.
.stc_class_colours <- c(similar = "#0072B2", x = "#E69F00", y = "#CC79A7", below = "#BDBDBD")
.stc_sig_colours <- c(positive = "#D55E00", negative = "#0072B2", ns = "#BDBDBD")

# The pixels that one sample lacks (drawn under every map, without a legend entry), and the pixels with
# missing values.
.stc_unshared_fill <- "grey85"
.stc_na_fill <- "grey60"

#' Standard figure of one gene in a compareSpatial() result
#'
#' @description One figure for one gene of a \code{\link{compareSpatial}()}
#'   result: the gene's expression in the two samples at native resolution
#'   (the cells or spots, when the objects before rasterization are given),
#'   on the shared grid of pixels, the similarity class of every pixel, and
#'   the expression at matched pixels as a scatter plot labelled with the
#'   gene's statistics. The panels are lettered and the figure is ready to
#'   print, to save with \code{ggplot2::ggsave()} or to adjust with patchwork
#'   and ggplot2.
#'
#' @details \strong{Panels.} In this order, each map with the sample name as
#'   its title:
#'   \describe{
#'   \item{\code{"native"}}{Two maps of the cells or spots of \code{native},
#'   coloured by the gene's expression in \code{nativeAssay}. Drawn only when
#'   \code{native} is given.}
#'   \item{\code{"raster"}}{Two maps of the shared pixels (those the result
#'   used), coloured by the gene's expression in \code{assay}, the values the
#'   tests compared. The pixels are drawn with their geometry when the objects
#'   were rasterized with SEraster (\code{colData()} column \code{geometry},
#'   which needs the sf package), otherwise as square tiles the size of the
#'   distance between neighbouring pixels. The pixels that one sample lacks
#'   (present in the other sample only, so not used by the tests) are drawn
#'   in light grey in every map, so that the outline of the tissue is
#'   complete and the maps are comparable.}
#'   \item{\code{"similarity"}}{The shared pixels coloured by their class in
#'   the similarity test, recomputed from the data with the result's settings
#'   (\code{foldChange}, \code{minPixels}, and the gene's thresholds
#'   \code{t1} and \code{t2}):
#'   similar (\eqn{|\log_2(y/x)| \le}{|log2(y / x)| <=} \code{foldChange}),
#'   higher in the first sample, higher in the second sample, or below the
#'   threshold in both. The legend gives the share of each class among the
#'   pixels above the threshold, so that the share of similar pixels is the
#'   gene's \code{percentSimilarity}. Drawn when the result has the
#'   similarity columns and the gene's \code{percentSimilarity} is not
#'   \code{NA}: not for a gene with negative values (a fold change needs
#'   non-negative values), nor for a gene with fewer pixels above the
#'   threshold than \code{minPixels} allows; the scatter plot then says why.}
#'   \item{\code{"scatter"}}{The expression at matched pixels, the first
#'   sample horizontally and the second vertically, on equal axes, coloured
#'   by the same classes (plain when they are not available), with the line
#'   \eqn{y = x} (solid), the fold-change band (dashed, slopes
#'   \eqn{2^{foldChange}}{2^foldChange} and \eqn{2^{-foldChange}}{2^-foldChange})
#'   and the thresholds (dotted); the figure's caption is the key to these
#'   lines. Its subtitle gives the gene's row of the
#'   result: \code{correlationCoef} (\code{r}), \code{pValuePermute}
#'   (\code{p}), \code{pValuePermuteAdj} (\code{adjusted p}),
#'   \code{pValueNaive} (\code{naive p}), the number of permutations and
#'   \code{percentSimilarity} (\code{similarity}), each to two significant
#'   digits; \code{NA}, with the reason from the result's \code{message}
#'   column, for a gene that was skipped or failed, whose status is then
#'   shown. \code{p} is marked \code{(minimum)} when it is
#'   \eqn{1 / (nPermutations + 1)}, the smallest p-value the permutations can
#'   give, and a p-value of 0 is printed as \code{< 2.2e-16}.}
#'   }
#'   A result without the correlation columns (\code{tests = "similarity"})
#'   gives a scatter subtitle with the similarity alone; a result without the
#'   similarity columns (\code{tests = "correlation"}) has no similarity panel
#'   and a plain scatter plot.
#'
#'   \strong{Layout and size.} Without native maps the panels are in one row;
#'   with them, the four expression maps are in the first row and the
#'   similarity map and the scatter plot, larger, in a second row. The maps
#'   keep the shape of the tissue, so the figure is as wide as the tissue
#'   makes it and the panels are centred in the canvas. With the default
#'   \code{baseSize = 10}, \code{ggplot2::ggsave(..., width = 11, height =
#'   3.5)} suits one row and \code{width = 11, height = 6.5} two rows;
#'   reduce the width for a tissue that is taller than wide. For a journal
#'   figure 180 mm wide, use \code{baseSize = 7} with \code{ggsave(...,
#'   width = 180, height = 110, units = "mm")} for two rows or
#'   \code{height = 60} for one row; for a single column of 90 mm, use
#'   \code{baseSize = 6} and \code{ncol = 2}. Every text size follows
#'   \code{baseSize}, and patchwork's \code{&} operator rescales the whole
#'   figure: \code{plotGene(...) & ggplot2::theme(text =
#'   ggplot2::element_text(size = 7))}. The layout can be changed with the
#'   arguments of \code{patchwork::wrap_plots()} passed through \code{...},
#'   for example \code{ncol = 2}, and the figure with patchwork's operators,
#'   for example \code{plotGene(...) + patchwork::plot_annotation(tag_levels =
#'   NULL)} to remove the letters.
#'
#'   \strong{Colour scales.} The expression maps use the viridis palette.
#'   With \code{sharedScale = TRUE} the two native maps share one scale (one
#'   legend) and so do the two raster maps, so that the colours of the two
#'   samples are comparable. When the samples' values are on different scales
#'   (for example raw counts in one and normalized values in the other, or
#'   very different sequencing depths), \code{sharedScale = FALSE} gives every
#'   map its own scale and legend. \code{transform = "log1p"} draws the
#'   expression maps and the scatter axes on a \eqn{\log(1 + v)} scale, which
#'   spreads out the low values of counts; the classes and the statistics are
#'   unchanged, and the fold-change band is then curved. It needs values that
#'   are not negative.
#'
#'   \strong{Checks.} The objects must be those that \code{compareSpatial()}
#'   compared: the function warns when they share another number of pixels
#'   than the result used (the classes are then recomputed on the pixels
#'   given), and when \code{x} and \code{y} appear to be swapped (the gene's
#'   thresholds \code{t1} and \code{t2} are those of the other sample). A
#'   named \code{native} list is put in the order of the result's samples
#'   when its names are the sample names, and taken in order, with a warning,
#'   when they are other names.
#'
#' @param result A \code{\link{compareSpatial}()} result.
#' @param gene The gene (one row name of \code{result}) to plot.
#' @param x,y The two rasterized \code{SpatialExperiment} objects that
#'   \code{compareSpatial()} compared, in the same order (or, as there, a list
#'   of the two as \code{x} with \code{y = NULL}). The sample names come from
#'   \code{result}.
#' @param native Optional: a list of the two \code{SpatialExperiment} objects
#'   before rasterization, in the same order and in the same coordinate frame
#'   as the rasterized objects (for example aligned with STalign), for the
#'   native-resolution maps. Default \code{NULL}: no native maps.
#' @param assay The assay of \code{x} and \code{y} to draw: a name or an
#'   index, or two (one per object). Default \code{NULL}: the assay that
#'   \code{compareSpatial()} used, from \code{attr(result, "params")}.
#' @param nativeAssay The assay of the \code{native} objects to draw, a name
#'   or an index, or two. Default \code{NULL}: the first assay.
#' @param panels The panels to draw, in the order of Details; the figure keeps
#'   that order. Default: all. \code{"native"} needs \code{native} and is
#'   dropped without it; \code{"similarity"} is dropped when the classes
#'   cannot be computed (Details).
#' @param transform The scale of the expression maps and of the scatter axes:
#'   \code{"none"} (default) or \code{"log1p"}.
#' @param sharedScale One colour scale for the two maps of a pair
#'   (default \code{TRUE}), or one per map (\code{FALSE}); see Details.
#' @param pointSize The size of the points of the native maps. Default
#'   \code{NULL}: from the number of cells or spots, between 0.1 and 2 at
#'   \code{baseSize = 10}, scaled with \code{baseSize}.
#' @param nativeUnit What the points of the native maps are, for their
#'   subtitles: \code{"cells"} (default), \code{"spots"} for Visium, or two
#'   (one per sample).
#' @param baseSize The base font size in points: the size of the axis titles,
#'   from which every other text, the legend keys and the points are scaled.
#'   Default 10, for a figure about 11 inches wide; 7 for a figure 180 mm
#'   wide (Details).
#' @param ... Passed to \code{patchwork::wrap_plots()}: \code{ncol},
#'   \code{nrow}, \code{design}, \code{widths}, \code{heights}, \code{guides}.
#'
#' @return A patchwork object (a ggplot2 figure of the panels, lettered,
#'   titled with the gene and subtitled with the samples). Printing it draws
#'   the figure; \code{ggplot2::ggsave()} saves it.
#'
#' @seealso \code{\link{plotOverview}()} for the figure of a whole result;
#'   \code{\link{compareSpatial}()} for the columns the label shows. The
#'   legacy plots of \code{spatialSimilarity()} results:
#'   \code{\link{linearRegression}()}, \code{\link{pixelClass}()},
#'   \code{\link{savePlots}()}; and \code{\link{plotCorrelationGeneExp}()}.
#'
#' @export
#'
#' @examples
#' data(speKidney)
#' rast <- SEraster::rasterizeGeneExpression(speKidney, assay_name = "counts",
#'                                           resolution = 0.2, fun = "mean",
#'                                           square = FALSE)
#' ab <- compareSpatial(list(A = rast$A, B = rast$B), nPermutations = 999,
#'                      nThreads = 2)
#'
#' # the pixels, the similarity classes and the scatter plot
#' plotGene(ab, "Gene", rast$A, rast$B)
#'
#' # with the cells at native resolution, on a log1p colour scale
#' plotGene(ab, "Gene", rast$A, rast$B, native = speKidney[c("A", "B")],
#'          transform = "log1p")
#'
#' # save it: without native maps an 11 x 3.5 inch figure, with them 11 x 6.5;
#' # for a journal figure 180 mm wide, baseSize = 7
#' \dontrun{
#' ggplot2::ggsave("Gene.pdf", plotGene(ab, "Gene", rast$A, rast$B),
#'                 width = 11, height = 3.5)
#' ggplot2::ggsave("Gene-180mm.pdf",
#'                 plotGene(ab, "Gene", rast$A, rast$B, baseSize = 7),
#'                 width = 180, height = 60, units = "mm")
#' }
plotGene <- function(result, gene, x, y = NULL, native = NULL, assay = NULL, nativeAssay = NULL,
                     panels = c("native", "raster", "similarity", "scatter"),
                     transform = c("none", "log1p"), sharedScale = TRUE, pointSize = NULL,
                     nativeUnit = "cells", baseSize = 10, ...) {
  # ---- arguments ----
  if (!is.data.frame(result) || !"gene" %in% names(result)) {
    stop("result must be a compareSpatial() result (a data frame with a gene column)")
  }
  if (!is.character(gene) || length(gene) != 1L || is.na(gene)) stop("gene must be one gene name")
  if (!gene %in% result$gene) {
    stop(sprintf("gene %s is not in result; its genes are %s", gene, .stc_some(result$gene, 5L)))
  }
  if (!is.character(panels) || !length(panels) || !all(panels %in% c("native", "raster", "similarity", "scatter"))) {
    stop("panels must be some of \"native\", \"raster\", \"similarity\" and \"scatter\"")
  }
  panels <- intersect(c("native", "raster", "similarity", "scatter"), panels)
  transform <- match.arg(transform)
  if (!is.logical(sharedScale) || length(sharedScale) != 1L || is.na(sharedScale)) stop("sharedScale must be TRUE or FALSE")
  if (!is.null(pointSize) && (!is.numeric(pointSize) || length(pointSize) != 1L || is.na(pointSize) || pointSize <= 0)) {
    stop("pointSize must be NULL or a positive number")
  }
  if (!is.character(nativeUnit) || !length(nativeUnit) %in% 1:2 || anyNA(nativeUnit)) {
    stop("nativeUnit must be one or two words, for example \"cells\" or \"spots\"")
  }
  nativeUnit <- rep_len(nativeUnit, 2L)
  .stc_check_base_size(baseSize)
  if (!is.null(native) && (!is.list(native) || inherits(native, "SpatialExperiment") || length(native) != 2L)) {
    stop("native must be NULL or a list of two SpatialExperiment objects")
  }
  if ("native" %in% panels && is.null(native)) panels <- setdiff(panels, "native")
  if (!length(panels)) stop("no panel to draw: give native for the native maps, or other panels")
  pr <- attr(result, "params")
  row <- as.data.frame(result)[match(gene, result$gene), , drop = FALSE]

  # ---- the rasterized samples, the gene's values on the shared pixels ----
  if (!is.null(y) && is.list(x) && !inherits(x, "SpatialExperiment")) {
    stop("give x and y as two SpatialExperiment objects, or a list of two as x with y = NULL")
  }
  objs <- if (is.null(y)) x else list(x, y)
  if (!is.list(objs) || inherits(objs, "SpatialExperiment") || length(objs) != 2L) {
    stop("give the two rasterized SpatialExperiment objects: plotGene(result, gene, x, y), or a list of two as x")
  }
  lab <- if (!is.null(pr$samples)) pr$samples else if (is.null(y) && !is.null(names(x))) names(x) else c("x", "y")
  for (i in 1:2) {
    if (!inherits(objs[[i]], "SpatialExperiment")) stop(sprintf("%s must be a SpatialExperiment object", lab[i]))
    if (!gene %in% rownames(objs[[i]])) stop(sprintf("gene %s is not in %s", gene, lab[i]))
  }
  if (is.null(assay)) assay <- .stc_params_assay(pr$assay, objs)
  d <- .stc_compare_input(objs[[1]], objs[[2]], assay, gene, verbose = FALSE)
  same_pixels <- is.null(pr$nPixels) || pr$nPixels == length(d$pixels)
  if (!same_pixels) {
    warning(sprintf(paste("x and y share %d pixels, but the result was computed on %d: are they the objects that",
                          "compareSpatial() compared? The similarity classes are recomputed on the pixels given"),
                    length(d$pixels), pr$nPixels), call. = FALSE)
  }
  vx <- d$X[, 1L]
  vy <- d$Y[, 1L]
  if (transform == "log1p" && any(c(vx, vy) < 0, na.rm = TRUE)) {
    stop("transform = \"log1p\" needs values that are not negative (counts, or normalized values on a linear scale)")
  }
  layers <- .stc_pixel_layers(objs, d$pixels)
  fc <- if (!is.null(pr$foldChange)) pr$foldChange else 1
  has_thresholds <- all(c("t1", "t2") %in% names(row)) && !is.na(row$t1) && !is.na(row$t2)
  # the gene's thresholds are quantiles of x and of y: swapped objects swap them
  if (has_thresholds && same_pixels && !is.null(pr$minQuantile)) {
    q <- c(.stc_col_quantile(d$X, pr$minQuantile), .stc_col_quantile(d$Y, pr$minQuantile))
    same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-8))
    if (!(same(q[1L], row$t1) && same(q[2L], row$t2)) && same(q[1L], row$t2) && same(q[2L], row$t1)) {
      warning(sprintf("x and y appear to be swapped: the thresholds of the gene in the result are those of %s for x and of %s for y",
                      lab[2L], lab[1L]), call. = FALSE)
    }
  }

  # ---- the similarity classes, recomputed with the result's settings and the gene's thresholds ----
  cls <- NULL
  sim_note <- NULL  # why the similarity is NA, for the statistics
  if ("percentSimilarity" %in% names(row)) {
    msg <- if (!is.null(row$message) && !is.na(row$message)) row$message else ""
    if (has_thresholds && all(c(vx, vy) >= 0, na.rm = TRUE)) {
      min_px <- if (!is.null(pr$minPixels)) pr$minPixels else 0
      s <- .stc_similarity(d$X, d$Y, minPixels = min_px, foldChange = fc, t1 = row$t1, t2 = row$t2, details = TRUE)
      det <- attr(s, "details")
      if (det$scored[1L] && !is.na(row$percentSimilarity)) {
        l <- det$log[[1L]]
        cls <- rep("below", length(vx))
        cls[det$kept[[1L]]] <- ifelse(is.na(l) | (l >= -fc & l <= fc), "similar", ifelse(l < -fc, "x", "y"))
        cls <- factor(cls, levels = names(.stc_class_colours))
      } else if (is.na(row$percentSimilarity) && !det$scored[1L]) {
        sim_note <- sprintf("%d of %d pixels above threshold, fewer than minPixels = %s", length(det$kept[[1L]]),
                            length(vx), format(min_px))
      }
    }
    if (is.null(cls) && is.null(sim_note) && grepl("negative values", msg, fixed = TRUE)) sim_note <- "negative values"
  }
  if (is.null(cls)) panels <- setdiff(panels, "similarity")
  cls_labels <- if (!is.null(cls)) .stc_class_labels(cls, lab)

  # ---- panels ----
  plots <- list()
  maps <- character(0)
  if ("native" %in% panels) {
    native <- .stc_native_order(native, lab)
    if (is.null(nativeAssay)) nativeAssay <- 1
    if (length(nativeAssay) == 1L) nativeAssay <- rep(nativeAssay, 2L)
    nat <- lapply(1:2, function(i) .stc_native_values(native[[i]], gene, nativeAssay[[i]], lab[i]))
    if (transform == "log1p" && any(unlist(lapply(nat, `[[`, "value")) < 0, na.rm = TRUE)) {
      stop("transform = \"log1p\" needs values that are not negative; the native objects have negative values")
    }
    lims <- if (sharedScale) range(unlist(lapply(nat, `[[`, "value")), na.rm = TRUE)
    box <- apply(do.call(rbind, lapply(nat, `[[`, "pos")), 2L, range)
    for (i in 1:2) {
      plots[[length(plots) + 1L]] <- .stc_native_map(nat[[i]], lab[i], nativeUnit[i], lims, box, transform, pointSize,
                                                     baseSize, legend = !sharedScale || i == 2L)
    }
    maps <- c(maps, "native", "native")
  }
  if ("raster" %in% panels) {
    lims <- if (sharedScale) range(c(vx, vy), na.rm = TRUE)
    for (i in 1:2) {
      p <- .stc_pixel_map(layers[[i]], if (i == 1L) vx else vy) +
        ggplot2::scale_fill_viridis_c(name = d$assay[i], limits = lims, transform = .stc_transform(transform),
                                      breaks = .stc_breaks(range(c(vx, vy), na.rm = TRUE), transform), na.value = .stc_na_fill) +
        ggplot2::labs(title = lab[i], subtitle = sprintf("%d pixels", length(d$pixels))) +
        .stc_map_theme(baseSize, legend = !sharedScale || i == 2L)
      plots[[length(plots) + 1L]] <- p
    }
    maps <- c(maps, "raster", "raster")
  }
  if ("similarity" %in% panels) {
    plots[[length(plots) + 1L]] <- .stc_pixel_map(layers[[1L]], cls, legend = TRUE) +
      ggplot2::scale_fill_manual(name = NULL, values = .stc_class_colours, labels = cls_labels, drop = FALSE) +
      ggplot2::labs(title = "Similarity classes", subtitle = sprintf("similarity %s", .stc_fmt(row$percentSimilarity))) +
      .stc_map_theme(baseSize, legend = TRUE)
    maps <- c(maps, "similarity")
  }
  if ("scatter" %in% panels) {
    plots[[length(plots) + 1L]] <- .stc_scatter_panel(vx, vy, cls, cls_labels, lab, d$assay, fc, row, transform,
                                                      baseSize, legend = !"similarity" %in% panels,
                                                      stats = .stc_stat_label(row, sim_note))
    maps <- c(maps, "scatter")
  }

  # ---- the figure ----
  # one row, or the expression maps in a first row and the similarity map and the scatter plot in a taller
  # second row (each row its own patchwork, whose cells fit its fixed-aspect panels, so that the letters and
  # titles sit next to the panels; a shorter row is centred); the arguments of wrap_plots() in ... replace
  # this layout
  dots <- list(...)
  top <- sum(maps %in% c("native", "raster"))
  bottom <- length(maps) - top
  fig <- if (any(c("ncol", "nrow", "design", "widths", "heights") %in% names(dots))) {
    do.call(patchwork::wrap_plots, c(plots, dots))
  } else if ("native" %in% maps && bottom > 0L) {
    patchwork::wrap_plots(
      do.call(patchwork::wrap_plots, c(plots[seq_len(top)], list(nrow = 1L))),
      do.call(patchwork::wrap_plots, c(plots[top + seq_len(bottom)], list(nrow = 1L))),
      ncol = 1L, heights = c(1, 1.25), guides = dots$guides)
  } else {
    do.call(patchwork::wrap_plots, c(plots, list(nrow = 1L), dots))
  }
  fig + patchwork::plot_annotation(
    title = gene,
    subtitle = sprintf("%s vs %s, %d shared pixels", lab[1], lab[2], length(d$pixels)),
    caption = if ("scatter" %in% panels) .stc_scatter_caption(fc, !is.null(cls)),
    tag_levels = "A",
    theme = .stc_figure_theme(baseSize, hjust = 0.5))
}

#' Overview figure of a compareSpatial() result
#'
#' @description One figure for all the genes of a \code{\link{compareSpatial}()}
#'   result: the correlation of every gene against its adjusted p-value, the
#'   distribution of the similarity, and the correlation against the
#'   similarity, with the significant genes coloured by the sign of their
#'   correlation.
#'
#' @details The panels, each drawn when the result has its columns:
#'   \describe{
#'   \item{A}{\code{correlationCoef} against
#'   \eqn{-\log_{10}}{-log10}(\code{pValuePermuteAdj}), with a dashed line at
#'   \code{alpha}. The genes at the top are those whose p-value reached the
#'   limit of the permutations (\code{nPermutations}): their p-values are not
#'   resolved below \eqn{1 / (nPermutations + 1)}, so their heights are tied.}
#'   \item{B}{A histogram of \code{percentSimilarity} (bins of 0.05 over 0 to
#'   1), stacked by significance when the correlation was tested. Genes with
#'   \code{NA} similarity are not drawn; the subtitle counts them.}
#'   \item{C}{\code{correlationCoef} against \code{percentSimilarity}.}
#'   }
#'   A gene is significant when \code{pValuePermuteAdj < alpha}, and is
#'   coloured by the sign of \code{correlationCoef}; the other genes
#'   (including those skipped or failed) are grey. The genes named in
#'   \code{label} are labelled in panels A and C, where their labels do not
#'   overlap each other.
#'
#'   With the default \code{baseSize = 10}, \code{ggplot2::ggsave(..., width
#'   = 11, height = 4)} suits the three panels; for a journal figure 180 mm
#'   wide use \code{baseSize = 7} and \code{height = 70, units = "mm"}, and
#'   for a single column of 90 mm \code{baseSize = 6} with \code{ncol = 1}
#'   (\code{height = 180}). The panels can be arranged with the arguments of
#'   \code{patchwork::wrap_plots()} passed through \code{...}.
#'
#' @param result A \code{\link{compareSpatial}()} result.
#' @param alpha The threshold on \code{pValuePermuteAdj} below which a gene is
#'   significant. Default 0.05.
#' @param label Genes to label in the scatter plots. Default \code{NULL}: none.
#' @param baseSize The base font size in points, from which every text, the
#'   points and the legend are scaled. Default 10.
#' @param ... Passed to \code{patchwork::wrap_plots()}: \code{ncol},
#'   \code{nrow}, \code{widths}, \code{heights}, \code{design}.
#'
#' @return A patchwork object (a ggplot2 figure of the panels, lettered, with
#'   one legend), to print or to save with \code{ggplot2::ggsave()}.
#'
#' @seealso \code{\link{plotGene}()} for the figure of one gene.
#'
#' @export
#'
#' @examples
#' data(simRanPatternRasts)
#' # two samples of 30 simulated genes on the pixels shared by 60 of the rasters: 25
#' # pairs of independent patterns, four shared patterns (genes 1 to 4) and one
#' # opposite pattern (gene 5)
#' px <- Reduce(intersect, lapply(simRanPatternRasts[1:60], colnames))
#' values <- function(k) {
#'   t(vapply(simRanPatternRasts[k], function(s) SummarizedExperiment::assay(s)[1, px],
#'            numeric(length(px)))) - 6
#' }
#' make <- function(m) {
#'   rownames(m) <- paste0("gene", seq_len(nrow(m)))
#'   SpatialExperiment::SpatialExperiment(
#'     assays = list(pixelval = m),
#'     spatialCoords = SpatialExperiment::spatialCoords(simRanPatternRasts[[1]])[px, ])
#' }
#' x <- values(1:30)
#' y <- values(31:60)
#' y[1:4, ] <- x[1:4, ] + y[1:4, ] / 2
#' y[5, ] <- 8 - x[5, ] + y[5, ] / 2
#' res <- compareSpatial(list(S1 = make(x), S2 = make(y)), nPermutations = 199,
#'                       nThreads = 2)
#' plotOverview(res, label = c("gene1", "gene5"))
#' plotGene(res, "gene5", make(x), make(y))
plotOverview <- function(result, alpha = 0.05, label = NULL, baseSize = 10, ...) {
  if (!is.data.frame(result) || !"gene" %in% names(result)) {
    stop("result must be a compareSpatial() result (a data frame with a gene column)")
  }
  .stc_check_unit_interval(alpha, "alpha", zero = FALSE)
  if (!is.null(label) && (!is.character(label) || anyNA(label))) stop("label must be NULL or a character vector of genes")
  .stc_check_base_size(baseSize)
  df <- as.data.frame(result)
  has_cor <- all(c("correlationCoef", "pValuePermuteAdj") %in% names(df))
  has_sim <- "percentSimilarity" %in% names(df)
  if (!has_cor && !has_sim) stop("result has neither the correlation nor the similarity columns")
  pr <- attr(result, "params")
  miss <- setdiff(label, df$gene)
  if (length(miss)) warning(sprintf("%d of the genes to label are not in result: %s", length(miss), .stc_some(miss, 5L)), call. = FALSE)
  df$label <- ifelse(df$gene %in% label, df$gene, NA_character_)

  # significance classes (grey without the correlation test)
  df$class <- "ns"
  if (has_cor) {
    sig <- !is.na(df$pValuePermuteAdj) & df$pValuePermuteAdj < alpha
    df$class[sig & df$correlationCoef > 0] <- "positive"
    df$class[sig & df$correlationCoef < 0] <- "negative"
  }
  df$class <- factor(df$class, levels = names(.stc_sig_colours))
  n_sig <- table(df$class)
  sig_labels <- c(positive = sprintf("r > 0 (%d)", n_sig[["positive"]]),
                  negative = sprintf("r < 0 (%d)", n_sig[["negative"]]),
                  ns = sprintf("not significant (%d)", n_sig[["ns"]]))
  # one fill scale for the points (filled circles) and the histogram; only the first panel shows its legend,
  # which patchwork collects, with a key for every class (show.legend = TRUE), present in the panel or not
  colours <- function() {
    ggplot2::scale_fill_manual(name = sprintf("adjusted p < %s", format(alpha)), values = .stc_sig_colours,
                               labels = sig_labels, drop = FALSE, guide = if (has_cor) "legend" else "none")
  }
  scale <- baseSize / 10
  base <- ggplot2::theme_classic(base_size = baseSize) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
                   legend.key.size = ggplot2::unit(1.2 * baseSize, "pt"),
                   plot.title = ggplot2::element_text(size = ggplot2::rel(1), face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85), colour = "grey30"))
  # the genes to label, drawn on top of the points
  label_layer <- function(mapping) {
    if (!any(!is.na(df$label))) return(NULL)
    ggplot2::geom_text(data = df[!is.na(df$label), , drop = FALSE], mapping = mapping, size = 3 * scale, vjust = -0.7,
                       check_overlap = TRUE, show.legend = FALSE, na.rm = TRUE)
  }
  plots <- list()
  if (has_cor) {
    ok <- !is.na(df$pValuePermuteAdj) & !is.na(df$correlationCoef)
    plots$volcano <- ggplot2::ggplot(df[ok, , drop = FALSE],
                                     ggplot2::aes(x = .data$correlationCoef, y = -log10(.data$pValuePermuteAdj))) +
      ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "grey50") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
      ggplot2::geom_point(ggplot2::aes(fill = .data$class), shape = 21, stroke = 0, size = 2 * scale, alpha = 0.8,
                          show.legend = TRUE) +
      label_layer(ggplot2::aes(label = .data$label)) +
      colours() +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = "correlation (correlationCoef)", y = "-log10(adjusted p)",
                    title = "Spatial correlation",
                    subtitle = sprintf("%d of %d genes tested, %d significant", sum(ok), nrow(df), sum(sig))) +
      base
  }
  if (has_sim) {
    ok <- !is.na(df$percentSimilarity)
    plots$hist <- ggplot2::ggplot(df[ok, , drop = FALSE], ggplot2::aes(x = .data$percentSimilarity)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = .data$class), binwidth = 0.05, boundary = 0, closed = "right",
                              colour = "white", linewidth = 0.2, show.legend = FALSE) +
      colours() +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      ggplot2::labs(x = "similarity (percentSimilarity)", y = "genes", title = "Spatial similarity",
                    subtitle = sprintf("%d genes%s", sum(ok),
                                       if (any(!ok)) sprintf(", NA for %d", sum(!ok)) else "")) +
      base
    if (has_cor) {
      plots$both <- ggplot2::ggplot(df[ok & !is.na(df$correlationCoef), , drop = FALSE],
                                    ggplot2::aes(x = .data$correlationCoef, y = .data$percentSimilarity)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
        ggplot2::geom_point(ggplot2::aes(fill = .data$class), shape = 21, stroke = 0, size = 2 * scale, alpha = 0.8,
                            show.legend = FALSE) +
        label_layer(ggplot2::aes(label = .data$label)) +
        colours() +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(x = "correlation (correlationCoef)", y = "similarity (percentSimilarity)",
                      title = "Correlation and similarity") +
        base
    }
  }
  dots <- list(...)
  if (!any(c("ncol", "nrow", "design") %in% names(dots))) dots$nrow <- 1L
  if (is.null(dots$guides)) dots$guides <- "collect"
  fig <- do.call(patchwork::wrap_plots, c(unname(plots), dots))
  fig + patchwork::plot_annotation(
    title = if (!is.null(pr$samples)) sprintf("%s vs %s", pr$samples[1], pr$samples[2]),
    subtitle = sprintf("%d gene%s%s", nrow(df), if (nrow(df) == 1L) "" else "s",
                       if (!is.null(pr$nPixels)) sprintf(" on %d shared pixels", pr$nPixels) else ""),
    tag_levels = "A",
    theme = .stc_figure_theme(baseSize, hjust = 0))
}

# ---------------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------------

.stc_check_base_size <- function(baseSize) {
  if (!is.numeric(baseSize) || length(baseSize) != 1L || is.na(baseSize) || baseSize <= 0) {
    stop("baseSize must be a positive number (a font size in points)")
  }
  invisible(baseSize)
}

# The assay to draw by default: the one the result records (by name, or by index as a character string
# when the assay had no name, which is turned back into an index), or the first assay.
.stc_params_assay <- function(recorded, objs) {
  if (is.null(recorded)) return(1)
  assay <- rep_len(recorded, 2L)
  if (is.character(assay)) {
    named <- vapply(1:2, function(i) assay[i] %in% SummarizedExperiment::assayNames(objs[[i]]), NA)
    if (!any(named) && all(grepl("^[0-9]+$", assay))) assay <- as.integer(assay)
  }
  assay
}

# The native objects in the order of the samples: a list named with the sample names is reordered, other
# names are taken in order with a warning.
.stc_native_order <- function(native, lab) {
  nm <- names(native)
  if (is.null(nm) || identical(nm, lab)) return(native)
  if (setequal(nm, lab)) return(native[lab])
  warning(sprintf("native is named %s, but the samples of result are %s: taken in this order as %s",
                  paste(dQuote(nm, FALSE), collapse = " and "), paste(lab, collapse = " and "),
                  paste(lab, collapse = " and ")), call. = FALSE)
  native
}

# ---------------------------------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------------------------------

# How to draw the pixels of the two objects: for each, the pixels named in `pixels` (the shared ones) and,
# as `other`, the pixels of either object that the other lacks, drawn under every map. Their sf geometry
# (SEraster's colData column "geometry"; needs the sf package), otherwise square tiles at the spatial
# coordinates, with sides the typical distance between neighbouring pixels.
.stc_pixel_layers <- function(objs, pixels) {
  layers <- lapply(objs, .stc_pixel_layer, pixels = pixels)
  tiles <- vapply(layers, function(l) is.null(l$geometry), NA)
  if (tiles[1L] == tiles[2L]) {
    other <- if (!tiles[1L]) {
      list(geometry = c(layers[[1L]]$other$geometry, layers[[2L]]$other$geometry))
    } else {
      list(x = c(layers[[1L]]$other$x, layers[[2L]]$other$x), y = c(layers[[1L]]$other$y, layers[[2L]]$other$y))
    }
    for (i in 1:2) layers[[i]]$other <- other
  }
  layers
}

.stc_pixel_layer <- function(spe, pixels) {
  idx <- match(pixels, colnames(spe))
  other <- setdiff(seq_len(ncol(spe)), idx)
  cd <- SummarizedExperiment::colData(spe)
  if ("geometry" %in% names(cd) && inherits(cd$geometry, "sfc") && requireNamespace("sf", quietly = TRUE)) {
    return(list(geometry = cd$geometry[idx], other = list(geometry = cd$geometry[other])))
  }
  pos <- SpatialExperiment::spatialCoords(spe)[, 1:2, drop = FALSE]
  size <- .stc_pixel_pitch(pos[idx, , drop = FALSE])
  list(x = as.double(pos[idx, 1L]), y = as.double(pos[idx, 2L]), size = size,
       other = list(x = as.double(pos[other, 1L]), y = as.double(pos[other, 2L])))
}

# The median distance from a pixel to its nearest neighbour, from at most 200 evenly spaced pixels (each
# compared with all the others; the session's RNG is not used); 1 for a single pixel.
.stc_pixel_pitch <- function(pos, n = 200L) {
  N <- nrow(pos)
  if (N < 2L) return(1)
  from <- if (N > n) unique(round(seq(1, N, length.out = n))) else seq_len(N)
  nn <- vapply(from, function(i) {
    d2 <- (pos[, 1L] - pos[i, 1L])^2 + (pos[, 2L] - pos[i, 2L])^2
    d2[i] <- Inf
    sqrt(min(d2))
  }, 0)
  nn <- nn[is.finite(nn) & nn > 0]
  if (!length(nn)) 1 else stats::median(nn)
}

# A map of the pixels of a layer (.stc_pixel_layers()) filled by `fill` (values or classes), without a
# scale, over the unshared pixels in light grey. legend = TRUE draws a key for every level of a class scale.
.stc_pixel_map <- function(layer, fill, legend = NA) {
  if (!is.null(layer$geometry)) {
    p <- ggplot2::ggplot(sf::st_sf(fill = fill, geometry = layer$geometry))
    if (length(layer$other$geometry)) {
      p <- p + ggplot2::geom_sf(data = sf::st_sf(geometry = layer$other$geometry), fill = .stc_unshared_fill,
                                colour = NA, show.legend = FALSE)
    }
    p + ggplot2::geom_sf(ggplot2::aes(fill = .data$fill), colour = NA, show.legend = legend) +
      ggplot2::coord_sf(expand = FALSE, clip = "off")
  } else {
    p <- ggplot2::ggplot(data.frame(x = layer$x, y = layer$y, fill = fill), ggplot2::aes(x = .data$x, y = .data$y))
    if (length(layer$other$x)) {
      p <- p + ggplot2::geom_tile(data = data.frame(x = layer$other$x, y = layer$other$y), fill = .stc_unshared_fill,
                                  width = layer$size, height = layer$size, show.legend = FALSE)
    }
    p + ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), width = layer$size, height = layer$size, show.legend = legend) +
      ggplot2::coord_fixed(expand = FALSE, clip = "off")
  }
}

# The gene's values and the coordinates of the cells of a native (not rasterized) object.
.stc_native_values <- function(spe, gene, assay, label) {
  if (!inherits(spe, "SpatialExperiment")) stop(sprintf("native: %s must be a SpatialExperiment object", label))
  if (!gene %in% rownames(spe)) stop(sprintf("gene %s is not in the native object of %s", gene, label))
  an <- SummarizedExperiment::assayNames(spe)
  na <- length(SummarizedExperiment::assays(spe))
  ok <- if (is.character(assay)) assay %in% an else is.numeric(assay) && assay >= 1 && assay <= na
  if (!ok) {
    stop(sprintf("the native object of %s has no assay %s; its assays are: %s", label,
                 if (is.character(assay)) dQuote(assay, FALSE) else assay,
                 if (length(an)) paste(an, collapse = ", ") else paste0(na, " unnamed")))
  }
  value <- as.double(as.vector(SummarizedExperiment::assay(spe, assay)[gene, ]))
  pos <- SpatialExperiment::spatialCoords(spe)
  if (is.null(pos) || ncol(pos) < 2L) stop(sprintf("the native object of %s needs two spatial coordinates", label))
  assay_label <- if (is.character(assay)) assay else if (length(an) >= assay && nzchar(an[assay])) an[assay] else "expression"
  list(value = value, pos = matrix(as.double(pos[, 1:2]), ncol = 2L), assay = assay_label)
}

# The map of one native object: the cells as points coloured by the gene's values, the highest drawn last.
.stc_native_map <- function(nat, label, unit, lims, box, transform, pointSize, baseSize, legend) {
  o <- order(nat$value, na.last = FALSE)
  df <- data.frame(x = nat$pos[o, 1L], y = nat$pos[o, 2L], value = nat$value[o])
  if (is.null(pointSize)) pointSize <- max(0.1, min(2, 60 / sqrt(nrow(df)))) * baseSize / 10
  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$value), size = pointSize, shape = 16, na.rm = TRUE) +
    ggplot2::scale_colour_viridis_c(name = nat$assay, limits = lims, transform = .stc_transform(transform),
                                    breaks = .stc_breaks(range(nat$value, na.rm = TRUE), transform), na.value = .stc_na_fill) +
    ggplot2::coord_fixed(xlim = box[, 1L], ylim = box[, 2L], clip = "off") +
    ggplot2::labs(title = label, subtitle = sprintf("%d %s", nrow(df), unit)) +
    .stc_map_theme(baseSize, legend = legend)
}

# The scatter plot of the values at matched pixels, coloured by class when there are classes, with the
# identity line, the fold-change band, the thresholds and the statistics as the subtitle (the key to the
# lines is the figure's caption).
.stc_scatter_panel <- function(vx, vy, cls, cls_labels, lab, assay, fc, row, transform, baseSize, legend, stats) {
  df <- data.frame(x = vx, y = vy)
  lims <- range(0, vx, vy, na.rm = TRUE)
  tr <- .stc_transform(transform)
  scale <- baseSize / 10
  # the band as curves, so that it is right under any axis transform
  grid <- seq(lims[1L], lims[2L], length.out = 200L)
  band <- data.frame(x = rep(grid, 3L), y = c(grid, grid * 2^fc, grid * 2^-fc),
                     line = rep(c("identity", "upper", "lower"), each = length(grid)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y))
  if (all(c("t1", "t2") %in% names(row)) && !is.na(row$t1) && !is.na(row$t2)) {
    p <- p + ggplot2::geom_vline(xintercept = row$t1, linetype = "dotted", colour = "grey60") +
      ggplot2::geom_hline(yintercept = row$t2, linetype = "dotted", colour = "grey60")
  }
  p <- p + ggplot2::geom_line(data = band, ggplot2::aes(group = .data$line, linetype = .data$line),
                              colour = "grey40", linewidth = 0.4 * scale, na.rm = TRUE) +
    ggplot2::scale_linetype_manual(values = c(identity = "solid", upper = "dashed", lower = "dashed"), guide = "none")
  if (!is.null(cls)) {
    df$class <- cls
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = .data$class), size = 1.2 * scale, alpha = 0.75,
                                 na.rm = TRUE, show.legend = TRUE) +
      ggplot2::scale_colour_manual(name = NULL, values = .stc_class_colours, labels = cls_labels, drop = FALSE)
  } else {
    p <- p + ggplot2::geom_point(size = 1.2 * scale, alpha = 0.6, colour = "grey30", na.rm = TRUE)
  }
  # the statistics are the subtitle: outside the panel, so that they never cover the points
  p + ggplot2::scale_x_continuous(transform = tr, breaks = .stc_breaks(lims, transform), n.breaks = 4L) +
    ggplot2::scale_y_continuous(transform = tr, breaks = .stc_breaks(lims, transform), n.breaks = 4L) +
    ggplot2::coord_fixed(xlim = lims, ylim = lims) +
    ggplot2::labs(x = sprintf("%s in %s", assay[1L], lab[1L]), y = sprintf("%s in %s", assay[2L], lab[2L]),
                  title = "Matched pixels", subtitle = stats) +
    ggplot2::theme_classic(base_size = baseSize) +
    ggplot2::theme(legend.position = if (legend && !is.null(cls)) "right" else "none",
                   legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
                   legend.key.size = ggplot2::unit(1.1 * baseSize, "pt"),
                   plot.title = ggplot2::element_text(size = ggplot2::rel(1), face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85), colour = "grey10", lineheight = 1.1))
}

# The key to the lines of the scatter plot, the caption of the figure (where it has the whole width).
.stc_scatter_caption <- function(fc, thresholds) {
  sprintf("Matched pixels: solid line y = x; dashed: %g-fold band%s", 2^fc, if (thresholds) "; dotted: thresholds" else "")
}

# The theme of the maps: no axes, a compact legend on the right (or none); every size relative to baseSize.
.stc_map_theme <- function(baseSize, legend = TRUE) {
  ggplot2::theme_void(base_size = baseSize) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(1), face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85), colour = "grey30", hjust = 0),
                   legend.position = if (legend) "right" else "none",
                   legend.title = ggplot2::element_text(size = ggplot2::rel(0.8)),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(0.75)),
                   legend.key.width = ggplot2::unit(baseSize, "pt"),
                   legend.key.height = ggplot2::unit(1.4 * baseSize, "pt"),
                   legend.key.spacing.y = ggplot2::unit(0.3 * baseSize, "pt"),
                   plot.margin = ggplot2::margin(4, 4, 4, 4))
}

# The theme of the figure's annotation (its title, subtitle, caption and panel letters, and the position
# of a collected legend), complete so that its sizes are relative to baseSize.
.stc_figure_theme <- function(baseSize, hjust) {
  ggplot2::theme_classic(base_size = baseSize) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3), hjust = hjust),
                   plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.9), colour = "grey30", hjust = hjust),
                   plot.tag = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2)),
                   plot.caption = ggplot2::element_text(size = ggplot2::rel(0.75), colour = "grey30", hjust = 1),
                   legend.position = "bottom")
}

# The ggplot2 transform of the colour scales and the scatter axes, and their breaks: ggplot2's own for the
# identity, and 0, 1, 3, 10, 30, ... up to the largest value for log1p.
.stc_transform <- function(transform) if (transform == "log1p") "log1p" else "identity"

.stc_breaks <- function(lims, transform) {
  if (transform != "log1p") return(ggplot2::waiver())
  b <- c(0, as.vector(outer(c(1, 3), 10^(0:12))))
  b <- b[b <= lims[2L]]
  if (length(b) > 7L) b <- b[seq(1L, length(b), by = 2L)]
  b
}

# The legend labels of the similarity classes: the share of each class among the pixels above the
# threshold, and the number of pixels below it.
.stc_class_labels <- function(cls, lab) {
  kept <- sum(cls != "below")
  share <- function(k) if (kept) sprintf("%.0f%%", 100 * sum(cls == k) / kept) else "-"
  c(similar = sprintf("similar (%s)", share("similar")),
    x = sprintf("higher in %s (%s)", lab[1L], share("x")),
    y = sprintf("higher in %s (%s)", lab[2L], share("y")),
    below = sprintf("below threshold (%d px)", sum(cls == "below")))
}

# The subtitle of the scatter plot: the gene's statistics, two significant digits each, on up to three
# short lines (r and p; the adjusted and the naive p; the permutations, the similarity and the status).
# `note` says why the similarity is NA when the status does not.
.stc_stat_label <- function(row, note = NULL) {
  lines <- character(0)
  last <- character(0)
  if ("correlationCoef" %in% names(row)) {
    p <- row$pValuePermute
    L <- row$nPermutations
    minimum <- !is.na(p) && !is.na(L) && L > 0 && identical(row$stop, "limit") && p <= (1 + 1e-9) / (L + 1)
    lines <- c(sprintf("r = %s, p %s%s", .stc_fmt(row$correlationCoef, 3L), .stc_fmt_p(p), if (minimum) " (minimum)" else ""),
               sprintf("adjusted p %s, naive p %s", .stc_fmt_p(row$pValuePermuteAdj), .stc_fmt_p(row$pValueNaive)))
    if (!is.na(L) && L > 0) last <- sprintf("%s permutations", format(L, big.mark = ","))
  }
  if ("percentSimilarity" %in% names(row)) {
    why <- if (is.na(row$percentSimilarity) && !is.null(note)) sprintf(" (%s)", note) else ""
    last <- c(last, sprintf("similarity = %s%s", .stc_fmt(row$percentSimilarity), why))
  }
  if (!is.null(row$status) && !identical(row$status, "ok")) {
    msg <- if (!is.null(row$message) && !is.na(row$message) && nzchar(row$message)) sprintf(" (%s)", row$message) else ""
    last <- c(last, sprintf("status: %s%s", row$status, msg))
  }
  paste(c(lines, if (length(last)) paste(last, collapse = ", ")), collapse = "\n")
}

# A number to `digits` significant digits, "NA" when missing.
.stc_fmt <- function(v, digits = 2L) {
  if (is.null(v) || length(v) != 1L || is.na(v)) "NA" else format(signif(v, digits))
}

# A p-value with its relation: "= 0.012", "= NA", or "< 2.2e-16" for a p-value that underflowed to 0.
.stc_fmt_p <- function(v) {
  if (!is.null(v) && length(v) == 1L && !is.na(v) && v == 0) "< 2.2e-16" else paste("=", .stc_fmt(v))
}
