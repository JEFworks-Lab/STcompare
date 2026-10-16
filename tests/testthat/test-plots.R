# Tests of the standard figures of compareSpatial() results (R/plotGene.R): plotGene() and plotOverview()
# return patchwork figures with the expected panels, every panel builds, the similarity classes agree with
# the result, and the inputs are checked. Labels are explained in helper-fixtures.R.

# The ggplot2 plots of a patchwork figure, in panel order: a patchwork holds its other plots in
# patches$plots and is itself the last plot; nested patchworks (the rows of a two-row figure) are walked, and
# spacers are left out.
pl_panels <- function(fig) {
  if (!inherits(fig, "patchwork")) return(if (inherits(fig, c("spacer", "plot_filler"))) list() else list(fig))
  last <- fig
  class(last) <- setdiff(class(last), "patchwork")
  last$patches <- NULL
  c(do.call(c, lapply(fig$patches$plots, pl_panels)), pl_panels(last))
}

# Build every panel of a figure (errors in a scale, a layer or a coord surface here).
pl_build <- function(fig) {
  for (p in pl_panels(fig)) ggplot2::ggplot_build(p)
  invisible(fig)
}

# The rasterized kidney samples A and B and a compareSpatial() result on them, once per session.
pl_kidney <- function() fx_memo("plots_kidney", {
  rk <- fx_speKidney_raster()
  list(x = rk$A, y = rk$B,
       res = compareSpatial(list(A = rk$A, B = rk$B), nPermutations = 99, nThreads = fx_threads,
                            progress = FALSE, verbose = FALSE))
})

# A rasterized object without pixel geometries (coordinates only), with the same assay.
pl_strip_geometry <- function(spe) {
  SpatialExperiment::SpatialExperiment(assays = SummarizedExperiment::assays(spe),
                                       spatialCoords = SpatialExperiment::spatialCoords(spe))
}

test_that("portable: plotGene() draws the four panels, or six with native maps, and every panel builds", {
  k <- pl_kidney()
  fig <- plotGene(k$res, "Gene", k$x, k$y)
  expect_s3_class(fig, "patchwork")
  panels <- pl_panels(fig)
  expect_length(panels, 4L)
  expect_identical(vapply(panels, function(p) p$labels$title, ""), c("A", "B", "Similarity classes", "Matched pixels"))
  pl_build(fig)
  # the figure is titled with the gene and the samples
  expect_identical(fig$patches$annotation$title, "Gene")
  expect_match(fig$patches$annotation$subtitle, "^A vs B, 273 shared pixels$")
  # the maps use the pixel geometry (sf) and the scatter the values of the shared pixels
  expect_s3_class(panels[[1]]$data, "sf")
  # the pixels that one sample lacks are drawn under the shared pixels, in every map
  unshared <- length(union(colnames(k$x), colnames(k$y))) - 273L
  expect_gt(unshared, 0L)
  expect_identical(nrow(panels[[1]]$layers[[1]]$data), unshared)
  expect_identical(nrow(panels[[3]]$layers[[1]]$data), unshared)
  d <- ggplot2::layer_data(panels[[4]], 4L)  # the layers: two thresholds, the band, the points
  pixels <- intersect(colnames(k$x), colnames(k$y))
  expect_equal(d$x, unname(SummarizedExperiment::assay(k$x)["Gene", pixels]))
  expect_equal(d$y, unname(SummarizedExperiment::assay(k$y)["Gene", pixels]))
  # with native maps: six panels, the native ones first, in a two-row design
  data(speKidney, package = "STcompare")
  fig6 <- plotGene(k$res, "Gene", list(k$x, k$y), native = speKidney[c("A", "B")], transform = "log1p")
  panels6 <- pl_panels(fig6)
  expect_length(panels6, 6L)
  expect_identical(panels6[[1]]$labels$subtitle, "1229 cells")
  expect_identical(panels6[[2]]$labels$subtitle, "1242 cells")
  rows <- fig6$patches$plots  # two rows, each a patchwork: the four maps, then the other two panels, taller
  expect_length(rows, 2L)
  expect_length(pl_panels(rows[[1]]), 4L)
  expect_length(pl_panels(rows[[2]]), 2L)
  expect_equal(fig6$patches$layout$heights, c(1, 1.25))
  pl_build(fig6)
  # the native maps name their points
  expect_identical(pl_panels(plotGene(k$res, "Gene", k$x, k$y, native = speKidney[c("A", "B")],
                                      nativeUnit = c("cells", "spots")))[[2]]$labels$subtitle, "1242 spots")
  # the panels can be chosen, and their order is fixed
  fig2 <- plotGene(k$res, "Gene", k$x, k$y, panels = c("scatter", "raster"))
  expect_identical(vapply(pl_panels(fig2), function(p) p$labels$title, ""), c("A", "B", "Matched pixels"))
  # a native panel asked for without native objects is dropped silently; nothing left is an error
  expect_length(pl_panels(plotGene(k$res, "Gene", k$x, k$y, panels = c("native", "scatter"))), 1L)
  expect_error(plotGene(k$res, "Gene", k$x, k$y, panels = "native"), "no panel to draw")
})

test_that("portable: plotGene() classifies the pixels as the result does and labels the statistics", {
  k <- pl_kidney()
  fig <- plotGene(k$res, "Gene", k$x, k$y)
  panels <- pl_panels(fig)
  row <- as.data.frame(k$res)["Gene", ]
  # the classes of the similarity map: shares among the pixels above the threshold = the result's columns
  cls <- panels[[3]]$data$fill
  kept <- cls != "below"
  expect_equal(mean(cls[kept] == "similar"), row$percentSimilarity)
  expect_equal(mean(cls[kept] == "x"), row$percentDissimilarityX)
  expect_equal(mean(cls[kept] == "y"), row$percentDissimilarityY)
  expect_identical(sum(kept), row$numPixelInThresh)
  # the legend labels carry the shares
  labels <- panels[[3]]$scales$get_scales("fill")$labels
  expect_identical(labels[["similar"]], sprintf("similar (%.0f%%)", 100 * row$percentSimilarity))
  expect_match(labels[["x"]], "^higher in A \\(")
  expect_match(labels[["y"]], "^higher in B \\(")
  expect_identical(panels[[3]]$labels$subtitle, sprintf("similarity %s", format(signif(row$percentSimilarity, 2))))
  # every class has a key, present in the figure or not (A vs C: no similar pixel, none higher in A)
  rk <- fx_speKidney_raster()
  ac <- compareSpatial(list(A = rk$A, C = rk$C), tests = "similarity", verbose = FALSE)
  keys <- ggplot2::get_guide_data(pl_panels(plotGene(ac, "Gene", rk$A, rk$C))[[3]], "fill")
  expect_identical(keys$fill, unname(STcompare:::.stc_class_colours))
  expect_identical(keys$.label[1:3], c("similar (0%)", "higher in A (0%)", "higher in C (100%)"))
  # the subtitle of the scatter plot: every statistic of the row, the p-value marked at the minimum
  lab <- panels[[4]]$labels$subtitle
  expect_match(lab, sprintf("^r = %s, p = %s%s\n", format(signif(row$correlationCoef, 3)),
                            format(signif(row$pValuePermute, 2)),
                            if (row$stop == "limit" && row$pValuePermute <= 1 / (row$nPermutations + 1)) " \\(minimum\\)" else ""))
  expect_match(lab, sprintf("\nadjusted p = %s, naive p = %s\n", format(signif(row$pValuePermuteAdj, 2)),
                            format(signif(row$pValueNaive, 2))))
  expect_match(lab, sprintf("\n%d permutations, similarity = %s$", row$nPermutations,
                            format(signif(row$percentSimilarity, 2))))
  # the key to the lines is the figure's caption, so that nothing covers the points and it never clips
  expect_null(panels[[4]]$labels$caption)
  expect_identical(fig$patches$annotation$caption, "Matched pixels: solid line y = x; dashed: 2-fold band; dotted: thresholds")
  # a p-value that underflowed to 0 is printed as < 2.2e-16
  expect_identical(STcompare:::.stc_fmt_p(0), "< 2.2e-16")
  expect_identical(STcompare:::.stc_fmt_p(0.0123), "= 0.012")
  expect_identical(STcompare:::.stc_fmt_p(NA), "= NA")
  # a gene whose similarity is NA because too few pixels pass the thresholds (minPixels): no similarity
  # panel, and the scatter says why
  sparse <- function(spe) {
    v <- SummarizedExperiment::assay(spe, "pixelval")["Gene", ]
    v[seq_along(v) %% 25 != 0L] <- 0
    SummarizedExperiment::assay(spe, "pixelval")["Gene", ] <- v
    spe
  }
  rare <- sparse(k$x)
  rare_y <- sparse(k$y)
  res_rare <- compareSpatial(list(A = rare, B = rare_y), tests = "similarity", verbose = FALSE)
  expect_true(is.na(res_rare$percentSimilarity))
  fig_rare <- plotGene(res_rare, "Gene", rare, rare_y)
  expect_length(pl_panels(fig_rare), 3L)
  expect_match(pl_panels(fig_rare)[[3]]$labels$subtitle,
               sprintf("^similarity = NA \\(%d of 273 pixels above threshold, fewer than minPixels = 0.1\\)$",
                       res_rare$numPixelInThresh))
  # a gene with negative values: no similarity, with the reason
  neg <- k$x
  SummarizedExperiment::assay(neg, "pixelval")["Gene", ] <- SummarizedExperiment::assay(neg, "pixelval")["Gene", ] - 5
  res_neg <- suppressWarnings(compareSpatial(list(A = neg, B = k$y), tests = "similarity", verbose = FALSE))
  expect_identical(pl_panels(plotGene(res_neg, "Gene", neg, k$y))[[3]]$labels$subtitle, "similarity = NA (negative values)")
  expect_error(plotGene(res_neg, "Gene", neg, k$y, transform = "log1p"), "needs values that are not negative")
})

test_that("portable: plotGene() works without geometry, with sparse assays, and for one test alone", {
  k <- pl_kidney()
  # pixels as tiles when the objects have no geometry
  x0 <- pl_strip_geometry(k$x)
  y0 <- pl_strip_geometry(k$y)
  fig <- plotGene(k$res, "Gene", x0, y0)
  panels <- pl_panels(fig)
  expect_length(panels, 4L)
  expect_false(inherits(panels[[1]]$data, "sf"))
  expect_s3_class(panels[[1]]$layers[[2]]$geom, "GeomTile")  # the unshared pixels first, in grey, then the shared
  expect_identical(panels[[1]]$layers[[1]]$aes_params$fill, STcompare:::.stc_unshared_fill)
  expect_equal(panels[[1]]$data$fill, unname(SummarizedExperiment::assay(x0)["Gene", intersect(colnames(x0), colnames(y0))]))
  pl_build(fig)
  # the tiles have the size of the pixel spacing (SEraster hexagons of resolution 0.2: 0.2 apart)
  expect_equal(panels[[1]]$layers[[2]]$aes_params$width, 0.2, tolerance = 1e-6)
  # sparse assays
  skip_if_not_installed("Matrix")
  sparse <- function(spe) {
    SummarizedExperiment::assay(spe, "pixelval") <- Matrix::Matrix(SummarizedExperiment::assay(spe, "pixelval"), sparse = TRUE)
    spe
  }
  xs <- sparse(k$x)
  ys <- sparse(k$y)
  expect_s4_class(SummarizedExperiment::assay(xs), "dgCMatrix")
  data(speKidney, package = "STcompare")
  nat <- lapply(speKidney[c("A", "B")], function(s) {
    SummarizedExperiment::assay(s, "counts") <- Matrix::Matrix(SummarizedExperiment::assay(s, "counts"), sparse = TRUE)
    s
  })
  figs <- plotGene(k$res, "Gene", xs, ys, native = nat)
  expect_length(pl_panels(figs), 6L)
  pl_build(figs)
  expect_equal(pl_panels(figs)[[3]]$data$fill, pl_panels(fig)[[1]]$data$fill)
  # similarity alone: no correlation lines in the subtitle; correlation alone: no similarity map, plain scatter
  sim <- compareSpatial(list(A = k$x, B = k$y), tests = "similarity", verbose = FALSE)
  fig_sim <- plotGene(sim, "Gene", k$x, k$y)
  expect_length(pl_panels(fig_sim), 4L)
  lab <- pl_panels(fig_sim)[[4]]$labels$subtitle
  expect_match(lab, "^similarity = ")
  expect_no_match(lab, "r = ")
  cor_only <- compareSpatial(list(A = k$x, B = k$y), tests = "correlation", nPermutations = 9, nThreads = fx_threads,
                             progress = FALSE, verbose = FALSE)
  fig_cor <- plotGene(cor_only, "Gene", k$x, k$y)
  panels <- pl_panels(fig_cor)
  expect_length(panels, 3L)
  expect_identical(panels[[3]]$labels$title, "Matched pixels")
  expect_length(panels[[3]]$layers, 2L)  # the band, the points: no thresholds, no classes
  expect_no_match(panels[[3]]$labels$subtitle, "similarity")
  expect_identical(fig_cor$patches$annotation$caption, "Matched pixels: solid line y = x; dashed: 2-fold band")
  pl_build(fig_cor)
  # a failed gene (missing values) still plots, with NA statistics and its status
  yna <- k$y
  SummarizedExperiment::assay(yna, "pixelval")["Gene", 3L] <- NA
  res_na <- suppressWarnings(compareSpatial(list(A = k$x, B = yna), nPermutations = 9, nThreads = fx_threads,
                                            progress = FALSE, verbose = FALSE))
  expect_identical(res_na$status, "failed")
  fig_na <- plotGene(res_na, "Gene", k$x, yna)
  panels <- pl_panels(fig_na)
  expect_length(panels, 3L)  # no similarity map: the thresholds are NA
  expect_identical(panels[[3]]$labels$subtitle,
                   "r = NA, p = NA\nadjusted p = NA, naive p = NA\nsimilarity = NA, status: failed (missing or infinite values in B)")
  pl_build(fig_na)
  # objects whose assays have no names: the result records the index, which is used again
  ua <- lapply(list(A = k$x, B = k$y), function(s) { names(SummarizedExperiment::assays(s)) <- NULL; s })
  rua <- compareSpatial(ua, tests = "similarity", verbose = FALSE)
  expect_identical(attr(rua, "params")$assay, c("1", "1"))
  expect_length(pl_panels(plotGene(rua, "Gene", ua$A, ua$B)), 4L)
})

test_that("portable: plotGene() scales every size from baseSize", {
  k <- pl_kidney()
  data(speKidney, package = "STcompare")
  sizes <- function(baseSize) {
    fig <- plotGene(k$res, "Gene", k$x, k$y, native = speKidney[c("A", "B")], baseSize = baseSize)
    panels <- pl_panels(fig)
    th <- lapply(panels, function(p) ggplot2::ggplot_build(p)$plot$theme)
    c(map_title = ggplot2::calc_element("plot.title", th[[1]])$size,
      map_legend = ggplot2::calc_element("legend.text", th[[2]])$size,
      scatter_axis = ggplot2::calc_element("axis.text", th[[6]])$size,
      scatter_subtitle = ggplot2::calc_element("plot.subtitle", th[[6]])$size,
      figure_title = ggplot2::calc_element("plot.title", fig$patches$annotation$theme)$size,
      tag = ggplot2::calc_element("plot.tag", fig$patches$annotation$theme)$size,
      point = panels[[6]]$layers[[4]]$aes_params$size,
      native_point = panels[[1]]$layers[[1]]$aes_params$size)
  }
  s10 <- sizes(10)
  s7 <- sizes(7)
  expect_equal(s10[["map_title"]], 10)
  expect_equal(s10[["figure_title"]], 13)
  expect_equal(s10[["tag"]], 12)
  expect_equal(unname(s7 / s10), rep(0.7, length(s10)), tolerance = 1e-6)
  # patchwork's & operator rescales the figure too (the sizes are relative)
  fig <- plotGene(k$res, "Gene", k$x, k$y, baseSize = 10) & ggplot2::theme(text = ggplot2::element_text(size = 5))
  expect_equal(ggplot2::calc_element("plot.title", ggplot2::ggplot_build(pl_panels(fig)[[1]])$plot$theme)$size, 5)
  expect_equal(ggplot2::calc_element("plot.tag", fig$patches$annotation$theme)$size, 6)
  expect_error(plotGene(k$res, "Gene", k$x, k$y, baseSize = 0), "baseSize must be a positive number")
  expect_error(plotOverview(k$res, baseSize = "big"), "baseSize must be a positive number")
})

test_that("portable: plotGene() checks its inputs", {
  k <- pl_kidney()
  expect_error(plotGene(k$res, "Other", k$x, k$y), "gene Other is not in result; its genes are Gene")
  expect_error(plotGene(k$res, c("a", "b"), k$x, k$y), "gene must be one gene name")
  expect_error(plotGene(as.data.frame(k$res)[, 1:2], "Gene", k$x), "give the two rasterized SpatialExperiment")
  other <- k$y
  rownames(other) <- "Z"
  expect_error(plotGene(k$res, "Gene", k$x, other), "gene Gene is not in B")
  expect_error(plotGene(k$res, "Gene", k$x, list()), "B must be a SpatialExperiment object")
  expect_error(plotGene(k$res, "Gene", list(k$x, k$y), k$y), "give x and y as two SpatialExperiment objects, or a list")
  expect_error(compareSpatial(list(k$x, k$y), k$y), "give x and y as two SpatialExperiment objects, or a list")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, nativeUnit = 1), "nativeUnit must be")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, assay = "nope"), "has no assay \"nope\"")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, panels = "maps"), "panels must be some of")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, transform = "log10"), "should be one of")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, pointSize = -1), "pointSize must be")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, native = k$x), "native must be NULL or a list of two")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, native = list(k$x, other)), "gene Gene is not in the native object of B")
  expect_error(plotGene(k$res, "Gene", k$x, k$y, native = list(k$x, k$y), nativeAssay = "counts"), "native object of A has no assay")
  # objects that are not those of the result: pixels with other coordinates, or another number of them
  moved <- k$y
  SpatialExperiment::spatialCoords(moved) <- SpatialExperiment::spatialCoords(moved) + 0.1
  expect_error(plotGene(k$res, "Gene", k$x, moved), "different coordinates in the two samples")
  pixels <- intersect(colnames(k$x), colnames(k$y))
  expect_warning(plotGene(k$res, "Gene", k$x[, pixels[1:200]], k$y),
                 "share 200 pixels, but the result was computed on 273.*classes are recomputed on the pixels given")
  # x and y in the other order: the gene's thresholds give it away
  expect_warning(plotGene(k$res, "Gene", k$y, k$x), "x and y appear to be swapped")
  # a named native list is put in the order of the samples; other names are taken in order, with a warning
  data(speKidney, package = "STcompare")
  fig <- plotGene(k$res, "Gene", k$x, k$y, native = speKidney[c("B", "A")])
  expect_identical(pl_panels(fig)[[1]]$labels$subtitle, "1229 cells")  # A
  expect_warning(fig <- plotGene(k$res, "Gene", k$x, k$y, native = list(k1 = speKidney$B, k2 = speKidney$A)),
                 "native is named \"k1\" and \"k2\", but the samples of result are A and B: taken in this order as A and B")
  expect_identical(pl_panels(fig)[[1]]$labels$subtitle, "1242 cells")  # B, as given
})

test_that("portable: plotOverview() draws the panels the result has, with labels, and checks its inputs", {
  # 12 simulated genes on one grid: 10 independent pairs and 2 shared patterns
  data(simRanPatternRasts, package = "STcompare")
  px <- Reduce(intersect, lapply(simRanPatternRasts[1:24], colnames))
  values <- function(k) t(vapply(simRanPatternRasts[k], function(s) SummarizedExperiment::assay(s)[1, px], numeric(length(px))))
  make <- function(m) {
    rownames(m) <- paste0("g", seq_len(nrow(m)))
    SpatialExperiment::SpatialExperiment(assays = list(v = m),
                                         spatialCoords = SpatialExperiment::spatialCoords(simRanPatternRasts[[1]])[px, ])
  }
  mx <- values(1:12)
  my <- values(13:24)
  my[1:2, ] <- mx[1:2, ]
  res <- fx_memo("plots_overview", compareSpatial(list(S1 = make(mx), S2 = make(my)), nPermutations = 99,
                                                  nThreads = fx_threads, progress = FALSE, verbose = FALSE))
  fig <- plotOverview(res, label = c("g1", "g2"))
  expect_s3_class(fig, "patchwork")
  panels <- pl_panels(fig)
  expect_length(panels, 3L)
  expect_identical(vapply(panels, function(p) p$labels$title, ""),
                   c("Spatial correlation", "Spatial similarity", "Correlation and similarity"))
  pl_build(fig)
  expect_identical(fig$patches$annotation$title, "S1 vs S2")
  expect_match(fig$patches$annotation$subtitle, sprintf("^12 genes on %d shared pixels$", length(px)))
  # the significance classes count the genes, and the labelled genes are drawn
  df <- as.data.frame(res)
  sig <- df$pValuePermuteAdj < 0.05
  labels <- panels[[1]]$scales$get_scales("fill")$labels
  expect_identical(labels[["positive"]], sprintf("r > 0 (%d)", sum(sig & df$correlationCoef > 0)))
  expect_identical(panels[[1]]$scales$get_scales("fill")$name, "adjusted p < 0.05")
  expect_identical(labels[["ns"]], sprintf("not significant (%d)", sum(!sig)))
  # every class has a key, present or not
  expect_identical(ggplot2::get_guide_data(panels[[1]], "fill")$fill, unname(STcompare:::.stc_sig_colours))
  # the histogram spans 0 to 1 whatever the data
  expect_equal(ggplot2::layer_scales(panels[[2]])$x$get_limits(), c(0, 1))
  expect_equal(range(ggplot2::layer_data(panels[[2]], 1L)$xmin), c(0, 0.95))
  expect_setequal(ggplot2::layer_data(panels[[1]], 4L)$label, c("g1", "g2"))
  expect_identical(ggplot2::layer_data(panels[[1]], 1L)$yintercept, -log10(0.05))
  # alpha changes the threshold line and the labels
  fig2 <- plotOverview(res, alpha = 0.01)
  expect_identical(ggplot2::layer_data(pl_panels(fig2)[[1]], 1L)$yintercept, -log10(0.01))
  expect_identical(pl_panels(fig2)[[1]]$scales$get_scales("fill")$name, "adjusted p < 0.01")
  # one test alone: the panels of the other are left out
  only_cor <- plotOverview(res[, setdiff(names(res), c("percentSimilarity", "percentDissimilarityX", "percentDissimilarityY",
                                                      "numPixelInThresh", "t1", "t2"))])
  expect_length(pl_panels(only_cor), 1L)
  only_sim <- plotOverview(res[, c("gene", "nPixels", "percentSimilarity", "status", "message")])
  expect_length(pl_panels(only_sim), 1L)
  expect_identical(pl_panels(only_sim)[[1]]$labels$title, "Spatial similarity")
  pl_build(only_sim)
  # inputs
  expect_error(plotOverview(res[, c("gene", "status")]), "neither the correlation nor the similarity columns")
  expect_error(plotOverview(res, alpha = 0), "alpha must be")
  expect_warning(plotOverview(res, label = c("g1", "nope")), "1 of the genes to label are not in result: nope")
})
