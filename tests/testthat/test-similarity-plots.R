# Tests of the fixes to spatialSimilarity(), the plotting functions and documented inputs
# (dev/investigation/09-bug-audit-verified.md, B12, B13, B16-B19, B21, B33, B34): the input checks of
# spatialSimilarity(), the kept pixels of genes below minPixels, the assay recorded in the result and used by the
# plots, savePlots() (the second sample in panel 2, no packages attached), plotCorrelationGeneExp() with negative
# values, NA p-values and compareSpatial() results, and one-row matrices and missing gene names in the
# correlation functions.
# spatialSimilarity() itself is checked against hand-computed values in test-reference-kernel.R and against
# compareSpatial() in test-compareSpatial.R. Labels are explained in helper-fixtures.R.

# Rasterized speKidney A and B with a second assay, lognorm = log10(pixelval + 1).
sp_kidney_ab <- function(geometry = TRUE) {
  rk <- fx_speKidney_raster()
  ab <- list(A = rk$A, B = rk$B)
  lapply(ab, function(s) {
    SummarizedExperiment::assay(s, "lognorm") <- log10(SummarizedExperiment::assay(s, "pixelval") + 1)
    if (!geometry) {
      s <- SpatialExperiment::SpatialExperiment(assays = SummarizedExperiment::assays(s),
                                                spatialCoords = SpatialExperiment::spatialCoords(s))
    }
    s
  })
}

test_that("portable: spatialSimilarity() reports the kept pixels of genes below minPixels and stops on invalid input", {
  px <- paste0("px", 1:20)
  x <- rbind(g1 = c(0, 0, 0, 1, 2, 4, 4, 2, 8, 3, 6, 1, 1, 5, 5, 2, 2, 7, 7, 9), z = rep(0, 20))
  y <- rbind(g1 = c(0, 1, 3, 2, 1, 2, 8, 4, 4, 6, 3, 0, 1, 5, 10, 1, 4, 7, 14, 9), z = rep(0, 20))
  colnames(x) <- colnames(y) <- px
  mk <- function(m) SpatialExperiment::SpatialExperiment(assays = list(counts = m), spatialCoords = cbind(x = seq_len(20), y = 0))
  s <- spatialSimilarity(list(mk(x), mk(y)))
  z <- s$similarityTable[s$similarityTable$gene == "z", ]
  # no pixel is above the thresholds: 0 kept pixels (reported as 1 before, audit B19), so no similarity
  expect_identical(c(z$numPixelInThresh, z$numPixelOutThresh), c(0L, 20L))
  expect_identical(z$pixelIDInThresh[[1]], character(0))
  expect_identical(z$pixelIDOutThresh[[1]], px)
  expect_true(is.na(z$percentSimilarity))
  expect_identical(s$pixelLogTransformation$gene, "g1")
  # the assay is recorded (audit B17)
  expect_identical(s$parameters$assayName, 1)
  expect_identical(spatialSimilarity(list(mk(x), mk(y)), assayName = "counts")$parameters$assayName, "counts")
  # invalid values stop with the gene named (audit B12 and B34)
  x_na <- x
  x_na["g1", 5] <- NA
  expect_error(spatialSimilarity(list(mk(x_na), mk(y))), "missing or infinite values in the first object.*g1")
  y_neg <- y
  y_neg["g1", 2] <- -1
  expect_error(spatialSimilarity(list(mk(x), mk(y_neg)), t1 = 1, t2 = 1),
               "negative values in the second object, for example g1")
  # no shared pixels (NaN without a warning before) and invalid arguments
  y_other <- y
  colnames(y_other) <- paste0("q", 1:20)
  expect_error(spatialSimilarity(list(mk(x), mk(y_other))), "share no pixel")
  expect_error(spatialSimilarity(list(mk(x))), "list of two SpatialExperiment objects")
  expect_error(spatialSimilarity(list(mk(x), mk(y)), t1 = c(1, 2)), "t1 must be NULL or a single number")
  expect_error(spatialSimilarity(list(mk(x), mk(y)), minQuantile = 1.5), "minQuantile")
  expect_error(spatialSimilarity(list(mk(x), mk(y)), foldChange = -1), "foldChange")
})

test_that("portable: the plots use the assay of spatialSimilarity(); savePlots() shows both samples and attaches nothing", {
  skip_if_not_installed("patchwork")
  ab <- sp_kidney_ab()
  pixels <- intersect(colnames(ab$A), colnames(ab$B))
  value <- function(s, assay, cols = pixels) unname(SummarizedExperiment::assay(s, assay)["Gene", cols])
  s <- spatialSimilarity(ab, assayName = "lognorm")
  # linearRegression() and pixelClass() default to the assay that was classified (audit B17)
  expect_identical(linearRegression(s, "Gene")$data$x, value(ab$A, "lognorm"))
  expect_identical(linearRegression(s, "Gene", assayName = "pixelval")$data$y, value(ab$B, "pixelval"))
  expect_s3_class(pixelClass(s, "Gene"), "ggplot")
  # savePlots(): no package is attached (audit B18); the expression panels show the classified assay of each
  # sample (B17), with the sample names as titles
  attached <- search()
  p <- savePlots("Gene", s, ab)$Gene
  expect_identical(search(), attached)
  expect_s3_class(p, "patchwork")
  panels <- p$patches$plots
  expect_length(panels, 3L)
  for (i in 1:2) {
    expect_equal(unname(panels[[i]]$layers[[1]]$data$fill), value(ab[[i]], "lognorm", colnames(ab[[i]])))
    expect_identical(panels[[i]]$labels$title, names(ab)[i])
  }
  # without pixel geometries: panel 2 is the second sample (it repeated the first, audit B16); a PDF is written
  # when filePath is given
  ng <- sp_kidney_ab(geometry = FALSE)
  out <- tempfile("savePlots-")
  dir.create(out)
  # (drawing the scatter panel warns about the points beyond the 95% quantiles, where linearRegression() cuts
  # the axes)
  pn <- suppressWarnings(savePlots("Gene", spatialSimilarity(ng), ng, filePath = out))$Gene$patches$plots
  for (i in 1:2) {
    expect_identical(pn[[i]]$data$value, value(ng[[i]], "pixelval"))
    expect_identical(pn[[i]]$labels$title, names(ng)[i])
  }
  expect_true(file.exists(file.path(out, "Gene.pdf")))
  unlink(out, recursive = TRUE)
})

test_that("portable: plotCorrelationGeneExp() draws negative values and shows a rounded p-value, NA if either is NA", {
  ab <- sp_kidney_ab()
  N <- length(intersect(colnames(ab$A), colnames(ab$B)))
  # centred values (half of them negative) are all drawn; values below 0 were dropped by the axis limits
  ab <- lapply(ab, function(s) {
    m <- SummarizedExperiment::assay(s, "pixelval")
    SummarizedExperiment::assay(s, "centred") <- m - mean(m)
    s
  })
  tab <- data.frame(correlationCoef = -0.94721, pValuePermuteX = 1 / 3, pValuePermuteY = 0.002, row.names = "Gene")
  p <- plotCorrelationGeneExp(ab, tab, "Gene", assayName = "centred")
  d <- ggplot2::layer_data(p)
  expect_identical(sum(stats::complete.cases(d[, c("x", "y")])), N)
  expect_lt(min(d$x), 0)
  expect_identical(p$labels$title, "Gene  r =  -0.947  p_E =  0.333")
  expect_identical(p$labels$x, "centred in A")
  # the greater of the two p-values, NA if either is NA
  for (col in c("pValuePermuteX", "pValuePermuteY")) {
    t2 <- tab
    t2[[col]] <- NA
    expect_match(plotCorrelationGeneExp(ab, t2, "Gene")$labels$title, "p_E =  NA$")
  }
  # a compareSpatial() result: r and padj
  res <- data.frame(gene = "Gene", r = 0.5, p = 0.001, padj = 0.0123456, row.names = "Gene")
  expect_identical(plotCorrelationGeneExp(ab, res, "Gene")$labels$title, "Gene  r =  0.5  padj =  0.0123")
  expect_error(plotCorrelationGeneExp(ab, tab, "Other"), "gene Other is not a row of spatialCorrelation")
})

test_that("portable: documented inputs of the correlation functions are accepted (audit B13, B21)", {
  local_default_rng()
  set.seed(1)
  x <- stats::rnorm(50)
  y <- stats::rnorm(50)
  pos <- cbind(stats::runif(50), stats::runif(50))
  # X and Y as vectors or as matrices with one row
  expect_identical(spatialCorrelation(matrix(x, 1), matrix(y, 1), pos, nPermutations = 5, deltaX = 0.5, deltaY = 0.5),
                   spatialCorrelation(x, y, pos, nPermutations = 5, deltaX = 0.5, deltaY = 0.5))
  # the within-sample function needs gene names
  spe <- SpatialExperiment::SpatialExperiment(assays = list(v = rbind(x, y)), spatialCoords = pos)
  rownames(spe) <- NULL
  expect_error(spatialCorrelationGeneExpWithinSample(spe, nPermutations = 5, verbose = FALSE),
               "input must have row names")
})
