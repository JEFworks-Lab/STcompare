# TODO: new block that implements iterative permutations


# iterative permutations with BH default MHT ###################################

# helper ---------------------------------------------------------------
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
      deltaX = deltaX,
      deltaY = deltaY,
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

.get_genes_to_repermute <- function(results_df, alpha, nPermutes) {
  t <- alpha / nPermutes
  keep <- results_df$pValuePermuteX < t & results_df$pValuePermuteY < t
  rownames(results_df)[keep]
}


# main ----------------------------------------------------------------

spatialCorrelationGeneExp_test <- function(
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

  #If lists of deltas to test are not supplied, try 0.1 to 0.9 for both datasets
  #for each gene
  if (is.null(deltaX)){
    deltaX <- rep(list(seq(0.1,0.9,0.1)), length(rownames(source)))
  }

  if (is.null(deltaY)){
    deltaY <- rep(list(seq(0.1,0.9,0.1)), length(rownames(source)))
  }


  shared_pixels <- intersect(rownames(SpatialExperiment::spatialCoords(source)),
                             rownames(SpatialExperiment::spatialCoords(target)))
  pos <- SpatialExperiment::spatialCoords(source)[shared_pixels,]

  # if name of assay to use in the SpatialExperiment object is not provided,
  # use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }

  genes_all <- rownames(source)
  n_genes <- length(genes_all)


  # creates an index for each gene name
  gene_index <- stats::setNames(seq_along(genes_all), genes_all)

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
      genes_current <- .get_genes_to_repermute(iter_res, alpha = alpha, nPermutes = nPermutations[k])
    }

  }

  # mht correct for pValuePermuteX and pValuePermuteY seperately
  final_results$pValuePermuteX <- stats::p.adjust(final_results$pValuePermuteX, method = adjustMethod)
  final_results$pValuePermuteY <- stats::p.adjust(final_results$pValuePermuteY, method = adjustMethod)

  # order rows back to original gene order
  final_results <- final_results[genes_all, , drop = FALSE]
  final_results
}

# testing  ---------------------------------------------------------------

## 100 random dataset   --------------------------------------------------
# code from file:
# system.file("extdata", "simRanPatternResults.RData", package = "STcompare")

data("simRanPatternRasts")

set.seed(0)
t <- Sys.time()

testOneIteration <- isTRUE(get0("testOneIteration", ifnotfound = FALSE))
nSimRanPatternRasts <- length(simRanPatternRasts)
simRanPatternRastPairs <- if (testOneIteration) {
  data.frame(i = 1L, j = min(2L, nSimRanPatternRasts))
} else {
  expand.grid(i = seq_len(nSimRanPatternRasts), j = seq_len(nSimRanPatternRasts))
}
pvalueCorrected <- matrix(NA_real_, nrow = nSimRanPatternRasts, ncol = nSimRanPatternRasts)

for (idx in seq_len(nrow(simRanPatternRastPairs))) {
  i <- simRanPatternRastPairs$i[idx]
  j <- simRanPatternRastPairs$j[idx]

  print(paste0(i,":", j))
  # shared pixels between rasters i and j
  vi <- intersect(rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[i]])),
                  rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[j]])))

  g1 <- SummarizedExperiment::assay(simRanPatternRasts[[i]], "pixelval")[1, vi]
  g2 <- SummarizedExperiment::assay(simRanPatternRasts[[j]], "pixelval")[1, vi]

  test <- cor.test(g1, g2)

  if (i != j) {
    rastGexpListAB <- list(i = simRanPatternRasts[[i]],
                           j = simRanPatternRasts[[j]])
    results <- spatialCorrelationGeneExp_test(
      rastGexpListAB,
      nPermutations = c(1e2, 1e3),
      # nThreads = parallel::detectCores() - 4, # should be 20 on franklin 
      nThreads = 1, 
      verbose = FALSE,
      BPPARAM = BiocParallel::MulticoreParam()
    )
    pvalueCorrected[i, j] <- results$pValuePermuteX
  } else {
    pvalueCorrected[i, j] <- test$p.value
  }
}

diag(pvalueCorrected) <- NA
tf <- Sys.time() - t
print(tf)

cors <- matrix(NA, nrow = length(simRanPatternRasts), ncol = length(simRanPatternRasts))
corspv <- cors

for (idx in seq_len(nrow(simRanPatternRastPairs))) {
  i <- simRanPatternRastPairs$i[idx]
  j <- simRanPatternRastPairs$j[idx]
  vi <- intersect(rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[i]])),
                  rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[j]])))

  g1 <- SummarizedExperiment::assay(simRanPatternRasts[[i]], "pixelval")[1, vi]
  g2 <- SummarizedExperiment::assay(simRanPatternRasts[[j]], "pixelval")[1, vi]

  ctest <- cor.test(g1, g2)
  cors[i, j] <- ctest$estimate
  corspv[i, j] <- ctest$p.value
}

library(ggplot2)
library(reshape2)

cors.df <- reshape2::melt(cors, value.name = "cors")
corspv.df <- reshape2::melt(corspv, value.name = "corspv")
pvalueCorrected.df <- reshape2::melt(pvalueCorrected, value.name = "corspv_corrected")

cors_df <- data.frame(cors.df,
                      corspv = corspv.df$corspv,
                      corspv_corrected = pvalueCorrected.df$corspv_corrected)



cors_df$corspv_corrected[cors_df$corspv_corrected == 0] <- 0.01

test_df <- na.omit(cors_df)

corrected <- ggplot(test_df) +
  geom_point(aes(x =cors, y = -log10(corspv)), alpha = 0.1, size = 0.5, color = "blue") +
  geom_point(aes(x =cors, y = -log10(corspv_corrected)), alpha = 0.1, size = 0.5, color = "green") +
  theme_classic() +
  labs(x = "Correlation", y = "-log10(p-value)",
       title = "Naïve vs Spatially Corrected p-values") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "black")  +
  ylim(min(-log10(test_df$corspv), na.rm = TRUE), max(-log10(test_df$corspv), na.rm = TRUE)) +
  labs(x = "Correlation" , y = "-log10(p-value)")
corrected

# save(cors_df, corrected,
#      file = file.path("inst", "extdata", "simRanPatternResultsTestingCorrection.RData"))



## kidney datasets    --------------------------------------------------

# TODO: this should be on the rasterized and aligned?
# TODO: before or after figuring out STalign?









