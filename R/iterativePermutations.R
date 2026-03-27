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

.get_genes_to_repermute <- function(results_df, alpha) {
  keep <- results_df$pValuePermuteX < alpha & results_df$pValuePermuteY < alpha
  rownames(results_df)[keep]
}


.apply_mht_correction <- function(results_df, total_hypotheses_tested,
                                  method = "BH") {
  results_df$pValuePermuteXAdj <- p.adjust(
    results_df$pValuePermuteX,
    method = method,
    n = total_hypotheses_tested
  )

  results_df$pValuePermuteYAdj <- p.adjust(
    results_df$pValuePermuteY,
    method = method,
    n = total_hypotheses_tested
  )

  results_df
}


# main ----------------------------------------------------------------

spatialCorrelationGeneExp <- function(
    input,
    alpha = 0.05,
    nPermutations = c(100, 1000),
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

  #If lists of deltas to test are not supplied, try 0.1 to 0.9 for both datasets
  #for each gene
  if (is.null(deltaX)){
    deltaX <- rep(list(seq(0.1,0.9,0.1)), length(rownames(source)))
  }

  if (is.null(deltaY)){
    deltaY <- rep(list(seq(0.1,0.9,0.1)), length(rownames(source)))
  }


  #Determine the positions of shared pixels between two rasterized spatial
  #experiments
  source <- input[[1]]
  target <- input[[2]]
  shared_pixels <- intersect(rownames(SpatialExperiment::spatialCoords(source)),
                             rownames(SpatialExperiment::spatialCoords(target)))
  pos <- SpatialExperiment::spatialCoords(source)[shared_pixels,]

  ## if name of assay to use in the SpatialExperiment object is not provided,
  ## use the first assay as a default
  if (is.null(assayName)) {
    assayName <- 1
  }

  genes_all <- rownames(source)
  n_genes <- length(genes_all)


  # creates an index for each gene name
  gene_index <- stats::setNames(seq_along(genes_all), genes_all)

  # storage for latest result per gene
  final_results <- NULL

  # for each iteration in nPermutations, store the num of genes
  # that are tested
  hypotheses_tested_per_iteration <- integer(length(nPermutations))

  # num of genes for next iteration, kept for debugging purposes
  n_repermuted_per_iteration <- integer(length(nPermutations))

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

      hypotheses_tested_per_iteration[k] <- 0L
      n_repermuted_per_iteration[k] <- 0L
      next
    }

    hypotheses_tested_per_iteration[k] <- length(genes_current)

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
      seed = seed + k + 1L # TODO: seed shouldn't change
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
      genes_next <- .get_genes_to_repermute(iter_res, alpha = alpha)
      n_repermuted_per_iteration[k] <- length(genes_next)
      genes_current <- genes_next
    } else {
      n_repermuted_per_iteration[k] <- 0L
    }
  }

  total_hypotheses_tested <- sum(hypotheses_tested_per_iteration)

  final_results <- .apply_mht_correction(
    results_df = final_results,
    total_hypotheses_tested = total_hypotheses_tested,
    method = adjustMethod
  )

  # order rows back to original gene order
  final_results <- final_results[genes_all, , drop = FALSE]
  final_results
}

# testing  ---------------------------------------------------------------

## 100 random dataset   --------------------------------------------------
data("simRanPatternRasts")

## kidney datasets    --------------------------------------------------













