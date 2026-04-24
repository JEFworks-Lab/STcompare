
data(simRanPatternRasts)

t <- Sys.time()

pvalueCorrected_iter <- do.call(rbind, lapply(1:length(simRanPatternRasts), function(i) {
  sapply(1:length(simRanPatternRasts), function(j) {

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
        nThreads = 22,
        verbose = FALSE,
        BPPARAM = BiocParallel::MulticoreParam()
      )
      return(results$pValuePermuteX)
    } else {
      return(test$p.value)
    }
  })
}))

diag(pvalueCorrected_iter) <- NA
tf <- Sys.time() - t
print(tf)

cors <- matrix(NA, nrow = length(simRanPatternRasts), ncol = length(simRanPatternRasts))
corspv <- cors

for (i in 1:length(simRanPatternRasts)) {
  for (j in 1:length(simRanPatternRasts)) {
    vi <- intersect(rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[i]])),
                    rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[j]])))

    g1 <- SummarizedExperiment::assay(simRanPatternRasts[[i]], "pixelval")[1, vi]
    g2 <- SummarizedExperiment::assay(simRanPatternRasts[[j]], "pixelval")[1, vi]

    ctest <- cor.test(g1, g2)
    cors[i, j] <- ctest$estimate
    corspv[i, j] <- ctest$p.value
  }
}

library(reshape2)
cors.df <- reshape2::melt(cors, value.name = "cors")
corspv.df <- reshape2::melt(corspv, value.name = "corspv")
pvalueCorrected.df <- reshape2::melt(pvalueCorrected_iter, value.name = "corspv_corrected")

cors_df <- data.frame(cors.df,
                      corspv = corspv.df$corspv,
                      corspv_corrected = pvalueCorrected.df$corspv_corrected)



cors_df$corspv_corrected[cors_df$corspv_corrected == 0] <- 0.01
cors_df <- na.omit(cors_df)

save(cors_df, file = file.path("inst", "extdata", "simRanPatternResults.RData"))

