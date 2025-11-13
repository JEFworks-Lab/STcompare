
set.seed(0)
t <- Sys.time()

pvalueCorrected <- do.call(rbind, lapply(1:length(simRanPatternRasts), function(i) {
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
      results <- spatialCorrelationGeneExp(
        rastGexpListAB,
        nPermutations = 100,
        nThreads = 5,
        verbose = FALSE,
        BPPARAM = BiocParallel::MulticoreParam()
      )
      return(results$pValuePermuteX)
    } else {
      return(test$p.value)
    }
  })
}))

diag(pvalueCorrected) <- NA
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
pvalueCorrected.df <- reshape2::melt(pvalueCorrected, value.name = "corspv_corrected")

cors_df <- data.frame(cors.df,
                      corspv = corspv.df$corspv,
                      corspv_corrected = pvalueCorrected.df$corspv_corrected)



cors_df$corspv_corrected[cors_df$corspv_corrected == 0] <- 0.01


# TODO: remove below

# corrected <- ggplot(cors_df) +
#   geom_point(aes(x = cors, y = -log10(corspv)),
#              alpha = 0.1, size = 0.5, color = "blue") +
#   geom_point(aes(x = cors, y = -log10(corspv_corrected)),
#              alpha = 0.1, size = 0.5, color = "green") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   theme_classic() +
#   labs(x = "Correlation", y = "-log10(p-value)",
#        title = "Naïve vs Spatially Corrected p-values") +
#   ylim(min(-log10(cors_df$corspv), na.rm = TRUE),
#        max(-log10(cors_df$corspv), na.rm = TRUE))
# corrected



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

save(cors_df, corrected,
     file = file.path("inst", "extdata", "simRanPatternResults.RData"))


