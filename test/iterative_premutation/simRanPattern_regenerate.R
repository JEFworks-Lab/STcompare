
data("simRanPatternRasts")
set.seed(0)
library(STcompare)

# TODO: print the start time next time 

print("Starting run ==============")
t <- Sys.time()
print("Start time: ")
print(t)

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

t_end <- Sys.time()

pvalueCorrected_orig_regen <- pvalueCorrected
save(pvalueCorrected_orig_regen, file = "simRanPatternRasts_regen.RData")

print("End time: ")
print(t_end)

tf <- t_end - t
print("===== Runtime =====")
print(tf)


# TODO: save a log file when running this

