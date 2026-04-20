# testing  ---------------------------------------------------------------

## 100 random dataset   --------------------------------------------------
# code from file:
# system.file("extdata", "simRanPatternResults.RData", package = "STcompare")

data("simRanPatternRasts")
set.seed(0)
library(STcompare)

print("Starting run ==============")

# simple 1 iteration tests 

rastGexpListAB <- list(simRanPatternRasts[[1]], simRanPatternRasts[[2]])

# test thread of 1 on test
test_1 <- spatialCorrelationGeneExp_test(
  rastGexpListAB,
  nPermutations = c(1e2, 1e3),
  nThreads = 1,
  verbose = FALSE,
  BPPARAM = BiocParallel::MulticoreParam()
)

# make sure test matches that of old (no premutation function)
test_2 <- STcompare::spatialCorrelationGeneExp(
  rastGexpListAB,
  nThreads = 1,
  verbose = FALSE,
  BPPARAM = BiocParallel::MulticoreParam()
)

# Franklin can handle nThreads = 22 
test_3 <- spatialCorrelationGeneExp_test(
  rastGexpListAB,
  nPermutations = c(1e2, 1e3),
  nThreads = 22,
  verbose = FALSE,
  BPPARAM = BiocParallel::MulticoreParam()
)



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

diag(pvalueCorrected) <- NA
tf <- Sys.time() - t
print(tf)
save(pvalueCorrected, file = "simRanPatternRasts_permuted_test.RData")

load("simRanPatternRasts_permuted_test.RData")

# cors <- matrix(NA, nrow = length(simRanPatternRasts), ncol = length(simRanPatternRasts))
cors_iter <- matrix(NA, nrow = length(pvalueCorrected), ncol = length(pvalueCorrected))

corspv <- cors_iter

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

cors_df_iter <- data.frame(cors.df,
                      corspv = corspv.df$corspv,
                      corspv_corrected = pvalueCorrected.df$corspv_corrected)



cors_df_iter$corspv_corrected[cors_df_iter$corspv_corrected == 0] <- 0.01

test_df_iter <- na.omit(cors_df_iter)
test_df_prev <- na.omit(cors_df)


hist(cors_df_iter$corspv)
hist(cors_df_iter$corspv_corrected)
hist(-log10(cors_df_iter$corspv_corrected))




sig_iter <- test_df[test_df_iter$corspv_corrected < 0.05, ]
sig_prev <- test_df_prev[test_df_prev$corspv_corrected < 0.05, ]

# TODO: confused about why these two have diff numbers 
dim(test_df_iter) 
dim(test_df_prev)

dim(sig_iter)
dim(sig_prev)

# Create a unique key for each Var1-Var2 pair
sig_iter$key <- paste(sig_iter$Var1, sig_iter$Var2, sep = "_")
sig_prev$key <- paste(sig_prev$Var1, sig_prev$Var2, sep = "_")

# Rows in sig_iter but NOT in sig_prev
only_in_iter <- sig_iter[!(sig_iter$key %in% sig_prev$key), ]

# Rows in sig_prev but NOT in sig_iter
only_in_prev <- sig_prev[!(sig_prev$key %in% sig_iter$key), ]

only_in_iter
only_in_prev

# PROBLEM: there is little overlap with the two
dim(only_in_iter)
dim(only_in_prev)



## kidney datasets    --------------------------------------------------

# TODO: this should be on the rasterized and aligned?
# TODO: before or after figuring out STalign?



