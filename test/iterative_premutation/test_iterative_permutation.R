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

# comparing after regenerating -----------------------------------------
# spatialCorrelation without iterative permutation, with MHT 

# pvalueCorrected 
load("test/iterative_premutation/data/simRanPatternRasts_permuted_test.RData")
pvalueCorrected_iter <- pvalueCorrected

# pvalueCorrected_orig_regen 
load("test/iterative_premutation/data/simRanPatternRasts_regen.RData")

data(simRanPatternRasts)


make_df <- function(pvalueCorrected_mat, simRanPatternRasts) {
  diag(pvalueCorrected_mat) <- NA
  
  cors <- matrix(NA, 
                 nrow = length(simRanPatternRasts), 
                 ncol = length(simRanPatternRasts))
  corspv <- cors
  
  for (i in seq_along(simRanPatternRasts)) {
    for (j in seq_along(simRanPatternRasts)) {
      vi <- intersect(
        rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[i]])),
        rownames(SpatialExperiment::spatialCoords(simRanPatternRasts[[j]]))
      )
      
      g1 <- SummarizedExperiment::assay(simRanPatternRasts[[i]], "pixelval")[1, vi]
      g2 <- SummarizedExperiment::assay(simRanPatternRasts[[j]], "pixelval")[1, vi]
      
      ctest <- cor.test(g1, g2)
      cors[i, j] <- ctest$estimate
      corspv[i, j] <- ctest$p.value
    }
  }
  
  cors.df <- reshape2::melt(cors, value.name = "cors")
  corspv.df <- reshape2::melt(corspv, value.name = "corspv")
  corrected.df <- reshape2::melt(pvalueCorrected_mat, value.name = "corspv_corrected")
  

  cors_df <- data.frame(
    cors.df,
    corspv = corspv.df$corspv,
    corspv_corrected = corrected.df$corspv_corrected
  )
  
  cors_df$corspv_corrected[cors_df$corspv_corrected == 0] <- 0.01
  
  test_df <- na.omit(cors_df)
  return(test_df)
}

library(reshape2)

iter_df <- make_df(
  pvalueCorrected_mat = pvalueCorrected_iter,
  simRanPatternRasts = simRanPatternRasts
)

orig_regen_df <- make_df(
  pvalueCorrected_mat = pvalueCorrected_orig_regen,
  simRanPatternRasts = simRanPatternRasts
)

head(iter_df)
head(orig_regen_df)

# overlap 

sig_iter <- iter_df[iter_df$corspv_corrected < 0.05, ]
sig_orig_regen <- orig_regen_df[orig_regen_df$corspv_corrected < 0.05, ]

dim(sig_iter)
dim(sig_orig_regen)

# Create a unique key for each Var1-Var2 pair
sig_iter$key <- paste(sig_iter$Var1, sig_iter$Var2, sep = "_")
sig_orig_regen$key <- paste(sig_orig_regen$Var1, sig_orig_regen$Var2, sep = "_")

# Rows in iter but NOT in orig_regen
only_in_iter <- sig_iter[!(sig_iter$key %in% sig_orig_regen$key), ]

# Rows in orig_regen but NOT in iter
only_in_orig_regen <- sig_orig_regen[!(sig_orig_regen$key %in% sig_iter$key), ]

only_in_iter
only_in_orig_regen

dim(only_in_iter)
dim(only_in_orig_regen)

# Pairs that were significant in both
common_keys <- intersect(sig_iter$key, sig_orig_regen$key)

both_in_iter <- sig_iter[match(common_keys, sig_iter$key), ]
both_in_orig_regen <- sig_orig_regen[match(common_keys, sig_orig_regen$key), ]

dim(both_in_iter)
dim(both_in_orig_regen)

# compare regen vs previous to see pairs removed as a result of MHT -----

# cors_df 
load("inst/extdata/simRanPatternResults.RData")

orig_not_regen_df <- na.omit(cors_df)

sig_orig_not_regen <- orig_not_regen_df[orig_not_regen_df$corspv_corrected < 0.05, ]

sig_orig_not_regen$key <- paste(
  sig_orig_not_regen$Var1,
  sig_orig_not_regen$Var2,
  sep = "_"
)

dim(sig_orig_not_regen)

# Rows in orig_regen but NOT in iter
only_in_orig_regen <- sig_orig_regen[!(sig_orig_regen$key %in% sig_orig_not_regen$key), ]
only_in_orig_not <- sig_orig_not_regen[!(sig_orig_not_regen$key %in% sig_orig_regen$key), ]

dim(only_in_orig_regen)
dim(only_in_orig_not)

common_keys <- intersect(sig_orig_not_regen$key, sig_orig_regen$key)
both_in_orig_regen <- sig_orig_regen[match(common_keys, sig_orig_regen$key), ]
dim(both_in_orig_regen)

# box plot 
library(ggplot2)
library(reshape2)

# 1) BH correction on corspv
iter_df$corspv_bh <- p.adjust(iter_df$corspv, method = "BH")

# 2) Avoid log(0)
iter_df$corspv_bh[iter_df$corspv_bh == 0] <- 1e-10
iter_df$corspv_corrected[iter_df$corspv_corrected == 0] <- 1e-10

# 3) Take -log10
iter_df$log_bh <- -log10(iter_df$corspv_bh)
iter_df$log_corrected <- -log10(iter_df$corspv_corrected)

# 4) Reshape for plotting
plot_df <- data.frame(
  niave_pval_bh_correct = iter_df$log_bh,
  emperical_pval = iter_df$log_corrected
)

plot_df_long <- melt(plot_df, variable.name = "type", value.name = "neglog10_pval")

# 5) Boxplot
plt <- ggplot(plot_df_long, aes(x = type, y = neglog10_pval)) +
  geom_boxplot() +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  theme_classic() +
  labs(
    y = "-log10(p-value)"  
  )
plt
ggsave("test_plot.png", plt, width = 7, height = 5, dpi = 300)

niave_sig <- plot_df[plot_df$niave_pval_bh_correct > -log(0.05), ]
niave_sig_per <- dim(niave_sig)[1] / dim(plot_df)[1]
niave_sig_per

emperical_sig <- plot_df[plot_df$emperical_pval > -log(0.05), ]
emperical_per <- dim(emperical_sig)[1] / dim(plot_df)[1]
emperical_per


## kidney datasets    --------------------------------------------------

# TODO: this should be on the rasterized and aligned?
# TODO: before or after figuring out STalign?



