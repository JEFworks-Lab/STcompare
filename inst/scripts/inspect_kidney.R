
library(dplyr)

# kidneyCorrelation
load("inst/extdata/kidneyCorrelation.RData")

head(kidneyCorrelation)
dim(kidneyCorrelation)

# 1. did the number of significant genes change? ##################

# saving names of the significantly positively correlated svg genes 
svgSigPos <- kidneyCorrelation %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%  
  dplyr::filter(correlationCoef > 0) %>% 
  rownames()

# previously 360 genes are SVGs that are significantly positively correlated (from the tutorial)
# now, 713 
length(svgSigPos)


# saving names of the significantly negatively correlated svg genes 
svgSigNeg <- kidneyCorrelation %>%
  dplyr::filter(pValuePermuteX < 0.05 & pValuePermuteY < 0.05) %>%  
  dplyr::filter(correlationCoef < 0) %>% 
  rownames()
# 5 genes are SVGs tha
# now 26 
length(svgSigNeg)


# 1b. previously (was this data stored somewhere?): 



# 1c. we use the previous example genes? (--> try finding in old repo, must be stored somewhere?)
# here is the MHT comparision and the nPermutations=1000 comparision, there are multiple things that changed 
# FIRST: just try to reproduce the plots again

# 1d. the number of genes that required nPermutation = 1000 
# from the deltaStarX length 

# Add nPermutations column = length of deltaStarX for each row
kidneyCorrelation$nPermutations <- sapply(
  kidneyCorrelation$deltaStarX,
  length
)

head(kidneyCorrelation)

# Filter rows where nPermutations is 100 or 1000
kidneyCorrelation_hund <- kidneyCorrelation[kidneyCorrelation$nPermutations == 100, ]
kidneyCorrelation_thou <- kidneyCorrelation[kidneyCorrelation$nPermutations == 1000, ]

dim(kidneyCorrelation_hund)
dim(kidneyCorrelation_thou)



# 2. what genes changed ##################


# kidney MHT, no iterative: ###########################

#visualize the correlation coefficient for the significantly positively correlated, significantly negatively correlated, and not significantly correlated svg genes as violin plots with all points overlaid
fig_2j <- kidneyCorrelation %>%
  dplyr::mutate(Sig = dplyr::case_when(rownames(kidneyCorrelation) %in% svgSigPos ~ "SigPos",
                                       rownames(kidneyCorrelation) %in% svgSigNeg ~ "SigNeg",
                                       .default = "")) %>%
  dplyr::mutate(pValueEmpirical = dplyr::case_when(pValuePermuteY > pValuePermuteX ~ pValuePermuteY,
                                                   .default = pValuePermuteX)) %>%
  dplyr::mutate(pValueEmpiricalRound = dplyr::case_when(pValueEmpirical == 0 ~ 0.001,
                                                        .default = pValueEmpirical)) %>%
  ggplot2::ggplot(ggplot2::aes(x= Sig, y = correlationCoef, color = Sig)) +
  ggplot2::geom_boxplot() +
  #ggplot2::geom_jitter(width = 0.2, alpha = 0.1) +
  ggplot2::ylim(-1,1) +
  ggplot2::scale_color_manual(values = c("SigPos" = "green", "SigNeg" = "blue", "NotSig" = "grey")) +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "Significance", y = "Correlation Coefficient", title = "Correlation Coefficient for SVGs")
fig_2j



ggsave(
  filename = "plots/kidneyBox.png", 
  plot = fig_2j, 
  height = 4, 
  width = 4
)

