

library(STcompare)

data(speKidney)

rastGexpListAB <- list(A = rastKidney$A, B = rastKidney$B)
scAB <- spatialCorrelationGeneExp(rastGexpListAB)
scAB
