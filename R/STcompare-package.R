#' @details The main function is \code{\link{compareSpatial}()}: for every gene
#'   measured in two samples rasterized onto one grid of pixels, it tests
#'   whether the spatial patterns are correlated beyond what their spatial
#'   autocorrelation alone would give, and measures how often the expression at
#'   matched pixels is within a fold change. The functions of the published
#'   analyses (\code{\link{spatialCorrelationGeneExp}()},
#'   \code{\link{spatialCorrelationGeneExpIterPermutations}()} and
#'   \code{\link{spatialSimilarity}()}) reproduce them exactly.
#'
#' @references Clifton K, Jiang V, Peixoto RdS, Singh S, Matsuura R, Rabb H,
#'   Fan J (2026). STcompare: comparative spatial transcriptomics data analysis
#'   of structurally matched tissues to characterize differentially spatially
#'   patterned genes. \emph{Bioinformatics} 42(9):btag644.
#'   \doi{10.1093/bioinformatics/btag644}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib STcompare, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom ggplot2 .data
## usethis namespace: end
NULL
