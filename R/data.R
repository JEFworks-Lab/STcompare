#' Simulated Kidney SpatialExperiment Data
#'
#' A dataset containing three simulated spatial transcriptomics experiments
#' of kidney tissue. Each element is a \code{SpatialExperiment} object with
#' one gene, simulated for \code{A} and \code{C} to be similar spatial expression
#' patterns and \code{B} to have an opposite spatial expression pattern.
#'
#' @format A named list with three elements: \code{A}, \code{B}, and \code{C}.
#' Each element is a \code{SpatialExperiment} object with:
#' \describe{
#'   \item{assays}{A single \code{counts} matrix (1 × N), where N is the number of cells.}
#'   \item{rownames}{One gene: \code{"Gene"}.}
#'   \item{colnames}{Unique spot/cell IDs (\code{cell1}, \code{cell2}, ...).}
#'   \item{colData}{Contains one variable: \code{sample_id}.}
#'   \item{spatialCoords}{Two spatial coordinates: \code{x}, \code{y}.}
#' }
#'
#' Dimensions:
#' \itemize{
#'   \item \code{A}: 1 gene × 1229 cells
#'   \item \code{B}: 1 gene × 1242 cells
#'   \item \code{C}: 1 gene × 1297 cells
#' }
#'
#' @usage data(speKidney)
#'
#' @source Simulated for demonstration purposes.
#'
#' @docType data
#'
"speKidney"
