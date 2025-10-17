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

#' Simulated Random Spatially Patterned Kidney Datasets
#'
#' A list of 100 simulated \code{SpatialExperiment} objects representing kidney-shaped
#' datasets, each containing one independently generated spatially patterned gene with
#' no correlation between datasets. Each dataset consists of \eqn{N = 5000} simulated
#' cells distributed within a kidney-shaped region, with spatial coordinates and
#' expression values generated from Gaussian random fields.
#'
#' @details
#' Each simulated dataset was generated as follows:
#' \itemize{
#'   \item Cell coordinates \eqn{(x, y)} were uniformly sampled within a bounding box
#'   from \eqn{x, y \in [0,1]}.
#'   \item A Matern covariance matrix \eqn{Cov_{ij}} of size \eqn{N \times N} was
#'   computed as:
#'   \deqn{
#'     Cov_{ij} = \frac{2^{1-\nu}}{\Gamma(\nu)}
#'       \left(\frac{h \sqrt{2\nu}}{\kappa}\right)^{\nu}
#'       K_\nu \left(\frac{h \sqrt{2\nu}}{\kappa}\right)
#'   }
#'   where \eqn{h = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}}, \eqn{\nu = 0.5},
#'   \eqn{\kappa = 0.1}, \eqn{\Gamma} is the gamma function, and
#'   \eqn{K_\nu} is the modified Bessel function of the second kind.
#'   For diagonal elements (\eqn{i = j}), \eqn{Cov_{ij} = 1}.
#'
#'   \item The covariance matrix was used to generate a multivariate normal vector
#'   \eqn{W_i} with mean 0 via \code{MASS::mvrnorm}.
#'
#'   \item Independently and identically distributed noise \eqn{Z_i} was sampled using
#'   \code{stats::rnorm} with mean 0 and standard deviation \eqn{\sqrt{0.3}}.
#'
#'   \item The simulated gene expression for each cell was computed as:
#'   \deqn{G_i = W_i + Z_i + 10}
#'
#'   \item To form a kidney-shaped spatial pattern, each cell position was transformed
#'   to \eqn{(x, y) = (6x - 3, 6y - 3)} and retained only if its Euclidean distance from
#'   the origin satisfied:
#'   \deqn{\sqrt{x^2 + y^2} \le 1 - 2\sin(\theta) + \cos(\theta^2)}
#'   where \eqn{\theta = \mathrm{atan2}(y, x)}.
#' }
#'
#' Each element of \code{simRanPatternRasts} is a \code{SpatialExperiment} object with:
#' \itemize{
#'   \item One simulated gene (1 row)
#'   \item 250–300 spatial pixels (columns)
#'   \item An assay named \code{"pixelval"} containing expression values
#'   \item \code{colData} with columns: \code{num_cell}, \code{cellID_list},
#'   \code{geometry}, and \code{sample_id}
#'   \item Spatial coordinates (\code{x}, \code{y})
#' }
#'
#' @format A list of length 100, where each element is a \code{SpatialExperiment}
#' object containing:
#' \describe{
#'   \item{\code{assays}}{Matrix of simulated expression values named \code{"pixelval"}.}
#'   \item{\code{colData}}{Data frame of per-cell metadata including coordinates and IDs.}
#'   \item{\code{spatialCoords}}{Matrix of spatial \eqn{(x, y)} coordinates.}
#' }
#'
#' @source Simulated by the STcompare development team.
#'
#' @examples
#' data(simRanPatternRasts)
#' simRanPatternRasts[[1]]
#' assays(simRanPatternRasts[[1]])$pixelval[1, 1:5]
#'
#' @seealso \linkS4class{SpatialExperiment}, \code{\link[MASS]{mvrnorm}},
#'   \code{\link[stats]{rnorm}}
#'
#' @name simRanPatternRasts
#' @docType data
#' @keywords datasets
"simRanPatternRasts"


