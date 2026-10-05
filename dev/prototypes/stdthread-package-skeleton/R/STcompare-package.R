#' @keywords internal
"_PACKAGE"

## Package-level roxygen tags for the compiled backend. NOTE: these only reach
## NAMESPACE if roxygen2 manages it (the current hand-written NAMESPACE with
## exportPattern() is skipped by roxygen2; edit it by hand or convert it).
#' @useDynLib STcompare, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
