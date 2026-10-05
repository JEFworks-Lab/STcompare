# Load the Bioconductor namespaces that the tests use before any test runs. Loading SparseArray (pulled in
# by SEraster -> SpatialExperiment -> SummarizedExperiment) leaves one unbalanced PROTECT in SparseArray 1.10,
# so R prints "Warning: stack imbalance in ..." when the namespace is first loaded. That message is not an R
# condition, so testthat does not count it. Loading here keeps the known message out of the tests: a stack
# imbalance printed while the tests run (for example by a compiled backend) is then new.
invisible(requireNamespace("SEraster", quietly = TRUE))
