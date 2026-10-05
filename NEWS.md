# STcompare (development version)

## Changes to results

* Empirical p-values are now `(b + 1) / (B + 1)`.
  * `b` is the number of null correlations whose absolute value is at least the absolute value of the
    observed correlation, and `B` is the number of permutations.
  * Previously they were `b / B` with a strict `>`, so they could be exactly 0. The smallest possible p-value
    is now `1 / (B + 1)`.
  * Affected functions: `viladomatCorrelation()`, `spatialCorrelation()`, `spatialCorrelationGeneExp()`,
    `spatialCorrelationGeneExpIterPermutations()` and `spatialCorrelationGeneExpWithinSample()`.
  * Null correlations and delta stars are unchanged. With the default `alpha = 0.05`,
    `spatialCorrelationGeneExpIterPermutations()` screens the same genes into each round as before.

* `spatialCorrelationGeneExp()` now applies `adjustMethod` across all genes, separately to
  `pValuePermuteX` and `pValuePermuteY`, as documented. Previously it was applied to one gene at a time,
  which left every p-value unadjusted.

## Bug fixes and minor improvements

* `spatialCorrelationGeneExp()` now uses a user-supplied `BPPARAM`. It also checks `adjustMethod` before
  computing anything.

* The documentation of `spatialCorrelationGeneExpIterPermutations()` now describes the screening rule the code
  applies: a gene is carried forward when both unadjusted p-values are below `100 * alpha / nPermutations[k]`.
