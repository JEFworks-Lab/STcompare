# STcompare 0.1.0.9000

## A compiled engine replaces the R implementation

* `viladomatCorrelation()`, `spatialCorrelation()`, `spatialCorrelationGeneExp()`,
  `spatialCorrelationGeneExpIterPermutations()` and `spatialCorrelationGeneExpWithinSample()` are now computed
  by compiled code (C++ through Rcpp). Their names, arguments and output structure are unchanged.
  * The method is unchanged: the code reimplements locfit's local-constant Gaussian kernel smoother (its
    adaptive k-d tree with vertex interpolation) and `geoR::variog()`'s binned variogram, and it draws the same
    permutations and noise as before.
  * It reproduces the published analyses in `inst/extdata`: on the test genes (35 AKI kidney genes on 311
    pixels, 30 brain genes on 2170 pixels, 100 permutations), delta\* is identical for every permutation, and
    the null correlations agree within 3e-14 relative. The empirical p-values are therefore the same.
  * Speed, per gene with 100 permutations in both directions (macOS arm64, M1 Ultra), against the R
    implementation on one core (19-20 s for AKI and 65-69 s for brain):
    * one thread: 0.029 s (AKI, 11 deltas) and 0.13 s (brain, 9 deltas), about 650 and 500 times faster;
    * 16 threads: 0.0035 s and 0.016 s per gene.
  * `nThreads` is now the number of threads of the compiled code. A non-`NULL` `BPPARAM` only sets that number
    (`BiocParallel::bpnworkers(BPPARAM)`); no BiocParallel back-end is started and nothing is forked. Results
    do not depend on the number of threads.
  * `spatialCorrelationGeneExp()` computes all genes in one call, and
    `spatialCorrelationGeneExpIterPermutations()` extends the genes it carries forward from
    `nPermutations[k] + 1` instead of recomputing their first permutations. The results are the same as a fresh
    run with the larger number of permutations.
  * `spatialCorrelationGeneExpWithinSample()` computes the permutations of every gene once and correlates them
    with all the other genes, instead of permuting both genes again for every pair, and it builds its result
    table without looping over the pairs, whose number grows with the square of the number of genes: 600 genes
    (179,700 pairs) with 20 permutations take about 2.6 s on 8 threads of an M1 Ultra.
  * The global random number generator state is left unchanged. The permutations are still drawn with
    `set.seed(seed)` under the session's RNG kinds, and the noise of permutation b with `set.seed(seed + b)`
    under L'Ecuyer-CMRG.
* Failures no longer print errors. Where the R implementation returned an NA row (for example a constant gene
  or a missing value), there is still an NA row, plus one `warning()` per gene that gives the reason. Inputs that
  crashed R or failed with an unrelated error also give an NA row and a warning:
  * `1 <= delta * N < 2`, where locfit killed the R session;
  * exactly duplicated coordinates with a small delta, where locfit overflowed the C stack;
  * fewer than 3 complete pairs, where `cor.test()` failed and the error handler raised "object 'corDF' not
    found".
  * `spatialCorrelationGeneExpIterPermutations()` no longer stops when a gene has an NA row: such genes are not
    carried forward to the next round.
* An NA row now has the same column types as other rows (numeric p-values, list columns with `NA` elements).
  `viladomatCorrelation()` returns `NA` elements with a warning instead of an error, and when Y is constant or
  has missing values it returns the permutations with `NA` nulls and p-value, as before, plus a warning.
* `verbose = TRUE` prints one message when the computation starts and one when it ends, instead of one per gene.
* `matchingVariograms()` is removed. It was the R helper of the old implementation, and its interface took a
  geoR variogram object.
* geoR and locfit are no longer needed at run time; they moved to Suggests (the tests use them as references).
* Installing STcompare from source needs a C++17 compiler. A source install stops with an explanatory error if
  the compiler flags include `-ffast-math`, `-Ofast` or similar unsafe floating-point options (for example from
  `~/.R/Makevars`): the compiled code reproduces R's arithmetic exactly.

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
    With other values of `alpha`, a gene at the screening boundary can be decided differently. When
    `100 * alpha` is a whole number, a gene is now carried forward exactly when it has fewer than
    `100 * alpha` exceedances in each direction; the earlier comparison of `b / B` with the threshold also
    kept some genes with exactly `100 * alpha`, depending on floating-point rounding (for example
    `alpha = 0.07` with 1000 permutations). Otherwise a gene with exactly `floor(100 * alpha)` exceedances
    may no longer be carried forward: with `alpha = 0.053` and 100 permutations, a gene can have at most 4
    exceedances instead of 5.

* `spatialCorrelationGeneExp()` now applies `adjustMethod` across all genes, separately to
  `pValuePermuteX` and `pValuePermuteY`, as documented. Previously it was applied to one gene at a time,
  which left every p-value unadjusted.

## Bug fixes and minor improvements

* `spatialCorrelationGeneExp()` checks `adjustMethod` before computing anything.

* The documentation of `spatialCorrelationGeneExpIterPermutations()` now describes the screening rule the code
  applies: a gene is carried forward when both unadjusted p-values are below `100 * alpha / nPermutations[k]`.

* The examples of `spatialCorrelationGeneExpWithinSample()` and `spatialCorrelationGeneExpIterPermutations()`
  run (they passed a list where one object was needed, and printed an undefined object).
