# Developer notes (not part of the package build)

This directory is excluded from `R CMD build` via `.Rbuildignore`. It records an investigation of STcompare (upstream commit `2983c99`) made on 2026-10-03. The aim was to make the package faster with a C++ backend, easier to use, and better documented.

- **`cpp-backend-plan.md`** — Start here. It covers where the time goes, what a C++ backend must reproduce, the proposed architecture, phases, and open decisions.
- **`engine-spec.md`** — The specification of the compiled engine behind the legacy functions: what it must reproduce bit for bit, its random streams, threads, batches and adaptive stopping.
- **`compare-spatial-spec.md`** — The specification of `compareSpatial()`, the new main function (phase 3): arguments, checks, the correlation and similarity tests, the output, progress and reproducibility.
- **`handoff-docs.md`** — The documentation phase's requests to the other phases (all done) and the questions it leaves for the maintainers.
- **`investigation/`** — Detailed reports, one per topic. Each claim cites source file:line or a reproduced measurement.
  - `01-locfit-smoother.md` — Exactly what `locfit(..., kern="gauss", nn=delta)` + `fitted()` computes: the adaptive `rbox` tree and bilinear interpolation, shown to be linear in the data and reproduced bit for bit.
  - `02-geoR-variogram.md` — Exact `geoR::variog` binning and estimator semantics, distance-formula and FMA subtleties, and lattice ties at `max.dist`.
  - `03-baseline-profiling.md` — Per-gene cost against pixel count and permutation count, a profile breakdown, BiocParallel overhead, and `spatialSimilarity` scaling.
  - `04-datasets-and-test-tiers.md` — Every dataset used or shipped, whether the `inst/extdata` results can be reproduced (they can, to about 1e-15), and the tiered test-data design implemented in `data-raw/` and `tests/testthat/fixtures/`.
  - `05-method-review-and-speedups.md` — Viladomat et al. (2014) compared with this implementation, statistical issues, speedups classed as exact, approximate or method-changing, and related tools.
  - `06-docs-and-ux-audit.md` — Places where the docs and code disagree, gaps a new user hits, API ergonomics, and a proposed user-facing design and docs plan.
  - `07-cpp-toolchain.md` — BLAS, threading (`std::thread`, RcppParallel or OpenMP), RNG, package plumbing, and portability measurements.
  - `08-bug-audit.md` and `09-bug-audit-verified.md` — The 37 findings and the baseline `R CMD check`. Each finding was re-run by an independent verifier: none was refuted, 23 were confirmed and 14 partly confirmed.
- **`prototypes/`** — Working proof-of-concept code. These are not package code, and their paths and APIs are rough.
  - `lfsmooth.cpp` — C++ replica of locfit's tree smoother and its factored linear operator. Its output is `identical()` to `fitted()` (verified on 24 coordinate/delta configurations up to N = 5041). Compile with `PKG_CXXFLAGS=-ffp-contract=off`.
  - `variog_fast.cpp` and `variog_fast.R` — Ports of `geoR::variog` that precompute the pair-to-bin table. They are bit-identical to geoR, including on lattice ties and on 300 random layouts.
  - `fastgene.cpp`, `fastvil.R`, `factorS.R` and `test_cpp.R` — An end-to-end engine for one gene and one direction. At N = 4984 it gives the same delta* and p-values as the package, with permutations within 5.6e-11. It takes 1.1 s per gene, against about 91 s.
  - `rng.cpp` — Per-permutation RNG streams that do not depend on the thread count.
  - `stdthread-package-skeleton/` — The DESCRIPTION, NAMESPACE, `R/` and `src/` changes that add compiled code to this package with `std::thread`. They passed the compiled-code checks of `R CMD check`.

The reports mention scratch paths under `/private/tmp/...`. Those were temporary and no longer exist. The relevant scripts are summarised inside each report.
