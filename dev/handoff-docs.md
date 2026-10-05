# Handoff from the documentation phase (phase 4)

The documentation phase of `cpp-backend-plan.md` (phase 4) rewrote the tutorials, added two articles,
and rewrote the README and `_pkgdown.yml`. This file lists the changes to files that other phases own and
the open questions for the maintainers.

**Status (2026-10-05).** Every request of section 2 is done, each marked **Done** below, except the two
optional items 7 (move the README figure to `man/figures/`) and 8 (continuous integration), which are left
to the maintainers. Of the questions of section 4, items 3 and 5 were updated when `surrogate = "remap"`
became the default of `compareSpatial()` (2026-10-05; `bench/calibration-results.md`); items 1, 2 and 4 are
still open.

**Files changed by the documentation phase:**

- `vignettes/`: three package vignettes and three pkgdown-only articles (see the layout below). The old
  knitr figure folders (`vignettes/acute-kidney-injury-10x-visium/`, `vignettes/brain-MERFISH-10x-visium/`
  and `vignettes/getting-started-with-STcompare-figures/`, about 13 MB of PNGs) are deleted; nothing
  produces them any more.
- `README.md` and `_pkgdown.yml`.
- `inst/extdata/vignette-aki-svg-genes.txt` and `inst/extdata/vignette-brain-svg-genes.txt` (9 KB): the
  spatially variable genes of the two case studies (see below).

## 1. Layout: vignettes and articles

Run times are for `rmarkdown::render()` on an M1 Ultra. The package vignettes use 2 threads and the case
studies 8; the first run of a case study also downloads its data.

| File | Kind | Data | Run time |
|---|---|---|---|
| `vignettes/getting-started-with-STcompare.Rmd` | package vignette | built-in (`speKidney`, `simRanPatternRasts`) | 13 s |
| `vignettes/how-STcompare-works.Rmd` | package vignette | built-in | 10 s |
| `vignettes/parameters-performance-reproducibility.Rmd` | package vignette | built-in | 6 s |
| `vignettes/articles/Install.Rmd` | pkgdown article | none | 1 s |
| `vignettes/articles/acute-kidney-injury-10x-visium-rasterized.Rmd` | pkgdown article | 30 MB from Zenodo, cached | 49 s (88 s with the first download) |
| `vignettes/articles/brain-MERFISH-10x-visium.Rmd` | pkgdown article | 75 MB from Zenodo and 10x Genomics, cached | 86-88 s (107 s with the first download) |

**Why this split.** The three package vignettes use only the data that ship with the package, run in well
under a minute in total with 2 threads, and therefore work under `R CMD check` and on the Bioconductor
builders; `browseVignettes("STcompare")` gives them offline. The two case studies download about 105 MB
and run for one to two minutes on 8 threads, so they are pkgdown articles (`vignettes/articles/`), which
`R CMD build` does not build. The installation page is only useful before the package is installed, so it is
an article too. All URLs of the website are unchanged (`articles/<name>.html`), so no redirects are needed.

**Case-study settings.**

- The articles use `threads <- min(8L, max(1L, parallel::detectCores() - 1L, na.rm = TRUE))` and
  `nPermutations = 1000`. On these data, 1000 gives the same significant genes as the default of 10,000:
  with the default rank-remapped surrogates (2026-10-05), 752 of the 1046 AKI genes (726 positive, 26
  negative) and 173 of the 325 brain genes, none skipped; with gaussian surrogates and their sqrt(N)
  detection filter it was 738 of the 1039 tested AKI genes and 132 of the 297 tested brain genes (7 and 28
  skipped). The articles explain the trade-off.
- **Downloads.** They go to `file.path(tools::R_user_dir("STcompare", "cache"), "downloads")`. They are
  checked against the same MD5 sums as `data-raw/download_data.R`, written to `<file>.part` first, and use a
  timeout of at least 30 minutes. `R_USER_CACHE_DIR` relocates the cache.
- **Inputs.** The articles rebuild exactly the published rasterized inputs: identical CPM and log values, the
  same 311 and 2170 shared pixels, and the same 325 analysed brain genes, checked against
  `data-raw/build_inputs_*.R`.

## 2. Required changes in files owned by the housekeeping phase

1. **`.Rbuildignore`.** Done by the housekeeping phase while this phase ran:
   - `^vignettes/articles$` is listed, so the case-study articles and their download code stay out of the
     tarball;
   - the entries of the deleted figure folders are gone.

   `R CMD build` builds the three vignettes in about 40 s and gives a 2.4 MB tarball. **Done.**
2. **`DESCRIPTION`:**
   - **Already right as of this writing:** `VignetteBuilder: knitr`, `Suggests: knitr, rmarkdown, patchwork`
     (used by the vignettes; ggplot2, SEraster and SpatialExperiment are imports), and
     `Config/Needs/website: jsonlite, Matrix, rhdf5` (used only by the articles).
   - **Keep:** reshape2, scatterbar, gridExtra and MERINGUE out of Suggests; nothing uses them any more.

   **Keep two housekeeping changes the tutorials rely on:**
   - `plotCorrelationGeneExp()` accepts a `compareSpatial()` result. "Getting started" and the brain case
     study pass one.
   - `linearRegression()` and `pixelClass()` accept an explicit `assayName`, which the kidney case study
     passes.

   **Done** (DESCRIPTION also requires R 4.5 now: SEraster is in Bioconductor from release 3.21).
3. **`NEWS.md`.** Add a "Documentation" section with these bullets:
   - The tutorials run the analyses live with `compareSpatial()` instead of loading precomputed results.
     "Getting started" uses the built-in data, including a null-calibration check on 50 independent pairs of
     `simRanPatternRasts` that replaces `simRanPatternResults.RData`. The case studies "Acute kidney injury
     (10x Visium)" and "Comparison of MERFISH and Visium for mouse brain" download and cache their inputs
     and run in one to two minutes on 8 threads.
   - New articles: "How STcompare works" (the tests step by step, with figures from the built-in data) and
     "Parameters, performance and reproducibility" (choosing the settings, run times, threads, seeds, and the
     legacy functions compared with `compareSpatial()`).
   - The case studies no longer need MERINGUE or scatterbar. The spatially variable genes of the published
     analyses ship in `inst/extdata/vignette-aki-svg-genes.txt` and `vignette-brain-svg-genes.txt`.
   - The README has the installation with BiocManager and its C++17 compiler requirement, a quick start, links
     to every tutorial, and corrected descriptions: alignment is done with STalign and rasterization with
     SEraster, and it is "Pearson's", not "Person", correlation.
   - The website has Tutorials and Articles menus, a grouped function reference and a changelog. The case
     studies and the installation page are pkgdown articles, not package vignettes.

   **Done:** section "Documentation" of `NEWS.md`.
4. **Roxygen wording.** Suggestions only; the topics render correctly as they are.
   - **Titles in the reference index.** These topics show the function name as their title:
     - `spatialCorrelationGeneExp`: "Spatial correlation test for every gene of two samples";
     - `spatialCorrelationGeneExpIterPermutations`: "Spatial correlation test with more permutations for
       promising genes";
     - `spatialCorrelationGeneExpWithinSample`: "Spatial correlation test between pairs of genes of one
       sample";
     - `spatialCorrelation`: "Spatial correlation test of two vectors of values at the same locations";
     - `viladomatCorrelation`: "Spatially autocorrelated surrogates and permutation p-value for one
       direction";
     - `plotCorrelationGeneExp`: "Plot a gene's values in two samples at the shared pixels".

     These titles also suit other topics:
     - `linearRegression` ("Generates linear regression plot for a given gene"): "Plot a gene's values in two
       samples, coloured by similarity class";
     - `pixelClass`: "Map the similarity class of every pixel";
     - `savePlots`: "Maps, pixel classes and scatter plot of genes, optionally saved as PDF".
   - **`compareSpatial()` examples.** `ac[, c("gene", "r", ...)]` keeps the class `"STcompareResult"` and
     prints a header with "Correlation:" and no settings. Either write `as.data.frame(ac)[, ...]`, as the
     vignettes do, or add a `[.STcompareResult` method that returns a plain data frame when columns are
     selected.
   - **`compareSpatial()` `@seealso`.** Point to `vignette("getting-started-with-STcompare", package =
     "STcompare")`, `vignette("how-STcompare-works", package = "STcompare")` and
     `vignette("parameters-performance-reproducibility", package = "STcompare")`.
   - **`compareSpatial()` `@param x,y`.** Say that the list returned by
     `SEraster::rasterizeGeneExpression(list(a = ..., b = ...))` can be passed directly as `x`. The tutorials
     do this (`compareSpatial(rast, ...)`).

   **Done:** the nine titles, the examples (`as.data.frame()`, and a `[` method that gives a plain data frame
   when columns are selected), the `@seealso` and the `@param x,y` wording.
5. **`inst/CITATION`.** It exists now. The README keeps the authors' citation text (the same paper). It
   could add one line: "`citation("STcompare")` gives the reference in BibTeX". **Done** (the README's
   Citation section).
6. **`data-raw/` (optional).**
   - **What to add.** A `data-raw/vignette_gene_lists.R` that regenerates the two gene lists.
   - **What it should do.** The derivation is summarised in section 3, with the essential code. It needs
     MERINGUE from GitHub and the inputs of `data-raw/build_inputs_*.R`.

   **Done:** `data-raw/vignette_gene_lists.R` writes both lists, identical to the shipped files.
7. **`images/`.** The README shows `images/overview_figure1.png` through its URL on the `main` branch of
   JEFworks-Lab/STcompare, so it works on GitHub and on the website. If the figure moves to
   `man/figures/`, the pkgdown convention, change the README to `man/figures/overview_figure1.png`.
   **Not done** (optional): the figure works through its URL.
8. **Continuous integration (if added).**
   - **What it needs.** A pkgdown workflow needs network access, the `Config/Needs/website` packages, and
     about 5 to 10 minutes on 4 cores, most of it the two case studies.
   - **Caching.** Cache `~/.cache/R/STcompare` (Linux) between runs, so that the downloads happen once.

   **Not done** (optional): no workflow was added.

## 3. The shipped gene lists

**`vignette-aki-svg-genes.txt`.**

- **What it holds.** The 1046 genes of the published AKI analysis, which are the row names of
  `kidneyCorrelation.RData`, now in `bench/published/`. They are given in the gene symbols of the 10x feature
  files.
- **How they were selected.** The published analysis named them with `make.names()`, which changed 19
  symbols, for example `mt.Nd1` for `mt-Nd1`. The genes were spatially variable in both sections by
  MERINGUE's Moran's I, with an adjusted p-value of exactly 0 in both.
- **Rerunning MERINGUE.** MERINGUE 1.0 on the rebuilt inputs gives 1044 genes. Three published genes are
  missing (Snhg3, Mrpl23 and Baiap2l2) and one gene is added (Por). The p.adj == 0 criterion depends on
  floating-point underflow, so the article ships the published list and says so.
- **Gene name format.** The 10x feature files have 40 duplicated symbols. The article uses `make.unique()`,
  which keeps the first occurrence's name. All 1046 published genes are first occurrences, so the list maps
  one to one.

**`vignette-brain-svg-genes.txt`.**

- **What it holds.** 230 of the 325 analysed genes: those spatially variable in both technologies.
- **How it was computed.** The original vignette's MERINGUE steps were rerun exactly: `getSpatialPatterns()`
  and `filterSpatialPatterns()` with adjustPv = TRUE, alpha = 0.05 and minPercentCells = 0.01, on the MERFISH
  pixels (lognorm, filterDist 21) and on the Visium spots (libnorm, filterDist 25).
- **Check.** It gives 230 and 95 genes, the numbers of the rendered published vignette.

**Essential code.** This is condensed from the script this phase used. Run it from the repository root, after
`Rscript data-raw/build_inputs_aki.R` and `build_inputs_brain.R`.

```r
source("data-raw/download_data.R")
library(SpatialExperiment)
moran <- function(mat, coords, filterDist, filter = FALSE) {
  w <- MERINGUE::getSpatialNeighbors(coords, filterDist = filterDist)
  I <- MERINGUE::getSpatialPatterns(mat, w)
  if (!filter) return(I)
  MERINGUE::filterSpatialPatterns(mat = mat, I = I, w = w, adjustPv = TRUE, alpha = 0.05,
                                  minPercentCells = 0.01, details = TRUE)
}
# AKI: the published genes (row names of kidneyCorrelation, in make.names() form), mapped back to 10x symbols
aki <- readRDS(stc_input_file("aki_rast.rds"))
symbols <- as.character(rhdf5::h5read(stc_download(group = "aki")[["NL3_filtered_feature_bc_matrix.h5"]],
                                      "matrix/features")$name)
akiGenes <- symbols[match(aki$genes, make.names(symbols, unique = TRUE))]
# check: MERINGUE (p.adj == 0 in both sections) gives the same list up to 4 genes
iC <- moran(assay(aki$rast$AKI_ctrl, "CPM"), spatialCoords(aki$rast$AKI_ctrl), 10)
iA <- moran(assay(aki$rast$AKI_aki, "CPM"), spatialCoords(aki$rast$AKI_aki), 10)
svgMeringue <- intersect(rownames(iC)[iC$p.adj == 0], rownames(iA)[iA$p.adj == 0])
setdiff(aki$genes, svgMeringue)   # Snhg3, Mrpl23, Baiap2l2
setdiff(svgMeringue, aki$genes)   # Por
# brain: the MERFISH pixels (lognorm of the 325 analysed genes, filterDist 21) and the Visium spots
# before rasterization (assay libnorm of the 466 shared genes, filterDist 25): Visium_SE is built by
# data-raw/build_inputs_brain.R, which does not save it, so run its lines up to Visium_SE first
brain <- readRDS(stc_input_file("brain_merfish_visium_rast.rds"))
iM <- moran(assay(brain$rast$MERFISH, "lognorm"), spatialCoords(brain$rast$MERFISH), 21, filter = TRUE)
iV <- moran(assay(Visium_SE, "libnorm"), spatialCoords(Visium_SE), 25, filter = TRUE)
brainGenes <- intersect(rownames(iM), rownames(iV))   # 230 genes
```

Write each list with a comment header, as in the shipped files.

## 4. Open questions for the maintainers

1. **AKI rotation.** The published preprocessing, which the article reproduces so that the published inputs
   are matched exactly, rotates each section with its own maximum: `y = max(x) - x`.
   - **The problem.** After the STalign alignment, the aligned AKI positions have a maximum of 128.95 and the
     control's is 127. The rotation therefore shifts the AKI section by 1.95 array units relative to the
     control, about 0.4 of a pixel at resolution 5.
   - **The likely fix.** One common offset for both sections is probably what was intended. It would change
     the rasterization slightly, and with it the gene list and the published numbers.
2. **Similarity scale in the brain case study.**
   - **What changed.** The published analysis computed the similarity on `log10(CPM + 1)`, where a ratio is
     not a fold change. The article computes it on the CPM (`pixelval`) and explains why.
   - **The effect.** The median similarity of the analysed genes is then about 0.05: the two technologies'
     levels rarely agree within two-fold.
   - **Why the correlation is unaffected.** It is still computed on the log values, as published.
3. **Published numbers quoted.** (Updated 2026-10-05 for the default `surrogate = "remap"`, with
   `minDetected = NULL` meaning no filter.)
   - **What the article quotes.** The AKI article quotes 707 positive and 24 negative genes for the
     published, legacy analysis.
   - **What the article computes live.** 726 positive and 26 negative genes (no gene is skipped). With
     gaussian surrogates and the sqrt(N) filter it was 711 and 27 (7 skipped).
   - **Where the agreement comes from.** The articles "Getting started" and "Parameters, performance and
     reproducibility" quote the agreement measured with the default settings of `compareSpatial()` (seed 0,
     10,000 permutations; `nPermutations = 1000` gives the same significant genes): 98% (AKI) and 86%
     (brain) of genes are classified the same, counting a gene as significant with the legacy function when
     both of its adjusted p-values are below 0.05.
     - AKI: 729 genes significant with both methods, 2 with the legacy method only (Rps3a1 and Rps15a, padj
       0.076 and 0.056), 23 with `compareSpatial()` only (padj 0.02 to 0.05).
     - Brain: 128 with both (every published gene), 0 legacy only, 45 `compareSpatial()` only, all positively
       correlated (r = 0.06 to 0.24). The calibration study (`bench/calibration-results.md`, item 5) attributes
       them to the remapped surrogates: gaussian surrogates are conservative for skewed, zero-inflated genes
       with a spatial pattern, as the log-normalized brain pixels are. With `surrogate = "gaussian"` and the
       sqrt(N) filter the brain agreement was 95% (122 with both, 6 legacy only of which 4 skipped, 10 new
       only, 28 skipped), and the AKI agreement 98% (722, 8 of which 1 skipped, 16, 7 skipped).
4. **Example genes in the brain case study.** The brain article takes the example genes from its own
   results, so they can differ from the published figure. They are the strongest positive SVG (Slc17a7) and
   the non-significant non-SVG with the lowest correlation (Gpr160 with the remapped surrogates); the gene
   shown before the test is Slc17a6 (Baiap2, which the published analysis found significant, has padj 0.04
   with the remapped surrogates; it had 0.06 with gaussian ones).
5. **A rare matched cell type.** Resolved by the remapped default (2026-10-05): CT-16 (in 28 of 2174 MERFISH
   pixels) is tested in the brain article and is significant (r = 0.76, padj 0.0036 at 1000 permutations),
   as in the published analysis, so all 9 matched cell types are significantly positively correlated. With
   `surrogate = "gaussian"` the sqrt(N) filter (47 pixels) still skips it.

## 5. Checks done

- **Rendering.** Every vignette and article rendered with `rmarkdown::render()` in a mirror of the
  repository, with the package installed from the mirror (`R CMD INSTALL`, `-O2`), and again within
  `pkgdown::build_site()`. The phase report lists the times.
- **Website.** `pkgdown::build_site(override = list(destination = <temporary directory>))` builds without
  errors. Its only notes are:
  - the "stack imbalance" lines printed by SparseArray, which are upstream (see `data-raw/README.md`);
  - pandoc's "Deprecated: --highlight-style" warning, which comes from rmarkdown with pandoc 3.9 (Homebrew).
  The missing alt text of the README figure is fixed.
- **Downloads.** The download code was run from an empty cache for both case studies.
- **Package build.** `R CMD build` (with `_R_CHECK_LIMIT_CORES_=TRUE`) builds the three vignettes and gives a
  2.4 MB tarball without `vignettes/articles/`.
