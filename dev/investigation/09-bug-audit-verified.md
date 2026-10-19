# Bug audit: findings with independent verification

Audit of upstream STcompare at commit 2983c99, 2026-10-03.

An auditor reported 37 findings. Each was then re-run by a separate skeptic agent that tried to refute it:

- 18 critical, high or medium findings were checked one at a time.
- The 19 low findings were checked in batches.

Results: none refuted, 23 confirmed, 14 partly confirmed. A partly confirmed finding is real but had some details wrong; its corrected description is shown under **Verified description**.

The `Severity` column is the verifier's assessment. For the full narrative, including the R CMD check output, see `08-bug-audit.md`.

| ID | Severity | Verdict | Location | Title |
|---|---|---|---|---|
| B01 | critical | confirmed | `R/spatialCorrelation.R:837` | spatialCorrelationGeneExp applies p.adjust to one gene at a time, so no multiple-testing correction happens |
| B02 | high | confirmed | `R/spatialCorrelation.R:326` | Empirical p-value is extreme/B with a strict '>', so it can be exactly 0; BH keeps 0, and the docs claim the minimum is 0.01 |
| B03 | high | confirmed | `R/iterativePermutations.R:65` | spatialCorrelationGeneExpIterPermutations crashes ('subscript out of bounds') when any gene has NA p-values |
| B04 | high | partially-confirmed | `R/spatialCorrelation.R:109` | locfit C-stack-overflow segfault kills the R session when delta*N is in [1,2); deltas are not checked against N |
| B05 | medium | partially-confirmed | `R/spatialCorrelation.R:591` | Error handler uses corDF, which does not exist when cor.test itself fails ('object corDF not found') |
| B06 | medium | confirmed | `R/spatialCorrelation.R:830` | User-supplied BPPARAM is ignored by spatialCorrelationGeneExp |
| B07 | medium | confirmed | `R/iterativePermutations.R:64` | Iterative rerun threshold is (alpha/nPermutes)*100, not the documented alpha/nPermutations[k] |
| B08 | medium | confirmed | `R/spatialCorrelation.R:242` | set.seed() inside package functions overwrites the user's global RNG stream |
| B09 | medium | confirmed | `R/spatialCorrelation.R:827` | deltaX/deltaY given as a numeric vector are silently used one element per gene |
| B10 | medium | confirmed | `R/spatialCorrelation.R:785` | Pixels are matched by name only: separately rasterized inputs are silently mis-paired, and stale spatialCoords rownames crash |
| B11 | medium | confirmed | `R/spatialCorrelation.R:822` | A gene missing from the second object aborts the whole run (genes are not intersected) |
| B12 | medium | partially-confirmed | `R/spatialCorrelation.R:534` | A single NA/NaN silently disables the permutation test; spatialSimilarity errors on NA |
| B13 | medium | confirmed | `R/spatialCorrelation.R:977` | spatialCorrelationGeneExpWithinSample documents deltaX/deltaY but its argument is delta; 1-gene input fails |
| B16 | medium | confirmed | `R/visualizationFunctions.R:395` | savePlots non-geometry branch plots sample 1's expression in both panels |
| B17 | medium | confirmed | `R/visualizationFunctions.R:422` | spatialSimilarity does not store its assay, and savePlots never forwards assayName, so a different assay is plotted than was classified |
| B14 | low | partially-confirmed | `R/data.R:97` | R CMD check ERROR: the simRanPatternRasts example calls assays() without a namespace |
| B15 | low | confirmed | `R/spatialCorrelation.R:969` | spatialCorrelationGeneExpWithinSample example passes a list of 3 objects and assay 'A', so it errors |
| B18 | low | partially-confirmed | `R/visualizationFunctions.R:351` | savePlots calls library() on undeclared gridExtra/patchwork, attaching them to the user's session and masking BiocGenerics::combine |
| B19 | low | confirmed | `R/packageFunction.R:249` | numPixelInThresh is dim(thresh)[1], which is always 1, for genes skipped by minPixels |
| B20 | low | confirmed | `R/packageFunction.R:17` | getGenePixelDF default assayName = assayName is self-referential |
| B21 | low | partially-confirmed | `R/spatialCorrelation.R:512` | spatialCorrelation documents '1 x N matrix' input for X and Y but errors on it |
| B22 | low | confirmed | `DESCRIPTION:26` | Undeclared dependencies and imports (sf, stats/utils functions, vignette packages); dplyr version floor missing |
| B23 | low | confirmed | `NAMESPACE:1` | exportPattern('^[[:alpha:]]+') exports internal helpers; NAMESPACE is not managed by roxygen |
| B24 | low | confirmed | `R/iterativePermutations.R:302` | Non-ASCII em dash in R code causes an R CMD check WARNING |
| B25 | low | confirmed | `R/packageFunction.R:162` | Rd problems: undocumented 'verbose' argument (WARNING) and 89 'Lost braces' NOTEs |
| B26 | low | partially-confirmed | `inst/extdata/brain-MERFISH-10x-visium/brainCorrelation_1.RData:1` | Package is 48.6 MB as a tarball (33.8 MB installed) and contains a byte-identical duplicate RData |
| B27 | low | partially-confirmed | `R/spatialCorrelation.R:760` | Other example defects: nThreads=5 fails under the CRAN/Bioc core limit, a 3.6-minute example, undefined negCorrelation, and the wrong sample in the savePlots example |
| B28 | low | confirmed | `vignettes/getting-started-with-STcompare.Rmd:332` | Getting-started vignette: precomputed simRanPattern empirical p-values are called 'BH-corrected' and 'max of pX/pY', but are neither |
| B29 | low | confirmed | `vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:419` | Kidney vignette: nThreads described as 'parallelize genes across threads', and spatialCorrelationGeneExp described as returning adjusted p-values |
| B30 | low | confirmed | `R/packageFunction.R:51` | Assorted documentation errors (pseudo-count, threshold semantics, percentages, copied description, dataset size) |
| B31 | low | partially-confirmed | `R/spatialCorrelation.R:831` | Same seed for every gene and both directions: shared permutations and noise, order-dependent p-values, and stage 2 recomputing stage 1 |
| B32 | low | partially-confirmed | `R/spatialCorrelation.R:97` | matchingVariograms ignores i and defaults to seed = 0, so every permutation gets identical noise when called as in its own example |
| B33 | low | partially-confirmed | `R/spatialCorrelation.R:1145` | plotCorrelationGeneExp drops negative values, mishandles NA in one direction, and prints an unrounded p-value |
| B34 | low | confirmed | `R/packageFunction.R:269` | spatialSimilarity with negative values gives proportions that sum above 1 and NA pixel IDs; no shared pixels gives NaN silently |
| B35 | low | partially-confirmed | `R/spatialCorrelation.R:586` | Input-validation gaps, and errors print()ed to stdout instead of raised as warnings |
| B36 | low | partially-confirmed | `R/packageFunction.R:25` | getGenePixelDF densifies the entire assay twice per gene (quadratic cost); output grown with rbind in a loop |
| B37 | low | partially-confirmed | `R/spatialCorrelation.R:298` | Redundant bplapply extraction passes, a new MulticoreParam per call, and a hard-coded, undocumented variogram subsample (N_s = 1000) |

## B01: spatialCorrelationGeneExp applies p.adjust to one gene at a time, so no multiple-testing correction happens

- **Location:** `R/spatialCorrelation.R:837`
- **Category:** statistical
- **Severity:** auditor critical, verifier critical
- **Verdict:** confirmed

**Description:** p.adjust() is called on output$pValuePermuteX and output$pValuePermuteY inside the per-gene lapply (lines 808-845). Each call sees a single p-value, so BH, Bonferroni and the other methods all return it unchanged. The returned p-values are raw, even though the adjustMethod argument (default "BH") is documented as correcting the final columns (lines 715-718), and the kidney vignette says this function 'returns two adjusted empirical p-values'. adjustMethod is also not validated up front: an invalid value only errors after the first gene's permutations have run. spatialCorrelationGeneExpIterPermutations does this correctly, adjusting across all genes at iterativePermutations.R:343-344.

**Evidence (auditor):** WD/02_padjust_bpparam.R (6 genes, B=20). The table shows pX_BH = pX_none = pX_bonferroni = 0.00, 0.55, 0.25, 1, 1, 0.75. The expected BH values are 0, 1, 0.75, 1, 1, 1. 'identical(BH, none): TRUE  identical(bonferroni, none): TRUE'. adjustMethod='BHH' errors only after gene 1 has finished.

**Verifier reasoning:** I tried to refute this finding and could not. Every part of it holds up.

(1) How the bug happens. The p.adjust calls at R/spatialCorrelation.R:836-842 sit inside the per-gene lapply (808-845). They act on `output`, which spatialCorrelation() builds as a one-row data.frame (comment at :550-551, construction at :553-577). Each call therefore sees a single p-value. stats::p.adjust returns any length-1 input unchanged: its source has `if (n <= 1) return(p0)`, and that check runs after `match.arg(method)`. I checked all 8 p.adjust.methods.

(2) Effect on results. With the default adjustMethod = "BH", spatialCorrelationGeneExp returns exactly the raw empirical p-values (extreme/B). The documentation says otherwise. @param adjustMethod at :715-718 says the correction is applied to the final pValuePermuteX/pValuePermuteY columns. vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:430 says "`STcompare::spatialCorrelationGeneExp` returns two adjusted empirical p-values". One nuance: that vignette's code chunk actually calls spatialCorrelationGeneExpIterPermutations, and its loaded kidneyCorrelation.RData is correctly adjusted. So the vignette's own numbers are fine, but its text tells readers this function adjusts.

(3) Validation. spatialCorrelationGeneExp does not check adjustMethod before starting. An invalid value only fails, through match.arg inside p.adjust, after gene 1's permutations have finished. spatialCorrelationGeneExpIterPermutations checks at iterativePermutations.R:242-244 and adjusts across all genes at :343-344, as the finding says.

(4) Shipped results. The package ships output from this function. inst/scripts/KidneyNoIter.R:208 runs it on 1046 genes with the default method. The saved inst/extdata/kidneyCorrelationNoIter.RData (added in commit 338ab08 on 2026-04-28, after the MHT commit 002b151 on 2026-04-22) holds p-values that are all exact multiples of 1/100, meaning they were never adjusted.

Severity: the defaults silently return p-values without multiple-testing correction, while the argument and docs promise BH. In normal multi-gene use this inflates the false discovery rate, which matches 'critical' on the given scale. How much the calls change depends on how many genes are truly null. It is small for gene sets pre-filtered to spatially variable genes (18 of 757 calls in the shipped kidney data) and much larger for null-heavy gene sets. The suggested fix is sound: adjust once across genes after rbind, and validate up front with match.arg(adjustMethod, p.adjust.methods). Note for whoever applies it: the fix will change results from scripts like KidneyNoIter.R, and users who currently run p.adjust themselves afterwards would end up adjusting twice. It deserves a NEWS entry.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE); library(SpatialExperiment)
data(simRanPatternRasts)
px <- Reduce(intersect, lapply(simRanPatternRasts[1:6], colnames))
mk <- function(ii){ m <- do.call(rbind, lapply(ii, function(i) as.numeric(assay(simRanPatternRasts[[i]])[1, px]))); dimnames(m) <- list(paste0('g',1:3), px); SpatialExperiment(assays=list(pixelval=m), spatialCoords=spatialCoords(simRanPatternRasts[[1]])[px,]) }
inp <- list(mk(1:3), mk(4:6)); dl <- rep(list(c(0.2,0.5)),3)
a <- spatialCorrelationGeneExp(inp, nPermutations=20, deltaX=dl, deltaY=dl, verbose=FALSE, adjustMethod='bonferroni')
b <- spatialCorrelationGeneExp(inp, nPermutations=20, deltaX=dl, deltaY=dl, verbose=FALSE, adjustMethod='none')
identical(a$pValuePermuteX, b$pValuePermuteX)  # TRUE
```

**Suggested fix:** Remove the per-row p.adjust calls. After do.call(rbind, ...), apply p.adjust once to each of pValuePermuteX and pValuePermuteY across all genes. Validate the method with match.arg(adjustMethod, p.adjust.methods) at the start of the function.

## B02: Empirical p-value is extreme/B with a strict '>', so it can be exactly 0; BH keeps 0, and the docs claim the minimum is 0.01

- **Location:** `R/spatialCorrelation.R:326`
- **Category:** statistical
- **Severity:** auditor high, verifier high
- **Verdict:** confirmed

**Description:** The code computes extreme <- sum(abs(cor.global) > abs(cor.global.obs)) and then p.value.global <- extreme / B (lines 326-327). It has no +1 correction and uses a strict inequality, so a gene with zero exceedances gets p = 0. Under an exchangeable null this happens with probability about 1/(B+1). p.adjust leaves 0 unchanged, so those null genes stay significant at any FDR, whatever the number of tests. The documentation says 'Default is 100, such that the smallest p-value is 0.01' (lines 169, 358, 647, 866). The valid estimator is (1 + #{|r_b| >= |r_obs|})/(B+1) (Phipson & Smyth 2010).

**Evidence (auditor):** The package examples return pValuePermuteX = pValuePermuteY = 0 for kidney A-B and A-C (STcompare-Ex.Rout). Null simulation (WD/04_pvalue_zero.R): 50 independent simRanPatternRasts pairs, B=19, default deltas. 'fraction pX == 0: 0.12  fraction pY == 0: 0.06  fraction both == 0: 0.06'; the expected rate is 0.05. 'BH-adjusted max(pX,pY) for pairs with p == 0: 0 0 0', and with the +1 rule these become 0.833. p.adjust(c(0, rep(1,9999)), 'BH')[1] = 0. In the shipped inst/extdata/kidneyCorrelation.RData, 444 of 1046 genes have adjusted pX == 0 and 414 have both adjusted p-values == 0.

**Verifier reasoning:** I tried to refute this and could not. Every factual claim checks out.

(1) The code: R/spatialCorrelation.R:326-327 computes `extreme <- sum(abs(cor.global) > abs(cor.global.obs)); p.value.global <- extreme / B`. This is the only p-value formula in the package (grep found no other). It feeds spatialCorrelation, spatialCorrelationGeneExp, spatialCorrelationGeneExpWithinSample and spatialCorrelationGeneExpIterPermutations.

(2) The docs: lines 169, 358, 647 and 866 claim the smallest p-value is 0.01. So do man/viladomatCorrelation.Rd:36, man/spatialCorrelation.Rd:34, man/spatialCorrelationGeneExp.Rd:37, man/spatialCorrelationGeneExpWithinSample.Rd:28 and the built docs/reference/*.html. The page docs/reference/spatialCorrelationGeneExp.html makes this claim and also prints example output with pValuePermuteX = pValuePermuteY = 0.

(3) p = 0 really happens. I reproduced it on the package's own kidney example.

(4) Under a calibrated null it happens at about 1/(B+1). This holds in theory and in my null simulations. In those simulations the exceedance counts look uniform, so the surrogate null is roughly calibrated and the missing +1 is the defect.

(5) BH and Bonferroni both leave 0 at 0.

(6) The shipped-data counts (444/1046 and 414) match exactly.

(7) Phipson & Smyth's (1+b)/(B+1) is the standard valid Monte Carlo p-value.

Context that sharpens the finding but does not overturn it:
- The strict '>' adds nothing in practice. I saw no ties in 950 null pairs or in the kidney runs. The real defect is the missing +1 in the numerator and denominator.
- In spatialCorrelationGeneExp, p.adjust runs inside the per-gene lapply on a length-1 vector (lines 837-842). It is therefore a no-op for every value, not just 0. This is a separate bug. In spatialCorrelationGeneExpIterPermutations BH runs across genes and does keep zeros at 0.
- For a single test at B=100 with a strict '< 0.05', ex/B and (1+ex)/(B+1) give identical decisions (both 5/101). The inflation shows with '<= 0.05' (6/101) and at alpha = 0.01 (2/101 vs 1/101).
- The impact on the authors' shipped real-data BH calls is small. Hundreds of genes are true positives, so p = 0 there is an invalid floor rather than a false positive. The impact is large under Bonferroni or FWER control, and in null-heavy analyses, where BH cannot remove the roughly 1/(B+1) of null genes that land on 0.
- The authors work around zeros by hand. inst/scripts/simRanPatternSpatialCorrelation.R:68 replaces 0 with 0.01 before saving the shipped null results. inst/scripts/biological-replicates-example.R:350, 687, 714, 727 and inspect_kidney.R:73 map 0 to 0.001 for plotting.

Severity: high rather than critical. The p-values are invalid in common use, and FDR/FWER control fails in null-heavy or Bonferroni settings. But typical single-test decisions at B=100 and alpha = 0.05 are unchanged, and BH calls on the shipped real data change by only 0-4%.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
data(speKidney)
rk <- SEraster::rasterizeGeneExpression(speKidney, assay_name='counts', resolution=0.2, fun='mean', square=FALSE, BPPARAM=BiocParallel::SerialParam())
r <- spatialCorrelationGeneExp(list(rk$A, rk$B), verbose=FALSE)
r$pValuePermuteX  # 0 (docs: min 0.01)
p.adjust(c(0, rep(1, 9999)), 'BH')[1]  # still 0
```

**Suggested fix:** Use p = (1 + sum(abs(null) >= abs(obs))) / (B + 1) and correct the documentation. The C++ port should use the same formula.

## B03: spatialCorrelationGeneExpIterPermutations crashes ('subscript out of bounds') when any gene has NA p-values

- **Location:** `R/iterativePermutations.R:65`
- **Category:** correctness
- **Severity:** auditor high, verifier high
- **Verdict:** confirmed

**Verified description:** The core claim is accurate. Corrections on scope: (1) The trigger is an exactly constant gene on the shared pixels (for example all zeros in one sample), any NA, NaN or Inf value, or any runtime error inside the permutation code, including a transient error in a non-significant gene. Merely sparse, non-constant genes do not trigger it: 0/36 runs gave NA with as little as 1 of 277 pixels non-zero. (2) A crash requires length(nPermutations) >= 2, which includes the default c(100, 1000), and an NA in a non-final round. One NA gene is enough on its own. (3) The vignettes partly guard against the constant-gene case by recommending gene filters, but those filters are computed on all pixels, not shared pixels, are not applied to the cell-type run, and do not cover NA/Inf values or worker errors.

**Evidence (auditor):** WD/06b_iter_na.R: .get_genes_to_repermute(...) returns [1] "g1" NA. The traced rounds are 'g1,g2,g3' then 'g1,NA', followed by 'ERROR: subscript out of bounds'. The same data run through spatialCorrelationGeneExp completes, with an NA row for g3.

**Verifier reasoning:** I read the cited code and reproduced the crash independently with my own 3-gene data. The repository was not modified; git status is clean.

Mechanism, checked against the code:
(1) At R/iterativePermutations.R:65-66, `keep <- pX < t & pY < t` gives NA when both p-values are NA, or when one is NA and the other is below t. `rownames(results_df)[keep]` then returns NA.
(2) Round 2 calls .run_spatial_correlation_iteration with an NA gene. `match()` returns NA, and line 38, `SummarizedExperiment::assays(source)[[assayName]][gene, shared_pixels]`, fails with 'subscript out of bounds'. conditionCall() points to that exact expression. The message is the same for a base matrix (what SEraster produces for speKidney) and for a dgCMatrix (simRanPatternRasts).
(3) No partial results are returned, so all round-1 work is lost. Any round-2 genes ordered before the NA are also wasted, and round 2 is the expensive one (1000 permutations by default).
(4) The NA p-values come from the catch-all tryCatch in spatialCorrelation (R/spatialCorrelation.R:531, handler at 585-618). For an all-zero target gene, the reverse direction raises 'NA/NaN/Inf in foreign function call (arg 4)'.

The finding's own evidence matched my run line for line: rounds 'g1,g2,g3' then 'g1,NA', then the error. spatialCorrelationGeneExp on the same data completes and gives an NA row for g3. Both suggested fixes work: patching the helper with `which()` or with `%in% TRUE` makes the run finish, with an NA row for g3 and p.adjust handling the NA.

Scope details the finding leaves out (none changes the verdict):
(a) Only exactly constant genes trigger it. A sweep with the default delta grid and 100 permutations gave 0/36 NA results even with a single non-zero pixel out of 277, placed randomly or clustered. So 'common in sparse data' is true for genes that are entirely zero on one sample's shared pixels, not for genes that are merely sparse.
(b) A single NA, NaN or Inf value in either sample triggers it, and so does any transient error in any permutation. A simulated error in one permutation of a non-significant gene (p about 0.7) crashed the run.
(c) One NA gene on its own is enough: a single-gene input crashes, and so does a set where no other gene passes the screen.
(d) It needs at least 2 rounds in nPermutations and an NA in a non-final round. The default c(100, 1000) qualifies; with nPermutations = 10 (one round) the run completes with an NA row.
(e) If one p-value is NA and the other is at or above t, there is no crash. In practice this does not happen, because the catch-all sets both to NA.

Severity: high. This is a hard crash in the main exported entry point. All vignettes and inst/scripts use it, and it fails with default arguments, losing hours of computation (the brain vignette reports 1.78 h). The documentation softens this: the kidney vignette warns that zero-expression genes 'do not work' and recommends filtering genes to more than 5% (kidney) or 1% (brain) non-zero pixels. However, that filter is computed on all pixels rather than shared pixels. It is not applied to the cell-type proportion run, and it does not cover NA/NaN/Inf values or worker errors. The non-iterative function handles the same input gracefully, so users would not expect a crash here. It is not critical, because nothing is silently wrong.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
STcompare:::.get_genes_to_repermute(data.frame(pValuePermuteX=c(0,NA), pValuePermuteY=c(0,NA), row.names=c('g1','g2')), alpha=0.05, nPermutes=100)  # 'g1' NA
# Full repro: build 3-gene src/tgt (see B01 helper), set assay(tgt)['g3',] <- 0, then
# spatialCorrelationGeneExpIterPermutations(list(src,tgt), nPermutations=c(10,20), deltaX=dl, deltaY=dl)  -> Error: subscript out of bounds
```

**Suggested fix:** Use keep <- which(pX < t & pY < t), or keep <- (pX < t & pY < t) %in% TRUE. Handle constant genes and NA values explicitly before running permutations, and report them with a status column.

## B04: locfit C-stack-overflow segfault kills the R session when delta*N is in [1,2); deltas are not checked against N

- **Location:** `R/spatialCorrelation.R:109`
- **Category:** correctness
- **Severity:** auditor high, verifier high
- **Verdict:** partially-confirmed

**Verified description:** matchingVariograms (R/spatialCorrelation.R:109-111) passes each delta to locfit::lp(nn=) without checking it against N, the number of shared pixels. locfit converts nn to an integer neighbour count k = (int)(N*nn + 1e-12). When k == 1 (1 <= N*delta < 2), unbounded recursion in locfit's adaptive-tree builder (atree_grow, ~128 KiB per stack frame) overflows the C stack. R prints 'Error: segfault from C stack overflow' and jumps to top level, bypassing every tryCatch error handler, including spatialCorrelation's (only finally/on.exit run). Under Rscript or R CMD BATCH the process halts with exit status 1. In an interactive session the session survives, but the entire top-level call (for example all genes in spatialCorrelationGeneExp) is aborted and returns nothing. With nThreads > 1, the forked worker dies, BiocParallel raises a cryptic error ('wrong args for environment subassignment', or later "cor(permutations, Y): 'x' must be numeric"), and spatialCorrelation prints it and returns NA p-values. In spatialCorrelationGeneExpIterPermutations (the kidney vignette's call with MulticoreParam), those NAs then make .get_genes_to_repermute return NA gene names, and round 2 fails with 'subscript out of bounds', losing all results. Shared pixels are computed once per dataset pair, so every gene of an affected pair fails. Affected N: default grid seq(0.1,0.9,0.1) for N = 2-19; vignette/inst-scripts grid c(0.01,0.05,seq(0.1,0.9,0.1)) for N = 2-39 and 100-199. For example, the built-in speKidney data rasterized at resolution 0.25 (N=188) with the vignette grid crashes, while the documented resolution 0.2 (N=277) works. For very small N, locfit occasionally runs out of vertex space first and raises a catchable 'newsplit: out of vertex space' error instead. Fix: reject or drop (with a warning) any delta with floor(N*delta) < 2, which also removes the k=0 deltas whose smoothing collapses to a constant.

**Evidence (auditor):** Each call ran in a fresh process (WD/03b_segfault.R, WD/03c_segfault_big.R). Crashes: (N, nn) = (8, 0.2), (10, 0.1), (12, 0.1), (15, 0.1), (100-199, 0.01), (30-39, 0.05), (1000-1500, 0.001). No crash: (8, 0.1), (20, 0.1), (200, 0.01), (40, 0.05), (2000, 0.001). Through the package (WD/03d_edge.R), spatialCorrelation with N=150 and the vignette deltas gives 'Error: segfault from C stack overflow / Execution halted' with nThreads=1. With nThreads=2 it prints 'wrong args for environment subassignment' and returns pValuePermuteX/Y = NA. N=210 runs fine. The original 03_errors.R run also died at N=10 with the default deltas.

**Verifier reasoning:** The core bug is real and reproduces exactly. Nothing validates delta against the number of shared pixels N; line 109 is the locfit call. locfit (src/locfit.c:355) turns nn into an integer neighbour count k = (int)(N*nn + 1e-12). When k == 1, the adaptive-tree bandwidth at each vertex is the distance to its nearest data point. That distance shrinks as cells close in on a point, so atree_split (ev_atree.c:55-78) never stops splitting and atree_grow recurses (ev_atree.c:83-124). Each frame holds a ~128 KiB array (Sint nce[1<<MXDIM], MXDIM=15), so R's ~7.6 MB C stack overflows after about 60 levels, before the vertex-space guard (ev_main.c:202) fires. With k == 0 every vertex has h = 0 and a size-based fallback score ends the recursion, so N*delta < 1 does not crash. Confirmed as claimed: the [1,2) band with exact boundaries, every listed (N, nn) point, the package-level N=150 vs N=210 results, the nThreads=2 message and NA p-values, and that no tryCatch error handler can intercept it. Details that are wrong or incomplete: (1) 'kills the R session / exits' is true only for non-interactive R (Rscript, R CMD BATCH). R's SIGSEGV handler jumps to top level, so in an interactive session the whole top-level call is aborted (all genes lost, only finally/on.exit run) but the session survives and stays usable. (2) With nThreads>1 the result is not quite silent: the package prints a cryptic simpleError that does not name the cause. In spatialCorrelationGeneExpIterPermutations, the function the kidney vignette actually calls with MulticoreParam, the NA p-values make .get_genes_to_repermute (R/iterativePermutations.R:62-66) return NA gene names, and round 2 fails with 'subscript out of bounds', losing all results. (3) The ranges are incomplete. The default grid crashes for N = 2-19, not only 10-19, because delta 0.2 and up also give k=1 when N < 10 (their own (8, 0.2) data point shows this). The vignette grid crashes for N = 2-39 and 100-199. (4) N is shared pixels per dataset pair (R/spatialCorrelation.R:785, R/iterativePermutations.R:266), not per gene, so every gene of an affected pair fails. (5) For very small N, locfit sometimes runs out of vertex space first and raises a catchable R error instead (2 of 234 in-band runs). Severity is high, not medium. The trigger is easy to hit with documented settings: the package's own speKidney data rasterized at resolution 0.25 instead of 0.2 (N=188), with the vignette/inst-scripts delta grid, kills an Rscript process. The shipped examples (N=273-288) and the vignette's own kidney run (0.01 was selected, 0 NA in 1046 genes) are unaffected, and the default grid only fails for N <= 19. It is not critical because the results are not silently wrong.

```r
# run in a throwaway R session; it will die
set.seed(150); n <- 150; x <- runif(n); y <- runif(n); v <- rnorm(n)
locfit::locfit(v ~ locfit::lp(x, y, nn = 0.01, deg = 0), kern = 'gauss', maxk = 300)
# Error: segfault from C stack overflow
# via package: spatialCorrelation(v, rnorm(n), cbind(x,y), nPermutations=2, deltaX=c(0.01,0.05,seq(0.1,0.9,0.1)), deltaY=c(0.01,0.05,seq(0.1,0.9,0.1)))
```

**Suggested fix:** Before smoothing, drop (with a warning) or clamp every delta with delta*N < 2 (or a safer minimum such as 3/N), and check 0 < delta <= 1. A native smoother in a C++ rewrite should enforce a minimum neighbourhood size.

## B05: Error handler uses corDF, which does not exist when cor.test itself fails ('object corDF not found')

- **Location:** `R/spatialCorrelation.R:591`
- **Category:** correctness
- **Severity:** auditor medium, verifier medium
- **Verdict:** partially-confirmed

**Verified description:** In spatialCorrelation() (R/spatialCorrelation.R:531-620), corDF is assigned only by the first statement of the tryCatch (line 534). If cor.test() itself errors, corDF is never bound. The handler's corDF$estimate and corDF$p.value (lines 591-592 and 605-606) then raise "object 'corDF' not found". The same tryCatch does not catch handler errors, so the error escapes. spatialCorrelationGeneExp() has no per-gene protection (lines 808-845), so the whole run aborts and all genes already completed are lost, instead of an NA row being returned for the bad gene. The real cause is printed to stdout by print(cond) at line 586, but the condition raised to the caller is the misleading 'corDF not found'.

cor.test fails when there are fewer than 3 complete (x, y) pairs:
- n < 3;
- a gene that is all NA or NaN, or nearly so;
- no shared pixels (character(0));
- both SPEs' spatialCoords() lack rownames. shared_pixels comes from rownames(spatialCoords()) at lines 785-786, and intersect(NULL, NULL) is NULL. This happens with SpatialExperiment(spatialCoords = <unnamed matrix>) even when colnames exist.

Missing colnames alone is NOT a trigger: with named coords, the call fails earlier with 'subscript out of bounds' at line 821. Constant or all-zero genes are also NOT triggers: cor.test only warns and returns NA, so the handler works and returns an NA row.

Extra consequence: corDF is resolved lexically up to .GlobalEnv. If the user has a global object named corDF, it is used silently. A list with estimate and p.value gives fabricated correlationCoef and pValueNaive values; an unrelated data.frame gives yet another misleading error ('arguments imply differing number of rows: 0, 1').

Fix: initialise corDF (for example, list(estimate = NA_real_, p.value = NA_real_)) before the tryCatch. Better, validate inputs up front: at least 3 complete shared pixels, and non-NULL pixel names.

**Evidence (auditor):** WD/03_errors.R. For n = 2, an all-NA X, empty vectors, and spatialCorrelationGeneExp on unnamed SpatialExperiments, the output is '<simpleError in cor.test.default(...): not enough finite observations>' followed by '"ERROR: object 'corDF' not found"'.

**Verifier reasoning:** The core bug is real and reproduces exactly as described. In spatialCorrelation(), corDF is assigned by the first statement of the tryCatch (R/spatialCorrelation.R:534). When cor.test() itself errors, corDF is never bound in the function frame. The handler then evaluates corDF$estimate and corDF$p.value (lines 591-592 and 605-606), which raises "object 'corDF' not found". A tryCatch does not catch errors raised inside its own handler, so this error escapes spatialCorrelation(); its call is value[[3L]](cond), i.e. the handler. spatialCorrelationGeneExp() has no per-gene protection (do.call(rbind, lapply(...)), lines 808-845), so the whole multi-gene run aborts and every gene already computed is discarded. I confirmed this with a 3-gene run.

The line numbers are correct, and the suggested fix works (a patched copy returns an NA row).

One detail in the description is wrong. It says "SpatialExperiments without colnames (so shared_pixels is NULL)", but shared_pixels is built from rownames(spatialCoords(.)) at lines 785-786 and does not depend on colnames:
- An SPE with no colnames but with named coords does NOT hit this bug. It fails earlier with 'subscript out of bounds' at line 821.
- An SPE WITH colnames but whose spatialCoords lack rownames DOES hit it: intersect(NULL, NULL) is NULL, which yields 0-length X and Y. This happens with the common constructor call SpatialExperiment(spatialCoords = as.matrix(df)).

Nuance on the 'misleading message': print(cond) at line 586 writes the real cause ('not enough finite observations') to stdout first, but the condition raised to the caller is the misleading one.

Additional consequence not in the finding: because corDF is resolved lexically (STcompare namespace -> imports -> base -> R_GlobalEnv), a user object named corDF in the global environment is used silently. A list containing estimate and p.value produces fabricated correlationCoef and pValueNaive values for an input that cor.test rejected.

Severity is medium, not high:
- The documented SEraster workflow is unaffected. The built-in data have no NAs and 277 shared pixels.
- The most common degenerate case in real data, a constant or all-zero gene, does not trigger the bug: cor.test only warns and returns NA, corDF is bound, and the handler returns an NA row.
- The triggers are edge cases: fewer than 3 complete pixel pairs (n < 3, or NA/NaN-laden genes), no shared pixels, or user-built SPEs with unnamed coordinate rows.
- The impact is a loud abort with a misleading message, which can throw away a long multi-gene run. Silently wrong values occur only when a global object named corDF exists.

Side observation, separate from B05: an n = 3 input kills the R session with 'segfault from C stack overflow' inside locfit (nn = 0.5). No tryCatch catches it, including the package's own.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
spatialCorrelation(c(1,2), c(2,1), cbind(1:2, 1:2), nPermutations=3, deltaX=0.5, deltaY=0.5)
# Error: object 'corDF' not found
```

**Suggested fix:** Initialise corDF <- NULL before the tryCatch and use NA_real_ in the handler when it is NULL. Better, validate up front that there are at least 3 complete, shared, non-constant pixels, and return a row with a status message.

## B06: User-supplied BPPARAM is ignored by spatialCorrelationGeneExp

- **Location:** `R/spatialCorrelation.R:830`
- **Category:** api
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** spatialCorrelationGeneExp passes BPPARAM = NULL to spatialCorrelation(), which therefore always builds MulticoreParam(workers = nThreads). The BPPARAM created at lines 777-779 is never used. A user's SnowParam (needed on Windows, where MulticoreParam does not work), SerialParam or tuned MulticoreParam is silently discarded, contradicting the documentation at lines 699-703. spatialCorrelationGeneExpWithinSample and the iterative function both pass BPPARAM through correctly.

**Evidence (auditor):** WD/02_padjust_bpparam.R, with trace() on viladomatCorrelation. 'User passed SnowParam(workers = 3); viladomatCorrelation received: "MulticoreParam workers = 1" "MulticoreParam workers = 1"'. For comparison, WithinSample with SerialParam() receives 'SerialParam workers = 1'.

**Verifier reasoning:** I tried to refute this and could not. The cited code and my own runs both show the bug.

What the source shows:
- R/spatialCorrelation.R:777-779 replaces a NULL BPPARAM with MulticoreParam(workers = nThreads). After that the variable is never read again.
- Line 830 calls spatialCorrelation(..., nThreads = nThreads, BPPARAM = NULL, ...).
- spatialCorrelation() therefore always takes its own NULL branch (lines 507-509), builds MulticoreParam(workers = nThreads), and passes that to both viladomatCorrelation() calls (lines 541, 547).
- The sibling entry points do pass the argument through: spatialCorrelationGeneExpWithinSample (line 1033), and spatialCorrelationGeneExpIterPermutations via .run_spatial_correlation_iteration (iterativePermutations.R:49 and :322).
- git blame shows the hard-coded NULL has been there since the first commit (6bc3af5, original line 330). In that same commit WithinSample already passed BPPARAM = BPPARAM, so this looks like an oversight, not a deliberate choice.

The documentation it contradicts:
- Lines 699-703, as the finding cites.
- More explicitly, the nThreads entry at lines 695-697 (and man/spatialCorrelationGeneExp.Rd:85-87): "If BPPARAM argument is not NULL, the BPPARAM argument would override nThreads argument."

The bug also works in the opposite direction: SerialParam() together with nThreads = 4 still forks worker processes.

On Windows, per the installed BiocParallel 1.44.0 source, MulticoreParam() warns "MulticoreParam() not supported on Windows, use SnowParam()" and forces workers = 1. So a Windows user's SnowParam is dropped, they get serial execution plus a misleading warning, and they have no way to parallelize this function. On Windows the discard is therefore not fully silent, only misleadingly reported; on macOS and Linux it is silent. These are refinements, not errors in the finding.

The author's own inst/scripts/KidneyNoIter.R:208-216 calls this path with nThreads = 22 and BPPARAM = MulticoreParam(). The bug stays hidden there because nThreads is also set.

The suggested fix (BPPARAM = BPPARAM at line 830) is correct.

Why medium severity: statistical output is unaffected. Null distributions, p-values and delta* are identical across every backend I tried, because X.randomized is drawn in the parent process and matchingVariograms calls set.seed(seed + i) itself. Nothing crashes. What remains is a documented API contract that is false for a key behavior (parallelism) in a package whose main pain point is runtime: user-requested parallelism is lost (2.4x slower with 4 workers in my test), and Windows users cannot parallelize this function at all. That matches "misleading docs about key behavior" = medium.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
trace('viladomatCorrelation', quote(print(class(BPPARAM))), where=asNamespace('STcompare'))
data(simRanPatternRasts)
spatialCorrelationGeneExp(list(simRanPatternRasts[[1]], simRanPatternRasts[[2]]), nPermutations=2, deltaX=list(0.5), deltaY=list(0.5), BPPARAM=BiocParallel::SnowParam(3), verbose=FALSE)
# prints MulticoreParam
```

**Suggested fix:** Pass BPPARAM = BPPARAM at line 830.

## B07: Iterative rerun threshold is (alpha/nPermutes)*100, not the documented alpha/nPermutations[k]

- **Location:** `R/iterativePermutations.R:64`
- **Category:** docs-code-mismatch
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** The code uses t <- (alpha / nPermutes) * 100. The roxygen (lines 90-94) says genes are carried forward when both p-values are below alpha / nPermutations[k], which is 100 times smaller. The vignette says genes with p < 0.05 are rerun. The code matches the vignette only when nPermutations[1] == 100. With a different first stage it behaves very differently: c(10, 20) gives a threshold of 0.5, so every gene with p < 0.5 is rerun.

**Evidence (auditor):** WD/06_iter.R. nPermutes=10: code threshold 0.5 (reruns p0, p004, p03, p2, p45) vs documented 0.005. nPermutes=100: 0.05 vs 0.0005. nPermutes=1000: 0.005 vs 5e-05.

**Verifier reasoning:** Every factual claim checks out.

(1) Code vs roxygen: R/iterativePermutations.R:64 is `t <- (alpha / nPermutes) * 100`. The roxygen for @param alpha (lines 90-94) says the threshold is `alpha / nPermutations[k]`, which is 100 times smaller. The same wrong text appears in man/spatialCorrelationGeneExpIterPermutations.Rd:36-40 and in the pkgdown page docs/reference/spatialCorrelationGeneExpIterPermutations.html:96-99.

(2) Vignette: vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:383 and :407 say genes with an empirical p-value below 0.05 are rerun. The code's threshold equals 0.05 only when nPermutations[k] == 100 and alpha is the default 0.05.

(3) Origin: git shows commit 0a27a4c ("fixed iterativePermutation bug", 2026-04-28) changed `alpha / nPermutes` to `(alpha / nPermutes) * 100`. The roxygen, Rd and pkgdown docs were never updated.

(4) What the rule really means: spatialCorrelation.R:326-327 computes p = extreme/B, which is k/B. So the code's rule p < 100*alpha/B is the same as k < 100*alpha, i.e. at most 4 exceedances at the default alpha, for any B. The documented rule p < alpha/B is the same as k < alpha, i.e. zero exceedances only. As a p-value threshold the code's cutoff is therefore 5/B: 0.5 at B=10, 0.25 at B=20, 0.05 at B=100, 0.005 at B=1000. This matches the finding's evidence.

(5) The consequence is real when run end to end: with c(10,20), genes with stage-1 p-values up to 0.4 were rerun. One wording nit: it should say both pValuePermuteX and pValuePermuteY must be < 0.5, which the finding states earlier.

No floating-point off-by-one: 5/B is identical to t, so a p-value exactly at the cutoff is not rerun.

Severity: medium is reasonable. The help page gets the only effect of `alpha` wrong by a factor of 100 at every setting, including the default. That counts as misleading docs about key behavior, and it matters for the runtime of a function whose second stage is 10 times more expensive. It does not produce biased or invalid p-values, though. Genes that are not rerun keep valid stage-1 p-values. When B <= 100, a gene that fails the screen has raw p >= 0.05 in at least one direction, so it could never pass the 'both adjusted p < 0.05' rule anyway. Every in-repo caller uses c(100, 1000), where the code and the vignette agree.

Side notes, not part of B07: the output has no column saying which B produced each gene's p-values (you can only infer it from length(nullCorrelationsX)). p-values of exactly 0 can occur (k/B). p-values from different B levels are BH-adjusted together.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
STcompare:::.get_genes_to_repermute(data.frame(pValuePermuteX=0.3, pValuePermuteY=0.3, row.names='g'), alpha=0.05, nPermutes=10)  # 'g' (threshold 0.5; docs say 0.005)
```

**Suggested fix:** Choose one rule (for example, rerun when p < alpha, or when the number of exceedances is below a constant), then implement it and document it consistently in the roxygen and the vignette.

## B08: set.seed() inside package functions overwrites the user's global RNG stream

- **Location:** `R/spatialCorrelation.R:242`
- **Category:** reproducibility
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** viladomatCorrelation (line 242) and matchingVariograms (line 99) call set.seed() on the global RNG. After any correlation call, all of the user's later random draws are determined by the package's seed, not the user's. Simulations or bootstraps that interleave calls to STcompare silently lose their randomness. (Results are reproducible across parallel back-ends, which I verified; that is not the problem.)

**Evidence (auditor):** WD/05_rng.R. 'runif(3) after call with user seed 123: 0.5116576 0.08492106 0.9947705' and 'runif(3) after call with user seed 999: 0.5116576 0.08492106 0.9947705', identical. Null correlations were identical for SerialParam, MulticoreParam(1), (2) and (4).

**Verifier reasoning:** I tried to refute this and could not. The core claim, the line numbers and the stated consequence are all correct.

What the code does:
- R/spatialCorrelation.R:242 runs `set.seed(seed)` (default 0) in the calling process.
- Right after that, still in the calling process, it draws from the reseeded stream: line 255 `sample(N, N_s)` (only when N > 1000) and lines 286-288 (B calls to `sample(X, ...)`).
- spatialCorrelation calls viladomatCorrelation twice (lines 537 and 543), both times with the same `seed`.
- The RNG state left after any call is therefore a fixed function of (seed, N, B). The user's own seed has no effect on it.
- The functions document `seed`, but nowhere say they change the global RNG.

One refinement, which does not contradict the finding:
- matchingVariograms:99 also calls `set.seed()` on the global RNG, but inside spatialCorrelation it runs within `BiocParallel::bplapply` (lines 290-295).
- BiocParallel 1.44 saves and restores the caller's `.Random.seed` around that call on the serial path. That path is also what the default `MulticoreParam(workers = 1)` falls back to. On the fork path the parent's state is never touched.
- So line 99 adds nothing to the state that spatialCorrelation leaves behind. That state equals exactly `set.seed(0)` followed by the main-process `sample()` draws.
- Line 99 does change the user's global RNG when someone calls the exported `matchingVariograms()` directly, as its own roxygen example does.

The consequence is real and worse than 'loses randomness' may suggest:
- In a Monte-Carlo loop that simulates new data and then calls spatialCorrelation, every iteration from the second one on gets the exact same 'random' dataset.
- In a bootstrap loop, every resample from the second one on uses the same index vector.

All the wrappers are affected. Besides viladomatCorrelation and spatialCorrelation, I confirmed it for spatialCorrelationGeneExp and spatialCorrelationGeneExpIterPermutations.

The side claim that results are reproducible across back-ends is also true.

Severity is medium:
- The package's own p-values are not wrong, and the repo's own scripts do not mix their own random draws with calls (the only `set.seed(111111)` in inst/scripts is for plot colours).
- But a user who runs a simulation or bootstrap around STcompare gets silently degenerate results, and that side effect is undocumented.

The suggested fix is sound: save and restore `.Random.seed` (or use withr, which would need to be added to Imports), or use independent streams in a C++ backend.

Related but separate issue: because the same seed is reused, the forward and reverse directions (and every gene) use identical shuffle index vectors at lines 286-288. That makes the X-permutation and Y-permutation nulls dependent on each other, which may be worth its own finding.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
set.seed(1); x <- rnorm(50); y <- rnorm(50); p <- cbind(runif(50), runif(50))
set.seed(123); invisible(spatialCorrelation(x, y, p, nPermutations=3, deltaX=.5, deltaY=.5)); a <- runif(1)
set.seed(999); invisible(spatialCorrelation(x, y, p, nPermutations=3, deltaX=.5, deltaY=.5)); b <- runif(1)
identical(a, b)  # TRUE
```

**Suggested fix:** Use withr::local_seed()/with_seed(), or save and restore .Random.seed. Preferably use independent per-permutation streams (dqrng or L'Ecuyer substreams) in the C++ backend, without touching global state.

## B09: deltaX/deltaY given as a numeric vector are silently used one element per gene

- **Location:** `R/spatialCorrelation.R:827`
- **Category:** correctness
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** spatialCorrelation() takes deltaX as a numeric vector (a grid of candidate deltas). spatialCorrelationGeneExp and the iterative helper (iterativePermutations.R:44-45) index it with deltaX[[i]], so passing the same kind of vector (for example c(0.1, 0.5, 0.9)) silently gives gene 1 only delta 0.1, gene 2 only 0.5, and so on. That changes the null distribution and the p-values. If there are more genes than elements, the run errors part-way. The list length is never checked against the number of genes.

**Evidence (auditor):** WD/03_errors.R. With deltaX = c(0.1, 0.5, 0.9) and 2 genes, the unique deltaStarX values per gene are 0.1 and 0.5 respectively. With deltaX = c(0.1, 0.5) and 3 genes: 'ERROR: subscript out of bounds'.

**Verifier reasoning:** Every factual claim reproduces, using my own input.

**Code.**
- spatialCorrelation.R:826-827 passes `deltaX = deltaX[[i]], deltaY = deltaY[[i]]`.
- iterativePermutations.R:29 sets `gene_idx <- match(gene, rownames(source))`, and lines 44-45 pass `deltaX[[gene_idx]]`. That is the same positional indexing.
- Only NULL gets the safe default `rep(list(seq(0.1,0.9,0.1)), nGenes)` (spatialCorrelation.R:791-797, iterativePermutations.R:272-278).
- Nothing in R/ checks the type or length of deltaX/deltaY against the number of genes. A grep for `length(delta`, `is.list(delta` and `stop(...delta` finds only the loop inside matchingVariograms.

**Why it is a trap.** The same argument name means different things in sibling exported functions:
- spatialCorrelation() documents deltaX as `numeric`, a single value or a vector grid (spatialCorrelation.R:360-371).
- The GeneExp and IterPermutations versions document a `list` of length nrow (spatialCorrelation.R:649-662, iterativePermutations.R:102-115). Yet they also say "a sequence of deltas can be inputted", and so does the kidney vignette prose (acute-kidney-injury-10x-visium-rasterized.Rmd:381).
- An atomic vector behaves exactly like a "list of single numerics": `list(0.1,0.5)` and `c(0.1,0.5)` gave identical results.

**Mechanism of the part-way error.** The `deltaX[[i]]` promise is forced by `is.null(deltaX)` at spatialCorrelation.R:523. That is before the `tryCatch` at line 531, so spatialCorrelation's per-gene NA fallback never sees the error. The whole lapply aborts, and the results for genes already computed are thrown away.

**Context that limits severity.**
- Every in-repo caller passes a correctly shaped list: the vignette and inst/scripts all use `deltaList`. So the trigger is a type the docs do not describe, not the documented path.
- Results are silently wrong only when nGenes <= length(vector). That covers any single-gene input, which is what the package's speKidney examples and simRanPattern validation use, and small gene panels. With more genes than elements, the run fails late with an uninformative 'subscript out of bounds'.
- The same problem exists in spatialCorrelationGeneExpWithinSample (`delta`, spatialCorrelation.R:996-998 and 1029). There, a 9-element vector with 4 genes silently gives the genes 0.1, 0.2, 0.3 and 0.4. A 3-element vector with 4 genes fails at `names<-`. Its roxygen also documents deltaX/deltaY, but the actual argument is `delta`.

**Size of the effect.** When it triggers, it can badly inflate false positives. The vignette's own grid `c(0.01,0.05,0.1..0.9)` passed as a vector silently fixes delta at 0.01 for a single-gene input. On 50 independent random-field pairs, the false-positive rate went from 0.08 to 0.42. The naive test, with no autocorrelation correction at all, gives 0.48.

**Severity.** I rate it medium: wrong results in an edge case where the input contradicts the documented type. It is at the high end of medium, because the result is silent, the effect is large, and the API is inconsistent across the exported functions.

**Suggested fix.** The proposed fix is reasonable. Treat an atomic vector as a grid shared by all genes. Otherwise require a list of length nGenes and stop with a clear message. Validate before the gene loop, so no compute is wasted. Applying the same fix to WithinSample's `delta` is advisable.

```r
# using inp from B01
r <- spatialCorrelationGeneExp(lapply(inp, function(s) s[1:2,]), nPermutations=5, deltaX=c(0.1,0.5,0.9), deltaY=c(0.1,0.5,0.9), verbose=FALSE)
lapply(r$deltaStarX, unique)  # 0.1 and 0.5: one delta per gene
```

**Suggested fix:** If deltaX is an atomic numeric vector, treat it as a shared grid (rep(list(deltaX), nGenes)). Otherwise require a list of length nGenes (or one named by gene) and stop with a clear message.

## B10: Pixels are matched by name only: separately rasterized inputs are silently mis-paired, and stale spatialCoords rownames crash

- **Location:** `R/spatialCorrelation.R:785`
- **Category:** correctness
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** shared_pixels is computed as intersect(rownames(spatialCoords(source)), rownames(spatialCoords(target))), and pos is taken from the source only. The coordinates of same-named pixels are never compared. The same pattern appears in iterativePermutations.R:266-268, plotCorrelationGeneExp, savePlots and pixelClass. SEraster names pixels by grid index, so two objects rasterized in separate calls with different extents reuse names for different locations, and the correlation is computed on mis-paired pixels with no warning. Separately, these rownames are used to index assay columns, but SpatialExperiment does not update rownames(spatialCoords) on colnames<-, which leads to a crash (spatialSimilarity, which uses colnames, still works).

**Evidence (auditor):** WD/08b_pixels.R. With speKidney$B restricted to y > 1.5 and rasterized separately: 166 shared names, all 166 with different coordinates (max offset 4.291), giving correlationCoef -0.2909. Rasterized jointly: 205 shared names, 0 coordinate mismatches, correlationCoef -0.9503. In the second case, colnames are 'px170...' while rownames(spatialCoords) are 'pixel170...', and spatialCorrelationGeneExp fails with 'subscript out of bounds'.

**Verifier reasoning:** I tried to refute this and could not. Both parts reproduce with the exact numbers given.

(1) Pairing by name only. spatialCorrelation.R:785-787 builds shared_pixels from intersect(rownames(spatialCoords(source)), rownames(spatialCoords(target))) and takes pos from the source only. The same three lines are at iterativePermutations.R:266-268 and plotCorrelationGeneExp (spatialCorrelation.R:1137-1141). Nothing ever compares the coordinates of same-named pixels. SEraster 1.2.0 explains why this goes wrong. rasterizeMatrix names pixels paste0('pixel', grid index) on sf::st_make_grid(bbox). rasterizeGeneExpression uses one shared bbox when given a list, but when given a single object it uses that object's own floor/ceiling bbox. So two objects rasterized in separate calls give the same name to different locations. The functions then correlate mis-paired pixels and emit no warning (I captured none).

(2) Stale rownames crash. colnames<- on a SpatialExperiment updates colnames and colData rownames but not rownames(spatialCoords), and validObject still passes. spatialCorrelationGeneExp then indexes the assay with the stale names (lines 821-822) and fails with 'subscript out of bounds'. spatialCorrelationGeneExpIterPermutations (iterativePermutations.R:38-39) and plotCorrelationGeneExp fail the same way. spatialSimilarity, which uses colnames (packageFunction.R:20), still runs.

The description is accurate. Three nuances do not change the verdict:
(a) savePlots and pixelClass use intersect(rownames(spatialCoords)) only in their no-'geometry' branches (visualizationFunctions.R:221-224, 389-395). SEraster output always has a geometry column, so for it they take the colnames path instead. That path still matches by name only. savePlots crashes on stale rownames only without geometry; pixelClass did not crash in either branch.
(b) Mis-pairing needs the floor/ceiling-rounded bboxes to differ in a way that shifts the grid indexing, not just any difference in extent. Unmodified speKidney A/B/C rasterized separately pair correctly, as do some xmax/ymax-only changes. With real micron-scale coordinates, though, separate rasterization will almost always give different bboxes.
(c) The finding leaves out that spatialSimilarity, linearRegression and pixelClass are mis-paired by separate rasterization just as badly, because they also match on names (via colnames).

Related failure (not claimed): if spatialCoords has no rownames, which is what the SpatialExperiment constructor gives for a plain matrix, intersect(NULL, NULL) is NULL and spatialCorrelationGeneExp fails with a confusing "object 'corDF' not found".

Severity: medium. The vignettes and examples always rasterize jointly with one list call. The roxygen for @param input states the same-name-means-same-location rule, but it does not warn against separate SEraster calls. The result is silently wrong output when a user makes this mistake, and a confusing crash after renaming columns. Both are edge cases or misuse rather than the documented normal workflow.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE); library(SpatialExperiment)
data(speKidney); B <- speKidney$B[, spatialCoords(speKidney$B)[,'y'] > 1.5]
f <- function(x) SEraster::rasterizeGeneExpression(x, assay_name='counts', resolution=0.2, fun='mean', square=FALSE, BPPARAM=BiocParallel::SerialParam())
rA <- f(speKidney$A); rB <- f(B); j <- f(list(A=speKidney$A, B=B))
spatialCorrelationGeneExp(list(rA, rB), nPermutations=3, verbose=FALSE)$correlationCoef  # -0.29 (mis-paired)
spatialCorrelationGeneExp(list(j$A, j$B), nPermutations=3, verbose=FALSE)$correlationCoef  # -0.95
```

**Suggested fix:** Match on colnames, check that spatialCoords agree (within a tolerance) for the shared names and stop if they do not. Document that the objects must be rasterized together in a single SEraster::rasterizeGeneExpression(list(...)) call.

## B11: A gene missing from the second object aborts the whole run (genes are not intersected)

- **Location:** `R/spatialCorrelation.R:822`
- **Category:** correctness
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** spatialCorrelationGeneExp loops over rownames(source) and indexes target by gene name (lines 821-822; the same at iterativePermutations.R:38-39). This happens outside spatialCorrelation's tryCatch, so a gene absent from target throws 'subscript out of bounds' after all earlier genes have been computed. spatialSimilarity, by contrast, uses intersect(rownames(x), rownames(y)). Mismatched gene sets are common when comparing technologies or panels.

**Evidence (auditor):** WD/03_errors.R. Progress shows '1: g1', '2: g2', then '"ERROR: subscript out of bounds"'. For the same inputs, spatialSimilarity returns genes 'g1' 'g3'.

**Verifier reasoning:** Code reading. spatialCorrelationGeneExp loops over rownames(source) (spatialCorrelation.R:809-813) and indexes each assay by gene name (lines 821-822) with no check that the gene exists in target. The only tryCatch is inside spatialCorrelation (lines 531-620), which covers that function's body and nothing else. Neither the lapply nor do.call(rbind, ...) in spatialCorrelationGeneExp has error handling. So a gene missing from input[[2]] aborts the whole call, and the rows already computed are discarded. spatialCorrelationGeneExpIterPermutations does the same thing: genes_all <- rownames(source) at iterativePermutations.R:286, with indexing at lines 38-39. spatialSimilarity, by contrast, intersects the gene sets (packageFunction.R:193).

I reproduced every claim independently. The error is raised by exactly the expression on line 822. The number of genes fully computed before the abort equals the number of genes ahead of the missing one. Nothing is returned. The iterative function fails the same way. spatialSimilarity returns g1 and g3 for the same input. The claim that mismatched gene sets are realistic holds: the brain vignette says "We will only be able to compare genes that are shared across the two technologies" (vignettes/brain-MERFISH-10x-visium.Rmd:154), and both vignettes subset genes before calling.

Extra details that do not contradict the finding:
(a) The behavior is asymmetric. Genes present only in input[[2]] are silently ignored with no warning, so swapping the argument order turns the crash into a silent subset.
(b) Sparse dgCMatrix assays fail with the same 'subscript out of bounds' message, raised from Matrix's .subscript.2ary.
(c) Genes in a different order in target are handled correctly, because indexing is by name. The problem is limited to missing genes.
(d) The man/*.Rd pages for spatialCorrelationGeneExp and spatialCorrelationGeneExpIterPermutations do not state the precondition. Only the vignettes hint at it, e.g. the AKI vignette line 234: 'dimensions of the assays must be same'.
(e) spatialSimilarity also drops genes silently, without a warning.
(f) The fix suggests indexing deltaX/deltaY by gene name, but the default and documented delta lists are unnamed and positional (length = nrow(source)). The fix should map with match(g, rownames(source)) or accept named lists.

Severity is medium, not high. The failure is loud, not silently wrong. The packaged workflows avoid it by pre-intersecting genes. The real cost is that the error comes late, has an uninformative message, and throws away all earlier work, which can be hours on real panels: the brain vignette run took 1.78 h on 22 threads.

```r
# using inp from B01
spatialCorrelationGeneExp(list(inp[[1]], inp[[2]][c(1,3),]), nPermutations=3, deltaX=rep(list(.5),3), deltaY=rep(list(.5),3))  # Error: subscript out of bounds
```

**Suggested fix:** Compute shared_genes <- intersect(rownames(source), rownames(target)) once, and warn about any dropped genes. Index deltaX/deltaY by gene name.

## B12: A single NA/NaN silently disables the permutation test; spatialSimilarity errors on NA

- **Location:** `R/spatialCorrelation.R:534`
- **Category:** statistical
- **Severity:** auditor medium, verifier medium
- **Verdict:** partially-confirmed

**Verified description:** A single NA or NaN in X or Y makes spatialCorrelation report correlationCoef and pValueNaive from cor.test's complete-case analysis. A single Inf makes both NaN. In every case, pValuePermuteX, pValuePermuteY, deltaStar* and nullCorrelations* come back NA, and the error is only print()ed (spatialCorrelation.R:585-619). Any non-finite pixel therefore discards that gene's permutation test without a warning or status column.

The failing call is not locfit. locfit silently drops NA rows through na.omit and returns N-1 fitted values that are misaligned with the coordinates. The 'NA/NaN/Inf in foreign function call' error comes from geoR::variog:
- N <= 1000: at the target variogram (line 272), before any permutation runs.
- N > 1000 with the NA outside the variogram subsample: at the permutation variogram (line 114), or as 'incompatible dimensions' from cor(permutations, Y) (line 323).
So a fix must remove non-finite pairs, and their coordinates, before both the correlation and the permutations. Guarding geoR::variog alone is not enough.

spatialSimilarity errors on NA or NaN at quantile() (packageFunction.R:221 and 226), before threshold()'s documented na.omit runs. Supplying the threshold for the dataset that holds the NA avoids the error. Even then, the NA pixel is counted in numPixelOutThresh, as if it were below threshold.

Constant genes give an all-NA row: cor.test warns and returns NA r and p, and the permutation step fails. In the combined output, deltaStarX, deltaStarY and nullCorrelations* hold logical NA for these rows and numeric vectors for normal rows.

NaN can come from the AKI vignette's unguarded CPM step (colSums of 0). That only happens when a pixel has zero total counts, and the bundled rasterized data has no such pixel.

**Evidence (auditor):** WD/03_errors.R, single NA in X: correlationCoef -0.1606388, pValueNaive 0.01473639, pValuePermuteX NA, pValuePermuteY NA. WD/07_similarity_plots.R: 'ERROR: missing values and NaN's not allowed if 'na.rm' is FALSE'. With t1 = t2 = 0 the same data works.

**Verifier reasoning:** The core behaviour is real, and the cited numbers reproduce exactly. One NA, NaN or Inf in X or Y gives NA for both empirical p-values. Only cor.test's numbers survive: complete-case r and p for NA/NaN, but NaN for Inf. The error handler at spatialCorrelation.R:585-619 only print()s the condition, so no warning or status is recorded.

The spatialSimilarity claim is also true. quantile() at packageFunction.R:221 and 226 errors on NA before threshold()'s documented na.omit (line 57) can run.

The constant-gene claim is true. These genes give an all-NA row whose deltaStarX is logical NA, while normal rows hold numeric vectors.

The mechanism detail is wrong: locfit does not fail on NA. Its model.frame uses na.action = na.omit, so it silently drops the NA row and returns N-1 fitted values, shifted relative to the coordinates. The quoted 'NA/NaN/Inf in foreign function call (arg 4)' comes from geoR::variog (lapply(as.data.frame(data), bin.f) -> .C('binit')):
- N <= 1000: the target variogram at line 272 fails before any permutation runs.
- N > 1000, NA pixel outside the 1000-point subsample: the permutation variogram at line 114 fails on the shifted locfit output, or cor(permutations, Y) at line 323 fails with 'incompatible dimensions'.
Every path ends in NA, never in a silently wrong p-value.

Smaller corrections:
- For spatialSimilarity, supplying only the threshold of the dataset that holds the NA is enough (t1 alone works when the NA is in x). When it does run, the NA pixel is silently counted in numPixelOutThresh, as if it were below threshold.
- NaN from CPM is plausible: the AKI vignette's CPM (vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:265-266) has no guard against zero totals. But it needs a pixel with zero total counts, and the rasterized speKidney has none, so 'arises naturally' is overstated.

Severity is medium, not critical. The p-values come back NA rather than wrong. NA input is an edge case, because SEraster output carries no NA unless the user's preprocessing creates it. Still, a gene's whole test is dropped with only a print, spatialSimilarity hard-errors, and threshold()'s docs promise NA removal.

Separate from B12: any all-NA row, from an NA pixel or a constant gene, is carried forward by .get_genes_to_repermute (iterativePermutations.R:65-66, where rownames are indexed by an NA logical) as an NA gene name. spatialCorrelationGeneExpIterPermutations then crashes in round 2 with 'subscript out of bounds' and the whole run is lost. That deserves its own finding, likely high.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE)
set.seed(1); n <- 200; p <- cbind(runif(n), runif(n)); x <- sin(5*p[,1]) + rnorm(n, sd=.2); y <- x + rnorm(n, sd=.2); x[5] <- NA
spatialCorrelation(x, y, p, nPermutations=5, deltaX=c(.2,.5), deltaY=c(.2,.5))[,1:4]  # r reported, empirical p NA
```

**Suggested fix:** Remove non-finite pairs once (and their coordinates) before both the correlation and the permutations, and report how many were removed. Use quantile(na.rm=TRUE) in spatialSimilarity. Return a status column instead of printing.

## B13: spatialCorrelationGeneExpWithinSample documents deltaX/deltaY but its argument is delta; 1-gene input fails

- **Location:** `R/spatialCorrelation.R:977`
- **Category:** docs-code-mismatch
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** The signature has 'delta' (line 977), but the roxygen (lines 868-889) documents deltaX and deltaY. R CMD check raises an Rd usage WARNING, and calling the function with the documented argument fails with 'unused argument'. With a single gene, combn(rownames(input), 2) errors with 'n < m', with no informative message. The @return section lists permutationsX but not permutationsY.

**Evidence (auditor):** R CMD check WARNING: 'Undocumented arguments ... delta; Documented arguments not in \usage ... deltaX deltaY'. WD/12_misc.R: deltaX=list(...) gives 'unused argument (deltaX = list(0.5, 0.5))'. With 1 gene: 'ERROR: n < m'.

**Verifier reasoning:** Every part of the claim holds, and the cited lines are correct.

1. The signature has `delta = NULL` at R/spatialCorrelation.R:977. The roxygen @param blocks at lines 868-889 document `deltaX` and `deltaY` instead. That text is a verbatim copy of spatialCorrelationGeneExp's docs (lines 649-670), and that function really does have separate deltaX/deltaY arguments (line 767). The generated man/spatialCorrelationGeneExpWithinSample.Rd and the pkgdown page (docs/reference/...html: arg-deltax and arg-deltay anchors, usage shows `delta = NULL`) both carry the mismatch.

2. R CMD check reports it as a WARNING under 'checking Rd \usage sections', with the exact text quoted in the finding.

3. The function has no `...`, so the documented arguments fail with 'unused argument'.

4. With 1 gene, `combn(rownames(input), 2)` at line 1007 stops with the uninformative 'n < m'.

5. spatialCorrelation() returns a `permutationsY` column when returnPermutations=TRUE (lines 563-564), and WithinSample passes it through unchanged. Its @return (lines 952-953) lists only permutationsX.

The real semantics of `delta`, verified at runtime: it is one list with an element per gene. Line 998 overwrites its names positionally with rownames(input). The X direction uses delta[[first gene]] and the Y direction uses delta[[second gene]] (lines 1029).

Related problems I found that the finding does not mention. They are not errors in it:
- (a) The function's own @examples (lines 962-973) cannot run. It passes the list `rastKidney`, which fails S4 dispatch in spatialCoords. Each element has only one gene, which gives 'n < m'. There is no assay named 'A' either; the rasterized assay is 'pixelval'.
- (b) A multi-gene SPE with NULL rownames also fails with 'n < m', so the suggested `nrow(input) >= 2` check alone is not enough.
- (c) Because `delta` is undocumented, plausible guesses at its format silently change the statistics. `delta = seq(0.1,0.9,0.1)` turns off the delta grid search: gene 1 gets a fixed 0.1 and gene 2 a fixed 0.2. A named list in a different row order, such as list(B=0.8, A=0.2), is relabelled by position, so the per-gene deltas get swapped. Other guesses (`delta=0.5`, or a length-1 list) fail with a cryptic names<- error.
- (d) The same R CMD check WARNING block also flags spatialSimilarity.Rd for an undocumented `verbose` argument. Fixing B13 alone will not clear the WARNING.

Severity: medium. This is misleading documentation about a key statistical parameter: the documented arguments error out loudly, and the silent cases need the user to guess the format. Nothing in vignettes/ or inst/scripts calls this exported function, so it is outside the main workflow, which argues against high.

```r
# using inp from B01
spatialCorrelationGeneExpWithinSample(inp[[1]][1:2,], nPermutations=3, deltaX=list(.5,.5))  # unused argument
spatialCorrelationGeneExpWithinSample(inp[[1]][1,], nPermutations=3)  # n < m
```

**Suggested fix:** Document 'delta' (a per-gene list) or rename it, check that nrow(input) >= 2, and add permutationsY to @return.

## B16: savePlots non-geometry branch plots sample 1's expression in both panels

- **Location:** `R/visualizationFunctions.R:395`
- **Category:** correctness
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** In the branch used when colData has no 'geometry' column (lines 383-420), dfb takes its coordinates from rastGexp[[2]] but its colour from assay(rastGexp[[1]], ...)[gene, sharedPixels] (line 395). Panel b therefore shows sample 1 again, mislabelled as sample 2. Neither panel has a title in this branch.

**Evidence (auditor):** WD/07_similarity_plots.R, with SEraster objects stripped of geometry: 'panel b colour values identical to sample A: TRUE ; identical to sample B: FALSE', even though cor(A,B) = -0.947 for this gene.

**Verifier reasoning:** I read R/visualizationFunctions.R. savePlots() checks only colData(rastGexp[[1]]) for a 'geometry' column (line 364). When it is absent, the else-branch at lines 383-420 runs. Line 394 builds dfb from spatialCoords(rastGexp[[2]]), but line 395 takes its colour from SummarizedExperiment::assay(rastGexp[[1]], assayName)[gene, sharedPixels]. So panel b is drawn with sample 1's values.

On the SEraster shared grid, shared pixels have identical coordinates in both objects. Panel b is therefore a pixel-for-pixel copy of panel a. Neither ggplot in this branch sets a title, whereas the geometry branch passes plotTitle = names(rastGexp)[[i]].

One wording point: 'mislabelled as sample 2' is loose, because panel b carries no label. It is presented as sample 2 only by its position and by the roxygen docs (lines 326-327: 'Panel 2: Spatial expression plot for the second experiment'). The finding already says the panels have no titles, so this is not a factual error. The mechanism, line number, consequence, evidence numbers and suggested fix all check out.

Positive control: the geometry branch (normal SEraster objects) is correct. The plotRaster layer fill for panel b equals sample B's values, so the defect is confined to the non-geometry branch.

Severity is medium, matching the finding. The output is a wrong figure, not a wrong statistic: the similarity table and panels 3-4 (pixelClass and linearRegression) are correct. That also makes the figure contradict itself, with two identical maps beside a classification showing large A+/B+ regions. The case is an edge case:
- The documented workflow uses SEraster::rasterizeGeneExpression output, which always has geometry, so it takes the correct branch.
- Both packaged datasets take the correct branch too: simRanPatternRasts carries geometry, and speKidney is rasterized with SEraster.
- The silent wrong figure needs inputs with no geometry column whose spatialCoords columns are named exactly x and y, such as SEraster output with geometry dropped or custom-binned objects.
- With other coordinate names (Visium-style pxl_col_in_fullres/pxl_row_in_fullres), the branch fails when the figure is drawn ('object x not found') rather than silently mis-plotting.

savePlots is user-facing: NAMESPACE is exportPattern("^[[:alpha:]]+").

Extra cosmetic issue, not part of the claim: in this branch guides(fill = guide_colorbar(...)) does nothing because the mapped aesthetic is colour. The 3-inch colourbar sizing used in the geometry branch is therefore not applied, which is visible in the rendered PNGs.

```r
devtools::load_all('/Volumes/Crucial SSD/Dropbox (Personal)/work/github.com/slowkow/STcompare', quiet=TRUE); library(SpatialExperiment); library(SummarizedExperiment)
data(speKidney); rk <- SEraster::rasterizeGeneExpression(speKidney, assay_name='counts', resolution=0.2, fun='mean', square=FALSE, BPPARAM=BiocParallel::SerialParam())
ng <- function(s){cd <- colData(s); cd$geometry <- NULL; SpatialExperiment(assays=assays(s), spatialCoords=spatialCoords(s), colData=cd)}
A <- ng(rk$A); B <- ng(rk$B); sh <- intersect(colnames(A), colnames(B)); A <- A[,sh]; B <- B[,sh]
p <- savePlots('Gene', spatialSimilarity(list(A=A,B=B)), list(A=A,B=B))
all.equal(unname(p$Gene[[2]]$data$color), unname(assay(A)['Gene',]))  # TRUE
```

**Suggested fix:** Use rastGexp[[2]] at line 395. Add titles, and factor out a single helper for both panels.

## B17: spatialSimilarity does not store its assay, and savePlots never forwards assayName, so a different assay is plotted than was classified

- **Location:** `R/visualizationFunctions.R:422`
- **Category:** api
- **Severity:** auditor medium, verifier medium
- **Verdict:** confirmed

**Description:** spatialSimilarity's $parameters (packageFunction.R:304-308) does not record assayName, so linearRegression and pixelClass default to assay 1. savePlots accepts assayName but calls pixelClass(spatialSimilarity, gene) at line 362 and linearRegression(input = spatialSimilarity, gene = gene) at line 422 without it, and calls SEraster::plotRaster with assay_name = NULL (lines 365 and 374), which plots the first assay. When similarity was computed on 'lognorm' or 'CPM', as in both real-data vignettes, the scatter shows 'counts'/'pixelval' values coloured by the lognorm classification, so the fold-change guide lines no longer match the points.

**Evidence (auditor):** WD/07_similarity_plots.R, with spatialSimilarity(..., assayName='lognorm'): 'linearRegression(sL, 'Gene') x-values come from assay 'pixelval': TRUE'; plotted x range 4.74-37.40 vs lognorm range 0.76-1.58. 'savePlots(..., assayName='lognorm') panel 4 x-values from 'pixelval': TRUE'. The plotRaster source contains 'if (is.null(assay_name)) mat <- SummarizedExperiment::assay(input)'.

**Verifier reasoning:** Each factual claim matches the source and a fresh run.

(1) R/packageFunction.R:304-308 returns parameters = list(foldChange, minPixels, input). It does not record assayName.

(2) linearRegression (R/visualizationFunctions.R:46-48) and pixelClass (lines 160-162) both default assayName to 1.

(3) savePlots calls pixelClass(spatialSimilarity, gene) at line 362 and linearRegression(input = spatialSimilarity, gene = gene) at line 422 without assayName. Its geometry branch calls SEraster::plotRaster at lines 365 and 374 without assay_name. SEraster 1.2.0's plotRaster does `if (is.null(assay_name)) mat <- SummarizedExperiment::assay(input)`, so it plots the first assay. rasterizeGeneExpression always names its output assay 'pixelval'. So whenever a normalized assay is added after rasterization, assay 1 is the raw 'pixelval'. Both real-data vignettes do this: AKI adds 'CPM' and MERFISH adds 'lognorm'.

The consequence reproduces. The scatter's x/y values come from pixelval while the colours come from the lognorm classification, so many 'similar' (blue) points fall outside the dashed fold-change lines. The worst case is A vs C: the title says S = 1 and every point is blue, yet every point lies outside the ±1 log2FC band. On pixelval itself S would be 0. Passing assayName='lognorm' to linearRegression gives zero mismatches. Inside savePlots this cannot be fixed by the user, because passing assayName does nothing for panels 1, 2 and 4 in the geometry branch. The only workaround is to reorder assays so the classified one is first; I checked that this makes the default plots correct.

Scope clarifications, which do not change the verdict:
(a) pixelClass's plot does not depend on assay values. Its layer data holds only 'fill' (from the stored classification) and 'geometry', and these are identical for the default call and assayName='lognorm'. So the missing assayName at line 362 is harmless and panel 3 is correct. The wrong panels are 1, 2 and 4.
(b) The vignettes' own figures are not affected. Both pass assayName explicitly to linearRegression (AKI Rmd lines 628-629, MERFISH Rmd lines 421 and 427). Neither calls savePlots, which appears only in its own roxygen example on single-assay speKidney data.
(c) In savePlots's non-geometry branch (lines 383-420), assayName is used for the two spatial panels. However, line 395 builds panel 2 from rastGexp[[1]], so that panel shows dataset 1's values. This is a separate bug; I reproduced it.
(d) In the vignette data, assay 1 is named 'pixelval', not 'counts'.

Severity is medium. The similarity scores and classifications are correct. The diagnostic scatter is silently inconsistent with its own title and guide lines whenever the classified assay is not assay 1 and assayName is omitted. With savePlots it is always wrong in that setting. The savePlots @param assayName documentation is misleading.

```r
# rk from B16
A <- rk$A; B <- rk$B; SummarizedExperiment::assay(A,'lognorm') <- log10(SummarizedExperiment::assay(A)+1); SummarizedExperiment::assay(B,'lognorm') <- log10(SummarizedExperiment::assay(B)+1)
s <- spatialSimilarity(list(A=A,B=B), assayName='lognorm'); range(linearRegression(s,'Gene')$data$x)  # pixelval scale, not lognorm
```

**Suggested fix:** Store assayName in spatialSimilarity()$parameters and use it as the default in linearRegression, pixelClass and savePlots. Forward assayName (as assay_name) to plotRaster.

## B14: R CMD check ERROR: the simRanPatternRasts example calls assays() without a namespace

- **Location:** `R/data.R:97`
- **Category:** packaging
- **Severity:** auditor medium, verifier low
- **Verdict:** partially-confirmed

**Verified description:** The roxygen example at R/data.R:97, shipped as man/simRanPatternRasts.Rd:75, calls the bare SummarizedExperiment function assays(). STcompare's NAMESPACE has no import directives (only exportPattern), and SummarizedExperiment is listed only in DESCRIPTION Imports, not Depends. Example code runs in the global environment with only STcompare attached, so in any clean session (R CMD check, example(simRanPatternRasts), or a script with only library(STcompare)) the line fails with 'could not find function "assays"'. The repo's built pkgdown page (docs/reference/simRanPatternRasts.html) already shows this error. The line works only when SummarizedExperiment is attached, for example via library(SpatialExperiment), as the vignettes do.

Under R CMD check this gives 'checking examples ... ERROR' (with --no-manual --ignore-vignettes: Status 1 ERROR, 3 WARNINGs, 2 NOTEs). The examples script halts at this example, so 6 of the 12 examples never run: spatialCorrelation, spatialCorrelationGeneExp, spatialCorrelationGeneExpIterPermutations, spatialCorrelationGeneExpWithinSample, spatialSimilarity and viladomatCorrelation. Among them, the spatialCorrelationGeneExpWithinSample example is itself broken (spatialCoords() is called on a list). R CMD check as a whole does not stop: with the default _R_CHECK_EXIT_ON_FIRST_ERROR_=false it continues to tests, vignettes and the manual.

The fix SummarizedExperiment::assay(simRanPatternRasts[[1]], 'pixelval')[1, 1:5] gives identical output and passes check. It must be made in R/data.R and the Rd regenerated with roxygen. On its own it leaves the check at 1 ERROR, because the error moves to the WithinSample example.

This is a packaging and documentation defect with no effect on computed results. Severity: low.

**Evidence (auditor):** check_full.log: 'Error in assays(simRanPatternRasts[[1]]) : could not find function "assays" ... Status: 1 ERROR, 3 WARNINGs, 2 NOTEs'.

**Verifier reasoning:** The core defect is real and reproduces four independent ways: devtools::load_all (with and without --vanilla), a clean session with only library(STcompare) loaded from an installed copy, example("simRanPatternRasts"), and my own R CMD build plus check. The cited line, error text, status line, masking effect and suggested fix all check out.

Three details are imprecise or overstated:
(1) Mechanism wording. "It is only imported" is loose. SummarizedExperiment is listed only under DESCRIPTION Imports. The NAMESPACE contains just exportPattern("^[[:alpha:]]+") with no import directives, so getNamespaceImports('STcompare') returns only base. This does not change the outcome: examples run in the global environment, so even a NAMESPACE import would not make assays() visible. The bare call only works when SummarizedExperiment is attached, for example via library(SpatialExperiment), which the vignettes do. That is presumably why the authors never noticed.
(2) Consequence wording. "R CMD check stops there" is overstated. Only the examples script STcompare-Ex.R halts. In R 4.5.2, run_examples() calls maybe_exit(1L), which quits only when _R_CHECK_EXIT_ON_FIRST_ERROR_ is true (default FALSE), so R CMD check goes on to demos, tests and vignettes. The cited run appeared to end because the package has no tests and the run used --ignore-vignettes --no-manual. The claim that all later examples go unchecked is correct: 6 of the 12 examples never run. One of them, spatialCorrelationGeneExpWithinSample, is itself broken, which confirms the masking point.
(3) Severity. Under the given scale this is low (packaging/minor), not medium. It is a documentation-example and packaging defect with no effect on any computed result, and users who attach SpatialExperiment never see it. It does block a clean R CMD check, BiocCheck or CI run, and the published pkgdown reference page already shows the error to readers.

The fix is easy to misread as complete. The suggested fix works (identical output; 'unstated dependencies in examples ... OK'). It must go into R/data.R and the Rd must be regenerated, because R CMD check runs man/simRanPatternRasts.Rd:75. Fixing B14 alone does not clear the check: the status stays at 1 ERROR, 3 WARNINGs, 2 NOTEs, with the error moving to the spatialCorrelationGeneExpWithinSample example.

```r
# in a session with only library(STcompare)
data(simRanPatternRasts); assays(simRanPatternRasts[[1]])  # could not find function 'assays'
```

**Suggested fix:** Use SummarizedExperiment::assay(simRanPatternRasts[[1]], 'pixelval')[1, 1:5].

## B15: spatialCorrelationGeneExpWithinSample example passes a list of 3 objects and assay 'A', so it errors

- **Location:** `R/spatialCorrelation.R:969`
- **Category:** packaging
- **Severity:** auditor medium, verifier low
- **Verdict:** confirmed

**Description:** The example calls spatialCorrelationGeneExpWithinSample(input = rastKidney, assayName = 'A'). But rastKidney is a list of 3 SpatialExperiments, and the function expects one; 'A' is a list element, not an assay name (the assay is 'pixelval'). Even with valid input, the kidney objects have only 1 gene, so combn() would fail. This failure is hidden in R CMD check behind B14.

**Evidence (auditor):** WD/examples/ex_spatialCorrelationGeneExpWithinSample.log: 'ERROR: unable to find an inherited method for function 'spatialCoords' for signature 'x = "list"''.

**Verifier reasoning:** I tried to refute each part of the finding. Every factual claim held up.

(1) The example at R/spatialCorrelation.R:960-973 (call at lines 969-972) is identical to the \examples section of man/spatialCorrelationGeneExpWithinSample.Rd. It is not wrapped in \dontrun, so both example() and R CMD check run it. It passes `input = rastKidney, assayName = "A"`.

(2) `rastKidney` is the output of SEraster::rasterizeGeneExpression() on the speKidney list. It is a plain list of 3 SpatialExperiments (A, C, B). The function documents `input` as a single SpatialExperiment (line 859), and its first data access is SpatialExperiment::spatialCoords(input) at line 992. On a list this raises exactly the quoted S4 dispatch error.

(3) 'A' is a list element name, not an assay. Each element's only assay is 'pixelval'. On a valid 2-gene SPE, assayName = 'A' fails at assay(input, assayName) (lines 1023-1024) with "'A' not in names(assays(<SpatialExperiment>))".

(4) Each kidney object has 1 gene ('Gene'), so even with valid input the call fails at combn(rownames(input), 2) (line 1007) with 'n < m'. This happens for assayName NULL, 'A' and 'pixelval' alike.

(5) R CMD check does hide this failure. In both --as-cran and default mode, the example run halts first in the simRanPatternRasts example: `assays()` is called without a namespace prefix, NAMESPACE has only exportPattern() and no imports, and cleanEx() detaches any packages attached by earlier examples. I cannot see the text of B14, but this earlier example error is a plausible match for it.

One precision the finding lacks: under --as-cran, which devtools::check() uses by default, there is a second blocker. Even after the assays() fix, the spatialCorrelationGeneExp example's nThreads = 5 makes BiocParallel stop with 'workers must be <= 2'. In default mode, fixing only assays() makes this example the next failure.

(6) The quoted error message and the line number are exact. The suggested fix works.

The broken example is also visible to users. The built pkgdown page (docs/reference/spatialCorrelationGeneExpWithinSample.html) renders the error, followed by "object 'sc_within_sample' not found". No vignette, script or README shows correct usage of this function, and no shipped dataset has more than 1 gene per object, so it cannot be run as-is.

Severity: I rate this low, one notch below the claimed medium; it sits on the low/medium boundary. It is a packaging and documentation defect (the finding's own category). It fails loudly and immediately, so it never produces wrong results, and the @param documentation is correct. It is still worth fixing alongside the other example failures: once those are fixed, it becomes a check-blocking ERROR.

A related but separate defect: the same Rd fails the R CMD check '\usage' test. It documents 'deltaX' and 'deltaY', but the function's argument is 'delta'.

```r
example('spatialCorrelationGeneExpWithinSample', package='STcompare')  # after installing
```

**Suggested fix:** Build a multi-gene SpatialExperiment for the example (for example by combining several simRanPatternRasts fields) and use assayName = 'pixelval'.

## B18: savePlots calls library() on undeclared gridExtra/patchwork, attaching them to the user's session and masking BiocGenerics::combine

- **Location:** `R/visualizationFunctions.R:351`
- **Category:** packaging
- **Severity:** auditor medium, verifier low
- **Verdict:** partially-confirmed

**Verified description:** savePlots() (R/visualizationFunctions.R:351-353) calls library(ggplot2), library(gridExtra) and library(patchwork). Every call attaches all three to the caller's search path, and they remain attached afterwards, even when the call errors.

The function also relies on that side effect:
- unit() (lines 370-416) and plot_layout() (line 424) are unqualified.
- The ggplot composition `a + b + pc + c` (line 424) needs patchwork's namespace loaded.
- Removing only the library() calls gives 'could not find function "unit"', and then "Can't add `b` to a <ggplot> object".

patchwork is undeclared and is not in STcompare's recursive Depends/Imports tree. A fresh remotes::install_github() therefore leaves it out, and savePlots then fails with "there is no package called 'patchwork'" after already attaching ggplot2 and gridExtra. gridExtra is also undeclared, but it is installed transitively because viridis (Imports) imports it.

The masking is conditional. gridExtra::combine masks BiocGenerics::combine and Biobase::combine only when the user has already attached those packages, e.g. via library(SpatialExperiment) as the vignettes do. In that case later combine() calls on data.frames or Bioc objects error. library(STcompare) alone attaches neither package.

R CMD check reports:
- a WARNING: 'library' or 'require' calls not declared from / in package code: ggplot2, gridExtra, patchwork;
- NOTEs for unit/plot_layout as undefined globals;
- the example harness detaching patchwork, gridExtra and ggplot2 after the savePlots example.

Fix:
- Drop the library() calls.
- Use ggplot2::unit (or grid::unit) and patchwork::wrap_plots(a, b, pc, c, ncol = 4).
- Declare patchwork, either in Imports with importFrom() in NAMESPACE, or in Suggests with a requireNamespace() check before composing.
- Remove the unused gridExtra.

In an installed package, keeping the `+` chain with only a trailing patchwork::plot_layout() and patchwork listed only in DESCRIPTION Imports would still fail: `a + b` is evaluated before patchwork's namespace is loaded.

**Evidence (auditor):** WD/07_similarity_plots.R: 'Attaching package: 'gridExtra' ... masked from 'package:BiocGenerics': combine'; 'newly attached: package:patchwork package:gridExtra'. The R CMD check WARNING says ''library' or 'require' calls in package code', and the check had to detach patchwork, gridExtra and ggplot2 after the savePlots example.

**Verifier reasoning:** The core claim holds. savePlots runs library(ggplot2), library(gridExtra) and library(patchwork) at R/visualizationFunctions.R:351-353. Those packages stay attached in the caller's session after the call. The function also relies on them being attached: unit() at lines 370-416 and plot_layout() at line 424 have no namespace prefix. patchwork is not declared, so savePlots crashes when patchwork is not installed. The R CMD check WARNING and the example-harness detach message both reproduce.

Three details are wrong or overstated.
(1) gridExtra IS installed by a fresh install_github(). viridis is in Imports and imports gridExtra, and remotes' default dependencies=NA covers Depends/Imports/LinkingTo. Only patchwork is missing from the dependency tree.
(2) The combine masking only happens when BiocGenerics/Biobase are already attached, e.g. after library(SpatialExperiment) as in the vignettes and inst/scripts. library(STcompare) alone attaches neither, so nothing is masked then.
(3) The mechanism is incomplete. The ggplot `a + b + pc + c` chain at line 424 also needs patchwork's namespace to be loaded, not just plot_layout. So the suggested fix `... + patchwork::plot_layout()` still fails if nothing has loaded patchwork before `a + b` is evaluated. In an installed package, listing patchwork in DESCRIPTION Imports does not load it while NAMESPACE has no import directives.

Severity: no statistical result is affected. The crash only happens when patchwork is absent, and the error names the missing package. The masking prints a message, and the downstream combine() calls error loudly rather than return wrong results. It is a packaging/hygiene problem that also blocks a clean R CMD check, so low. A case for medium rests only on fresh installs that lack patchwork.

```r
before <- search(); invisible(savePlots('Gene', spatialSimilarity(list(rk$A, rk$B)), list(rk$A, rk$B))); setdiff(search(), before)
```

**Suggested fix:** Drop the library() calls. Use patchwork::plot_layout / patchwork::wrap_plots and grid::unit (or ggplot2::unit). Put patchwork in Imports (or in Suggests guarded by requireNamespace). Remove gridExtra, which is unused.

## B19: numPixelInThresh is dim(thresh)[1], which is always 1, for genes skipped by minPixels

- **Location:** `R/packageFunction.R:249`
- **Category:** correctness
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** In the skip branch, thresh is the 1-row data.frame returned by threshold(), so dim(thresh)[1] is always 1. The intended value is nrow(threshDF). The non-skip branch (line 286) is correct.

**Evidence (auditor):** WD/07_similarity_plots.R, a gene with 5% non-zero pixels: numPixelInThresh 1, numPixelOutThresh 247. The true number passing the threshold is 26, and 1 + 247 does not equal the 273 shared pixels.

**Verifier reasoning:** threshold() (R/packageFunction.R:73-79) always returns a 1-row data.frame with list-columns, so dim(thresh)[1] at line 249 is always 1. The intended count is nrow(threshDF), the data frame built at lines 232-236. The skip branch therefore reports numPixelInThresh = 1 whatever the real count is, and that holds even when 0 pixels pass. As a result, numPixelInThresh + numPixelOutThresh does not add up to the number of shared pixels. The non-skip branch (line 286, dim(logTrans)[1]) is correct. The skip branch is reached in normal use for sparse genes, because with the default minQuantile = 0.05, t1 and t2 become 0. Nothing in R/ reads numPixelInThresh; it is only an output column, and the similarity score is NA for skipped genes anyway. So the harm is a misleading diagnostic field, and low severity fits.

```r
# add a gene with 95% zeros to both rk$A and rk$B, then spatialSimilarity(...)$similarityTable$numPixelInThresh  # 1 for that gene
```

**Suggested fix:** Use numPixelInThresh = nrow(threshDF) (and pixelIDInThresh = thresh$pixel) in the skip branch.

## B20: getGenePixelDF default assayName = assayName is self-referential

- **Location:** `R/packageFunction.R:17`
- **Category:** api
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** The default argument refers to itself, so calling the exported getGenePixelDF without assayName fails with a recursive-default-argument error.

**Evidence (auditor):** WD/07_similarity_plots.R: 'promise already under evaluation: recursive default argument reference or earlier problems?'.

**Verifier reasoning:** The signature at R/packageFunction.R:17 is getGenePixelDF(x, y, gene, assayName = assayName). A default is evaluated in the function's own frame, where 'assayName' is bound to that same promise, so leaving the argument out always gives a recursive-promise error. A global variable named assayName does not help. The helper is exported through exportPattern, and man/getGenePixelDF.Rd shows the usage line with this default. Every internal caller passes assayName explicitly (packageFunction.R:215, visualizationFunctions.R:50 and :164), so only users who call the helper directly are affected.

```r
getGenePixelDF(rk$A, rk$B, 'Gene')  # recursive default argument reference
```

**Suggested fix:** Use assayName = 1, or NULL mapped to 1.

## B21: spatialCorrelation documents '1 x N matrix' input for X and Y but errors on it

- **Location:** `R/spatialCorrelation.R:512`
- **Category:** api
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** The @param docs (R/spatialCorrelation.R:346-350) say X and Y may be 'a 1 x N numeric vector or matrix', but a 1 x N matrix errors. The failure is not at lines 512-515. data.frame(X = X, Y = Y, x = pos[,1], y = pos[,2]) succeeds: it recycles the 1-row matrices and silently builds an N x (2N+2) data frame with columns X.1..X.N, Y.1..Y.N, x, y. dataForward$X is then NULL, because no column is named exactly 'X' and the partial match is ambiguous. The dataReverse data.frame() at lines 517-520 then fails with 'arguments imply differing number of rows: 0, N'. Both statements run before the tryCatch at line 531, so the user gets the raw error. Numeric vectors and N x 1 matrices work. Internal callers (spatialCorrelation.R:822-823, assay[g, shared_pixels]) pass dropped vectors. Fix: coerce with as.numeric() and check lengths against nrow(pos), or correct the docs.

**Evidence (auditor):** WD/03_errors.R: 'ERROR: arguments imply differing number of rows: 0, 231'.

**Verifier reasoning:** The main claim is real: documented 1 x N matrix input errors before the tryCatch. The mechanism and line are wrong, though. The dataForward call at 512-515 does not fail. The error comes from dataReverse at 517-520, because dataForward$X is NULL after the matrix is expanded into many columns. That explains the '0' in the reported message ('0, 231'). The function errors loudly rather than returning wrong numbers, and only direct calls are affected, so severity is low.

```r
spatialCorrelation(matrix(rnorm(50),1), matrix(rnorm(50),1), cbind(runif(50),runif(50)), nPermutations=2, deltaX=.5, deltaY=.5)
```

**Suggested fix:** Coerce with as.numeric(X) and as.numeric(Y), and validate the lengths against nrow(pos).

## B22: Undeclared dependencies and imports (sf, stats/utils functions, vignette packages); dplyr version floor missing

- **Location:** `DESCRIPTION:26`
- **Category:** packaging
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** sf:: is used in pixelClass (visualizationFunctions.R:182) but sf is not declared. stats functions (cor, cor.test, dist, fitted, lm, median, na.omit, quantile, rnorm, p.adjust.methods) and utils::combn are used unqualified with no importFrom. dplyr::case_when(.default=) requires dplyr >= 1.1.0, but no version is given. The vignettes use MERINGUE, scatterbar, rhdf5, rjson, Matrix, BiocGenerics, patchwork and gridExtra, none of which are in Suggests; MERINGUE (GitHub-only) and scatterbar are not installed here, so those vignettes cannot be built. 'class' is in Suggests but is never used.

**Evidence (auditor):** R CMD check WARNING: ''::' or ':::' import not declared from: 'sf''. NOTE: 'Undefined global functions or variables: ... combn cor cor.test dist ... Consider adding importFrom("stats", ...)'. dplyr NEWS.md line 656 (section 'dplyr 1.1.0') introduces the case_when .default argument. requireNamespace('MERINGUE') and requireNamespace('scatterbar') both return FALSE.

**Verifier reasoning:** I verified every listed item. (1) sf:: is called at visualizationFunctions.R:182, and sf is in neither Imports nor Suggests. R CMD check gives a WARNING. The practical impact is small because SEraster imports sf, so sf is always installed. (2) The stats and utils functions are used without qualification and NAMESPACE has no imports. R CMD check suggests exactly the importFrom list in the finding. (3) The .default argument of case_when is used at spatialCorrelation.R:1131 (plotCorrelationGeneExp) and only exists from dplyr 1.1.0. In dplyr 1.0.10, case_when is function(...) and validate_formula aborts on any argument that is not a formula, while DESCRIPTION lists 'dplyr' with no version. (4) All 8 vignette packages are undeclared. MERINGUE and scatterbar are not installed, and their library() calls sit in evaluated setup chunks (no eval=FALSE) of the AKI and MERFISH vignettes, so those vignettes cannot be knitted here. MERINGUE is on neither CRAN nor Bioc 3.22, which supports 'GitHub-only'. scatterbar is on CRAN. (5) The 'class' package is never used anywhere in the repo, docs/ excluded. The finding leaves some things out but nothing it says is wrong. savePlots (visualizationFunctions.R:351-353) calls library(ggplot2), library(gridExtra) and library(patchwork) inside package code, so gridExtra and patchwork are runtime dependencies of the package, not just of the vignettes, and R CMD check flags this in the same WARNING. The vignettes also call devtools::load_all(), and devtools is undeclared. Also, 'checking for unstated dependencies in vignettes ... OK' in an R CMD check without built vignettes means nothing, because the checker only reads the built inst/doc.

**Suggested fix:** Add sf (and patchwork) to Imports; add @importFrom stats ... and utils combn; require dplyr (>= 1.1.0); list the vignette packages in Suggests (with Remotes: for MERINGUE) or make those chunks conditional; drop class.

## B23: exportPattern('^[[:alpha:]]+') exports internal helpers; NAMESPACE is not managed by roxygen

- **Location:** `NAMESPACE:1`
- **Category:** packaging
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** NAMESPACE contains only exportPattern("^[[:alpha:]]+") and has no roxygen header, so roxygen2 never regenerates it and the @export tags have no effect. Helpers without an @export tag (getGenePixelDF, threshold, assignFill) are exported, including the very generic name 'threshold', which can mask other packages. No imports are declared at all.

**Evidence (auditor):** getNamespaceExports('STcompare') lists 14 functions, including getGenePixelDF, threshold and assignFill. roxygen2::roxygenise() on a copy left NAMESPACE unchanged (WD/pkgcopy_roxy).

**Verifier reasoning:** NAMESPACE is the single line exportPattern("^[[:alpha:]]+") with no roxygen header. roxygen2 7.3.3 refuses to overwrite it, so the @export tags never reach NAMESPACE. Every name that starts with a letter is exported. That includes 3 helpers with no @export tag (assignFill, getGenePixelDF, threshold), alongside the 11 functions that do have @export. No imports are declared: the only import is base. The masking concern is plausible because 'threshold' is a generic name, but none of the packages installed here export 'threshold', so the conflict is possible rather than observed. Side note: man/spatialCorrelationGeneExpIterPermutations.Rd is also not roxygen-generated and was skipped as well.

```r
sort(getNamespaceExports('STcompare'))
```

**Suggested fix:** Delete NAMESPACE and let roxygen2 generate it from the @export and @importFrom tags. Mark helpers @noRd or @keywords internal.

## B24: Non-ASCII em dash in R code causes an R CMD check WARNING

- **Location:** `R/iterativePermutations.R:302`
- **Category:** packaging
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** The verbose message 'no genes left to test — skipping' contains a UTF-8 em dash (e2 80 94) in R code. Portable packages must use ASCII in code.

**Evidence (auditor):** check_full.log: 'checking code files for non-ASCII characters ... WARNING Found the following file with non-ASCII characters: R/iterativePermutations.R'. tools::showNonASCIIfile reports line 302.

**Verifier reasoning:** Line 302 of R/iterativePermutations.R has a UTF-8 em dash (bytes e2 80 94) inside a string literal passed to sprintf/message. R CMD check flags non-ASCII characters in R code even though DESCRIPTION declares Encoding: UTF-8. The other non-ASCII characters, in R/data.R and R/spatialCorrelation.R, sit in roxygen comments, which the check allows, so this is the only file flagged.

```r
tools::showNonASCIIfile('R/iterativePermutations.R')
```

**Suggested fix:** Replace it with '-' or the — escape.

## B25: Rd problems: undocumented 'verbose' argument (WARNING) and 89 'Lost braces' NOTEs

- **Location:** `R/packageFunction.R:162`
- **Category:** packaging
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** spatialSimilarity has a 'verbose' argument with no @param, which triggers an Rd usage WARNING. The @return sections use \itemize{\item{a}{b}}, which checkRd reports as 'Lost braces in \itemize; \value handles \item{}{} directly': 89 instances across 11 Rd files. savePlots.Rd:47 has unescaped braces in '{gene_name}.pdf'.

**Evidence (auditor):** check_full.log: 'Undocumented arguments in Rd file 'spatialSimilarity.Rd' 'verbose''; 'checking Rd files ... NOTE' with 89 'Lost braces' lines.

**Verifier reasoning:** Reproduced with the same functions R CMD check uses, and with a real `R CMD check --as-cran`. spatialSimilarity has `verbose = FALSE` (R/packageFunction.R:162) but no @param verbose. That gives a WARNING under 'checking Rd \usage sections': tools:::.check_packages calls warningLog for checkDocFiles output on non-internal Rd. The Rd-files check runs checkRd at minlevel = -1 (the default for non-base packages). That yields exactly 89 'Lost braces' lines in 11 Rd files, and since they are all negative-level checkRd lines R classes them as a NOTE. The 89 break down as 88 'Lost braces in \itemize; \value handles \item{}{} directly' plus 1 'savePlots.Rd:47: Lost braces; missing escapes or markup?' for "{gene_name}.pdf". The description slightly conflates these by attributing all 89 to the itemize message, but the evidence line ('89 Lost braces lines') is exact. Called directly at all levels, checkRd reports 177, because each \item also emits a level -3 message. man/ is in sync with roxygen: re-running roxygenise on a copy produced no diff, so the fix has to go in the roxygen source. The finding misses one thing: the same WARNING also lists spatialCorrelationGeneExpWithinSample.Rd, which has an undocumented 'delta' and documents 'deltaX'/'deltaY', which are not in \usage.

**Suggested fix:** Add @param verbose. Use \describe{\item{...}{...}} (or a plain \item list in \value) in @return. Escape the braces.

## B26: Package is 48.6 MB as a tarball (33.8 MB installed) and contains a byte-identical duplicate RData

- **Location:** `inst/extdata/brain-MERFISH-10x-visium/brainCorrelation_1.RData:1`
- **Category:** packaging
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** inst/extdata holds 32.3 MB of precomputed .RData results. inst/extdata/brain-MERFISH-10x-visium/brainCorrelation_1.RData is a byte-identical copy of brainCorrelation.RData (MD5 5661fd72aa53da9ebf788764326808f0, 2,759,214 bytes each), and nothing in R/, vignettes/ or inst/scripts references it. .Rbuildignore excludes only Rproj/doc/Meta/_pkgdown/docs/pkgdown, so images/ and the knitr fig.path output folders under vignettes/ (about 13 MB of PNGs) ship in the tarball. R CMD check also NOTEs "Non-standard file/directory found at top level: 'images'". `R CMD build --no-build-vignettes` gives a 48.6 MB tarball, and the installed size is 33.8 MB (extdata 32.3 MB, reported as INFO by R 4.5.2). The applicable Bioconductor rules are a source tarball under 10 MB (BiocCheck::checkPackageSize errors above 10 MB, hard limit 100 MB) and individual files of at most 5 MB (BiocCheck warning). The tarball is about 4.9 times the 10 MB limit, and three files each exceed the 5 MB per-file limit: kidneyCorrelation.RData (13.0 MB), merfishCorrelation.RData (6.7 MB) and merfishCorrelation_affine.RData (6.5 MB). CRAN's rule of thumb (data/docs at most 5 MB) is also exceeded.

**Evidence (auditor):** STcompare_0.1.0.tar.gz is 48,623,451 bytes. R CMD check INFO: 'installed size is 33.8Mb ... extdata 32.3Mb'. Both RData files have MD5 5661fd72aa53da9ebf788764326808f0 (2,759,214 bytes each).

**Verifier reasoning:** Every measurable claim reproduced: the duplicate file, the tarball size (within 3 bytes, a gzip-header difference), the installed size, the extdata size, and .Rbuildignore not excluding images/ or the figure folders. Those folders are knitr fig.path outputs, and README links images/ by absolute GitHub URL, so excluding them is safe. The one wrong detail is the limit: the current Bioconductor guidelines (contributions.bioconductor.org/general.html) say the R CMD build tarball 'should occupy less than 10 MB', with 'individual files must be <= 5MB', and BiocCheck's checkPackageSize defaults to size = 10L. The conclusion ('far above' the limit) still holds.

**Suggested fix:** Remove the duplicate. Host the large precomputed results on Zenodo or ExperimentHub and download them in the vignettes. Add images/ and the figure folders to .Rbuildignore (or reference them by URL).

## B27: Other example defects: nThreads=5 fails under the CRAN/Bioc core limit, a 3.6-minute example, undefined negCorrelation, and the wrong sample in the savePlots example

- **Location:** `R/spatialCorrelation.R:760`
- **Category:** packaging
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** All four example defects are real. (1) The spatialCorrelationGeneExp example (R/spatialCorrelation.R:760-761) uses nThreads = 5. Wherever _R_CHECK_LIMIT_CORES_ is set (R CMD check --as-cran sets it to TRUE if unset), BiocParallel's .enforceWorkers stops with "_R_CHECK_LIMIT_CORES_' environment variable detected, BiocParallel workers must be <= 2 was (5)". On Bioconductor builders the documented mechanism is IS_BIOC_BUILD_MACHINE=true, which caps workers at 4 with a warning and does not fail. I found no evidence that Bioconductor's builders set _R_CHECK_LIMIT_CORES_, so the 'Bioc' part of the title is unsupported. Today the failure is also masked in a real --as-cran check, because the example run halts earlier at the simRanPatternRasts example ("could not find function 'assays'"). (2) The spatialCorrelation example has no \dontrun or \donttest. It runs 100 permutations with the default 9-delta grids in both directions, then 100 permutations with 18 (deltaX) and 25 (deltaY) deltas on 998 quakes points. It took 4.13 min in my run (3.64 min claimed; shared machine). (3) The \dontrun IterPermutations example ends by printing negCorrelation (iterativePermutations.R:224), which it never assigns, so it fails with "object 'negCorrelation' not found". (4) The savePlots example (visualizationFunctions.R:344-345) computes the similarity for A vs B but passes the 3-element rastKidney, whose order is A, C, B. Panel 2 is therefore kidney C (titled 'C'), while panels 3-4 and S = 0.535 come from A vs B. Line 340 has a stray "#' #'", which renders as a harmless comment.

**Evidence (auditor):** WD/12_misc.R: '_R_CHECK_LIMIT_CORES_' environment variable detected, BiocParallel workers must be <= 2 was (5)'. The ex_spatialCorrelation.log result line reports 'elapsed: 3.642768 mins'. WD/06_iter.R: 'object 'negCorrelation' not found'. WD/07_similarity_plots.R: 'panel 2 title: C'.

**Verifier reasoning:** I ran each example verbatim (extracted with tools::Rd2ex) under devtools::load_all on a copy of the repo. The only inaccuracy is the consequence on Bioconductor. BiocParallel 1.44.0 errors only under _R_CHECK_LIMIT_CORES_. Under IS_BIOC_BUILD_MACHINE it warns and sets 4 workers: I verified this directly, and BiocParallel's own docs and vignette show the 4-worker warning. The description body itself attributes the error correctly to _R_CHECK_LIMIT_CORES_ set by --as-cran. The suggested fixes remain appropriate.

```r
Sys.setenv('_R_CHECK_LIMIT_CORES_'='TRUE'); BiocParallel::MulticoreParam(workers=5)  # error
```

**Suggested fix:** Use nThreads = 2 (or the default); cut nPermutations in the examples or wrap them in \donttest; print 'corr'; pass list(A = rastKidney$A, B = rastKidney$B) to savePlots.

## B28: Getting-started vignette: precomputed simRanPattern empirical p-values are called 'BH-corrected' and 'max of pX/pY', but are neither

- **Location:** `vignettes/getting-started-with-STcompare.Rmd:332`
- **Category:** docs-code-mismatch
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** The vignette says corspv_corrected is 'chosen to be the higher of pValuePermuteY and pValuePermuteX' (line 332) and titles the plot 'BH-corrected empirical p-value (pE)' (line 366). The generating script returns results$pValuePermuteX only (inst/scripts/simRanPatternSpatialCorrelation.R:29), and each pair is a separate 1-gene call, so no BH is applied across the 9,900 pairs. The script also silently replaces zeros with 0.01 (line 68). The vignette therefore compares BH-corrected naive p-values (about 43%) with uncorrected empirical ones (about 4%). The text also names the function spatialCorrelationGeneExpIterPermutation (missing 's', line 317).

**Evidence (auditor):** WD/13_vignette_simran.R: naive p<0.05 uncorrected 0.498, after BH 0.428; stored empirical p<0.05 0.039 (uncorrected); after actually applying BH across the 9,900 pairs: 0.

**Verifier reasoning:** Every element checks out, and recomputation goes beyond reading the script. inst/scripts/simRanPatternSpatialCorrelation.R:29 returns results$pValuePermuteX. Each ordered pair is a separate 1-gene call, so the p.adjust at iterativePermutations.R:343-344 acts on a single value and changes nothing, and no BH is ever applied across pairs. The null is deterministic (set.seed(seed), plus seed + i per permutation), so I re-ran two pairs with the repo code. Both reproduced the stored values exactly, and the stored value is pX, not max(pX, pY). The vignette's numbers reproduce exactly, so its plot title ('BH-corrected empirical p-value') and its line-332 comment ('higher of pValuePermuteY and pValuePermuteX') are both false. Line 68 does replace zeros with 0.01; zeros are possible because p = extreme/B (spatialCorrelation.R:326-327). Line 317 does say spatialCorrelationGeneExpIterPermutation, without the 's'. Two details the finding doesn't mention: the 9,900 rows are 4,950 unordered pairs times 2 orientations, with naive p-values duplicated across orientations; and under the max rule the vignette describes, the empirical rate would be 2.1%, not 3.9%. I rate it low because the qualitative conclusion survives any consistent comparison: 49.8% vs 3.9% uncorrected, or 42.8% vs 0% with BH on both.

```r
load(system.file('extdata','simRanPatternResults.RData', package='STcompare')); mean(p.adjust(na.omit(cors_df)$corspv_corrected,'BH') < 0.05)  # 0
```

**Suggested fix:** Regenerate the results with max(pX, pY), apply BH across pairs (or present both uncorrected rates as a calibration check), and fix the labels and the function name.

## B29: Kidney vignette: nThreads described as 'parallelize genes across threads', and spatialCorrelationGeneExp described as returning adjusted p-values

- **Location:** `vignettes/acute-kidney-injury-10x-visium-rasterized.Rmd:419`
- **Category:** docs-code-mismatch
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** Genes are processed sequentially (lapply at iterativePermutations.R:25 and spatialCorrelation.R:808); nThreads parallelises the B permutations within a single gene (bplapply(1:B, ...) at spatialCorrelation.R:290). Line 430 says spatialCorrelationGeneExp returns two adjusted empirical p-values, which B01 shows is false.

**Evidence (auditor):** Code inspection at the cited lines; B01 reproduction.

**Verifier reasoning:** Both stated mismatches are real. (1) Genes run sequentially: lapply over genes at iterativePermutations.R:25 and spatialCorrelation.R:808-809. Parallelism happens only inside viladomatCorrelation, via bplapply(1:B) at spatialCorrelation.R:290 (also 298 and 306), so line 419's 'parallelize genes across threads' is wrong. There is an extra nuance: in that very call the vignette also passes BPPARAM = MulticoreParam(), so nThreads = 22 has no effect at all. BPPARAM is built from nThreads only when it is NULL (iterativePermutations.R:258-260, spatialCorrelation.R:507-509 and 238-240), and the worker count comes from BiocParallel's default (detectCores() - 2, which is 18 here). (2) Line 430 says spatialCorrelationGeneExp returns adjusted p-values. In fact spatialCorrelationGeneExp calls p.adjust inside the per-gene lapply (lines 836-842), one value at a time, so its p-values are unadjusted (B01 is correct). Context the finding omits: the vignette code and inst/scripts/visiumKidneySpatialCorrelation.R:180 actually call spatialCorrelationGeneExpIterPermutations. That function adjusts across genes, and the stored kidneyCorrelation p-values are adjusted. So the vignette's significance calls are on adjusted p-values, and the line-430 error is naming the wrong function. Separately, spatialCorrelationGeneExp passes BPPARAM = NULL to spatialCorrelation (line 830), so it ignores any user-supplied BPPARAM.

**Suggested fix:** Correct the comments. Consider parallelising over genes, which would scale much better than parallelising over B.

## B30: Assorted documentation errors (pseudo-count, threshold semantics, percentages, copied description, dataset size)

- **Location:** `R/packageFunction.R:51`
- **Category:** docs-code-mismatch
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** threshold() doc says zeros become 0.001 (line 51), but the code uses 0.0001 (lines 69-71). numPixelInThresh is documented as 'above the threshold in both experiments' (line 132), but the code keeps x > t1 | y > t2 (line 60). percentSimilarity and the dissimilarity columns are called percentages but are proportions in [0,1]. The vignette defines similar as |log2(y/x)| < b, but the code uses <= (line 269). The plotCorrelationGeneExp description (spatialCorrelation.R:1049-1054) is copied from WithinSample. viladomatCorrelation's nullCorGlobal doc says the correlations are with X, but they are with Y. The simRanPatternRasts doc says N = 5000 cells (R/data.R:37), but the datasets contain 1,201-1,381 cells, and its colData list omits 'type' and 'resolution'. The 1,000-point variogram subsample is undocumented.

**Evidence (auditor):** WD/07_similarity_plots.R: threshold() output gives x '1e-04 2 3'. Sum of num_cell per simRanPatternRasts dataset ranges from 1201 to 1381. colData names: num_cell cellID_list type resolution geometry sample_id.

**Verifier reasoning:** All eight sub-claims reproduced, and man/ matches the roxygen source. (a) threshold() docs say zeros become 0.001 (packageFunction.R:51), but the code sets 0.0001 (lines 69-71). (b) numPixelInThresh is documented as 'above the threshold in both experiments' (line 132), but the code keeps x > t1 | y > t2 (line 60). (c) percentSimilarity and the two dissimilarity columns are documented as percentages (lines 126-128) but are proportions that sum to 1 (line 270). (d) Vignette line 406 defines similar as |log2(y/x)| < b, but the code uses >= -b & <= b (line 269), so exact two-fold pixels count as similar; the roxygen @param input text (-b <= log2(y/x) <= b) matches the code. (e) plotCorrelationGeneExp's @description (spatialCorrelation.R:1049-1054) is a verbatim copy of WithinSample's (852-857). (f) nullCorGlobal is documented as correlations 'between the permutations and X' (lines 196-198), but the code computes cor(permutations, Y) (line 323). (g) R/data.R:37 says each dataset consists of N = 5000 cells within the kidney region, but datasets hold 1,201-1,381 cells. The 5000 appears to be the count sampled in the bounding box before the kidney mask, so the sentence is wrong rather than the number being invented. The colData list (lines 79-80) also omits 'type' and 'resolution'. (h) The 1000-point variogram subsample (spatialCorrelation.R:251-255) appears only in a code comment and is mentioned in no roxygen, man, vignette or README text.

**Suggested fix:** Update the roxygen text to match the code (or the code to match the intent), and regenerate man/.

## B31: Same seed for every gene and both directions: shared permutations and noise, order-dependent p-values, and stage 2 recomputing stage 1

- **Location:** `R/spatialCorrelation.R:831`
- **Category:** reproducibility
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** All genes and both directions use one seed: spatialCorrelation.R:542/:548/:831, iterativePermutations.R:50. viladomatCorrelation runs set.seed(seed) and then sample() (:242, :286-288), so every gene and both the permute-X and permute-Y runs reuse the same B permutation index vectors. matchingVariograms then reseeds with seed + i (:294, :99), so every gene and direction also gets the same standard-normal noise draws, scaled by per-permutation coefficients.

Monte Carlo errors are therefore correlated across genes. The correlation is strong for genes with similar spatial patterns (p-value correlation across seeds 0.975 vs 0.03 with independent seeds) and weak for unrelated genes. This positive dependence is undocumented.

In spatialCorrelationGeneExpIterPermutations, each stage restarts from the same seed, so stage k recomputes all nPermutations[k-1] earlier permutations exactly: 10% wasted work per carried-forward gene with c(100, 1000). The final p-value also reuses the nulls that were used for screening.

Separately, p-values change when the source's column (pixel) order changes, because shared_pixels follows source order (:785). This holds for any seeded permutation scheme, not because of seed sharing. The shared seed does make results invariant to gene order and gene subsetting. A fix should key streams on gene name and direction (not row index) and reuse the stage-1 nulls, or document the shared permutations.

**Evidence (auditor):** WD/05_rng.R: permutation #1 indices are identical for g2 and g3, and for permute-X and permute-Y. WD/08_pixels.R: shuffling the source columns leaves r unchanged but moves pValuePermuteX from 0.6 to 0.4 (B=5); shuffling the target columns gives identical results. WD/06_iter.R: the first 10 null correlations of the B=20 run are identical to the B=10 run.

**Verifier reasoning:** Code check. Both directions get the same seed (R/spatialCorrelation.R:542 and :548). Every gene gets the same seed (:831 in spatialCorrelationGeneExp, and R/iterativePermutations.R:50 and :324 in the iterative path). viladomatCorrelation calls set.seed(seed) (:242) and draws all B permutations with sample() (:286-288) before any gene-specific RNG use. It then calls matchingVariograms with seed + i (:294), which reseeds at :99. Because N (the number of shared pixels) is the same for every gene, the permutation index vectors and the standard-normal noise come out identical for every gene and both directions. I confirmed this directly.

Shared Monte Carlo error is real, but it is strong only for genes with similar patterns. For two near-identical genes, the p-values across 30 seeds correlate 0.975 with the shared seed versus 0.033 with independent seeds. For unrelated genes the figures are 0.29 versus -0.12, which is not significant at n=30.

Stage-2 recomputation is confirmed end to end. With the default c(100, 1000), the first 100 nulls of stage 2 are exactly the stage-1 nulls, so about 10% of stage-2 work is repeated for the genes carried forward. The final p-value also reuses the same nulls that did the screening.

What is wrong is the mechanism behind the column-order effect. Results do change when the source columns are shuffled, but this is not a result of seed sharing. It reproduces with a single gene, and with viladomatCorrelation alone in one direction. Any seeded permutation scheme behaves this way, including the per-gene substreams the finding suggests.

The shared seed also has an unmentioned upside: per-gene results do not depend on gene (row) order or on which genes are included (verified). A fix keyed on gene index would lose that property; streams keyed on gene name would keep it.

**Suggested fix:** Derive per-gene and per-direction streams (for example seed, gene index and direction mapped to L'Ecuyer substreams), or document that permutations are shared on purpose. In the iterative scheme, reuse stage-1 nulls and generate only the extra permutations.

## B32: matchingVariograms ignores i and defaults to seed = 0, so every permutation gets identical noise when called as in its own example

- **Location:** `R/spatialCorrelation.R:97`
- **Category:** api
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** The exported matchingVariograms(..., i, seed = 0) never uses i and always calls set.seed(seed) (R/spatialCorrelation.R:96-99). Called as in its own roxygen example (:89-92, no seed), every permutation reuses the same standard-normal noise vector for each delta index, so the B outputs are not independent null draws.

In tests this did not reduce the null's spread: SD ratio fixed/varying seed was 1.02 on the quakes example and 1.14 on white noise, not significant. It did distort delta* selection, which concentrated on a few grid values. viladomatCorrelation passes seed + i (:294), so the main pipeline is unaffected.

Fix: drop i or use it, stop reseeding inside the helper, and fix the example.

**Evidence (auditor):** WD/14_misc2.R: 'recovered noise vectors for permutation i=1 and i=2 identical: TRUE'.

**Verifier reasoning:** The core claim is true. matchingVariograms (R/spatialCorrelation.R:96-97) takes i but never references it. Its seed defaults to 0 and it calls set.seed(seed) at :99. The roxygen example (:89-92) calls it without a seed, so every permutation gets the same standard-normal noise vector for each delta index k, scaled by that permutation's sqrt(|beta1|). I confirmed this by tracing rnorm and by recovering the noise algebraically. viladomatCorrelation avoids the problem by passing seed + i (:294), so the main pipeline is unaffected.

The stated consequence, that this 'reduces the variability of the null', did not reproduce:
- Quakes depth, the example's data: the null-correlation SD barely changes.
- Nugget-dominated data: the SD with the fixed seed is equal or larger, and the difference is not significant at B=80.

What I did see is that the permutations are not independent draws, and that the selected delta* collapses onto a few grid values. The shared per-k noise biases the variogram-fit residuals the same way in every permutation.

**Suggested fix:** Remove i or use it (seed + i), or do not reseed inside this helper. Fix the example.

## B33: plotCorrelationGeneExp drops negative values, mishandles NA in one direction, and prints an unrounded p-value

- **Location:** `R/spatialCorrelation.R:1145`
- **Category:** api
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** plotCorrelationGeneExp (R/spatialCorrelation.R:1107-1160) has several problems:
- xlim(0, max) and ylim(0, max) (:1145-1146) drop every point with a negative value. ggplot2 reports this only through its generic 'Removed N rows ... outside the scale range' warning at print time.
- case_when(pY > pX ~ pY, .default = pX) returns pX when pY is NA but NA when pX is NA, which contradicts the documented 'greater of the two'.
- p_E is printed unrounded (e.g. 0.333333333333333).
- A gene missing from the results table silently gives 'r = NA p_E = NA'.
- labs(fill = 'Data') triggers an 'Ignoring unknown labels' message on every call. guides(color = ...) refers to an unmapped aesthetic and is ignored without any message.

The dropped points only occur with assays that contain negative values.

**Evidence (auditor):** WD/09_plotcor.R: with z-scored data, 'points supplied: 231 ; points drawn: 113'. With pY NA and pX 0.03 the title shows 'p_E = 0.03'. Another title reads 'p_E = 0.333333333333333'.

**Verifier reasoning:** Every listed behavior reproduces:
- xlim/ylim(0, max) at R/spatialCorrelation.R:1145-1146 turns negative values into NA, and those points are not drawn.
- case_when(pY > pX ~ pY, .default = pX) at :1128-1131 returns pX when pY is NA, but NA when pX is NA. That is asymmetric and does not match the documented 'greater of the two'.
- The p-value is pasted unrounded (:1157).
- A gene missing from the results table gives an 'r = NA p_E = NA' title with no error.
- labs(fill = 'Data') prints an 'Ignoring unknown labels' message on every call.

Two details are wrong:
- Removal is not silent. Printing the plot emits ggplot2's generic warning 'Removed 164 rows containing missing values or values outside the scale range (`geom_point()`)', though it does not say why.
- The 'Ignoring unknown labels' message comes only from labs(fill). guides(color = ...) is ignored without any message.

With normal non-negative assays all points are drawn, so the dropped points are an edge case for scaled or residual assays.

**Suggested fix:** Use the data range (or coord_equal without forced limits), take pmax(pX, pY) with NA propagated, use signif() for the p-value, check that geneName exists, and remove the unused guide and label.

## B34: spatialSimilarity with negative values gives proportions that sum above 1 and NA pixel IDs; no shared pixels gives NaN silently

- **Location:** `R/packageFunction.R:269`
- **Category:** correctness
- **Severity:** auditor low, verifier low
- **Verdict:** confirmed

**Description:** log2(y/x) is NaN when x and y have opposite signs. Subsetting a data.frame with a logical index that contains NA returns NA rows, so those pixels are counted as similar, dissimilarX and dissimilarY at once (lines 269-276). With no shared pixel names, the function returns percentSimilarity = NaN without a warning. Negative values are never checked.

**Evidence (auditor):** WD/16_similarity_negative.R: 'percentSimilarity + percentDissimilarityX + percentDissimilarityY = 1.052 (should be 1)'; 'NA pixel IDs inside similarPixelID: 6 of 225'. WD/12_misc.R: percentSimilarity NaN, numPixelInThresh 0.

**Verifier reasoning:** Code (R/packageFunction.R): log2(y/x) at :265 is NaN when x and y have opposite signs. Zeros are replaced by +0.0001 at :70-71, so a zero paired with a negative value also gives NaN. The comparisons at :269 and :272-273 then yield NA, and logical row-subsetting a data.frame with NA returns an all-NA row. Each NaN pixel is therefore counted in similarPixels, dissimilarPixelsX and dissimilarPixelsY at once, so the three proportions (:270, :275-276) sum to more than 1 and the ID lists contain NA. No negative-value check exists.

With no shared pixel names, getGenePixelDF returns 0 rows and quantile() gives NA thresholds. The minPixels check (0 < 0.1*0) is FALSE, so the code goes on to compute 0/0 = NaN with no warning or error.

One nuance, which does not contradict the finding: the negative-value case does emit R's generic 'NaNs produced' warning from log2. Severity stays low because negative values are outside the domain of a fold-change metric. Non-negative data behaves correctly.

**Suggested fix:** Require non-negative input (or document the behaviour), subset with which(), and stop when no pixels are shared.

## B35: Input-validation gaps, and errors print()ed to stdout instead of raised as warnings

- **Location:** `R/spatialCorrelation.R:586`
- **Category:** api
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** Entry points do not validate their inputs:
- Lists longer than 2 are accepted and the extra elements ignored (spatialCorrelation.R:783-784, packageFunction.R:190-191).
- The default assayName = 1 silently pairs different assays when the two objects order their assays differently.
- nPermutations = 0 fails internally (1:0 = c(1, 0)) and gives a printed error plus NA. nPermutations = 1 gives p of 0 or 1, and p = 0 is common. The iterative function accepts 1.
- delta is never range-checked. A negative delta segfaults R in locfit, delta = 0 gives a printed error plus NA, and delta > 1 is silently accepted.

spatialCorrelation catches every error and print()s it to stdout (:586). Callers cannot catch it as a condition, suppressWarnings() and suppressMessages() do not hide it, and only capture.output() or sink() can. The NA row it returns has logical NA in the list columns (deltaStarX, nullCorrelationsX), so rbind'ed results have mixed element types.

**Evidence (auditor):** WD/12_misc.R: '3-element input accepted'; 'nPermutations = 0 -> pValuePermuteX = NA'; 'nPermutations = 1 -> pValuePermuteX = 1'. WD/14_misc2.R: assay orders 'pixelval log' vs 'log pixelval' are compared with r = 0.9974 and no warning. WD/03_errors.R: class of deltaStarX per row is 'numeric' 'numeric' 'logical'.

**Verifier reasoning:** All the validation gaps reproduce:
- A 3-element list is silently accepted by spatialCorrelationGeneExp (input[[1]]/input[[2]] at R/spatialCorrelation.R:783-784) and by spatialSimilarity (R/packageFunction.R:190-191).
- The default assayName = 1 (:801-803) silently compares different assays when the two objects store them in different orders. This matches the documentation ('first assay') but is a trap.
- nPermutations = 0: 1:B becomes c(1, 0), so X.randomized[[0]] fails. The error is caught and printed, and the result is NA.
- delta is never range-checked.
- The error handler print()s the condition to stdout (:586) and returns an NA row whose deltaStarX and nullCorrelationsX entries are logical NA instead of numeric or matrix values.

Three details are wrong or understated:
1. 'nPermutations = 1 gives p = 1' is not general. p = extreme/1 is 0 or 1, and was 0 for 3 of 4 test genes. A p of 0 is the more worrying outcome. The iterative function also accepts c(1, 2), because its check only requires values > 0.
2. 'Cannot be silenced' is overstated. capture.output() or sink() silence the printout, but suppressWarnings() and suppressMessages() do not, and neither tryCatch nor withCallingHandlers sees a condition.
3. The delta gap is worse than described. delta = -0.2 segfaults R inside locfit and the process aborts. delta = 0 gives a printed error and NA. delta = 1.5 or 5 is silently accepted.

**Suggested fix:** Validate inputs at the top of each entry point (length(input) == 2, assay present by name in both objects, nPermutations >= 1, 0 < delta <= 1). Replace print(cond) with warning() and a status column, and keep column types consistent.

## B36: getGenePixelDF densifies the entire assay twice per gene (quadratic cost); output grown with rbind in a loop

- **Location:** `R/packageFunction.R:25`
- **Category:** performance
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** getGenePixelDF (R/packageFunction.R:25-26) calls as.matrix() on the entire assay before extracting one row, once per sample per gene.

For Matrix-class assays this copies the full G x P matrix on every call: G*P*8 bytes per side and O(G^2 P) total time for spatialSimilarity. At G=12,800, P=4,000 that is about 63 ms and 391 MB per side for a dgCMatrix, or about 26 minutes in total versus about 2.8 with direct row extraction. This applies to both dgCMatrix and dgeMatrix, which is what SEraster and the vignette CPM/lognorm recipes produce.

For base R matrices as.matrix() is a no-op, so nothing extra is allocated and cost is linear in G.

The result tables are also grown with rbind inside the loop (:241, :278, :294). That is quadratic but a minor cost at realistic gene counts.

Fix: extract assay[gene, pixels] directly, or better, subset both assays to the shared pixels once and vectorise over genes.

**Evidence (auditor):** WD/11b_densify.R: with G=12,800 genes and P=4,000 pixels, densifying takes 0.117 s per call (391 MB) vs 0.007 s for direct row extraction, i.e. 49.9 min vs 2.99 min per run. WD/11_similarity_perf.R: G=800 takes 27.5 s (34.4 ms/gene, up from 23 ms/gene at G=100).

**Verifier reasoning:** The code is as described. R/packageFunction.R:25-26 calls as.matrix(assay(...)) on the whole assay before subsetting [gene, pixels], once per sample per gene. The output tables are grown with rbind inside the loop (:241, :278, :294).

The claim 'sparse or not' is wrong. For a base R matrix, as.matrix() is a no-op that copies and allocates nothing. Per-call cost then stays flat as G grows, so total time is linear in G; speKidney rasterizes to a base matrix and is unaffected.

For Matrix-class assays the whole G x P matrix is copied on every call, giving O(G^2 P) total time and G*P*8 bytes allocated per side per gene. This covers both dgCMatrix and the dense dgeMatrix class. These classes are what the real-data vignettes produce:
- SEraster output from a dgCMatrix input is a dgCMatrix.
- The AKI vignette's 'CPM' assay is a dgCMatrix.
- The MERFISH vignette's log10(x + 1) 'lognorm' assay is a dgeMatrix.

So the problem matters in practice for real data, but not 'regardless of class'. My magnitude at G=12,800, P=4,000 is about 26 minutes of getGenePixelDF time for sparse input versus 2.8 for base. The finding's 49.9 versus 2.99 minutes is about 2x my sparse figure but the same order.

The rbind growth is real but minor at these sizes.

**Suggested fix:** Extract assay[gene, pixels] directly, or better, subset both matrices to the shared pixels once and vectorise over genes. Collect results in a list and rbind once at the end.

## B37: Redundant bplapply extraction passes, a new MulticoreParam per call, and a hard-coded, undocumented variogram subsample (N_s = 1000)

- **Location:** `R/spatialCorrelation.R:298`
- **Category:** performance
- **Severity:** auditor low, verifier low
- **Verdict:** partially-confirmed

**Verified description:** viladomatCorrelation's two bplapply calls at R/spatialCorrelation.R:298-313 only extract list elements. vapply gives identical results in about 0.0002 s. Measured inside the real function (N=273, B=100, 9 deltas), the two calls cost:
- about 0.02 s per direction at the default nThreads=1 (about 0.4% of the time);
- 0.14-0.25 s per direction at nThreads=4 (about 10-13% of the time), because each call forks 4 children.
The claimed 0.053 s/0.104 s figures wrongly include the cost of building a MulticoreParam.

A MulticoreParam object is built once per gene, not per gene and direction. This happens at 507-509 only because spatialCorrelationGeneExp passes BPPARAM = NULL at line 830. That also discards its own object (777-779) and any BPPARAM the user supplies. viladomatCorrelation's construction at 238-240 runs only when it is called directly. Building the object costs one full gc(), through BiocParallel:::.snowCoresMax and showConnections(). That is about 0.05 s in a bare R session and about 0.3 s with SpatialExperiment loaded, roughly 3-10% of per-gene time.

The repeated fork/teardown is a separate issue from building the object. No caller ever calls bpstart, so each of the 6 bplapply calls per gene forks nThreads children (TransientMulticoreParam). This also happens when one object is reused (IterPermutations), and never at the default nThreads=1 (SerialParam fallback).

N_s = 1000 (line 253) is hard-coded and undocumented. One subsample per direction is used consistently for the target variogram and all permuted variograms, while locfit fits all N points. It runs correctly at N=1500.

Calling bpstart once does not help by itself: tested, it slowed the extraction calls and did not reduce total time. Replacing the bplapply calls with vapply, and passing BPPARAM through at line 830, are the effective fixes.

**Evidence (auditor):** WD/17_bplapply_overhead.R: the 2 extraction bplapply calls cost 0.053 s (1 worker) and 0.104 s (4 workers) per viladomatCorrelation call, versus about 0 for vapply. That is 17.7-34.7 min per 10,000 genes. WD/14_misc2.R: N=1500 run OK.

**Verifier reasoning:** Three of the claims hold. (1) The two bplapply calls at R/spatialCorrelation.R:298-313 only pull elements out of `output`, and vapply gives identical results. (2) `N_s <- 1000` at line 253 is hard-coded. It appears nowhere in man/, the vignettes or docs/; only the inline comment at lines 251-252 mentions it. locfit (line 109) fits all N points. The subsample is used consistently: at N=1500, all 13 variog calls received the same 1000 coordinates and the run finished cleanly. (3) Workers are forked repeatedly when nThreads > 1.

The mechanism and the numbers are partly wrong. (a) A MulticoreParam is not built 'for every gene and direction'. viladomatCorrelation's own construction at 238-240 is never reached from spatialCorrelation, because spatialCorrelation passes its BPPARAM down. One object per gene is built at 507-509 only because spatialCorrelationGeneExp passes `BPPARAM = NULL` at line 830. That throws away its own object (777-779) and also any BPPARAM the user supplies. The IterPermutations and WithinSample paths build a single object in total. (b) The repeated forking does not come from building new param objects. Any bplapply on an un-started MulticoreParam with workers > 1 runs as a TransientMulticoreParam and forks one child per task (`mcparallel` in `.send_to`). That is nThreads forks per call and 6 calls per gene. It happens the same way when one object is reused. With the default nThreads = 1, BiocParallel:::.bpinit falls back to SerialParam and nothing is forked. (c) The building cost that does exist is hidden: MulticoreParam() calls .snowCoresMax, then showConnections(), then a full gc(). (d) The benchmark timed MulticoreParam construction together with the two calls, so the cost is assigned to the wrong step. Measured inside the real function, the extraction costs differ from the claimed figures. (e) The suggested bpstart-once fix does not help on its own. With a started backend, the extraction calls get slower because the closures and their environments are serialised to the workers; the vapply change is what removes the cost. This is a performance issue only and does not change results, so severity is low.

**Suggested fix:** Use vapply for the extraction, start the BiocParallel backend once per top-level call (bpstart/bpstop), and expose and document the subsample size.

## R CMD check summary (baseline, before any changes)

R CMD build --no-build-vignettes succeeded: STcompare_0.1.0.tar.gz is 48.6 MB. R CMD check --no-manual --ignore-vignettes --no-build-vignettes (R 4.5.2, aarch64-apple-darwin, _R_CHECK_FORCE_SUGGESTS_=false, 52 s) gave Status: 1 ERROR, 3 WARNINGs, 2 NOTEs, plus INFO: installed size 33.8 Mb (extdata 32.3 Mb, data 1.2 Mb).

ERROR (examples): the simRanPatternRasts example (R/data.R:97) calls assays() without a namespace: 'could not find function "assays"'. The check stops at the first failing example, so I ran the rest individually. spatialCorrelationGeneExpWithinSample also ERRORS (spatialCoords called on a list; R/spatialCorrelation.R:969-972). The following run OK: spatialCorrelation (3.64 min), spatialCorrelationGeneExp (9.9 s), spatialSimilarity, viladomatCorrelation and speKidney, plus linearRegression, matchingVariograms, pixelClass, plotCorrelationGeneExp and savePlots inside the check. savePlots attaches patchwork, gridExtra and ggplot2.

WARNINGs:
(1) Non-ASCII characters in R/iterativePermutations.R (an em dash at line 302).
(2) Dependencies in R code: '::' import not declared from 'sf'; library() calls to undeclared 'ggplot2' 'gridExtra' 'patchwork'; library() calls in package code.
(3) Rd usage: spatialCorrelationGeneExpWithinSample documents deltaX/deltaY but its argument is 'delta'; spatialSimilarity's 'verbose' is undocumented.

NOTEs:
(1) R code: undefined globals cor, cor.test, dist, fitted, lm, median, na.omit, quantile, rnorm, p.adjust.methods, combn, unit and plot_layout, plus NSE variables x, y, X, Y, XGexp, YGexp, color, fill and pValuePermuteX. Fix with importFrom(stats/utils).
(2) checkRd: 89 'Lost braces in \itemize; \value handles \item{}{} directly' across 11 Rd files, plus lost braces at savePlots.Rd:47.

Not exercised by this check: under --as-cran (_R_CHECK_LIMIT_CORES_=TRUE), the nThreads = 5 examples would also error ('BiocParallel workers must be <= 2'). Vignettes were not built: they need the undeclared and uninstalled MERINGUE and scatterbar, and download data from Zenodo/10x. Logs: WD/check/check_full.log and WD/check/STcompare.Rcheck/00check.log, where WD = /private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/bugs
