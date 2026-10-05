# STcompare: review of the statistical method and algorithmic speedups

Scope: investigation only. Nothing in the repository was modified. All scripts, data, and outputs are under
`/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/method/`
(called `$W` below). Environment: R 4.5.2, locfit 1.5-9.12, geoR 1.9-6, BiocParallel 1.44.0, Apple M1 Ultra. Timings are single-threaded (`VECLIB_MAXIMUM_THREADS=1`) unless noted.

---

## 0. Key findings (short version)

**Statistics: problems with evidence**

1. **The multiple-testing correction in `spatialCorrelationGeneExp()` does nothing.** `p.adjust()` runs inside the per-gene `lapply` on a one-row data frame (`R/spatialCorrelation.R:836-843`), and `p.adjust(p, "BH")` on a single value returns `p`. The shipped results confirm it: in `kidneyCorrelationNoIter.RData`, the stored `pValuePermuteX` equals the raw `extreme/B` exactly, and BH-adjusted values would differ by up to 0.044. On that dataset (AKI Visium, 1046 genes), the vignette rule "both p < 0.05" gives **757** genes with the unadjusted values but **709** with BH on correctly computed p-values. The iterative function does apply BH across genes (`R/iterativePermutations.R:343-344`). So the two functions return columns with the same names but different meanings (raw versus BH-adjusted), and both overwrite the raw p-values.
2. **The p-value can be zero.** The code uses `extreme / B` with a strict `>` (`R/spatialCorrelation.R:326-327`). This copies the paper and its prototype code, but it is anti-conservative: under exchangeability, P(p ≤ 0.05) = 6/101 = 0.059 when B = 100. In the shipped results, 573/1046 genes (kidney, B = 100), 444/1046 genes (kidney, iterative, B = 1000) and 366/483 genes (MERFISH) have p = 0. The docs say "the smallest p-value is 0.01" (`:169, :358, :647`). MERINGUE, from the same lab, already uses (b+1)/(n+1).
3. **The iterative screening rule does not match its documentation.** The code threshold is `t = 100·alpha/nPermutations[k]` (`R/iterativePermutations.R:64`). The docs say `alpha/nPermutations[k]` (`:90-94`), which is 100 times stricter. Because every round reuses the same seed, the B = 1000 round contains the B = 100 round as its first 100 surrogates. I verified this bit-for-bit (`$W/nesting.R`). The scheme is therefore a two-stage sequential Monte Carlo procedure. That makes it roughly valid, apart from the p = 0 problem, but it recomputes the first 100 surrogates and saves little when most genes are significant: 757/1046 kidney genes and 397/483 MERFISH genes went to stage 2.
4. **Every gene uses the same random inputs.** Every call starts with `set.seed(seed)` (default 0). All genes and both directions therefore share the same 1000-point variogram subsample and the same B index permutations: `sample(X)` permutes indices the same way for any X of the same length (verified). Each permutation's noise comes from `set.seed(seed + i)`, which is also shared. The Monte Carlo errors are therefore correlated across genes. Also, `BiocParallel::bplapply()` runs workers with `RNGkind("L'Ecuyer-CMRG")`, so `set.seed(seed+i)` inside `matchingVariograms()` seeds L'Ecuyer, not Mersenne-Twister (measured). Reproducing the results exactly depends on this BiocParallel behaviour.
5. **The "both p < 0.05" rule and the plotting function's "larger of the two" rule are conservative.** On the package's own null simulation (`simRanPatternResults.RData`, 4950 independent field pairs, both directions available), the type-I error at 0.05 is **2.1%** with "both < 0.05" and 5.7% with "either < 0.05". The two directions' p-values have Spearman ρ = 0.94. Running BH separately on X and Y and then intersecting the two rejection sets has no formal FDR guarantee. BH on max(pX, pY) does.
6. **The default δ grid `seq(0.1, 0.9, 0.1)` is often too coarse at the small end.**
   - MERFISH stored results (vignette grid starting at 0.01): 23.5% of all permutations chose δ* = 0.01, the smallest value, and 87/483 genes had more than 50% of permutations at that boundary.
   - My synthetic random fields (same recipe as `simRanPatternRasts`): 84-91% of permutations chose δ* = 0.1, the lower bound of the default grid.
   - The paper itself adapts Δ to the autocorrelation range.
7. **Calibration under the null.** Exact replica, (b+1)/(B+1), B = 100, N = 273, 500 pairs per scenario:
   - Gaussian fields: type-I error 0.030 at 0.05 (conservative).
   - Heavy-tailed marginals: 0.048 at 0.05 but 0.016-0.018 at 0.01. This is suggestive tail inflation (SE about 0.0045). Smoothing plus Gaussian noise does not preserve the marginal distribution.
   - brainsmash offers a `resample=True` option that rank-matches surrogates to the observed values.
8. **The analytic modified t-test is an excellent fast reference.** I used `SpatialPack::modified.ttest` (Clifford, Richardson & Hémon 1989; Dutilleul 1993). On 4950 null pairs, its p-values have Spearman ρ = **0.966** with STcompare's stored permutation p-values. It takes 19 ms per pair at N ≈ 273 and 1.4 s at N ≈ 5000, versus about 90 s per gene for STcompare.
9. **Spatial similarity caveats** (demonstrated in `$W/similarity_demo.R`):
   - Identical patterns that differ by a global factor of 2.5 get **0%** similarity.
   - Two independent Poisson draws from the *same* spatial rate (mean 2 counts per pixel) get only **56%** similarity.
   - Replacing zeros with 1e-4 is not scale-invariant: the docs say 0.001, and the replacement exceeds the data when values are below about 1e-4.
   - Genes skipped by `minPixels` report `numPixelInThresh = 1`, which is a bug (`R/packageFunction.R:249`).

**Speed (the main enabler for C++)**

- **What the smoother really is.** I read the locfit C source and verified numerically. `locfit(lp(nn = δ, deg = 0), kern = "gauss")` followed by `fitted()` works in two steps:
  - It computes untruncated Gaussian Nadaraya-Watson means at the vertices of an adaptive tree. Each vertex uses bandwidth = distance to its k-th nearest data point, with k = ⌊Nδ⌋.
  - It then **bilinearly interpolates** those vertex values to the data points, because deg = 0 means no derivatives are stored.

  So **S_δ = W_δ V_δ exactly**: V_δ is m × N (vertex weights) and W_δ is N × m and sparse (about 4 non-zeros per row). The vertex count m depends on δ and the tissue shape but **not on N**: m ≈ 13 at δ = 0.9, 120 at δ = 0.1, about 1460 at δ = 0.01. The reconstructed operator matches `fitted()` to **2e-16**. It depends only on the coordinates and δ, so one operator per δ serves every gene, both directions, and every permutation.
- **The binned variogram is a fixed pair-to-bin accumulation.** I precomputed the pair list once per coordinate set. A C++ loop over that list is **bitwise identical** to `geoR::variog` and takes 0.18 ms versus 11.3 ms per call (63×).
- **Measured end to end.** I built an exact replica that combines (a) the precomputed operators, (b) the pair-list variogram, (c) closed-form OLS, (d) subsample-only smoothing during δ selection and (e) reuse of all of these across directions and genes.
  - Against the package with the same seed, it gives **identical δ\* sequences and p-values** and surrogates equal to about 1e-11 absolute (1e-14 relative) at N = 273, 1003, 2005 and 4984. The last three exercise the 1000-point subsample path.
  - The pure-R version is 12-16× faster than the package.
  - A small RcppArmadillo core (`$W/fastgene.cpp`) takes **1.19 s per gene** (2 directions × 100 permutations × 9 δ, N = 4984), of which 0.39 s is the R random-number generation kept for exact equality. The package takes about **91 s** for the same gene, so the C++ core is **about 76× faster single-threaded**, and gene-level threading comes on top of that.
- **Within-sample comparisons can reuse surrogates.** A gene's surrogates depend only on that gene and the coordinates, never on its partner, and the seed is the same everywhere. `spatialCorrelationGeneExpWithinSample()` therefore regenerates identical surrogates for every pair. Generating them once per gene is **exact** and saves a factor of (G − 1).
- **Stop permutations once the BH decision is settled.** I simulated a simplified sequential rule on the stored B = 1000 null sequences, in the spirit of MMCTest/SIMCTEST: stop once 99% Clopper-Pearson bounds clear the BH threshold. The 757 stage-2 kidney genes would need on average **308 draws instead of 1000**. Classic Besag-Clifford stopping saves only about 9% here, because most genes are strongly significant.

---

## 1. The method in Viladomat et al. (2014)

Viladomat, Mazumder, McInturff, McCauley & Hastie, *Assessing the significance of global and local correlations under spatial autocorrelation: a nonparametric approach*, Biometrics 70(2):409-418, doi:10.1111/biom.12139. The open-access copy is at [PMC4108159](https://pmc.ncbi.nlm.nih.gov/articles/PMC4108159/). The prototype R code is at [hastie.su.domains/Papers/biodiversity/viladomat_code.R](https://hastie.su.domains/Papers/biodiversity/viladomat_code.R), saved as `$W/viladomat_code.R`.

**Problem.** The paper tests the Pearson correlation between two spatially autocorrelated fields X_s and Y_s observed at N locations. The naive t-test is badly anti-conservative. Clifford et al.'s modified t-test corrects it through an effective sample size.

**Algorithm.** Each of B = 1000 draws runs these steps:

1. Randomly permute X over the locations. This destroys both the autocorrelation and the association with Y.
2. Smooth the permuted field with a local-constant kernel regression. The bandwidth at each location is the distance to its k = ⌊Nδ⌋ nearest points (the paper's eq. A.3), with a Gaussian kernel, for each δ in a grid Δ.
3. Compute the empirical variogram of each smoothed field. Regress the target variogram γ̂(X) on γ̂(X_δ) to get α̂_δ and β̂_δ.
4. Choose δ* as the value that minimises "the sum of squares of the residuals of the fit" (paper wording).
5. Output X̂ = |β̂|^{1/2} X_δ* + |α̂|^{1/2} Z, where Z is iid N(0,1).
6. Compute r(X̂, Y).

The p-value is p̂ = (1/B) Σ I(|r_b| > |r̂|). The paper uses a **kernel-smoothed variogram** (eq. A.1: Gaussian kernel exp{−(2.68x)²/2h²}, evaluated at 100 uniformly spaced distances), truncated at the 25th percentile of pairwise distances. The paper does not mention subsampling locations.

**Grids and results.** Δ is data-dependent:
- Biodiversity application (N = 19,926): Δ ∈ (0.005, 0.785).
- Simulations: (0.1,…,0.9) for φ = 0.3, (0.03,…,0.074) for φ = 0.1, (0.013,…,0.027) for φ = 0.05.

Type-I error in the simulations is about 0.048-0.052 versus 0.093 for Clifford's method in one setting (Table 1). The paper permutes only one field ("we are free to pick the most convenient one") and reports no runtimes. Local correlations use kernel-weighted Pearson correlations with the same surrogates.

## 2. Implementation versus paper (with code citations)

STcompare's `matchingVariograms()` / `viladomatCorrelation()` are a near-verbatim port of the prototype's `matching()` / `main()`. They add seeding, BiocParallel, a `B` parameter, both directions, and SpatialExperiment wrappers.

| Aspect | Paper | Prototype code | STcompare |
|---|---|---|---|
| Variogram estimator | Kernel-smoothed variogram (A.1), 100 distances | `geoR::variog(option="bin")`: 13 equal-width bins on [0, umax], classical estimator, pairs.min = 2 (geoR `.define.bins`, `binit` in `geoR.c:370-422`) | Same as prototype (`spatialCorrelation.R:114-117, 125-128, 272-274`) |
| Points used for variograms | All | Random subsample N_s = 1000 drawn once | Same; drawn after `set.seed(seed)`, so **identical for all genes and both directions** (`:242, :253-264`) |
| Max lag | 25th percentile of pair distances | 25th percentile of subsample distances | `maxDistPrctile = 0.25` (`:268-269`) |
| Smoother | Local constant, kNN bandwidth, Gaussian kernel, at each s_j | `locfit(lp(nn=δ, deg=0), kern="gauss", maxk=300)`, then `fitted()` | Same (`:109-112`). The default `ev` = adaptive tree plus bilinear interpolation (§5.1) |
| δ grid | Tailored to the autocorrelation range | `seq(0.1,0.9,0.1)` "user may have to experiment" | Default `seq(0.1,0.9,0.1)` (`:523-529`); vignettes use `c(0.01,0.05,0.1..0.9)` |
| δ* criterion | RSS of target ~ α + β γ(X_δ) | RSS between the variogram of the **noised** surrogate and the target (one extra variogram per δ; random) | Same as prototype (`:124-134`) |
| Surrogate | \|β\|^{1/2}X_δ* + \|α\|^{1/2}Z | Same | Same (`:124`) |
| B | 1000 | 1000 | Default 100; iterative 100 → 1000 |
| p-value | (1/B) Σ I(\|r_b\| > \|r̂\|) | `extreme/B` | Same (`:326-327`) |
| Direction | One field | X | Both, run separately with the same seed (`:537-548`); vignettes require both < 0.05; the plot shows the max (`:1127-1131`) |
| RNG | — | foreach/doMC, unseeded | `set.seed(seed)` in the parent; `set.seed(seed+i)` per permutation inside `bplapply`; same seed for every gene |

Comments on the specific choices:

- **1000-point subsample.**
  - Effect: variogram cost becomes O(10⁶) per field instead of O(N²). The same pairs are used for the target and the surrogates, which keeps the matching consistent.
  - Cost: extra variance in the target variogram, plus one random subsample shared by all genes and directions. On regular grids the full-data variogram can be computed exactly by FFT at similar cost (§5, item (i)), which would remove this variance (a method change).
  - geoR recomputes `dist()` on the subsample at every call (`variogram.R:103`): 3600 times per gene.
- **max.dist = 25th percentile.** This matches the paper. A quirk: geoR sets the last bin edge to `umax = max(u[u < max.dist])` and excludes pairs exactly at `umax` (`binit`: `dist < lims[ind]`).
- **δ grid.** The default reproduces the paper's grid for strong, long-range autocorrelation (φ = 0.3). For finer-scale patterns the optimum sits at or below 0.1:
  - Stored MERFISH results: 23.5% of permutations at δ = 0.01.
  - Kidney: δ* is mostly 0.05-0.1.
  - Synthetic random fields: 84-91% of permutations at δ = 0.1 with the default grid (`$W/beta_signs.R`).

  Choosing the boundary value means the surrogates are probably over-smoothed. That is consistent with the conservative calibration in §3.6.
- **locfit kernel and `maxk`.**
  - `kern="gauss"` is W(u) = exp(−(2.5u)²/2) with **no truncation** (`weight.c:31`, `GFACT 2.5` in `local.h:137`). Every data point gets weight, and the bandwidth is the k-th order statistic of the distances (`lf_nbhd.c` `nbhd()`/`compbandwid()`), with k = `(int)(n*nn+1e-12)` (`locfit.c:355`).
  - `maxk=300` does **not** cap the vertices at 300. It multiplies locfit's vertex allocation guess by maxk/100 (`ev_atree.c` `atree_guessnv`). I observed 1996 vertices at δ = 0.01 with N = 273.
  - The paper prints its kernel as exp(−2.5x²/2λ²) according to my transcription; locfit's actual kernel is exp(−(2.5x/λ)²/2). Treat locfit as authoritative.
- **Noise and scale terms.** For β, α ≥ 0, aX_δ + bZ has variogram a²γ_δ + b² at h > 0, so a = √β and b = √α match the regression. The `abs()` is a heuristic for negative coefficients. In my runs (4 synthetic genes, kidney A/B, 100 permutations each), β0 and β1 were never negative at δ* (0%). The post-noise RSS criterion probably steers away from such fits.
- **p-value.** Use (1 + #{|r_b| ≥ |r̂|})/(1 + B) (Davison & Hinkley 1997; North et al. 2002; Phipson & Smyth 2010). Ties are measure-zero for continuous surrogates, so `>` versus `≥` hardly matters. The +1 is what restores validity and prevents p = 0.
- **Two-sided test via `abs()`.** This is fine when the null is symmetric around 0, which holds here because surrogates are independent of Y. A one-sided option would be natural for replicate comparisons where r > 0 is expected.

## 3. P-values, screening, multiple testing, seeds (detailed)

### 3.1 The `p.adjust` no-op

`spatialCorrelationGeneExp()` calls `stats::p.adjust(output$pValuePermuteX, method = adjustMethod)` *inside* the per-gene `lapply` (`R/spatialCorrelation.R:836-843`), and `output` has one row. Verified: `p.adjust(0.03, "BH") == 0.03` and `p.adjust(0.03, "bonferroni") == 0.03` (`$W/nesting.R`). In the shipped `kidneyCorrelationNoIter.RData`, stored p equals raw `extreme/B` exactly (max |diff| = 0), while BH would change values by up to 0.044 (`$W/inspect_results.R`).

Decisions under the vignette rule "both < 0.05", for the same 1046 genes (`$W/analyze_stored.R`):

| p-values used | genes called significant |
|---|---|
| Raw, unadjusted (what the function returns) | 757 |
| BH per direction on raw `extreme/B` | 739 |
| BH per direction on (b+1)/(B+1) | 709 |
| BH on max of (b+1)/(B+1) | 709 |
| Bonferroni on max | **0** (1046 × 1/101 > 0.05) |

`spatialCorrelationGeneExpWithinSample()` applies no correction at all over its choose(G, 2) pairs.

The iterative function applies BH correctly across genes but writes the result into the same `pValuePermuteX/Y` columns. The two main entry points therefore return columns with identical names and different meanings. Recommendation: always return raw `pValuePermuteX/Y`, a combined `pValuePermute = max(pX, pY)`, and separate adjusted columns.

### 3.2 Zero p-values and resolution limits

`extreme/B` produces p = 0 with probability 1/(B+1) under the null, and it is the most frequent value in practice. Kidney iterative: 444 genes with pX = 0 at B = 1000. MERFISH: 366/483. BH treats p = 0 as always significant.

With (b+1)/(B+1) the floor is 1/(B+1), so the number of permutations needed depends on the multiplicity:
- Bonferroni is impossible at m ≈ 1000 even with B = 1000 (1046/1001 > 0.05).
- BH works only when many genes sit at the floor.
- For m = 20,000 genes and about 10 true discoveries, BH would need p ≈ 2.5e-5, so B ≥ 4e4. That is where sequential and tail methods matter (§5 (f)).

Fitting a Gaussian to Fisher-z of 100 surrogates gives p-values loosely consistent with the empirical B = 1000 values. For genes with b ≥ 5: Spearman 0.75-0.83, median ratio 1.07-1.15, IQR of the ratio 0.68-1.42. But some genes with b = 0 at B = 1000 get a Gaussian p up to 0.02, i.e. tails heavier than Gaussian for some genes (Shapiro p < 0.01 for 21-22% of genes; median |skew| 0.06-0.08). A parametric tail is therefore an option for ranking, not a drop-in replacement.

### 3.3 Iterative screening (`spatialCorrelationGeneExpIterPermutations`)

- **Threshold.** The code uses `t = (alpha/nPermutes)*100` (`:64`): 0.05 at B = 100, 0.1 at B = 50, 0.005 at B = 1000. The docs say `alpha/nPermutations[k]` (`:90-94`).
- **Selection rule.** A gene is carried forward only if `pX < t & pY < t` (`:65`).
- **Nesting.** Every round calls `set.seed(seed)`, so the B = 1000 run's first 100 surrogates are exactly the B = 100 run's surrogates. I verified this for B = 10 versus B = 20: the first 10 surrogates matched with max |diff| = 0. Stage 2 extends stage 1 rather than resampling, which makes the reported p-value a two-stage sequential Monte Carlo p-value.
- **Validity.** For u < t, P(reported ≤ u) ≤ P(stage-1 pX ≤ u, stopped) + P(continued, p₁₀₀₀ ≤ u) ≈ u. The procedure is roughly valid apart from the p = 0 issue. The final BH mixes resolutions (0.01 for stopped genes, 0.001 for continued ones). Even a stopped gene can carry pX = 0 when it was stopped because pY ≥ t.
- **Efficiency.** Stage 2 recomputes the first 100 surrogates (10% waste; caching them is exact). Pruning helps only when most genes are null: 757/1046 kidney genes and 397/483 MERFISH genes went to stage 2. Total draws per direction were 861,600 versus 1,046,000 for B = 1000 everywhere (an 18% saving).
- **Better replacements** (§5 (f)):
  - Besag-Clifford stopping (stop at h exceedances) makes null genes cheap, but saved only 9% on these strongly correlated datasets (83-96% of stage-2 genes reach B_max).
  - A decision-oriented rule that stops when the gene's status relative to the BH threshold is settled saves about 3× on exactly these genes (mean 308 versus 1000 draws; `$W/mmc_potential.R`). MMCTest (Gandy & Hahn 2014) does this with guaranteed resampling risk.

### 3.4 Two directions and the "larger of the two" rule

The paper permutes one field. STcompare computes both, which doubles the cost.

Under the null (package simulation, 4950 pairs × 2 directions):
- Spearman(pX, pY) = 0.94.
- Type-I error at 0.05: 5.7% with "either < 0.05", **2.1%** with "both < 0.05".

max(pX, pY) is a valid p-value. But "BH on X, BH on Y, then intersect" (the vignette rule on the iterative output) is a superset of "BH on max(pX, pY)" and has no FDR guarantee. On kidney the difference is small: 730 versus 728 genes, both using (b+1)/(B+1). Recommendation: either define one gene-level p-value, max(pX, pY), and run BH once, or permute one designated field as in the paper.

### 3.5 Seeds and RNG

- **Shared random inputs.** `set.seed(seed)` (default 0) runs in `viladomatCorrelation()` (`:242`) for every gene and both directions. Every gene and direction therefore shares:
  - the 1000-point subsample;
  - the B index permutations (sample of an equal-length vector);
  - the noise streams `set.seed(seed+i)`.

  These are common random numbers: each gene's p-value is valid on its own, but Monte Carlo errors are correlated across genes, and genes with near-identical patterns get near-identical nulls. Either use per-gene or per-direction streams (for example L'Ecuyer substreams via `parallel::nextRNGStream`, or a counter-based RNG keyed by gene, direction and permutation), or document the coupling as intentional. The subsample can and should stay shared; it is what makes the precomputation exact.
- **BiocParallel switches the RNG kind.** Inside `bplapply` (Serial and Multicore alike), `RNGkind()` is `"L'Ecuyer-CMRG"`, so `set.seed(seed+i)` seeds L'Ecuyer (measured, BiocParallel 1.44.0). Calling `matchingVariograms()` directly, as in its example, draws different noise. Exact reproduction therefore depends on BiocParallel's RNG policy.
- **Global side effects.** `set.seed()` inside package functions mutates the user's global RNG state.
- **Parallel settings.**
  - `spatialCorrelationGeneExp()` passes `BPPARAM = NULL` down (`:830`), so a user-supplied `BPPARAM` is ignored.
  - Two of the three `bplapply` calls per direction only extract list elements (`:298-313`), yet each one forks.
  - `MulticoreParam(workers=1)` costs about 0.29 s per call even with one worker (`$W/bp_overhead.R`).

### 3.6 Calibration checks (null data)

| Setting | Type-I @0.05 | Type-I @0.01 | Source |
|---|---|---|---|
| Package's simRanPattern (N ≈ 273, iterative, stored, zeros set to 0.01) | 0.039 (≤: 0.0485) | 0.0012 (≤: 0.0135) | `$W/inspect_sim.R` |
| Gaussian GRF, exact replica, (b+1)/(B+1), B = 100, 500 pairs | 0.030 | 0.004 | `$W/null_marginals.R` |
| Lognormal (exp(1.5 GRF)) | 0.048 | **0.018** | same |
| Zero-inflated Poisson(exp(GRF)) | 0.040 | **0.016** | same |
| Modified t-test, same pairs (G / LN / ZI) | 0.052 / 0.054 / 0.046 | 0.002 / 0.026 / 0.022 | same |
| Naive `cor.test` (G / LN / ZI) | 0.378 / 0.202 / 0.084 | — | same |

Interpretation:
- The method is roughly calibrated to conservative at 0.05; the mean null p is 0.55-0.58.
- There is suggestive inflation at 0.01 for heavy-tailed or zero-inflated marginals. Smoothing plus Gaussian noise does not reproduce outlier-driven structure.
- Worth a larger study. Candidate fixes: an option to rank-map surrogates onto the observed values (brainsmash `resample=True`), or a log/rank transform before testing.

## 4. Spatial similarity (fold change; `R/packageFunction.R`)

What the code does:
1. Keep pixels with `x > t1 | y > t2` (`:60`). The thresholds default to the 5th percentile of each gene (`:221, :226`).
2. Replace zeros by 1e-4 (`:69-71`). The docs at `:51` say 0.001.
3. Compute log2(y/x) and call a pixel similar if |log2(y/x)| ≤ `foldChange` (`:260-276`).
4. Return `NA` if fewer than `minPixels` × N pixels pass (`:239`).

Caveats, demonstrated in `$W/similarity_demo.R`:

| Synthetic case | percentSimilarity |
|---|---|
| y = x | 1.00 |
| y = 2.5·x (identical pattern, global scale) | **0.00** |
| x, y independent Poisson with the same spatial rate, mean 2 counts/pixel | **0.56** |
| Same, mean 100 counts/pixel | 1.00 |
| Values ~1e-5 with 50 zeros in y | Zeros become 1e-4, about equal to the data, so most zero pixels are called *similar* |
| 10 expressed pixels (skipped by `minPixels`) | NA, but `numPixelInThresh` is reported as **1** (bug: `dim(thresh)[1]` at `:249`) |

What this means:
- **Pattern and magnitude are conflated.** Any platform or library-size scale difference moves every log ratio, so normalization choices dominate the score. Consider a per-gene centring option, e.g. subtract the median log ratio.
- **Sampling noise biases similarity downward at low counts.** The score therefore depends on expression level and raster resolution. A reference such as "expected similarity if both share one rate" (Poisson or NB) would make values comparable across genes.
- **Thresholds behave inconsistently across genes.** The quantile thresholds become 0 for sparse genes, so their meaning varies by sparsity, and the per-gene denominator changes. `numPixelInThresh` is documented as "above threshold in both" but uses OR logic. `percentSimilarity` is a proportion, not a percentage. `foldChange` is a log2 threshold.
- **No uncertainty is attached and pixels are treated as independent.** There is no interval or test; pixel-level calls ignore spatial autocorrelation.
- **Efficiency.** `getGenePixelDF()` converts the *entire* assay to dense for every gene (`:25-26`), and results grow by `rbind` in a loop. This is O(G²·N) work. A vectorized version is exact and near-instant.

## 5. Algorithmic speedups

### 5.1 The enabling observation: what the smoother computes

Sources: locfit C code (`weight.c`, `lf_nbhd.c`, `ev_atree.c`, `ev_interp.c:64-80`, `startlf.c:147`). Checks: `$W/test_exact.R`, `$W/factorS.R`, `$W/test_factor.R`, `$W/vertex_counts.R`.

- **Linearity.** `fitted(locfit(y ~ lp(..., nn=δ, deg=0), kern="gauss"))` is linear in y for fixed coordinates. |f(a y₁ + b y₂) − (a f(y₁) + b f(y₂))| ≤ 1.3e-15. Rows of the operator sum to 1.
- **Structure.** The default evaluation structure is an adaptive tree. Cells split until side/h < `cut` = 0.8. Gaussian NW fits are computed at the vertices, each with its own kNN bandwidth. Because deg = 0, locfit stores no derivatives (`hasd = deg>0 | dc`), so fitted values are **multilinear (bilinear) interpolations** of the 4 corner values of each point's cell. Pseudo-vertices are interpolated, not fitted. Vertex values are stored relative to a parametric component (the global mean for deg = 0); this cancels because interpolation weights sum to 1.
- **Exact factorization.** S_δ = W_δ V_δ:
  - V_δ: Gaussian NW weights of the m fitted vertices (m × N, dense, untruncated).
  - W_δ: interpolation weights (N × m, about 4 non-zeros per row), extracted exactly by calling locfit's `sfitted` with unit coefficient vectors.
  - max |W V − S_unit-vector| ≤ 1.8e-16 for all 9 δ.
  - Building all 9 operators takes 0.5 s (N = 1003), 0.7 s (N = 2005), 2.1 s (N = 4984).
- **Vertex counts (fitted) are independent of N.**

| δ | 0.01 | 0.02 | 0.05 | 0.1 | 0.3 | 0.9 |
|---|---|---|---|---|---|---|
| m (N = 1003) | 1490 | 457 | 226 | 120 | 35 | 13 |
| m (N = 4984) | 1463 | 453 | 226 | 120 | 35 | 13 |

  Smoothing B permutations therefore costs O(m·N·B) as a GEMM.
- **Direct NW is a different smoother.** Direct NW at the data points (= `locfit(..., ev = dat())`, verified to 9e-14) differs from the default tree-plus-interpolation output by a median RMS of **11-25% of the smoothed field's SD** (correlation 0.97-0.99; `$W/approx_smoother.R`). An implementation that skips locfit's tree is therefore a *method change*, not a numerically equivalent one.

### 5.2 Package cost profile (current code)

- **N = 273 (speKidney A vs B):** 8.3 s per gene, serial. Profile: `geoR::variog` 35%, `locfit` 29% (plus `fitted`/`locfit.matrix` about 10%), `lm`/`model.frame` about 22% (`$W/time_baseline.R`). With `nThreads = 4`: 2.4 s.
- **N ≈ 5000 per permutation and direction:**
  - locfit+fitted: 10 ms (δ = 0.9) to 46 ms (δ = 0.1), about 185 ms over the default 9 δ. δ = 0.01 alone costs 425 ms, so the vignette grid roughly triples smoothing cost.
  - 18 `variog` calls × 11.3 ms ≈ 200 ms.
  - 9 `lm` calls × 0.84 ms.
  - Measured total ≈ 0.46 s, i.e. **about 91 s per gene** (2 × 100). About 60 s per gene at N = 1000-2000.

### 5.3 Exact replica and measured gains

`$W/fastvil.R` (R, factorized operators) and `$W/fastgene.cpp` (RcppArmadillo core) consume exactly the package's random inputs:
- permutations: `sample()` in the parent after `set.seed(seed)`;
- noise: N·K normals per permutation from L'Ecuyer-CMRG after `set.seed(seed+i)`.

Exactness and speed:

| N | Check vs package (same seed) | Package time | R replica (a-e) | C++ core (a-e) |
|---|---|---|---|---|
| 273 (all points) | δ\* identical, p identical, max\|Δsurrogate\| 6e-12 (scale 1400) | 4.0-4.2 s per direction (B = 100) | 0.16-0.50 s | — |
| 1003 (subsample path) | identical, 7e-12 | 61 s per gene (projected) | 4.75 s per gene (13×) | — |
| 2005 | identical, 3.5e-12 | 59 s per gene | 5.0 s (12×) | — |
| 4984 | identical, 1.4e-11 (vs package B = 10); C++ vs R replica 4.6e-13 | **91 s per gene** | 5.9 s (16×) | **1.19 s per gene** (0.80 s C++ + 0.39 s R RNG) ≈ **76×** |

Where the C++ time goes: about 80% is the pair-list variograms (2 per δ per permutation, P ≈ 125k pairs). That is a floor for exact equality, because the noise differs per δ. The other large item is generating N·K·B normals in R to stay bit-compatible.

### 5.4 Classification of candidate speedups

**EXACT** means identical output up to floating point for the same random draws. **APPROXIMATE** means the same null distribution but different numbers. **METHOD-CHANGING** means a different null construction.

| # | Idea | Class | Evidence | Expected speedup | Risk |
|---|---|---|---|---|---|
| (a) | Linear smoother as a precomputed per-δ operator shared by all genes, directions and permutations; use the factorized S_δ = W_δ V_δ (GEMM V·X_perm, sparse W) | **EXACT** | linearity ≤ 1.3e-15; factorization ≤ 1.8e-16; replica identical δ\*, p | Removes locfit (~40% of runtime); smoothing becomes < 5% of the C++ core | W extraction relies on locfit internals (`eva$coef`, `cell$s`, parametric component). Guard with a runtime self-check (S y vs `fitted()`), or port locfit's tree (GPL ≥ 2, compatible with GPL-3) |
| (b) | Variogram as a precomputed pair→bin list (geoR loop order, `hypot`, edges, nugget bin, pairs.min = 2) | **EXACT** (bitwise) | `$W/variog_cpp.R`: identical to `geoR::variog`; 0.18 ms vs 11.3 ms | 63× per variogram; variograms are ~45-50% of package runtime | Must replicate edge semantics (`dist ≤ max.dist`, last edge `umax` excluded) |
| (b') | Smoothed-field variogram as quadratic forms c'M_{δ,k}c with M = W_ids' L_k W_ids (m × m per bin) | **EXACT** (fp) | algebra | ~1.5× on variograms for δ ≥ 0.3 (m ≤ 35); none for small δ | Noise variogram still needs the pair pass |
| (c) | Closed-form OLS instead of `lm()` | **EXACT** (fp) | identical δ\* in all tests | `lm` ≈ 1.5 s per gene in R | none |
| (d) | δ selection needs only the subsample rows; full-length smoothing only for δ\* | **EXACT** | used in replica | Small with factorization; N/1000 for direct-NW implementations | Must still draw the same N·K noise to keep the RNG stream in sync |
| (e) | Reuse subsample, `dist`, `prctile`, bins and operators across directions, genes and iterative rounds; cache stage-1 surrogates in the iterative scheme | **EXACT** (with current seeding) | same ids verified; nesting verified | Removes 3600 internal `dist()` calls per gene; 10% of stage-2 work | If seeding becomes per-gene, keep the subsample shared by design |
| (e') | Within-sample all-pairs: generate each gene's surrogates once and reuse across partners | **EXACT** | surrogates depend only on (gene, coordinates, seed) | **(G − 1)×** fewer surrogate sets | Memory for G × N × B |
| (f1) | Besag-Clifford sequential p-values (stop at h exceedances) | APPROXIMATE (same null; a different, valid p estimator) | kidney stage-2 genes: 914/1000 draws (h = 10), 972 (h = 20) | Large for null-heavy screens (~h/p draws per null gene); small here | Needs the BC formula; resolution only where needed |
| (f2) | Decision-based sequential testing for BH (MMCTest / SIMCTEST) | APPROXIMATE (decisions with bounded resampling risk) | kidney: mean 308 vs 1000 draws for stage-2 genes | ~3× on these data; more for large m | Reports decisions or bounds rather than point p-values; more code |
| (f3) | Parametric tail (Gaussian on Fisher-z, or GPD, Knijnenburg et al. 2009) from about 100-200 surrogates | METHOD-CHANGING | Spearman 0.75-0.83 vs empirical; heavier tails for some genes | 5-50× fewer surrogates for strong genes | Tail mis-specification in exactly the region that matters |
| (g) | Restrict or adapt the δ search (pilot of 20 then ±1 step window; coarse-to-fine; warm starts) | METHOD-CHANGING | δ\* within ±1 step of the pilot median: median 87-90% of permutations per gene, 10th percentile 30-61%; ±2 steps: 93-96% (51-79%) | 2-4× fewer δ evaluations | Changes surrogate selection for 5-15% of permutations in a typical gene, much more for diffuse genes. Not needed once (a)-(e) make the full grid cheap; better to spend compute on a wider grid |
| (h1) | kNN-truncated kernel (u = d/h ≤ 3) | **EXACT** to ~1e-13 | max relative error 2.6e-13; non-zero share 8% / 33% / 58% at δ = 0.01 / 0.05 / 0.1 | Only for direct-NW implementations, or to bound V memory at very large N | negligible |
| (h2) | float32 GEMM/variograms | APPROXIMATE | not measured | ≤ 2× on parts that are not the bottleneck | δ\* flips in near-ties; not recommended |
| (i) | FFT/convolution on regular pixel grids (SEraster square or hex lattice in axial coordinates) | METHOD-CHANGING | — | Smoothing: little gain at N ≤ 10⁴ given (a). Variogram: full-data (no subsample) at O(G log G), a variance *improvement* at similar cost | kNN bandwidth grows near edges and tree interpolation differs, so the smoother differs; hex-grid bookkeeping |
| (j) | Batch genes: one setup, multi-gene GEMMs, shared pairs and operators | **EXACT** | follows from (a), (b), (e) | Amortizes setup; better BLAS efficiency | Memory (N × B per gene in flight) |
| (k) | Threads inside C++ instead of forking per call | **EXACT** with fixed RNG streams | fork overhead 0.29-0.45 s per `bplapply`, 6 calls per gene | Near-linear across genes (~1.2 s per gene per thread) | BLAS oversubscription (pin vecLib/OpenBLAS to 1 thread inside workers); RNG reproducibility needs per-(gene, direction, perm) streams. Moving the RNG into C++ (dqrng/Philox) is APPROXIMATE versus today's numbers but removes the 0.39 s per gene R-RNG cost and lets you draw noise only where needed |
| (m) | Paper-faithful δ\* criterion (regression RSS, no noise variogram); equivalent to the code's criterion in expectation when α, β > 0 | METHOD-CHANGING (small) | algebra: E[γ(aX + bZ)] = a²γ(X) + b² | Removes half the variograms and all per-δ noise generation (about 9× fewer normals). Estimate: ~3× beyond the exact C++ core (1.19 → about 0.35 s per gene at N ≈ 5000) | Changes δ\* where the noise realization currently decides |

### 5.5 Suggested C++ architecture (exact mode first)

1. **Per coordinate set (once per dataset):** subsample ids; pair list (i, j, bin); per δ: W_δ (sparse) and V_δ (dense m × N, or kNN-truncated beyond about 10⁵ points). Build these with locfit as in `factorS.R` and validate with a self-check, or port `ev_atree.c` with its interpolation.
2. **Per gene and direction:**
   - target variogram;
   - per δ: VX = V Xr (BLAS), X_ids = W_ids VX, γ(X_ids), closed-form OLS, H = a X_ids + b E_ids,k, γ(H), RSS;
   - argmin over δ;
   - surrogate = W VX + noise;
   - null correlations by one crossprod of standardized matrices.
3. **Random numbers:** "compat" mode draws in R exactly as now and passes the draws to C++. "fast" mode uses a counter-based C++ RNG keyed by (seed, gene, direction, permutation), so results are reproducible regardless of thread count; numbers differ from today's.
4. **Parallelism:** RcppParallel/TBB across genes, BLAS single-threaded inside workers; drop per-call forking.
5. **Statistics fixes alongside:** (b+1)/(B+1); raw plus adjusted columns; BH on max(pX, pY) (or a single direction); within-sample surrogate reuse; δ-grid diagnostics (boundary-hit rate) and an adaptive default grid (for example k from about 4 neighbours to N/2, log-spaced); optional sequential or decision-based B.

## 6. Related implementations (design ideas and validation references)

- **brainsmash** (Burt et al. 2020, *NeuroImage* 220:117038; [paper](https://www.sciencedirect.com/science/article/pii/S1053811920305243), [docs](https://brainsmash.readthedocs.io/en/latest/approach.html), [code](https://github.com/murraylab/brainsmash), read in `$W/src_pkgs/brainsmash/`). This is the same variogram-matching null.
  - **`Base`:**
    - Precomputes, once, the argsort of the full distance matrix and the per-δ k-nearest-neighbour indices and distances. Smoothing uses exactly k = int(δN) neighbours, excluding self, with kernels `exp` (default), `gaussian`, `invdist` or `uniform`. That is a kNN-truncated sparse kernel.
    - Variogram: smoothed variogram (Gaussian kernel with the 2.68 factor, `nh=25` distances, `pv=25` percentile), as in the paper.
    - δ selection: by **regression residual** (no noise), as in the paper text.
    - Surrogates are generated in vectorized batches (`batch_size` up to 500) with joblib parallelism.
    - `resample=True` rank-maps surrogates onto the empirical values to preserve the marginal distribution.
  - **`Sampled`, for large N:** memory-mapped distance and index matrices with rows pre-sorted and truncated to `knn=1000` neighbours; a fresh random subsample of `ns=500` rows for each surrogate's variogram; δ as a fraction of `knn`; `pv=70`.

  Takeaways for STcompare: precompute kNN structures once; batch surrogates; offer `resample`; for very large N, truncate neighbourhoods and subsample variogram rows.
- **neuromapr** ([CRAN vignette](https://cran.r-project.org/web/packages/neuromapr/vignettes/null-models.html)): an R port of the brain null models (`burt2020` variogram matching, `moran`, spin tests). It is a ready R cross-check for the null generator.
- **GET::GET.localcor** ([docs](https://rdrr.io/cran/GET/man/GET.localcor.html)): another R implementation of Viladomat et al.'s procedure (geoR + locfit, `Delta`, `nsim = 1000`, `N_s = 1000`, `maxk = 300`) with FDR/FWER global envelopes for local correlations (Mrkvička & Myllymäki 2023, *Stat. Comput.* 33:109). It is the closest independent implementation for validation.
- **MERINGUE** (Miller et al. 2021, *Genome Research*; [code](https://github.com/JEFworks-Lab/MERINGUE), read in `$W/src_pkgs/meringue/`). From the same lab as STcompare.
  - RcppArmadillo C++ for Moran's I with analytical moments and normal-approximation p-values (`moranTest_C`, `getSpatialPatterns_C` looping genes in C++), and for the spatial cross-correlation index (`spatialCrossCor_C`, `spatialCrossCorMatrix_C`). Weights are dense N × N, which limits N.
  - Permutation tests use **(b+1)/(n+1)** (`spatialCrossCorTest`, `moranPermutationTest`).
  - A toroidal-shift null (`spatialCrossCorTorTest`) preserves autocorrelation cheaply.

  Takeaway: adopt the same p-value convention; use sparse weights for scale.
- **SpatialPack::modified.ttest** ([refman](https://search.r-project.org/CRAN/refmans/SpatialPack/html/modified.ttest.html)). Implements Clifford, Richardson & Hémon (1989) and Dutilleul (1993): an effective sample size from Moran's I in `nclass = 13` distance classes, in C, O(N²).
  - Measured here: Spearman **0.966** with STcompare's stored p on 4950 null pairs; type-I 0.047 at 0.05. Time: 19 ms per pair (N ≈ 273), 0.06 s (N = 1003), 1.4 s (N = 4984).
  - Recommended as a validation reference and as a fast analytic screen or companion statistic. Its distance-class structure can also be precomputed per coordinate set.
- **Moran's-I / bivariate tools:**
  - spdep, wrapped by Bioconductor's Voyager/SpatialFeatureExperiment ([Voyager](https://bioconductor.org/packages/release/bioc/html/Voyager.html)): sparse neighbour lists, `moran.mc` permutations, Lee's L and other bivariate statistics.
  - squidpy `gr.spatial_autocorr` ([docs](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_autocorr.html)): numba kernels, analytic `pval_norm` plus optional permutations, joblib with numba threads pinned to 1 per job to avoid oversubscription.
  - SpatialDM (Li et al. 2023, *Nat. Commun.* 14:3995; [paper](https://www.nature.com/articles/s41467-023-39608-w)): an **analytical null** (first and second moments) for a bivariate Moran's R statistic. Reported >100× speedups over permutation and scaling to about 10⁶ spots on one CPU.

  Lesson: analytic moments replace permutations wherever the statistic allows it. A Viladomat-type null does not, but the modified t-test does.
- **Moran spectral randomization (MSR)** (Wagner & Dray 2015, *Methods Ecol. Evol.* 6:1169-1178; `adespatial::msr` [refman](https://search.r-project.org/CRAN/refmans/adespatial/html/msr.html); BrainSpace; neuromaps). It eigendecomposes the spatial weight matrix *once per coordinate set*: O(N³), or O(N·r²) with r eigenvectors. Each surrogate then randomizes coefficients in the Moran eigenvector basis at O(N·r) while preserving Moran's I spectrum. Same pattern as (a): heavy per-geometry setup, cheap per-surrogate work, shared across genes. It is a candidate fast alternative null (METHOD-CHANGING). Markello & Misic (2021, *NeuroImage* 236:118052; [paper](https://www.sciencedirect.com/science/article/pii/S1053811921003293)) benchmark these null families against each other.

## 7. Recommended test datasets and validation harness

1. **Exact-equivalence harness** (any new implementation in compat mode versus the current package at a fixed seed): δ\* sequences identical, p identical, surrogates within 1e-10.
   - `speKidney` A vs B (N = 273, all points): `$W/kidney_AB.rds`.
   - Synthetic kidney-shaped hexagonal grids with independent GRF genes (exponential covariance κ = 0.1, nugget sd √0.3, mean 10, as in `simRanPatternRasts`) at N = 1003, 2005, 4984; these exercise the 1000-point subsample: `$W/sim_N1000.rds`, `$W/sim_N2000.rds`, `$W/sim_N5000.rds` (4 X genes + 4 independent Y genes each; generator `$W/sim_data.R`).
   - Reference scripts: `$W/test_exact2.R`, `$W/test_scale.R`, `$W/test_cpp.R`.
2. **Null calibration:** `simRanPatternRasts` (9900 ordered pairs with stored STcompare p as a regression baseline) plus the modified t-test (`$W/modttest_null.R`, output `$W/simRanPattern_with_mtt.rds`); marginal-robustness scenarios in `$W/null_marginals.R` (outputs `$W/null_*.rds`).
3. **Real-data decision checks:** the shipped `inst/extdata/kidneyCorrelation*.RData` and `merfishCorrelation*.RData` (null correlations and δ\* for 1046 and 483 genes) for p-value, BH and δ-grid diagnostics (`$W/analyze_stored.R`, `$W/mmc_potential.R`).

## 8. Files produced (all under `$W`)

- `fastvil.R`: exact vectorized replica (variogram precompute via Rcpp `pairBins`, `variogCols`, `fastViladomat`, `fastViladomatF`, `prepareCoords`).
- `factorS.R`: exact factorization of locfit's smoother, S_δ = W_δ V_δ.
- `fastgene.cpp`: exact RcppArmadillo core for one gene direction.
- `test_exact.R`, `test_exact2.R`, `test_factor.R`, `test_scale.R`, `test_cpp.R`, `nesting.R`: exactness and timing tests.
- `time_baseline.R`, `vertex_counts.R`, `variog_cpp.R`, `bp_overhead.R`, `approx_smoother.R`, `beta_signs.R`, `check_warnings.R`: profiling and approximation checks.
- `inspect_results.R`, `inspect_sim.R`, `analyze_stored.R`, `mmc_potential.R`, `modttest_null.R`, `mtt_timing.R`, `null_marginals.R`, `similarity_demo.R`: statistical analyses.
- `sim_data.R` → `sim_N1000.rds`, `sim_N2000.rds`, `sim_N5000.rds`; `setup_data.R` → `kidney_AB.rds`; `Slist_kidney.rds` (unit-vector operators for N = 273).
- `viladomat_code.R` (original prototype code); `src_pkgs/` (locfit, geoR, brainsmash, MERINGUE sources read for this report); `Rlib/` (SpatialPack installed for validation).

## References

- Viladomat J, Mazumder R, McInturff A, McCauley DJ, Hastie T (2014). Biometrics 70(2):409-418. doi:10.1111/biom.12139. https://pmc.ncbi.nlm.nih.gov/articles/PMC4108159/ ; code: https://hastie.su.domains/Papers/biodiversity/viladomat_code.R
- Burt JB, Helmer M, Shinn M, Anticevic A, Murray JD (2020). Generative modeling of brain maps with spatial autocorrelation. NeuroImage 220:117038. https://www.sciencedirect.com/science/article/pii/S1053811920305243 ; https://github.com/murraylab/brainsmash
- Miller BF, Bambah-Mukku D, Dulac C, Zhuang X, Fan J (2021). MERINGUE. Genome Research. https://github.com/JEFworks-Lab/MERINGUE
- Clifford P, Richardson S, Hémon D (1989). Biometrics 45:123-134. Dutilleul P (1993). Biometrics 49:305-314. SpatialPack: https://search.r-project.org/CRAN/refmans/SpatialPack/html/modified.ttest.html
- Wagner HH, Dray S (2015). Methods Ecol Evol 6:1169-1178 (MSR). https://search.r-project.org/CRAN/refmans/adespatial/html/msr.html
- Markello RD, Misic B (2021). Comparing spatial null models for brain maps. NeuroImage 236:118052. https://www.sciencedirect.com/science/article/pii/S1053811921003293
- Li Z, Wang T, Liu P, Huang Y (2023). SpatialDM. Nat Commun 14:3995. https://www.nature.com/articles/s41467-023-39608-w
- Mrkvička T, Myllymäki M (2023). Stat Comput 33:109; GET::GET.localcor https://rdrr.io/cran/GET/man/GET.localcor.html
- Besag J, Clifford P (1991). Sequential Monte Carlo p-values. Biometrika 78(2):301-304.
- Gandy A (2009). JASA 104:1504-1511 (SIMCTEST). Gandy A, Hahn G (2014). Scand J Stat 41:1083-1101 (MMCTest).
- Sandve GK, Ferkingstad E, Nygård S (2011). Sequential Monte Carlo multiple testing. Bioinformatics 27(23):3235-3241.
- Phipson B, Smyth GK (2010). Permutation p-values should never be zero. Stat Appl Genet Mol Biol 9:39. North BV, Curtis D, Sham PC (2002). AJHG 71:439-441. Davison AC, Hinkley DV (1997). Bootstrap Methods and their Application.
- Knijnenburg TA et al. (2009). Fewer permutations, more accurate P-values. Bioinformatics 25(12):i161-i168.
- Voyager: https://bioconductor.org/packages/release/bioc/html/Voyager.html ; squidpy: https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_autocorr.html ; neuromapr: https://cran.r-project.org/web/packages/neuromapr/vignettes/null-models.html
