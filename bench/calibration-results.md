# Calibration, power and cost of the surrogate modes of `compareSpatial()`

`compareSpatial()` has two ways to build the surrogates of its correlation test (`surrogate`): `"gaussian"`, the
Viladomat surrogates as the smoothing and the added noise leave them (the surrogates of the legacy functions), and
`"remap"`, the same surrogates with their values replaced by the gene's own values in the surrogate's rank order
(the amplitude adjustment of AAFT surrogates), so that each surrogate has exactly the gene's distribution of
values. `dev/cpp-backend-plan.md`, section 8, records why: for genes detected in few pixels the null correlations
have heavy tails that Gaussian-valued surrogates do not reproduce, and their p-values are far too small.

This file holds the results of `bench/calibrate-surrogates.R`, the study that `dev/surrogate-remap-spec.md`,
section 4, asks for. The section between the two HTML comments is written by the script; the discussion above it
was written from that run. Every run used the default settings of `compareSpatial()` (adaptive p-values with
`exceedances = 10` and `nPermutations = 10000`, the extended delta grid) unless stated. At the time of the run
`"gaussian"` with the sqrt(N) detection filter was the default; on the recommendation below, `"remap"` with no
filter is the default now (`minDetected = NULL` means 0 pixels with `"remap"` and sqrt(N) pixels with
`"gaussian"`), and the script passes the mode and the filter explicitly so that the run reproduces.

## Summary

1. **Dense null** (all 4950 pairs of the 100 independent simulated fields of `data(simRanPatternRasts)`). Both
   modes are conservative and nearly the same: P(p <= 0.05) is 0.030 (gaussian) and 0.034 (remap), P(p <= 0.001)
   is 0.0006 and 0.0008, no BH discovery in either. Remapping a Gaussian field onto its own values is close to a
   linear transformation, so the two modes should agree here, and they do.

2. **Sparse null without spatial structure** (the acceptance review's setting: 1800 independent genes detected in
   1, 2, 3, 5, 10 or 20 of the 311 AKI pixels, `minDetected = 0`). Gaussian: P(p <= 0.001) reaches 0.033 for
   genes detected in 10 pixels and 0.013 for 20 pixels, P(p <= 0.01) 0.043 and 0.027, and BH finds 9 false
   discoveries among the 1800 null genes. Remap: P(p <= 0.001) is at most 0.0033 (one gene of 300) and
   P(p <= 0.01) at most 0.013 in every class, P(p <= 0.1) is 0.03 to 0.10, and BH finds no false discovery. One
   remapped p-value sits at the floor 1 / 10001, which is what uniform p-values give with probability 0.16 for 1800
   genes. The remapped run also used 3.7 times fewer permutations (71,396 against 260,835), because fewer null
   genes get small p-values and run long.

3. **Zero-inflated spatial fields** (Poisson counts of Gaussian random fields on the AKI and brain coordinates,
   1% to 100% of the pixels detected, `minDetected = 0`). Gaussian is anti-conservative wherever fewer than about
   10% of the pixels are detected: on the AKI grid P(p <= 0.001) is 0.028 at 1% detected (about 2 pixels) and
   0.013 at 3% (8 pixels), with P(p <= 0.01) of 0.044 and 0.038; on the brain grid 0.010 and 0.015 at 1% (19
   pixels). From 10% detected upwards it is calibrated or conservative. Remap keeps P(p <= alpha) at or below
   alpha at every fraction, on both grids and at every alpha: its largest rates are 0.11 (standard error 0.02) at
   alpha = 0.1 and 0.06 (0.017) at alpha = 0.05 for the dense Gaussian fields on the brain grid, within one
   standard error of nominal; P(p <= 0.01) never exceeds 0.01 and no remapped p-value is below 0.001 in any
   class. (37 of the 400 AKI genes at 1% detected are constant in one sample and skipped in both modes.)

4. **Power** (fields mixed to a population correlation of 0.2, 0.4 and 0.6; 300 pairs each, the same pairs and
   seeds for both modes). The two modes have the same power within the standard errors; remap is slightly
   higher at alpha = 0.01 (0.487 against 0.463 at rho = 0.4; 0.937 against 0.920 at rho = 0.6).

5. **The published inputs.** AKI (1046 genes, CPM): 738 genes significant with gaussian and 750 with remap, 737
   in both; of the 731 published genes, 723 and 728 are recovered (8 and 3 are missed). Brain (325 genes,
   lognorm): 132 against 171, with the 132 gaussian genes all among the remapped ones; of the 128 published
   genes, 122 and 124 are recovered. The remapped p-values of the brain genes are systematically smaller (the
   median ratio of the gaussian to the remapped p-value over the tested genes is 1.5, the upper quartile 3): the
   39 genes significant with remap only all have a positive correlation of 0.06 to 0.35, and they are not
   especially sparse. This agrees with the null results on the brain grid, where the gaussian mode is
   conservative for skewed, zero-inflated fields with a spatial pattern (P(p <= 0.05) of 0.015 to 0.02 at 3% to
   30% detected, against 0.015 to 0.04 for remap and 0.05 nominal), so the additional genes are consistent with
   power the gaussian surrogates lose on log-normalized zero-inflated data; a null experiment on these genes
   themselves is not possible, because unrelated genes of a tissue share its anatomy.

6. **The whole AKI raster** (32,285 genes, `minDetected = 0`, `nPermutations = 1000`). 14,383 genes are never
   detected in one sample and are skipped. Among the 4402 tested genes detected in at most 17 pixels of a sample
   (which the default `minDetected` excludes), gaussian calls 102 significant, with P(p <= 0.001) of 0.009 to
   0.018 in these classes, and one of them is a published spatially variable gene; remap calls 1, with
   P(p <= 0.001) at most 0.0006. For genes detected in 18 to 100 pixels gaussian calls 102 and remap 51; above
   100 pixels, 944 and 840 (507 and 484 of them published genes). At this resolution (p-values down to 1 / 1001,
   BH across about 17,900 tests) the modes differ for the genes near the BH threshold; on the 1046 published genes
   at the default 10,000 permutations remap found more genes than gaussian (item 5).

7. **Cost.** The remapping adds one sort of N values per permutation and direction: 160 against 171 microseconds
   times threads per permutation and direction on the AKI input, 1120 against 1210 on the brain input (7 to 8%).
   The wall times differ more where remap gives smaller p-values and more genes run to `nPermutations`: see the
   script's table below (AKI about 10% longer, brain about 40% longer, the whole AKI raster about 10%). On null
   genes remap is faster, because fewer of them get small p-values.

## Recommendation

Adopted on 2026-10-05: `surrogate = "remap"` is the default of `compareSpatial()`, and `minDetected = NULL` means 0
pixels with `"remap"` and sqrt(N) pixels with `"gaussian"`.

- **Make `surrogate = "remap"` the default.** It is calibrated in every null setting of this study, including the
  genes detected in one or two pixels where the gaussian surrogates give p-values ten to thirty times too small
  in the far tail and BH false discoveries, it has the same power on correlated Gaussian fields, it recovers the
  published genes at least as well (728 of 731 AKI genes, 124 of 128 brain genes), and it costs about 7% per
  permutation. For a gene without spatial structure it is the exact permutation test. `"gaussian"` stays
  available for comparability with the legacy functions and the published analyses.
- **`minDetected`.** With remapped surrogates the p-values were calibrated down to genes detected in a single
  pixel, so the filter is no longer needed for validity: `minDetected = 0` can be the default when
  `surrogate = "remap"` is. Such genes are cheap (they stop after about ten permutations; the whole AKI raster
  with 17,902 tested genes took about 70 s on 16 threads at 1000 permutations), and a gene detected in k pixels
  cannot get a p-value below about k / N anyway (the exact test is discrete). With `surrogate = "gaussian"` the
  filter is still needed: the sqrt(N) default removes the genes detected in fewer than 18 of 311 or 47 of 2170
  pixels, where the gaussian p-values are worst, and even above it the far tail is too heavy (P(p <= 0.001) of
  0.013 for 20 of 311 pixels). A clean way to express this is `minDetected = NULL` meaning 0 for `"remap"` and
  sqrt(N) pixels for `"gaussian"`, with the help page saying so.
- **Two observations worth knowing.** (a) Ties: for a gene detected in few pixels many rearrangements give
  exactly the observed correlation, and `cor()` rounds each of them differently, so without a tolerance about half
  of the ties would be lost and the remapped p-value of such a gene could be several times too small; remapped
  tasks therefore count a null within a relative 1e-9 of |r| as an exceedance (`kRemapTieRel` in
  `src/stc_engine.h`), as the exact test does. (b) For a gene without spatial structure correlated with a smooth
  gene, the Viladomat surrogate keeps some smoothness (the slope of the flat target variogram on the candidate
  variogram is noise, and the method uses its absolute value), so in that case the null correlations of both modes
  are about 30% wider than the plain permutation null: conservative, and the remapping cannot change it because
  it keeps the ranks.

<!-- calibrate-surrogates:begin -->
Run of `bench/calibrate-surrogates.R` on 2026-10-05 09:23 with 16 threads.
Machine: Apple M1 Ultra (arm64, 20 logical cores), R version 4.5.2 (2025-10-31). Package: STcompare 0.1.0.9000.

Wall time of each part and the load averages (1, 5, 15 min) after it:

| part | wall | load_after |
|---|---|---|
| dense | 56.5 s | 12.49 14.01 12.97 |
| sparse (a) gaussian | 6.4 s | 13.41 14.16 13.03 |
| sparse (a) remap | 2.7 s | 13.41 14.16 13.03 |
| sparse (b) aki gaussian | 7.1 s | 14.20 14.30 13.10 |
| sparse (b) aki remap | 3.1 s | 14.20 14.30 13.10 |
| sparse (b) brain gaussian | 15.2 s | 15.47 14.62 13.25 |
| sparse (b) brain remap | 11.7 s | 15.47 14.62 13.25 |
| power rho = 0.2 | 11.6 s | 14.66 14.48 13.22 |
| power rho = 0.4 | 72.9 s | 11.42 13.87 13.12 |
| power rho = 0.6 | 104.1 s | 15.36 14.77 13.56 |
| real aki gaussian | 105.3 s | 18.13 16.69 14.68 |
| real aki remap | 115.5 s | 18.13 16.69 14.68 |
| real brain gaussian | 2.1 min | 19.28 17.74 15.72 |
| real brain remap | 3.0 min | 19.28 17.74 15.72 |
| real AKI all genes gaussian | 66.0 s | 17.87 17.61 15.97 |
| real AKI all genes remap | 72.3 s | 17.87 17.61 15.97 |

## Null, dense fields

All 4950 pairs of the 100 independent simulated fields of `data(simRanPatternRasts)` (one gene per pair, 259 to 283 shared pixels), default settings, both modes on the same pairs and seeds. P(p <= alpha) with its binomial standard error; the pairs share fields, so the standard errors understate the uncertainty somewhat. Wall time 56.5 s on 16 processes.

| mode | pairs | failed | median L | BH < 0.05 | P(p <= 0.1) | P(p <= 0.05) | P(p <= 0.01) | P(p <= 0.001) |
|---|---|---|---|---|---|---|---|---|
| gaussian | 4950 | 0 | 17 | 0 | 0.0699 (0.0036) | 0.0303 (0.0024) | 0.0057 (0.0011) | 0.0006 (0.0003) |
| remap | 4950 | 0 | 17 | 0 | 0.0754 (0.0038) | 0.0343 (0.0026) | 0.0065 (0.0011) | 0.0008 (0.0004) |

## Null, sparse genes without spatial structure (the acceptance review's setting)

1800 independent genes on the AKI coordinates (311 pixels): for each number k of detected pixels, 300 genes in each sample with values round(rexp(k) * 20) + 1 at k random pixels (set.seed(123), as in the review). Default settings, minDetected = 0, seed 4, 16 threads. The default minDetected would skip every gene detected in fewer than 18 pixels.

| detected pixels | mode | genes | failed | median L | P(p <= 0.1) | P(p <= 0.05) | P(p <= 0.01) | P(p <= 0.001) |
|---|---|---|---|---|---|---|---|---|
| 1 | gaussian | 300 | 0 | 10.0 | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| 1 | remap | 300 | 0 | 10.0 | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| 2 | gaussian | 300 | 0 | 10.0 | 0.0067 (0.0047) | 0.0067 (0.0047) | 0.0067 (0.0047) | 0.0033 (0.0033) |
| 2 | remap | 300 | 0 | 10.0 | 0.0067 (0.0047) | 0.0067 (0.0047) | 0.0033 (0.0033) | 0.0000 (0.0000) |
| 3 | gaussian | 300 | 0 | 10.5 | 0.0233 (0.0087) | 0.0200 (0.0081) | 0.0133 (0.0066) | 0.0100 (0.0057) |
| 3 | remap | 300 | 0 | 10.0 | 0.0300 (0.0098) | 0.0267 (0.0093) | 0.0067 (0.0047) | 0.0000 (0.0000) |
| 5 | gaussian | 300 | 0 | 11.0 | 0.0200 (0.0081) | 0.0200 (0.0081) | 0.0133 (0.0066) | 0.0000 (0.0000) |
| 5 | remap | 300 | 0 | 10.0 | 0.0367 (0.0109) | 0.0333 (0.0104) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| 10 | gaussian | 300 | 0 | 13.0 | 0.0600 (0.0137) | 0.0567 (0.0133) | 0.0433 (0.0118) | 0.0333 (0.0104) |
| 10 | remap | 300 | 0 | 11.0 | 0.0867 (0.0162) | 0.0467 (0.0122) | 0.0133 (0.0066) | 0.0033 (0.0033) |
| 20 | gaussian | 300 | 0 | 15.0 | 0.0633 (0.0141) | 0.0533 (0.0130) | 0.0267 (0.0093) | 0.0133 (0.0066) |
| 20 | remap | 300 | 0 | 16.0 | 0.0967 (0.0171) | 0.0467 (0.0122) | 0.0067 (0.0047) | 0.0033 (0.0033) |

- gaussian: 9 BH false discoveries (padj < 0.05) among 1800 genes, 2 p-values at the floor 1 / 10001, 260,835 permutations in 6.4 s
- remap: 0 BH false discoveries (padj < 0.05) among 1800 genes, 1 p-values at the floor 1 / 10001, 71,396 permutations in 2.7 s

## Null, zero-inflated spatial fields by detection fraction

Independent pairs of Gaussian random fields (exponential covariance, range 10% of the coordinate span) turned into Poisson counts with the mean chosen so that the given share of pixels is detected (nonzero); 100% is the Gaussian field itself. 400 genes per fraction on the AKI coordinates and 200 on the brain coordinates, both modes on the same genes, default settings, minDetected = 0, seed 5, 16 threads. padj is adjusted within each run (all fractions together).

| coordinates | detection fraction | detected pixels (min of x, y), median | mode | genes | skipped (constant) | failed | median L | BH < 0.05 | P(p <= 0.1) | P(p <= 0.05) | P(p <= 0.01) | P(p <= 0.001) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aki (N = 311) | 0.01 | 2 | gaussian | 400 | 37 | 0 | 11.0 | 0 | 0.0523 (0.0117) | 0.0523 (0.0117) | 0.0441 (0.0108) | 0.0275 (0.0086) |
| aki (N = 311) | 0.01 | 2 | remap | 400 | 37 | 0 | 10.0 | 0 | 0.0386 (0.0101) | 0.0193 (0.0072) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| aki (N = 311) | 0.03 | 8 | gaussian | 400 | 0 | 0 | 15.0 | 0 | 0.0875 (0.0141) | 0.0675 (0.0125) | 0.0375 (0.0095) | 0.0125 (0.0056) |
| aki (N = 311) | 0.03 | 8 | remap | 400 | 0 | 0 | 10.0 | 0 | 0.0475 (0.0106) | 0.0350 (0.0092) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| aki (N = 311) | 0.10 | 28 | gaussian | 400 | 0 | 0 | 16.0 | 0 | 0.0575 (0.0116) | 0.0300 (0.0085) | 0.0100 (0.0050) | 0.0000 (0.0000) |
| aki (N = 311) | 0.10 | 28 | remap | 400 | 0 | 0 | 15.0 | 0 | 0.0650 (0.0123) | 0.0200 (0.0070) | 0.0050 (0.0035) | 0.0000 (0.0000) |
| aki (N = 311) | 0.30 | 90 | gaussian | 400 | 0 | 0 | 16.0 | 0 | 0.0575 (0.0116) | 0.0325 (0.0089) | 0.0075 (0.0043) | 0.0025 (0.0025) |
| aki (N = 311) | 0.30 | 90 | remap | 400 | 0 | 0 | 17.0 | 0 | 0.0850 (0.0139) | 0.0400 (0.0098) | 0.0075 (0.0043) | 0.0000 (0.0000) |
| aki (N = 311) | 1.00 | 131 | gaussian | 400 | 0 | 0 | 17.0 | 0 | 0.0700 (0.0128) | 0.0300 (0.0085) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| aki (N = 311) | 1.00 | 131 | remap | 400 | 0 | 0 | 17.0 | 0 | 0.0700 (0.0128) | 0.0300 (0.0085) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.01 | 19 | gaussian | 200 | 0 | 0 | 14.0 | 0 | 0.0850 (0.0197) | 0.0450 (0.0147) | 0.0150 (0.0086) | 0.0100 (0.0070) |
| brain (N = 2170) | 0.01 | 19 | remap | 200 | 0 | 0 | 10.0 | 0 | 0.0250 (0.0110) | 0.0150 (0.0086) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.03 | 61 | gaussian | 200 | 0 | 0 | 17.0 | 0 | 0.0350 (0.0130) | 0.0150 (0.0086) | 0.0050 (0.0050) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.03 | 61 | remap | 200 | 0 | 0 | 15.0 | 0 | 0.0350 (0.0130) | 0.0150 (0.0086) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.10 | 210 | gaussian | 200 | 0 | 0 | 15.0 | 0 | 0.0350 (0.0130) | 0.0150 (0.0086) | 0.0100 (0.0070) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.10 | 210 | remap | 200 | 0 | 0 | 17.0 | 0 | 0.0600 (0.0168) | 0.0200 (0.0099) | 0.0050 (0.0050) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.30 | 639 | gaussian | 200 | 0 | 0 | 16.0 | 0 | 0.0500 (0.0154) | 0.0200 (0.0099) | 0.0050 (0.0050) | 0.0000 (0.0000) |
| brain (N = 2170) | 0.30 | 639 | remap | 200 | 0 | 0 | 19.0 | 0 | 0.0700 (0.0180) | 0.0400 (0.0139) | 0.0050 (0.0050) | 0.0000 (0.0000) |
| brain (N = 2170) | 1.00 | 969 | gaussian | 200 | 0 | 0 | 17.0 | 0 | 0.1000 (0.0212) | 0.0600 (0.0168) | 0.0050 (0.0050) | 0.0000 (0.0000) |
| brain (N = 2170) | 1.00 | 969 | remap | 200 | 0 | 0 | 16.5 | 0 | 0.1100 (0.0221) | 0.0600 (0.0168) | 0.0100 (0.0070) | 0.0000 (0.0000) |

## Power, mixed simulated fields

300 random pairs (i, j) of the simulated fields per rho, Y = stc_mix(f_i, f_j, rho) (the recipe of the calibration fixture: a field with the same covariance model and population correlation rho with X = f_i), default settings, the same pairs and seeds for both modes.

| rho | mode | pairs | mean r | failed | median L | P(p <= 0.05) | P(p <= 0.01) | P(p <= 0.001) |
|---|---|---|---|---|---|---|---|---|
| 0.2 | gaussian | 300 | 0.208 | 0 | 38.0 | 0.2000 (0.0231) | 0.0667 (0.0144) | 0.0067 (0.0047) |
| 0.2 | remap | 300 | 0.208 | 0 | 40.0 | 0.2067 (0.0234) | 0.0767 (0.0154) | 0.0100 (0.0057) |
| 0.4 | gaussian | 300 | 0.402 | 0 | 704.5 | 0.6767 (0.0270) | 0.4633 (0.0288) | 0.2867 (0.0261) |
| 0.4 | remap | 300 | 0.402 | 0 | 802.0 | 0.6733 (0.0271) | 0.4867 (0.0289) | 0.3100 (0.0267) |
| 0.6 | gaussian | 300 | 0.587 | 0 | 10000.0 | 0.9667 (0.0104) | 0.9200 (0.0157) | 0.8100 (0.0226) |
| 0.6 | remap | 300 | 0.587 | 0 | 10000.0 | 0.9700 (0.0098) | 0.9367 (0.0141) | 0.8167 (0.0223) |

## The published inputs: significant genes and cost

compareSpatial() with the default settings on the inputs of the published analyses (seed 0, 16 threads). Significant: padj < 0.05. Published: both BH-adjusted direction p-values below 0.05 in bench/published (100 then 1000 permutations). The last column is the wall time per permutation and direction, times the threads, over the permutations kept.

| dataset | mode | genes | tested | skipped | failed | padj < 0.05 | in both modes | also published | published only | published total | median L | permutations | wall | us per task-permutation x threads |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aki | gaussian | 1046 | 1039 | 7 | 0 | 738 | 737 | 723 | 8 | 731 | 3764 | 5,272,989 | 105.3 s | 160 |
| aki | remap | 1046 | 1039 | 7 | 0 | 750 | 737 | 728 | 3 | 731 | 4858 | 5,417,273 | 115.5 s | 171 |
| brain | gaussian | 325 | 297 | 28 | 0 | 132 | 132 | 122 | 6 | 128 | 275 | 883,451 | 2.1 min | 1120 |
| brain | remap | 325 | 297 | 28 | 0 | 171 | 132 | 124 | 4 | 128 | 635 | 1,189,212 | 3.0 min | 1210 |

## The whole AKI raster with minDetected = 0

All 32285 genes of the AKI raster (CPM, 311 shared pixels), minDetected = 0, nPermutations = 1000, exceedances = 10, seed 0, 16 threads: wall 66.0 s (gaussian) and 72.3 s (remap). Genes by the number of pixels where both samples detect them; the default minDetected tests only the genes detected in at least 18 pixels of each sample. The 1046 published genes are spatially variable genes.

| detected pixels (min of x, y) | mode | genes | tested | padj < 0.05 | of them published genes | P(p <= 0.001) | median L |
|---|---|---|---|---|---|---|---|
| 0 | gaussian | 14383 | 0 | 0 | 0 |  | NA |
| 0 | remap | 14383 | 0 | 0 | 0 |  | NA |
| 1-2 | gaussian | 1731 | 1731 | 33 | 0 | 0.0087 (0.0022) | 10 |
| 1-2 | remap | 1731 | 1731 | 0 | 0 | 0.0000 (0.0000) | 10 |
| 3-5 | gaussian | 1066 | 1066 | 30 | 0 | 0.0178 (0.0041) | 12 |
| 3-5 | remap | 1066 | 1066 | 0 | 0 | 0.0000 (0.0000) | 10 |
| 6-17 | gaussian | 1605 | 1605 | 39 | 1 | 0.0100 (0.0025) | 15 |
| 6-17 | remap | 1605 | 1605 | 1 | 0 | 0.0006 (0.0006) | 12 |
| 18-50 | gaussian | 1791 | 1791 | 42 | 0 | 0.0128 (0.0027) | 17 |
| 18-50 | remap | 1791 | 1791 | 15 | 0 | 0.0039 (0.0015) | 19 |
| 51-100 | gaussian | 1541 | 1541 | 60 | 10 | 0.0234 (0.0038) | 21 |
| 51-100 | remap | 1541 | 1541 | 36 | 9 | 0.0130 (0.0029) | 22 |
| > 100 | gaussian | 10168 | 10168 | 944 | 507 | 0.0647 (0.0024) | 32 |
| > 100 | remap | 10168 | 10168 | 840 | 484 | 0.0623 (0.0024) | 33 |

<!-- calibrate-surrogates:end -->
