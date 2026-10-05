# What `geoR::variog` computes inside STcompare, and how to port it exactly to C++

Scope: the three variogram calls in `R/spatialCorrelation.R` (target at lines 272-274, candidates at
114-117 and 125-128), the `lm()` step at line 120, and the RSS at line 131. geoR version 1.9-6. The R code
of the installed CRAN binary (built "R 4.5.0; aarch64-apple-darwin20") is identical to the CRAN source
tarball (checked with `identical(deparse(...))`). Environment: R 4.5.2 on macOS arm64 (M1 Ultra), R linked
to Accelerate/vecLib BLAS, Apple clang 17.

All paths below are relative to the work directory
`/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/variog/`
(WD). geoR source: `WD/geoR/` (from `WD/geoR_1.9-6.tar.gz`). Nothing in the repository was modified.

---------------------------------------------------------------------------------------------------

## 0. Summary

* With STcompare's arguments, `geoR::variog` returns **at most 13 bins**. The edges are `k*(umax/13)`, where
  **`umax` is the largest inter-point distance strictly below `max.dist`**. It is not `max.dist`. A pair goes
  into a bin only if **`hypot(dx,dy) < umax`**, using `[lo, hi)` bins. Pairs exactly at `umax`, and pairs
  between `umax` and `max.dist`, are **dropped**. A pair that lies exactly on an inner edge goes to the upper
  bin. The value in each bin is `v = sum((z_i - z_j)^2 / 2) / n`, summed sequentially in double. Bins with
  `n < 2` are **left out** of the result rather than set to NA, so `length(v)` can be less than 13. On the
  kidney data it is 10.
* Bins, `u` and `n` depend **only on the coordinates and `max.dist`**. They were identical across 69
  different data vectors per coordinate set. `trend = "cte"` does **not** centre the data, and centring
  would change about 1% of `v` values in the last bit.
* geoR combines **three different distance formulas**:
  1. STcompare's `prctile` uses R `dist(cbind(lat, long))`.
  2. geoR's `umax` uses R `dist(cbind(long, lat))`. The columns are swapped relative to (1).
  3. The C binning uses `hypot`.

  On this machine, R's `dist()` is compiled with **FMA contraction**. It equals `sqrt(fma(d2, d2, d1*d1))`
  exactly (0 mismatches over 1.04 M pairs), so formulas (1) and (2) differ for 9-11% of pairs, and
  `hypot` differs from `dist()` for 5-6% of pairs, by 1 ulp. All three have to be reproduced.
* On rasterized (hexagonal-lattice) coordinates, **`max.dist` always lands on a cluster of distances that
  are tied mathematically**, and `umax` is a member of that same cluster 1 ulp lower. Which members of the
  cluster make it into the last bin therefore depends on ulp-level noise in the coordinates. For the
  kidney data: of 481 pairs at a distance of about 1.03923, 328 are kept, 134 are dropped and 19 fall
  outside `max.dist`. Using `dist()` instead of `hypot` for binning would change 2 pairs on the kidney data.
* I implemented a vectorized pure-R port (`WD/variog_fast.R`) and an Rcpp port (`WD/variog_fast.cpp`). Both
  precompute the pair-to-bin index once. They are **bit-identical to geoR** (`identical()` on `u`, `n`, `v`
  and `bins.lim`):
  * 7 real coordinate setups × 68-69 data vectors each;
  * 5 tie-heavy stress scenarios (exact edge ties, pairs at `max.dist`, co-located points);
  * 300 randomized configurations (hex, square, uniform and clustered layouts; N from 30 to 1000; spacing
    from 0.055 to 50; offsets up to ±1000; percentiles from 0.05 to 0.9).
* Speed per variogram:

  | N | geoR | pure R (batched) | Rcpp, one column | Rcpp, batched (transposed) |
  |---|---|---|---|---|
  | 273-300 | 0.76-1.19 ms | 86-112 µs | 11-13 µs | 2.3-2.7 µs |
  | 1000 | 9.8-10.4 ms | 1.2 ms | 122-126 µs | 31-33 µs |

  That is a **70-93x** speedup for single columns and **310-430x** batched. One gene with the defaults
  needs 3602 variograms: 2.7 s (N=273) to 38 s (N=1000) with geoR, versus 0.04-0.46 s with Rcpp.
  Precomputation, done once: 0.3-0.6 ms (N≈300) and 6.6 ms (N=1000) in C++.
* `lm()` is an ordinary least-squares fit with intercept and slope. The closed form agrees to within 7e-16
  relative on the slope and 3e-13 relative on the intercept. In 10 of 10 permutations it picks the same
  `delta*`. **`lm()` cannot be matched bit-for-bit reliably.** With Accelerate BLAS, the LINPACK QR
  (`dqrls`) gives different last bits depending on the 16-byte alignment of its buffers. If the candidate
  variogram is constant or nearly constant (`||x - mean(x)|| < 1e-7 ||x||`), the slope is NA. geoR then
  stops with an error on the NA data, and STcompare reports NA p-values for that gene. An all-zero gene
  hits exactly this path.

---------------------------------------------------------------------------------------------------

## 1. How STcompare calls variog

`R/spatialCorrelation.R`:

* `spatialCorrelation()` builds `dataForward <- data.frame(X, Y, x = pos[,1], y = pos[,2])` (512-515).
  `viladomatCorrelation()` then sets `lat <- data[,3]` and `long <- data[,4]` (247-248). So **lat = pos[,1]
  (x)** and **long = pos[,2] (y)**.
* `set.seed(seed)` (242). If N > 1000, `ids <- sample(N, 1000)` (253-255); otherwise `ids <- 1:N`.
  **Every gene and both directions use the same seed** (`seed = seed` at 542, 548 and 831), and N is fixed
  for a dataset pair. So `ids`, `prctile`, the coordinates and therefore the whole bin structure are
  **identical for every variog call in a `spatialCorrelationGeneExp` run**. One precomputation per dataset
  pair is enough.
* `prctile <- quantile(dist(cbind(lat_s, long_s)), probs = maxDistPrctile)` (268-269). The column order is
  (lat, long).
* `target_variog <- geoR::variog(data = X_s, coords = cbind(long_s, lat_s), max.dist = prctile, ...)`
  (272-274). The column order is (long, lat).
* For each permutation and each delta, `matchingVariograms()` runs:
  * `variog(X.delta[ids], cbind(long[ids], lat[ids]))` (114-117);
  * `lm(target$v ~ 1 + cand$v)` (120);
  * `hat = X.delta*sqrt(|b1|) + rnorm(N)*sqrt(|b0|)` (124). `rnorm` draws N_full values, not just the
    subsample.
  * `variog(hat[ids], ...)` (125-128);
  * `RSS = sum((v_hat - v_target)^2)` (131);
  * `which.min` (134).

  With the defaults, that is 2 × (1 + 100 × 9 × 2) = **3602 variog calls per gene**.

## 2. geoR 1.9-6 source walk-through (with defaults)

`WD/geoR/R/variogram.R` (`variog` is lines 1-277, `.define.bins` is lines 1019-1047) and
`WD/geoR/src/geoR.c` (`binit` is lines 370-422):

| step | source | behaviour |
|---|---|---|
| coercion | variogram.R:57-58 | `coords <- as.matrix(coords)`, `data <- as.matrix(data)`; no rescaling |
| `data.var` | :60 | `apply(data, 2, var)`: returned, not used |
| Box-Cox | :72-75 | skipped (`abs(lambda - 1) > 1e-4` is FALSE for lambda = 1) |
| trend | :79-99, auxiliar.R:322-325 | `trend.spatial("cte")` is a column of 1s; for "cte" only `beta.ols <- colMeans(data)` (:99). **The data are not centred or detrended.** |
| R-side distances | :103 | `u <- as.vector(dist(as.matrix(coords)))`: **R's `dist()`**, column order (long, lat) |
| nugget tolerance | :104-113 | `nugget.tolerance` missing, so `1e-12` and `nt.ind <- FALSE`; then `nt.ind <- TRUE` if `min(u) < 1e-12` (co-located points) |
| umax | :131-132 | `umax <- max(u[u < max.dist])`: **strict `<`**, and it is a distance, not `max.dist` |
| bins | :133, 1023-1028 | `uvec = 13`; `bins.lim <- seq(0, umax, l = 14)`; `bins.lim <- c(0, 1e-12, bins.lim[bins.lim > 1e-12])` (15 values); `uvec <- 0.5*(bins.lim[-1] + bins.lim[-15])` (14 centres) |
| nbins | :135 | `length(bins.lim) - 1` = 14 (1 nugget bin + 13 regular bins) |
| first edge | :137 | `if (bins.lim[1] < 1e-16) bins.lim[1] <- -1`; the centres were already computed with 0 |
| binning | :138-149, geoR.c:380-405 | `.C("binit", n, coords[,1], coords[,2], data, nbins, bins.lim, modulus = FALSE, max.dist, ...)` |
| pair loop | geoR.c:380-386 | `for j in 0..n-1, for i in j+1..n-1: dx = xc[i]-xc[j]; dy = yc[i]-yc[j]; dist = hypot(dx, dy);` |
| range test | geoR.c:388 | `if (dist <= *maxdist)`: **inclusive**, against `max.dist` (the prctile) |
| estimator | geoR.c:390-392 | `v = sim[i]-sim[j]; v = (v*v)/2.0;` (classical: half the squared difference) |
| bin search | geoR.c:393-402 | `ind = 0; while (ind < nbins && dist >= lims[ind]) ind++; if (dist < lims[ind]) { vbin[ind-1] += v; cbin[ind-1]++; }`. A comment at :394-395 shows the older version, which could read past the end of `lims` |
| mean | geoR.c:407-414 | `if (cbin[j]) vbin[j] = vbin[j]/cbin[j];` |
| pairs.min | variogram.R:151-152 | `indp <- n >= 2`; `v[!indp] <- NA` |
| nugget removal | :153-159 | if `!nt.ind`, the first bin (nugget) is removed from `uvec`, `indp`, `bins.lim` and the result |
| output | :165-168 | `keep.NA = FALSE`, so `u = uvec[indp]`, `v = v[indp]`, `n = n[indp]` (as double), `bins.lim` (14 values after nugget removal, first one 1e-12) |
| message | :258-263 | if `nt.ind`, prints `"variog: co-locatted data found, adding one bin at the origin"` with `cat()`, **even when `messages = FALSE`** |

R's `seq.default`, `length.out` branch: `as.vector(c(from, from + seq_len(length.out - 2L) * (del/n1), to))`
with `del = to - from`. So the edges are `0 + k*(umax/13)` for k = 1..12, and the last edge is **exactly
`umax`**. (The alternative `(k*umax)/13` happened to agree on one test set, while `umax*(k/13)` differed in
1 edge. Use the `seq()` form.)

### Answers to question 1 (defaults, STcompare arguments)

* **Number of bins.** Internally 14: one nugget bin `[-1, 1e-12)` and 13 regular bins. The nugget bin is
  removed unless some R `dist()` value is below 1e-12. Then every bin with `n < 2` is **omitted**. So
  `length(u) = length(v) = length(n)` is at most 13 and depends on the coordinates:
  * kidney A/B: 10;
  * simRanPatternRasts[[1]]: 9;
  * synthetic 1000-point subsample: 13;
  * across all 100 simRanPatternRasts: 9 bins (93 rasters) or 10 bins (7 rasters).
* **Edges.** `e_0 = 1e-12`, `e_k = k*(umax/13)` for k = 1..12, `e_13 = umax`. Regular bin k is
  `[e_{k-1}, e_k)`, and bin 1 is `[1e-12, umax/13)`.
* **Centres `u`.** `0.5*(e_{k-1} + e_k)`. Note that **bin 1's centre is `0.5*(1e-12 + umax/13)`, not
  `umax/26`**. The nugget centre, when kept, is `5e-13`. The centres are computed before the first limit is
  set to -1 (variogram.R:133-137).
* **How `max.dist` is applied.** A pair is binned only if `hypot <= max.dist` **and** `hypot < umax`.
  Because `umax < max.dist` by construction, this is the same as **`hypot(dx,dy) < umax`**. Checked:
  `sum(n) == #(hypot < umax)` for all three coordinate sets (9178, 9695, 124815), while
  `#(hypot <= max.dist)` is 9312, 10053, 124864. Pairs at `umax`, and pairs between `umax` and `max.dist`,
  are silently dropped. On an integer grid with `max.dist = 13 + 1e-9`: `umax = 13`, `sum(n) = 154,494 =
  #(d < 13)`, and the 2,820 pairs at exactly 13 are dropped. With `max.dist = 13` exactly: `umax = √164`, and
  4,580 pairs with `√164 <= d <= 13` are dropped.
* **Pairs exactly on an edge.** The intervals are closed on the left, so the pair goes to the **upper** bin.
  On the integer grid with `umax = 13`, the edges are exactly 1..12 and `u = 1.5, 2.5, ..., 12.5`. Pairs at
  distance k are counted in `[k, k+1)`. Bin 1, `[1e-12, 1)`, is empty and therefore omitted.
* **The estimator.** The classical one: `v_k = (1/n_k) * Σ_{pairs in bin k} (z_i - z_j)^2 / 2`. Each term
  is computed as `t = z_i - z_j; t = t*t; t = t/2.0`. The terms are added to a double accumulator that
  starts at 0, in loop order (j outer ascending, i > j inner ascending, which is also the order of
  `as.vector(dist())`). Then the sum is divided by `(double) n_k`.
* **How `n` is counted.** Each unordered pair is counted once. The count does not depend on the data and is
  returned as double.
* **Dropped bins, and NA.** Bins with `n < pairs.min = 2` are dropped. **`v` never contains NA.** If the
  data contain NA, NaN or Inf, `.C` stops with `NA/NaN/Inf in foreign function call (arg 4)` (checked).

### Question 2: does anything depend on the data values?

No. `umax`, the edges, the pair-to-bin assignment, `n`, `indp`, `u` and `bins.lim` depend only on the
coordinates (and their column order) and on `max.dist`. `indp` uses only the counts (variogram.R:151). In
`03_verify.R`, all 69 data vectors per set produced identical `u` and `n`; these included Gaussian, uniform,
Poisson with many ties, heavy-tailed, 1e6-offset, the real kidney X, 9 locfit-smoothed permutations and
rescaled-plus-noise versions. **Target and candidate variograms in STcompare always have the same bins**,
because they use the same `ids` and the same coordinates, so `lm()` always gets vectors of equal length.

`trend = "cte"` only records `colMeans` (variogram.R:99). The data are used **raw**. Mathematically,
centring would not change `v`. In floating point it changes 1% of the entries, by up to 4e-16 relative
(`09_misc.R` (a)). A bit-exact port must not centre.

### Question 3: distances and tolerances

* **Binning uses C `hypot(dx, dy)`** (geoR.c:386; changelog `inst/doc/CHANGES:125`, "replacement of pythag
  by hypot", geoR 1.6-33). It is symmetric in its arguments (0 mismatches). It differs from
  `sqrt(dx*dx+dy*dy)` for 4.3-4.6% of pairs, and from R `dist()` for 4.6-5.7% of pairs, by at most 1 ulp. In R, `Mod(complex(real = dx,
  imaginary = dy))` is **bit-identical to C `hypot`** (0 mismatches over 1.04 M pairs), which is how the
  pure-R port reproduces it.
* **`umax` and the nugget test use R `dist()`** on `coords` in the (long, lat) order. **STcompare's
  `prctile` uses R `dist()` in the (lat, long) order.** R's `R_euclidean` runs `dist += dev*dev` once per
  column. On this CRAN arm64 build the compiler fused it into an FMA:
  `dist(cbind(a, b)) == sqrt(fma(db, db, da*da))` for every pair tested (0 mismatches over 1.04 M pairs),
  and `dist(cbind(a,b)) != dist(cbind(b,a))` for 9.4-11.3% of pairs (`dist_arith_output.txt`). This
  depends on the platform. An x86-64 build without `-mfma` gives `sqrt(da*da + db*db)`, and arm64 GCC
  (`-ffp-contract=fast` by default) probably uses FMA.
* **No jitter, scaling or handling of duplicate coordinates**, apart from the nugget bin for distances
  below 1e-12. With co-located points (`nt.ind`), the nugget bin `[-1, 1e-12)` is kept if `n >= 2`,
  `u[1] = 5e-13`, and the message above is printed on every call. Checked: kidney coordinates plus 3
  duplicated points give `u[1] = 5e-13, n = 3`, and both ports match.
* **Tolerances:**
  * `nugget.tolerance = 1e-12`;
  * the `bins.lim[1] < 1e-16` test;
  * the `bl > 1e-12` filter on edges;
  * `u < max.dist` (strict) for `umax`;
  * `dist <= max.dist` (inclusive) and `dist >= lims[ind]` / `dist < lims[ind]` in C;
  * `pairs.min = 2`;
  * the `u[1:2] < 1e-11` adjustment at :262, which only runs when co-located data exist and never triggers
    in practice.

  There are no others: no fuzz in `seq()`, and `quantile` type 7 uses `fuzz = 0`.

## 3. Ports and verification

* `WD/variog_fast.R`:
  * `variog_prep_r(coords, max.dist)` uses `dist()` for `umax`, `seq()` for the edges,
    `Mod(complex())` for `hypot`, and `findInterval(d, lims)`, which reproduces the C while loop exactly.
  * `variog_compute_r(prep, Z)` takes Z as N × B and uses `rowsum()`. `rowsum` adds rows in order into a
    double accumulator, so the sums are bit-identical.
* `WD/variog_fast.cpp` (`Rcpp::sourceCpp`). It sets `#pragma clang fp contract(off)`; explicit `std::fma`
  is used only to emulate R's `dist()`.
  * `variog_prep_cpp(coords, max_dist, umax = NA, rdist_fma = 1)`: computes `umax` by emulating R's
    `dist()`, builds the edges exactly as `seq()` does, runs the `binit` loop with `std::hypot`, and stores
    the pairs per output bin in CSR form, in loop order.
  * `variog_compute_cpp(prep, Z)`: one column at a time.
  * `variog_compute_t_cpp(prep, t(Z))`: walks the pair list once for all columns; the inner loop over
    columns vectorizes, and the order of each column's sum is unchanged.
  * `variog_compute_il_cpp`: a lock-step variant across bins. Same results but slower, so it is not
    recommended.
  * `rdist_quantile_cpp(pos, p)`: STcompare's `prctile`, type 7 using `nth_element`.
  * `pair_distances_cpp`: the three distance formulas, for diagnostics.
* Test coordinates (`WD/coords.rds`, built by `02_make_coords.R`):
  * kidney A/B shared pixels: 273 pixels, from
    `SEraster::rasterizeGeneExpression(speKidney, 'counts', resolution = 0.2, fun = 'mean', square = FALSE)`;
    hexagonal lattice, spacing 0.2;
  * simRanPatternRasts[[1]]: 280 pixels;
  * a 1000-point subsample (`set.seed(0); sample(16800, 1000)`) of a 120 × 140 hex grid with spacing 0.2.

`verify_output.txt` (`03_verify.R`): with `maxDistPrctile = 0.25`, plus 0.05, 0.1, 0.5 and 0.9 on the kidney
coordinates, both ports give maximum absolute differences of **0** for `u`, `n` and `v`, `identical()` is
TRUE for `u`, `n`, `v` and `bins.lim`, the C++ `umax` and `prctile` are identical to R's, and so is the
transposed kernel.

| set | N | `prctile` (= the quantile) | pairs total | pairs ≤ `max.dist` | pairs used | bins out |
|---|---|---|---|---|---|---|
| kidney A/B | 273 | 1.0392304845413269 (0.2√27) | 37,128 | 9,312 | 9,178 | 10 |
| simRan[[1]] | 280 | 1.0583005244258361 (0.2√28) | 39,060 | 10,053 | 9,695 | 9 |
| synthetic subsample | 1000 | 8 | 499,500 | 124,864 | 124,815 | 13 |
| kidney, p = 0.05 / 0.1 / 0.5 / 0.9 | 273 | 0.4 / 0.6 / 1.6 / 3.079 | 37,128 | 1,905 / 3,772 / 18,568 / 33,454 | 1,812 / 3,612 / 18,428 / 33,401 | 3 / 5 / 12 / 13 |

Kidney geoR output: `u` = 0.19985 0.35973 0.43967 0.51962 0.59956 0.67950 0.75944 0.83938 0.91932 0.99926
and `n` = 740 693 676 1270 611 584 1147 1609 1030 818. Regular bins 1, 2 and 4 are empty: the smallest
distance is 0.2, and there is no lattice distance in `[0.2398, 0.3198)`.

`fuzz_output.txt` (`10_fuzz.R`, 300 random configurations): `umax`, `u`, `n`, both `v` kernels and
`prctile` were identical in **300 of 300**.

### Timing (`bench_output.txt`, `05_bench.R`; bench::mark medians, 1 thread)

| set | N | pairs used | geoR / call | R prep (once) | C++ prep (once) | R, 1 column | R, B=100, per variogram | Rcpp, 1 column | Rcpp, B=100, per variogram | Rcpp transposed, B=100, per variogram |
|---|---|---|---|---|---|---|---|---|---|---|
| kidney A/B | 273 | 9,178 | 756 µs | 2.5 ms | 0.30 ms | 270 µs | 86 µs | 11.0 µs | 8.9 µs | **2.3 µs** |
| hex patch | 300 | 10,999 | 1,186 µs | 2.5 ms | 0.61 ms | 272 µs | 112 µs | 12.8 µs | 10.6 µs | **2.7 µs** |
| hex patch | 1000 | 120,933 | 9,776 µs | 29 ms | 6.5 ms | 2,258 µs | 1,155 µs | 122 µs | 119 µs | **31 µs** |
| 1000-pt subsample | 1000 | 124,815 | 10,436 µs | 31 ms | 6.7 ms | 2,482 µs | 1,271 µs | 126 µs | 116 µs | **33 µs** |

* Speedup of Rcpp over geoR: 69-93x for a single column, 312-432x batched.
* All variog calls for one gene (3602 calls): 2.7 s (N=273), 4.3 s (N=300) and 35-38 s (N=1000) with geoR.
  With Rcpp (one prep plus single-column calls): 0.04, 0.05 and 0.45 s.
* Where geoR's time goes at N=1000 (`extra_output.txt`): 10.9 ms in total, of which R `dist()` is 1.6 ms,
  `dist()` plus the `u < max.dist` subset plus `max` is 4.4 ms, and the `.C("binit")` call alone is 6.0 ms
  (mostly 500k `hypot` calls).
* `prctile`: R `quantile(dist())` takes 24 ms at N=1000, against 4.0 ms with the C++ `nth_element`.
* End to end (`pipeline_output.txt`): `matchingVariograms` on kidney, 10 permutations × 9 deltas, takes
  0.36-0.50 s with geoR and 0.18-0.26 s with the precomputed kernel. What remains is locfit.

## 4. The `lm()` step (`06_lm.R`, `lm_output.txt`, `07_pipeline.R`)

* `lm(target$v ~ 1 + cand$v)` is an ordinary least-squares fit of the target's v on the candidate's v, with
  intercept and slope, solved by LINPACK QR (`dqrls`, `tol = 1e-7`). `lm()`, `.lm.fit(cbind(1, x), y)` and
  `qr.coef(qr(X, tol = 1e-7), y)` gave identical coefficients.
* Closed form, with `x` = candidate v and `y` = target v over the K ≤ 13 shared bins: `xbar`, `ybar`,
  `Sxx = Σ(x - xbar)^2`, `Sxy = Σ(x - xbar)(y - ybar)`, `b1 = Sxy/Sxx`, `b0 = ybar - b1*xbar`. STcompare uses
  `sqrt(|b1|)` and `sqrt(|b0|)`. Over 450 real (target, candidate) pairs from kidney locfit permutations:
  maximum relative difference 7.1e-16 on `b1` and 3.4e-13 on `b0`. Over 10 permutations: same `delta*` in
  10 of 10, maximum relative difference in RSS 4.2e-14, maximum absolute difference in `hat` 1.1e-13.
* **Bit-exact `lm()` is not achievable reliably.** Calling R's own `F77_CALL(dqrls)` from C++
  (`WD/ols_dqrls.cpp`) on identical inputs gave **two different results depending on whether the buffers
  are 16-byte aligned**, and neither matched `lm()` exactly (`07e_align.R`). The BLAS level-1 routines in
  Accelerate depend on alignment. `lm()` itself is stable across 200 repeated calls. With the reference
  BLAS the results would not depend on alignment.
* Degenerate cases:
  * **Constant candidate variogram**, for example an all-zero or constant gene, where locfit's fitted
    values are constant and v = 0. `lm()` returns `b0 = mean(y)` and `b1 = NA` (rank deficient). The
    closed form gives NaN for both.
  * **Nearly constant candidate.** `dqrdc2` treats the column as negligible when its norm after the
    intercept reflection, which is about `||x - xbar||`, falls below `1e-7 * ||x||`. Observed: the slope
    was finite at a relative variation of 1.2e-7 and NA at 9e-8. The closed form returns a finite number in
    both cases.
  * **What happens after an NA slope.** `hat` becomes all NA, and the next `geoR::variog` call stops with
    `NA/NaN/Inf in foreign function call (arg 4)`. Running `matchingVariograms` on an all-zero gene hits
    exactly this error. The `tryCatch` in `spatialCorrelation` (585-619) then returns NA p-values.
  * **K < 2 bins** also gives an NA slope. This happens with a small `maxDistPrctile`: kidney with p = 0.01
    gives 1 bin (n = 354, made of the members of the 0.2 nearest-neighbour cluster whose `hypot` falls
    below `umax`). With p = 0.02, the bins have n = 740 and 2.
  * **NA bins** cannot occur (they are omitted). If one were passed, `lm`'s `na.omit` would drop that row,
    which is the same as the closed form on complete cases (checked).
* RSS: `sum()` on this platform equals a sequential double sum (`sizeof(long double) = 8`; checked over
  2000 random cases). On x86-64, R's `sum` uses 80-bit long double, so a port would need `long double` to
  match there.

## 5. What makes exact equivalence hard (`04_ties.R`, `ties_output.txt`)

1. **Summation order.** Each bin's sum has to run sequentially over its own pairs, in (j, i) loop order, in
   double, with no FMA and no tree or SIMD reduction across pairs. Interleaving across bins or across data
   columns is fine. BLAS-based indicator-matrix products would break exactness; `rowsum()` and the C++
   kernels do not. Accumulating `d*d` and halving at the end is exact except for subnormals, but the C++
   port keeps the per-term `/2.0` to stay safe.
2. **Three distance formulas** (section 2). Binning with the wrong formula changes pair membership:
   * kidney defaults: R `dist()` instead of `hypot` moves 2 pairs into or out of the range;
   * kidney with `max.dist = 2.6 + 1e-9`: `dist()` (variog column order) moves 8 pairs to another bin;
     `dist()` (quantile column order) changes 24 pairs (3 into or out of range, 21 to another bin);
     `sqrt(dx*dx+dy*dy)` moves 27 pairs to another bin.
3. **Tie clusters at `max.dist` and `umax`, on every rasterized dataset.** `prctile` is a sorted distance
   (type 7 with h = 0.75 interpolates only when `x[lo] != x[hi]`; in all three sets they were equal). On a
   lattice, that value belongs to a cluster of mathematically equal distances that differ by a few ulps
   because of noise in the coordinates. `umax` is the next lower member of the cluster, 1 ulp below
   `prctile` in all three sets. So the cluster is split by float noise:

   | set | cluster size | kept (`hypot < umax`) | dropped (= `umax`) | dropped (in `(umax, max.dist]`) | outside `max.dist` | distinct doubles in cluster |
   |---|---|---|---|---|---|---|
   | kidney (≈1.03923) | 481 | 328 | 59 | 75 | 19 | 14 (span about 18 ulp) |
   | simRan (≈1.05830) | 990 | 8 | 43 | 315 | 624 | 8 |
   | synthetic (≈8) | 107 | 7 | 17 | 32 | 51 | 7 |

   An implementation must therefore use the **identical coordinate doubles**: no re-centering, no scaling,
   no float32 storage, and no recomputation of pixel centroids. Any change in how SEraster computes the
   centroids changes `n` and `v`.
4. **Inner edge ties.** With the default `maxDistPrctile` on the three test sets, no pair is within 4 ulp of
   edges 1-12 (the edges fall between lattice distances). Ties do occur for other values of `max.dist`. With
   `max.dist = 2.6 + 1e-9`, the edges are about `0.2k`; on the kidney coordinates, 11 to 216 pairs lie
   *exactly* on each edge and 125 to 1,026 lie within 4 ulp. geoR puts exact ties in the upper bin, and
   both ports reproduce this.
5. **`lm()`** is not reproducible to the last bit across BLAS implementations or memory alignment. Compare
   RSS with a relative tolerance of about 1e-12, and require `delta*` to match exactly.
6. **Platform dependence of R itself:** FMA in `dist()` (arm64 against x86-64), `long double` in `sum()`,
   and libm's `hypot`. The C++ port must call the same libm `hypot` that geoR uses (`std::hypot` does).
7. **RNG streams** (`sample`, `rnorm` after `set.seed(seed + i)`) must be consumed in the same order. That is
   outside the variogram, but it is needed for parity of the whole pipeline.

---------------------------------------------------------------------------------------------------

## 6. Spec for a C++ implementer

**Inputs.** `pos` (N_full × 2, columns `x = pos[,1]` and `y = pos[,2]` exactly as stored), `ids`
(`1:N_full` if N_full ≤ 1000, otherwise `set.seed(seed); sample(N_full, 1000)`), `p = maxDistPrctile`
(default 0.25). Let `M = length(ids)`, `L[m] = x[ids[m]]` ("lat") and `G[m] = y[ids[m]]` ("long").

`RDIST(a, b)` = `sqrt(fma(b, b, a*a))` if R's `dist()` is FMA-contracted on this platform, otherwise
`sqrt(a*a + b*b)`, evaluated as separate operations. Detect this once at runtime by comparing R's `dist()`
on about 1000 random probe points with both formulas. Alternatively, compute `prctile` and `umax` in R with
`dist()`, since they are needed only once per dataset pair.

**A. Precompute (once per dataset pair; shared by all genes, both directions and all permutations).**

1. `prctile`. Enumerate the pairs `j < i` (any order). Let `d_q = RDIST(L_i - L_j, G_i - G_j)` and
   `P = M(M-1)/2`. Set `index = 1 + (P-1)*p`, `lo = floor(index)` and `hi = ceil(index)`; let `x_(k)` be the
   k-th smallest `d_q`. Then `q = x_(lo)`; if `index > lo` and `x_(hi) != x_(lo)`, set `h = index - lo` and
   `q = (1-h)*x_(lo) + h*x_(hi)` (two rounded products, then one rounded add, no FMA).
2. Variogram coordinates: `X1 = G` (long) and `X2 = L` (lat), because STcompare passes
   `coords = cbind(long, lat)`.
3. `umax = max{ RDIST(X1_i - X1_j, X2_i - X2_j) : that value < q }`. **The arguments are swapped relative to
   step A1.** If the set is empty, report an error, as geoR does (`-Inf` followed by a failing `seq()`). Set
   `nt = (min over all pairs of the same RDIST) < 1e-12`.
4. `s = umax / 13.0`. The edges are `e_k = 0.0 + (double)k * s` for k = 1..12, and `e_13 = umax`. Keep only
   the edges with `e_k > 1e-12`; normally all 13 are kept. Limits for the bin search: `lims = [-1, 1e-12,
   e_1, ..., e_13]`, giving 14 bins numbered 0..13. Centres: `c_0 = 0.5*(1e-12 + 0)`,
   `c_1 = 0.5*(e_1 + 1e-12)`, `c_k = 0.5*(e_k + e_{k-1})` for k = 2..13.
5. For `j = 0..M-1`, then `i = j+1..M-1`, in **this order**: `dx = X1_i - X1_j`, `dy = X2_i - X2_j`,
   `d = hypot(dx, dy)` (libm `hypot`, **not** `sqrt`). If `d <= q`: `ind = 0`; `while (ind < 14 &&
   d >= lims[ind]) ind++`; if `d < lims[ind]`, assign the pair to bin `ind - 1`. This is equivalent to:
   include the pair iff `d < umax`; its bin is 0 if `d < 1e-12`, otherwise `1 + #{k in 1..12 : e_k <= d}`.
6. `n_b` is the number of pairs in bin b. Output bins are the b in 1..13 with `n_b >= 2`, in ascending
   order, plus bin 0 if `nt` is true and `n_0 >= 2`. Output `u_b = c_b` and `n_b` (as double). For each
   output bin, store its pairs `(i, j)` in the order of step A5 (CSR), and discard all other pairs. If `nt`
   is true, geoR prints the "co-locatted data" message on every call.

**B. Variogram of a data vector `z` (length M, the values at `ids`; raw, not centred).** Reject NA, NaN and
±Inf, as geoR does. For each output bin b: `acc = 0.0`; for each stored pair in order: `t = z_i - z_j`;
`t = t*t`; `t = t/2.0`; `acc = acc + t`. Then `v_b = acc / (double) n_b`. Compile with FP contraction off
(or keep these as separate statements), and do not reassociate within a bin. Many data vectors can be
processed together: lay them out row-major (B values contiguous per point), loop over pairs, and inside
that loop over columns with SIMD. This is bit-identical and about 4x faster than one column at a time.

**C. Matching step.** `y = v(target)` and `x = v(candidate)` have the same length K. Then `b1 = Sxy/Sxx`
and `b0 = ybar - b1*xbar`. If `K < 2`, or `||x - xbar|| < 1e-7 * ||x||` (including `x` all zero), treat
`b1` as NA. In STcompare this leads to an error and an NA p-value row; a port may choose to handle it more
gracefully, but should document that choice. Then `hat = X.delta*sqrt(|b1|) + rnorm(N_full)*sqrt(|b0|)`,
using R's RNG in the same order. Compute `RSS = Σ_b (v_b(hat[ids]) - y_b)^2` as a sequential double sum
(long double on x86-64 to match R), and take `delta* = first argmin`. Expected agreement with STcompare/R:
`u`, `n` and `v` bit-identical; `b0` and `b1` within about 1e-15 relative (3e-13 for `b0` near zero); RSS
within about 1e-13 relative; `delta*` identical.

---------------------------------------------------------------------------------------------------

## 7. Files produced (all under WD)

* `variog_fast.R`: pure-R prep and compute.
* `variog_fast.cpp`: Rcpp prep, compute (3 kernels), C++ quantile and distance diagnostics.
* `ols_dqrls.cpp`: OLS through R's `dqrls`.
* `coords.rds`: test coordinate sets (kidney `pos`, X and Y; simRan[[1]] `pos` and X; synthetic 1000-point
  `pos`).
* `scripts/`:
  * `01_dist_arith.R` and `dist_arith.cpp`: which distance arithmetic R and geoR use;
  * `02_make_coords.R`: builds `coords.rds`;
  * `03_verify.R`: verification;
  * `04_ties.R`: ties and stress scenarios;
  * `05_bench.R`: timing;
  * `06_lm.R`: the `lm` step;
  * `07_pipeline.R`: end-to-end `matchingVariograms`;
  * `07b`-`07e` and `align_test.cpp`: debugging and the BLAS alignment demonstration;
  * `08_extra.R`: geoR time breakdown, lock-step kernel, edge formula, C++ quantile;
  * `09_misc.R`: centring, inclusion criterion, small percentiles, bin counts over 100 rasters;
  * `10_fuzz.R`: 300 random configurations;
  * `sum_test.cpp`: R `sum` accumulation.
* Outputs: `verify_output.txt`, `ties_output.txt`, `bench_output.txt`, `lm_output.txt`,
  `pipeline_output.txt`, `extra_output.txt`, `misc_output.txt`, `fuzz_output.txt`,
  `dist_arith_output.txt`, `verify_table.rds`, `bench_table.rds`.
* `geoR_1.9-6.tar.gz` and the extracted `geoR/` source.

Side note: in this R 4.5.2 / Bioconductor 3.22 setup, the first `SpatialExperiment::spatialCoords()` call in
a fresh session prints `Warning: stack imbalance in '::'` and similar warnings. This comes from loading the
namespace and is unrelated to this code (reproduced with a one-liner).
