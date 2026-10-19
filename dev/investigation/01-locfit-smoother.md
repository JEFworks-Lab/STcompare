# What the locfit smoothing step in STcompare computes, and how to port it to C++

Scope: the smoothing call in `matchingVariograms()` (R/spatialCorrelation.R:109-112):

```r
fit <- locfit::locfit(X.randomized ~ locfit::lp(long, lat, nn = delta[k], deg = 0),
                      kern = "gauss", maxk = 300)
X.delta <- fitted(fit)
```

All statements below were checked against the locfit **1.5-9.12** sources (the installed version; downloaded with
`download.packages()` and unpacked into `scratchpad/locfit/locfit/`). Citations use `file:line` relative to
`locfit/R/` and `locfit/src/` in that tarball. Every quantitative claim comes from a script in the work directory.
§11 lists the scripts and logs.

Work dir: `/private/tmp/claude-502/-Volumes-Crucial-SSD-Dropbox--Personal--work-github-com-slowkow-STcompare/d28687f8-19db-4d98-9d9d-f548807d4dba/scratchpad/locfit/` (called `$WD` below).

---

## TL;DR

* **Kernel:** `W(u) = exp(-(2.5 u)^2 / 2) = exp(-3.125 u^2)`, with `u = ||x_i - x_0|| / h` (Euclidean distance in raw
  coordinate units). Weights are **not truncated**. Every data point whose weight does not underflow to 0
  (u < 15.44) enters every local fit. h is the distance to the **k-th nearest data point, k = (int)(n*nn + 1e-12)**,
  counted from the fit point; a data point at distance 0 counts. The bandwidth combines as h = max(fixed h = 0, kNN distance).
* **Local fit:** at a fit point the estimate is the Nadaraya-Watson weighted mean (followed by one Newton step that
  changes only the last few bits).
* **What `fitted()` returns:** not that estimate at each data point. The default `ev = rbox()` builds locfit's
  **adaptive tree** (`"tree"`, cut = 0.8; this is *not* the `"kdtree"` structure). Local fits are computed only at
  tree vertices (13–235 vertices in all tests, independent of N). `fitted()` **bilinearly interpolates** the vertex
  values inside each terminal cell. Some vertices are "pseudo-vertices" whose value is the mean of two parent vertices.
  There are no derivatives and no cubic Hermite blending, because `deg = 0`.
* **`maxk = 300`** only triples the vertex capacity. It never changes the result. If the capacity is exceeded, locfit
  raises a hard **error** (`"newsplit: out of vertex space"`): no warning, no silent truncation. All 216 test fits
  (12 coordinate sets × 9 deltas × 2 responses) ran with no warning or error. The largest `maxk` any of them needed
  was 82, so even the default `maxk = 100` would have worked.
* **Linearity:** `fitted()` is exactly linear in y, with relative error ≤ 6e-16. The operator factors as
  **L = M W**:
  * W (n_real × N, n_real = 9–174) holds the normalised Gaussian weights at the real vertices.
  * M (N × n_real, ≤ 5 non-zeros per row) holds the interpolation weights.
  * Both depend only on the coordinates and delta.
* **Size of the gap to exact evaluation:** default `fitted()` differs a lot from exact (direct) evaluation at the
  data points. For permuted, spatially white input:
  * relative L2 difference: median 0.23 (range 0.09–0.44)
  * max absolute difference: up to 1.5 sd of the smooth
  * var(default) / var(direct): median 0.73 (range 0.42–0.90)
* **Self-contained C++ (`$WD/lfsmooth.cpp`):** the exact-direct smoother and the default tree smoother both match
  locfit **bit for bit** (max |diff| = 0 on every tested case), when compiled with `-ffp-contract=off`.
  * Plugged into a mirror of `viladomatCorrelation()`, the replica gives **identical** permutations, delta* and
    p-values to the package.
  * Per call it is 5.6–18× faster than `locfit()+fitted()`.
  * The precomputed factored operator makes the smoothing for one gene (1800 smooths) **1400–2700× cheaper**:
    34 s → 25 ms at N = 5000.
* **RNG:** `locfit()`/`fitted()` never touch the RNG.
* **Recommendation:** replicate locfit's tree and interpolation exactly, and implement it as a per-delta factored
  operator, computed once per coordinate set and applied to all permutations of all genes. Offer exact direct
  evaluation only as an explicit method change. It changes delta* for about 35 % of permutations, it is O(N²) per
  smooth, and a dense operator for it needs 1.7 GB at N = 5000.

---

## 0. The call and how its arguments reach C

* **Coordinate order.** `viladomatCorrelation()` sets `lat <- data[,3]` and `long <- data[,4]`
  (R/spatialCorrelation.R:247-248). `spatialCorrelation()` builds `data` as `(X, Y, x = pos[,1], y = pos[,2])`
  (R/spatialCorrelation.R:514-515). The first locfit coordinate is therefore `pos[,2]` and the second is `pos[,1]`.
  The order matters only for tie-breaking (split dimension) and for the floating-point order of the bilinear
  interpolation. An exact port must keep it.
* **Defaults that enter.** `lp(..., nn = delta, h = 0, adpen = 0, deg = 0, scale = FALSE)` (locfit.r:314-332)
  stores `alpha = c(nn, h, adpen)`. `locfit()` passes `kern = "gauss", maxk = 300` through `...` to `locfit.raw()`,
  whose other defaults are `kt = "sph"`, `ev = rbox()`, `family = "qgaussian"` (y not logical) and `link = "default"`
  (locfit.r:48-54, 77-78).
* **String matching.** `"gauss"` partially matches `"gaussian"` and becomes WGAUS. `"sph"` becomes KSPH, and
  `rbox()$type = "tree"` becomes ETREE (lfstr.c:80-115).
* **Scale.** `scale = FALSE` becomes `1 - as.numeric(FALSE) = 1` (locfit.r:114-117). The C side rescales only when
  scale ≤ 0 (startlf.c:110-119), so coordinates are used **unscaled**. `lp(scale = TRUE)` would divide each
  coordinate by its sd.
* **R to C.** `guessnv` sizes the work arrays (locfit.r:125-131; S_enter.c:535-625). Then `slocfit` runs
  `startlf(procv)` followed by `ressumm` (S_enter.c:257-363, 319-325).
* **`fitted()`.** `fitted.locfit` re-evaluates the model frame in the calling frame (locfit.r:348-387,
  locfit.matrix 1246-1288) and calls C `sfitted`, which calls `fitted()` (S_enter.c:466-500; fitted.c:68-98).

---

## 1. Gaussian kernel

* **Formula.** `W(u) = exp(-SQR(GFACT*u)/2.0)` with `GFACT = 2.5` (weight.c:31; local.h:137). With
  `u = ||x_i - x_0|| / h` this is `exp(-3.125 u²)`, i.e. a Gaussian with **σ = h / 2.5**. The weight at the k-th
  neighbour (u = 1) is e^-3.125 = 0.0439.
* **Spherical kernel.** `weightsph()` evaluates `W(di/h)` on the spherical distance `di` (weight.c:77-90). If h = 0
  the weight is 1 for di == 0 and 0 otherwise (weight.c:87).
* **No truncation.** The `u > 1` cut-off exists only for the compact kernels (weight.c:17-30). `iscompact(WGAUS) = 0`
  (weight.c:42-46).
* **Neighbourhood of a fit point.** `nbhd()` (lf_nbhd.c:209-276) does the following:
  1. Computes the distance from the fit point to **all n data points** (lf_nbhd.c:246-252).
  2. Computes h from those distances (lf_nbhd.c:261-264).
  3. Evaluates the weight of every point inside `xlim` (no `xlim` here) and keeps every point with `w > 0`
     (lf_nbhd.c:265-275).

  For the Gaussian kernel the neighbourhood is therefore every point except where `exp()` underflows to exactly 0:
  `exp(-745) = 4.9e-324` and `exp(-746) = 0`, so the cut-off is u > 15.44 ($WD/out/06_operator_stats.log). With
  nn ≥ 0.1 no point in any tested domain is that far away, so **every point enters every local fit**. The kNN
  distance only sets the scale; it is not a neighbourhood cut-off.

## 2. Nearest-neighbour bandwidth

* **k.** `nbhd(lfd, des, (int)(lfd->n*nn(sp)+1e-12), 0, sp)` (locfit.c:354-355). This is truncation (floor for
  positive values) with a 1e-12 guard against round-off.
  * Example: n = 273, nn = 0.1 gives 27.3, so k = 27.
  * For the default grid `seq(0.1, 0.9, 0.1)` (whose 3rd and 7th values are 0.30000000000000004 and
    0.70000000000000007) the guard never changes k for any n ≤ 20000. It can matter for user-supplied deltas, so a
    port should compute k the same way.
* **h from k.** `compbandwid()` (lf_nbhd.c:94-111) works as follows:
  * k = 0: h = fixh (lf_nbhd.c:101).
  * k < n: `nnh = kordstat(di, k, n, ind)`, the **k-th smallest (1-based) of all n distances**
    (lf_nbhd.c:50-75, 103-104).
  * k ≥ n: `nnh = max(di) * (k/n)^(1/d)` (lf_nbhd.c:105-108). This is the only "adjustment factor", and it applies
    only for nn ≥ 1.
  * Result: **h = max(fixh, nnh)** (lf_nbhd.c:110). With `lp(h = 0)` this is just the kNN distance. There is no
    other inflation: adaptive penalty `adpen = 0`, `acri = "none"`.
* **Does the fit point count?** The distance vector includes every data point, so a data point at the fit location
  (distance 0) is the first order statistic.
  * With `ev = dat()` (fit at the data points) h is the distance to the (k-1)-th *other* point. The C++ direct
    evaluator includes self and matches `ev = dat()` bit for bit (§7).
  * With the default tree, fit points are vertices. A vertex counts a coincident data point only by coincidence,
    e.g. bounding-box corners of a full rectangular grid.
* **Distance metric.** Spherical, i.e. Euclidean: `rho = sqrt(Σ_j (u_j/s_j)²)` with s_j = 1 here
  (lf_nbhd.c:14-44, KSPH branch 40-44). h differs at every fit point (adaptive), and is larger near the domain edge
  and in sparse regions.

## 3. Local-constant (deg = 0) estimate at a fit point x0

* **Starting value.** `reginit()` (locfit.c:220-239). For the Gaussian family with `LINIT`, `res[ZDLL] = w*y`
  (family.c:93-95), giving `s1 = Σ w_i y_i`, `s0 = Σ w_i` (prior weights 1) and `cf = (s1 - 0)/s0` (locfit.c:226-238).
* **One Newton step.** `max_nr()` evaluates the likelihood (m_max.c:183), solves for the step, then sets
  `coef = old + 1.0*delta` (m_max.c:203). It re-evaluates, and `likereg()` returns NR_BREAK for Gaussian/identity
  ("prevent iterations", locfit.c:160-162), so `max_nr` returns with the updated coefficient (m_max.c:205-206).
* **Result.** `f̂(x0) = m0 + Σ w_i (y_i − m0) / Σ w_i`, where `m0 = Σ w_i y_i / Σ w_i`. In exact arithmetic this is
  the **Nadaraya-Watson estimate**; the Newton step only refines the last bits. For p = 1 the solve is
  `delta = ((f1·dg)/(Z·dg²))·dg` with `dg = 1/sqrt(Z)` (m_jacob.c:43-52, 70-73; m_eigen.c:71-100).
* **Parametric component (matters only for bit-exactness).**
  * `compparcomp()` fits the same local model once with unit weights at x̄: this is the global mean ȳ, plus a Newton
    step (pcomp.c:55-123, unit weights at 71-81).
  * `procvraw()` stores *vertex value − ȳ* (procv.c:36-37; subparcomp pcomp.c:139).
  * Evaluation adds ȳ back (addparcomp pcomp.c:185-196, called at ev_interp.c:252).
  * Interpolation weights sum to 1, so this cancels mathematically. It must be mimicked only to reproduce locfit to
    the last bit.
* **Degree.** For `deg = 0`, `makecfn()` returns `ncoef = 1` (lf_fitfun.c:74-82), so no vertex derivatives exist and
  `hasd = 0` (startlf.c:147; S_enter.c:384).

## 4. Evaluation structure and what `fitted()` does

### 4.1 It is the adaptive tree, not the kd-tree

* `ev = rbox()` means `rbox(cut = 0.8, type = "tree")` (locfit.r:499-508), and `"tree"` is ETREE (lfstr.c:107-115).
  The construction is in ev_atree.c ("the default evaluation structure used by Locfit", ev_atree.c:5-7).
* The `"kdtree"` structure (ev_kdtre.c) splits cells at data medians, and its capacity ignores `maxk`
  (`kdtre_guessnv`, ev_kdtre.c:15-35). It is not used here.

### 4.2 Tree construction (`atree_start`, ev_atree.c:126-160)

1. **Root cell.** The root is the bounding box of the data, `[min, max]` of each coordinate, with no padding
   (set_flim, startlf.c:60-88). The 4 corner vertices are fitted: vertex i has coordinate k at `ur[k]` if bit k of i
   is set and at `ll[k]` otherwise (ev_atree.c:144-154).
2. **Split rule.** `atree_grow()` (ev_atree.c:83-124) recursively splits cells. `atree_split()` (ev_atree.c:55-78)
   splits a cell when `max_j side_j / h_min > cut`:
   * `h_min` is the smallest positive bandwidth among the cell's 4 corner vertices.
   * The split dimension is the argmax of `side_j / h_min` (first index on ties).
   * Splitting stops when every side is ≤ 0.8 × h_min.
3. **New vertices.** The split is at the midpoint. A new vertex is created at the midpoint of each of the two edges
   parallel to the split dimension (newsplit, ev_main.c:191-225). Adjacent cells share midpoints by looking up the
   parent pair (lo, hi) (findpt, ev_main.c:178-185). The tree has no explicit cell list (`nce = 1`): the vertices,
   their parents, their h values and the split rule define it implicitly.
4. **Pseudo-vertices.** If the edge being split is already short relative to *both* endpoints
   (`side < cut · min(h_i0, h_i1)`, ev_atree.c:109-110), the new vertex is a **pseudo-vertex**. It gets no local fit,
   `s = 1`, and `h = (h_i0 + h_i1)/2` (ev_main.c:213-216). Otherwise it is a real vertex fitted by `procv`
   (ev_main.c:217-221).
5. **Unused work.** `procv` also runs `comp_vari` (procv.c:119), and after the tree is built `ressumm` interpolates
   three quantities at every data point (frend.c:30-86). Neither affects `fitted()`; both are pure overhead for
   STcompare.
6. **The tree ignores y.** The tree depends only on the vertex bandwidths, which depend only on the coordinates and
   nn. Verified: identical `xev`, `s` and `h` for two different y ($WD/out/06_operator_stats.log).

### 4.3 Capacity and `maxk`

* **Capacity formula.** `atree_guessnv()` (ev_atree.c:16-48) gives, for d = 2 and nn ≤ 1:
  `nvm = floor(maxk/100 · floor((5/(nn·cut²) + 1)·4))`. N does not enter.
* **Capacity values.** Measured from `fit$nvc` ($WD/out/05_maxk.log):

  | delta | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
  |---|---|---|---|---|---|---|---|---|---|
  | nvm, maxk = 100 (locfit default, locfit.r:53) | 316 | 160 | 108 | 82 | 66 | 56 | 48 | 43 | 38 |
  | nvm, maxk = 300 (STcompare) | 948 | 480 | 324 | 246 | 198 | 168 | 144 | 129 | 114 |

* **When capacity runs out.** `newsplit` calls `ERROR(("newsplit: out of vertex space"))` (ev_main.c:202-206).
  `ERROR` is R's `error()` (local.h:128), so **locfit() fails with an R error**. There is no warning and refinement
  never stops silently. Demonstrated on quakes, delta = 0.1, which needs 235 vertices
  ($WD/out/05_maxk.log):
  * `maxk = 74` → `ERROR: newsplit: out of vertex space`
  * `maxk = 75` → ok (nvm = 237)
  * `fitted()` is `identical()` for maxk = 75, 300 and 1000.
* **Effect inside STcompare.** Such an error would be caught by `spatialCorrelation`'s `tryCatch`
  (R/spatialCorrelation.R:531, 585-586), printed, and turned into NA p-values.
* **Vertices actually needed.** Real-vertex counts are in the second table; minimum maxk is from
  $WD/out/default_vs_direct.csv:

  | set (N) | d=0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | min maxk needed |
  |---|---|---|---|---|---|---|---|---|---|---|
  | kidneyAB (273) | 151 | 51 | 45 | 45 | 41 | 27 | 15 | 15 | 15 | 63 |
  | simRan1 (280) | 150 | 45 | 45 | 45 | 39 | 27 | 15 | 15 | 15 | 60 |
  | quakes (998, clustered) | 235 | 113 | 85 | 67 | 53 | 45 | 35 | 23 | 23 | 82 |
  | kidneyCellsA (1229, raw cells) | 148 | 53 | 45 | 45 | 41 | 27 | 15 | 15 | 15 | 63 |
  | hexdisk 1000/2000/5000 | 81 | 81 | 65 | 25 | 25 | 25 | 25 | 25 | 21 | 61 |
  | square 1024/2025/5041 | 81 | 81 | 25 | 25 | 25 | 25 | 25 | 21 | 21 | 56 |
  | hex2000 / hex5000 | 137/149 | 59/45 | 45/45 | 37/41 | 25/27 | 25/15 | 21/15 | 15/13 | 15/13 | 46/50 |

  Real (non-pseudo) vertices n_real range from 9 to 174. The vertex count does **not** grow with N for a given
  shape: the hexdisk sets at N = 1000, 2000 and 5000 give exactly the same counts. **maxk = 300 never binds and
  never warns here** (216 fits, 0 warnings, 0 errors). It only buys a 3× safety margin for strongly clustered or
  irregular point patterns. Outside the default grid, delta = 0.02 needs up to 1354 vertices (quakes), against a
  capacity of 4698 at maxk = 300.

### 4.4 How `fitted()` gets values at the data points

The call chain is `fitted()` (fitted.c:68-98) → `dointpoint(..., PCOEF, ETREE, i)` (ev_interp.c:227-254) →
`atree_int()` (ev_atree.c:162-205). There is **no direct evaluation**. For each data point x:

1. **Descend the tree.** Start from the root cell with its 4 stored vertex values. At each level, re-derive the split
   with the same `atree_split` rule; this is possible because h is stored for every vertex, including pseudo-vertices.
   Go to the lower half if `2(x_ns − ll_ns) < (ur_ns − ll_ns)`; a point exactly on the midpoint goes to the upper
   half.
2. **Update corner values.** Each midpoint vertex is found with `findpt`. A real vertex supplies its stored
   coefficient. A pseudo-vertex supplies **(v_a + v_b)/2**, where v_a and v_b are the current values of the two
   corners on its edge (exvvalpv, nc = 1, ev_interp.c:163-166). Pseudo-vertices can chain.
3. **Interpolate in the leaf.** Use **bilinear** interpolation (rectcell_interp, `nc == 1`, ev_interp.c:73-80):
   first along coordinate 2, then along coordinate 1, with `linear_interp(h, d, f0, f1) = ((d−h)f0 + h f1)/d`
   (ev_interp.c:8-12). This is multilinear and uses **no derivatives**:
   * `exvval` returns nc = 1 because `hasd = 0` (ev_interp.c:148; S_enter.c:384).
   * The cubic-Hermite branch (ev_interp.c:83-97) would be used only for deg ≥ 1 or `dc = TRUE`.
4. **Finish.** Add ȳ, the parametric component (ev_interp.c:252). Add base = 0 (fitted.c:82). `resid(..., RFIT)`
   returns the fit (fitted.c:93, 40). The R `trans` is the identity for the identity link (locfit.r:190-200).
5. **Curiosity.** If any corner value equals exactly the sentinel `NOSLN = 0.1278433`, the interpolator returns that
   sentinel (ev_interp.c:70; local.h:136). Stored values are vertex fit − ȳ, so this is practically impossible.

**Check.** A pure-R replica of this descent, working from the fitted object's stored tree, reproduces `fitted()`
with max |diff| = **0** for delta = 0.1, 0.3, 0.5 and 0.9 ($WD/explore2.R, $WD/out/explore2_interp_replica.log).
The global term `fit$eva$pc[3]` equals `mean(y)` to 15 digits.

## 5. Linearity and the linear operator

**Linearity.** For coordinates and delta held fixed, the tree, the vertex bandwidths, the Gaussian weights and the
interpolation weights are all independent of y. The vertex values are weighted means, and the interpolation and
pseudo-vertex averaging are affine combinations whose weights sum to 1. The empirical test in
$WD/out/02_linearity_operator.log measured
`max|f(a·y1 + b·y2) − a·f(y1) − b·f(y2)| / max|f(a·y1 + b·y2)|` with y1 ~ N(0,1), y2 ~ Exp(1) and random a, b:

| set | nn = 0.1 | 0.3 | 0.5 | 0.9 |
|---|---|---|---|---|
| kidneyAB (273) | 2.4e-16 | 2.3e-16 | 2.7e-16 | 3.8e-16 |
| simRan1 (280) | 3.0e-16 | 1.7e-16 | 3.0e-16 | 3.3e-16 |
| quakes (998) | 3.3e-16 | 3.2e-16 | 1.9e-16 | 5.2e-16 |
| hexdisk1000 | 3.1e-16 | 3.0e-16 | 4.0e-16 | 2.4e-16 |
| hexdisk2000 | 3.9e-16 | 5.3e-16 | 4.2e-16 | 3.8e-16 |
| square5041 | 4.0e-16 | 3.5e-16 | 6.0e-16 | 3.9e-16 |

Shift equivariance also holds: `f(y + 1000) − 1000 = f(y)` to 2e-13 – 9e-12 relative to max|f(y)|. That is
round-off on values of size 1000 (1000 × 1.1e-16 ≈ 1e-13).

**Structure: X.delta = L y with L = M W.**
* W is n_real × N. Row v holds `w_vi / Σ_i w_vi`, the normalised Gaussian weights at real vertex v.
* M is N × n_real. It holds the bilinear and pseudo-vertex weights, with **at most 5 non-zeros per row** across all
  sets (06_operator_stats.log) and rows summing to 1.
* So L has rank ≤ n_real (9–174), far below N.
* All 9 deltas' factors take 1.4 MB (N = 273) to 24 MB (N = 5000) when stored dense. Dense N×N operators for the
  same 9 deltas would take 5 MB to 1.7 GB.

**Three exact extraction routes.** Seconds; the error is max |L·y − fitted| on a fresh y
($WD/out/operator_extraction.csv).

| N (set), nn | (a) N unit vectors through locfit | (b) locfit internals: `geth = 1` + `fitted()` with unit vertex coefs | (c) C++ `lf_tree_operator` (+ dense M W) | error a / b / c |
|---|---|---|---|---|
| 273 (kidneyAB), 0.1 | 1.02 | 0.049 | 0.001 (+0.001) | 2.2e-16 / 1.3e-15 / 1.2e-15 |
| 273, 0.5 | 0.40 | 0.011 | 0.001 | 1.4e-16 / 1.9e-16 / 1.9e-16 |
| 1000 (hexdisk), 0.1 | 8.35 | 0.082 | 0.002 (+0.001) | 4.4e-16 / 8.9e-16 / 8.3e-16 |
| 1000, 0.5 | 3.88 | 0.020 | 0.002 | 1.4e-16 / 2.7e-16 / 3.1e-16 |
| 2000 (hexdisk), 0.1 | 30.6 | 0.137 | 0.004 (+0.005) | 4.2e-16 / 6.7e-16 / 6.4e-16 |
| 2000, 0.5 | 13.5 | 0.040 | 0.002 | 1.5e-16 / 1.6e-16 / 1.7e-16 |
| 5000 (hexdisk), 0.1 | not run (≈ 3–4 min expected) | 0.379 | 0.010 (+0.023) | – / 5.7e-16 / 6.0e-16 |

The brute-force L and the C++ L differ by ≤ 1.8e-16 elementwise.

**How route (b) works:**
* `locfit(..., geth = 1)` returns the N × nv matrix of vertex hat weights, i.e. Wᵀ. Pseudo-vertex columns are 0, and
  real columns sum to 1 (locfit.r:167-168; S_enter.c:331-333; procvhatm procv.c:213-227).
* Column j of M comes from `fitted()` on the fit object after setting `eva$coef[,1] <- e_j` and `eva$pc[3] <- 0`.
* $WD/explore3.R shows the method; the log confirms `t(Lv) %*% y = coef + ȳ` at real vertices (2.2e-16).

## 6. How far default `fitted()` is from exact direct evaluation

**What was compared.** Script `$WD/01_default_vs_direct.R`; full table in `$WD/out/default_vs_direct.csv`.
* "Direct" means the exact Gaussian NW estimate at every data point, with the same kNN bandwidth computed at that
  point (self included). This is exactly `locfit(..., ev = dat())`, reproduced bit for bit by `lf_direct()`.
* Two inputs were used:
  * **white:** a random permutation of real expression values (kidney X, simRan values, quakes depth, kidney cell
    counts) or N(0,1) for synthetic grids. This is the realistic input, because STcompare smooths *permuted* X.
  * **signal:** spatially structured y (unpermuted real values, or a sinusoidal surface plus noise on the grids).
* Coordinate sets:
  * rasterized speKidney A∩B shared pixels (N = 273);
  * `simRanPatternRasts[[1]]` (280);
  * quakes from the `spatialCorrelation` example, duplicates removed (998);
  * raw speKidney A cells (1229);
  * square grids (1024, 2025, 5041), hex grids (2000, 5000) and hex disks (1000, 2000, 5000).

**Ranges over the 9 deltas** (max_over_sd = max|default − direct| / sd(direct); relL2 = ‖default − direct‖ / ‖direct − mean‖):

| set (N) | y | max_over_sd | relL2 | cor | var(default)/var(direct) |
|---|---|---|---|---|---|
| kidneyAB (273) | white | 0.41–0.98 | 0.13–0.25 | 0.973–0.994 | 0.73–0.87 |
| kidneyAB | signal | 0.16–0.80 | 0.05–0.28 | 0.992–0.999 | 0.59–0.94 |
| simRan1 (280) | white | 0.48–0.98 | 0.18–0.30 | 0.957–0.992 | 0.68–0.76 |
| quakes (998) | white | 0.61–0.82 | 0.18–0.24 | 0.984–0.993 | 0.67–0.76 |
| kidneyCellsA (1229) | white | 0.63–1.53 | 0.20–0.38 | 0.939–0.987 | 0.60–0.78 |
| square1024 | white | 0.63–1.26 | 0.21–0.42 | 0.949–0.987 | 0.45–0.74 |
| hexdisk1000 | white | 0.34–1.43 | 0.12–0.37 | 0.941–0.996 | 0.64–0.83 |
| hex2000 | white | 0.35–0.84 | 0.11–0.26 | 0.982–0.996 | 0.64–0.88 |
| hexdisk2000 | white | 0.40–1.05 | 0.16–0.30 | 0.964–0.991 | 0.66–0.81 |
| square2025 | white | 0.50–0.97 | 0.16–0.33 | 0.964–0.990 | 0.58–0.84 |
| hex5000 | white | 0.64–1.19 | 0.18–0.42 | 0.923–0.988 | 0.57–0.82 |
| hexdisk5000 | white | 0.26–0.94 | 0.10–0.25 | 0.973–0.998 | 0.67–0.90 |
| square5041 | white | 0.57–1.14 | 0.18–0.44 | 0.945–0.990 | 0.42–0.81 |

Over all white-noise runs, relL2 has median **0.23** (max 0.44) and the variance ratio has median **0.73**
(range 0.42–0.90). For signal runs, relL2 has median 0.18 and the variance ratio median 0.79. In absolute terms the
max differences are of the same order as the smooth's own sd. The default (tree-interpolated) output is
systematically **smoother and lower-variance** than exact evaluation. The gap does not shrink with N, because the
vertex count is set by geometry and delta, not by N.

Exactness of the comparison itself: `lf_direct()` equals `fitted(locfit(..., ev = dat()))` with max |diff| = 0 for
all 18 tested (set, delta) combinations ($WD/out/test_cpp1_exactness.log).

**Downstream effect, method change vs locfit.** Script `$WD/04_downstream.R`; table in
`$WD/out/downstream_locfit_vs_direct.csv`. The direct smoother was plugged into an exact mirror of
`viladomatCorrelation()`. Runs: kidney A–B, kidney A–C and 20 independent simRan pairs (i vs i+50), both
directions, B = 100, seed 0 (44 runs).
* Only **65 %** of permutations on average (range 43–95 %) pick the same delta*. The median delta* moves **up by
  0.1 in 24/44 runs** (kidney: 0.2 → 0.3), because the rougher direct smoother needs more smoothing to match the
  variogram.
* The variogram-matching step compensates for most of this. Null-correlation sd ratio (direct/locfit) is 0.98 on
  average. p-values correlate at 0.996, with mean |Δp| = 0.022 and max 0.07 (p resolution is 0.01).
* For the 40 independent-field runs, p < 0.05 occurred 2/40 times with locfit and 1/40 with direct.

So direct evaluation is a genuine (moderate) method change: delta* changes substantially, p-values change modestly.

## 7. Self-contained reimplementation: `$WD/lfsmooth.cpp` (Rcpp)

**Build:** `Sys.setenv(PKG_CXXFLAGS = "-ffp-contract=off"); Rcpp::sourceCpp("lfsmooth.cpp")`. The code mirrors the
locfit C (comments cite file:line).

| function | what it computes |
|---|---|
| `lf_direct(xy, y, nn)` | exact direct evaluation at the data points (`ev = dat()` replica) |
| `lf_direct_multi(xy, Y, nn)` | direct evaluation of B columns at once |
| `lf_direct_operator(xy, nn)` | dense N×N direct operator |
| `lf_tree(xy, nn, maxk, cut)` | the adaptive tree only: vertices, h, pseudo flags, parents, capacity |
| `lf_tree_fitted(xy, y, nn, maxk)` | **exact replica of `fitted(locfit(...))`**; same overflow error |
| `lf_tree_operator(xy, nn)` | factored operator: `list(M = N×n_real, W = n_real×N, …)` |
| `lf_apply_factored(M, W, Y)` | apply the factored operator to the columns of Y |

**Agreement with locfit.**
* `lf_tree_fitted` vs `fitted()`: max |diff| = **0**, with an identical tree (same xev, s, lo, hi, h). Tested on
  kidneyAB and simRan1 for all 9 deltas, and on all benchmark sets (quakes, hexdisk 1000/2000/5000; `err_cpp_tree`
  = 0 in benchmark.csv).
* `lf_direct` vs `ev = dat()`: **0**.
* Factored operator `M %*% (W %*% y)` vs `fitted()`: ≤ 1e-15 (different summation order).

**FMA caveat.** The first build used clang's default FP contraction, which fuses `s += u*u` into an FMA. That
changed a few vertex h values by 1 ulp: trees were still identical, and fitted values differed by ≤ 2.5e-16. With
`-ffp-contract=off` both paths became bit-identical to the CRAN locfit binary, so the installed locfit binary does not use
FMA here. A port that wants bit-exact regression tests against locfit should compile its distance code without FP
contraction (or write it so contraction cannot occur).

**Per-call timing.** Medians, single thread, `VECLIB_MAXIMUM_THREADS=1`, Apple M1 Ultra
($WD/out/benchmark.csv). "Op" means the factored tree operator; "batched" means B = 100 vectors per call.

| set (N), nn | n_real | locfit+fitted (ms) | C++ exact tree (ms) | speed-up | C++ direct (ms) | op setup (ms) | op apply, 1 vec (µs) | op apply batched (µs/vec) | dense direct L apply batched (µs/vec) |
|---|---|---|---|---|---|---|---|---|---|
| kidneyAB (273), 0.1 | 127 | 3.61 | 0.46 | 7.9× | 0.79 | 0.77 | 26.7 | 1.2 | 1.1 |
| 273, 0.5 | 29 | 1.46 | 0.13 | 10.9× | 0.84 | 0.36 | 6.9 | 0.65 | 1.1 |
| 273, 0.9 | 13 | 1.03 | 0.057 | 18× | 0.70 | 0.24 | 3.9 | 0.51 | 1.1 |
| hexdisk1000, 0.1 | 81 | 7.99 | 1.21 | 6.6× | 11.5 | 2.1 | 61 | 3.0 | 16.8 |
| 1000, 0.5 | 25 | 3.43 | 0.47 | 7.2× | 13.3 | 1.1 | 21 | 1.8 | 17.3 |
| 1000, 0.9 | 13 | 2.48 | 0.27 | 9.3× | 11.2 | 0.97 | 11 | 1.5 | 16.9 |
| quakes (998), 0.1 | 174 | 14.8 | 2.64 | 5.6× | 13.0 | 4.0 | 136 | 5.1 | 19.4 |
| hexdisk2000, 0.1 | 81 | 14.4 | 2.22 | 6.5× | 42.1 | 3.8 | 126 | 6.1 | 55 |
| 2000, 0.5 | 25 | 6.31 | 0.88 | 7.2× | 48.8 | 2.1 | 45 | 5.2 | 62 |
| 2000, 0.9 | 13 | 4.62 | 0.53 | 8.8× | 44.7 | 1.9 | 25 | 3.6 | 57 |
| hexdisk5000, 0.1 | 81 | 34.8 | 5.79 | 6.0× | 274 | 9.8 | 346 | 19.8 | (200 MB/delta, not run) |
| 5000, 0.5 | 25 | 16.7 | 2.37 | 7.1× | 312 | 5.7 | 120 | 11.9 | – |
| 5000, 0.9 | 13 | 10.4 | 1.33 | 7.8× | 258 | 4.9 | 70 | 10.8 | – |

Reading the table:
* The exact tree replica is 5.6–18× faster per call, even without precomputation.
* Most of what remains is the tree build (kNN distance at each vertex). For N = 5000 and nn = 0.1 the build alone is
  2.6 ms, and it can be cached because it does not depend on y.
* Direct evaluation is O(N²): it is already *slower* than locfit's default from N ≈ 1000 (11–13 ms vs 2.5–8 ms), and
  reaches 260–310 ms per call at N = 5000.
* Where locfit's own time goes ($WD/out/07_overhead.log): `.C` takes about 75 % of a call, and model-frame/fitted
  overhead about 20 %. `locfit.raw` with `ev = none()` (no vertices) costs only 0.09–0.26 ms. The cost is the vertex
  fits plus `comp_vari`/`ressumm`, which STcompare does not need.

**Smoothing cost for one gene** under the default call pattern: 2 directions × 100 permutations × 9 deltas = 1800
smooths ($WD/out/per_gene_cost.csv).

| N | locfit (s) | C++ exact, per vector (s) | operator setup, once per coordinate set (s) | operator apply, per gene (s) | locfit / operator-apply |
|---|---|---|---|---|---|
| 273 | 2.73 | 0.27 (10×) | 0.0035 | 0.0010 | 2700× |
| 1000 | 7.74 | 1.07 (7×) | 0.013 | 0.0036 | 2100× |
| 2000 | 14.3 | 2.07 (7×) | 0.025 | 0.0087 | 1600× |
| 5000 | 34.3 | 4.99 (7×) | 0.061 | 0.025 | 1400× |

The max error of the operator path vs locfit over all 9 deltas is ≤ 5e-16. The operators depend only on pixel
coordinates and delta, so for an SE pair they are built once and reused for every gene, both directions and all
permutations.

**End-to-end check** ($WD/out/04_downstream.log). Kidney A vs B, B = 100, run through a mirror of
`viladomatCorrelation()` with each smoother:
* mirror + locfit vs `STcompare::viladomatCorrelation`: identical p, delta* and permutations (max diff 0);
* mirror + C++ exact tree vs package: **identical** (max diff 0);
* mirror + precomputed operator vs package: identical p and delta*, permutations within 2.9e-11 (rounding amplified
  by the variogram rescaling).

Elapsed time for one direction: package 3.8 s, C++ tree 2.9 s, operator 2.1 s. At N = 273 the remaining time is
dominated by the 1800 `geoR::variog` calls (0.78 ms each, $WD/explore4.R). This is outside the scope of this report,
but it will be the next bottleneck once the smoother is fixed. Because X.delta = L y with L fixed, each binned
semivariance is a quadratic form in y, which a port could exploit.

## 8. Random number generator state

* `.Random.seed` is `identical()` before and after each of the following: `locfit()`, `fitted()`, an `ev = dat()`
  fit plus `fitted()`, and a `geth = 1` fit. The next `runif()` is identical whether or not a locfit call ran in
  between ($WD/out/02_linearity_operator.log).
* Source: the only RNG use in locfit is Monte-Carlo integration in `monte()` (m_imont.c:22-43, GetRNGstate/unif_rand),
  which is called only for simultaneous-confidence-band constants (scb_cons.c:454-455). It is never on the
  regression/`fitted()` path.
* Side finding relevant to a port (not about locfit): `matchingVariograms()` calls `set.seed(seed + i)` (R/spatialCorrelation.R:99, 294),
  but it runs inside `BiocParallel::bplapply`, whose workers use **L'Ecuyer-CMRG**. This holds for both
  `MulticoreParam` and `SerialParam`. A serial call in a default Mersenne-Twister session gives a different
  `hat.X.delta.star`. Re-running serially under `RNGkind("L'Ecuyer-CMRG")` reproduces the worker result exactly
  ($WD/08_rng_workers.R, $WD/out/08_rng_workers.log). A port that wants identical output must reproduce
  `set.seed(seed+i)` + `rnorm()` under L'Ecuyer-CMRG, or call back into R's RNG.

## 9. Pseudo-code of the exact smoother (for the port)

```
input: coords (N×2, column 1 = pos[,2], column 2 = pos[,1]), y (N), nn, cut = 0.8, maxk = 300
k     = (int)(N*nn + 1e-12)
bandw(p) = k-th smallest of { sqrt((x_i1-p1)^2 + (x_i2-p2)^2) : i = 1..N }        # self included
w(d,h)   = h == 0 ? (d == 0) : exp(-((2.5*|d/h|)^2)/2)
ybar  = mean(y)  [+ Newton step]                                                  # parametric component
fit(p)= NW mean of y with weights w(d_i, bandw(p)) [+ one Newton step] - ybar

tree:  root = bounding box; vertices 0..3 = corners (bit k of id -> upper bound in coordinate k)
       h[v] = bandw(vertex v) for the 4 corners
grow(cell):
       hmin = min positive h over the 4 corners; score_j = side_j / hmin
       ns = argmax score (first index); if !(cut < score_ns) return
       for each of the 2 edges parallel to ns with endpoints (a, b):
           m = lookup(a, b) or create at midpoint:
               pseudo if side_ns < cut*min(h[a], h[b]):  h[m] = (h[a] + h[b])/2, no fit
               else h[m] = bandw(m), coef[m] = fit(m)
           (error if #vertices would exceed nvm = floor(maxk/100*floor((5/(nn*cut^2)+1)*4)))
       grow(lower half); grow(upper half)

fitted(x):
       cell = root, vals = coef[0..3]
       while split(cell) = ns != -1:
           lower = 2*(x_ns - ll_ns) < (ur_ns - ll_ns)
           for each pair (i, i+2^ns): m = lookup(ce[i], ce[i+2^ns])
               v_m = pseudo[m] ? (vals[i] + vals[i+2^ns])/2 : coef[m]
               replace corner i+2^ns (if lower) or corner i (if upper) by m / v_m
       bilinear interpolation in the leaf (coordinate 2 first, then 1) + ybar
```

The operator form computes the same thing as `M[i, ·]` (the interpolation weights expanded through
pseudo-vertices, ≤ 5 non-zeros per row) times `W[v, ·]` (normalised weights at each real vertex v).

## 10. Recommendation for the C++ port

Three candidates:
1. Replicate locfit's adaptive-tree interpolation exactly.
2. Use exact direct evaluation.
3. Precompute the linear operator.

**Recommended: (1) and (3) together.** Replicate the tree and bilinear interpolation exactly, and implement them as a
precomputed, factored, per-delta linear operator. Keep exact direct evaluation only as an opt-in alternative.

1. **Exact replication is cheap, proven and keeps results unchanged.**
   * About 300 lines of C++ reproduce `fitted()` bit for bit (with `-ffp-contract=off`).
   * Inside the full `viladomatCorrelation()` algorithm they give identical permutations, delta* and p-values to
     the current package.
   * This allows regression tests against the existing package (and existing results in `inst/extdata/`) at zero
     tolerance.
2. **Direct evaluation changes the method and is slower.**
   * Its output differs from locfit's by relL2 ≈ 0.23 (up to 0.44) and has about 1.4× the variance.
   * It changes delta* for about 35 % of permutations and p-values by up to 0.07 at B = 100.
   * Cost is O(N²) per smooth: 11 ms at N = 1000 and about 300 ms at N = 5000, i.e. slower than locfit itself at
     these sizes. nn = 0.1–0.9 means every bandwidth covers 10–90 % of the points, so neighbour truncation cannot make
     it sparse.
   * Dense precomputed direct operators cost 8·N² bytes per delta: 1.7 GB for 9 deltas at N = 5000.
   * If the authors consider direct evaluation more faithful to Viladomat et al., expose it as an explicit option
     (e.g. `smoother = "direct"`) and document it as a change.
3. **The factored operator is the fastest route and is exact to round-off.**
   * The tree, W (n_real × N, n_real ≤ 174 in all tests) and M (sparse, ≤ 5 nnz/row) depend only on the shared
     pixel coordinates and delta.
   * All 9 operators for a coordinate set cost 4–61 ms for N = 273–5000 and ≤ 24 MB dense; M could also be sparse.
   * Smoothing all 200 permutations of a gene for one delta is then two small dense products, which BLAS handles or
     which parallelise trivially.
   * Per-gene smoothing time falls from 2.7–34 s to 1–25 ms (1400–2700×), with errors ≤ 1e-15. In the end-to-end
     test, p-values and delta* were identical.
   * For bit-exact regression mode, also keep the per-vector exact path (`lf_tree_fitted`, 5.6–18× faster than locfit).
4. **Port details that matter:**
   * keep the coordinate order (pos[,2], pos[,1]);
   * k = (int)(N·nn + 1e-12);
   * bounding box = data range;
   * pseudo-vertex h = mean of its parents;
   * find vertices by (lo, hi) parent pair;
   * descend with `2(x−ll) < (ur−ll)`;
   * bilinear interpolation in the order coordinate 2 then 1;
   * for bit-exactness, subtract and add back ȳ and apply the Newton step;
   * either reproduce the `maxk` capacity error or drop the cap: results are identical whenever locfit succeeds, so
     dropping it only removes failures that STcompare currently turns into NA p-values.
5. **Next bottleneck.** After this change, `geoR::variog` dominates (3600 calls per gene, about 0.8 ms each at
   N ≈ 280). Since X.delta = L·y_perm with L fixed, the binned semivariances are quadratic forms, which opens further
   speed-ups.

## 11. Files produced (all under `$WD`)

**Sources and data**
* `locfit_1.5-9.12.tar.gz` and `locfit/`: the locfit sources that were read.
* `00_make_coords.R`, `00b_more_coords.R` → `data/coords.rds` (12 coordinate sets in STcompare's (long, lat)
  order) and `data/values.rds` (kidney X/Y, simRan values, quakes depth/mag, raw kidney cell counts).
* `lfsmooth.cpp`: the C++ reimplementation (§7).

**Exploration and exactness scripts**
* `explore1.R`: formula checks (vertex h = kNN distance, NW at the data points).
* `explore2.R`: pure-R replica of `fitted()` from the stored tree.
* `explore3.R`: the `geth = 1` operator route.
* `explore4.R`: variog timing and BiocParallel RNG kind.
* `test_cpp1.R`, `test_cpp2.R`: bit-exactness tests (`out/test_cpp1_exactness.log`).

**Experiments**

| script | question | outputs |
|---|---|---|
| `01_default_vs_direct.R` | §4.3 and §6 | `out/default_vs_direct.csv`, `out/01_default_vs_direct.log` |
| `02_linearity_operator.R` | §5 and §8 | `out/operator_extraction.csv`, `out/02_linearity_operator.log` |
| `03_benchmark.R` | §7 per-call timing | `out/benchmark.csv`, `out/03_benchmark.log` |
| `04_downstream.R` | end-to-end validation and method-change impact | `out/downstream_locfit_vs_direct.csv`, `out/04_downstream.log` |
| `05_maxk.R` | capacity and overflow error | `out/05_maxk.log` |
| `06_operator_stats.R` | operator sparsity, size, y-independence | `out/06_operator_stats.log` |
| `07_overhead.R` | where locfit spends time | `out/07_overhead.log`, `out/rprof_locfit.out` |
| `08_rng_workers.R` | RNG kind in BiocParallel workers | `out/08_rng_workers.log` |
| `09_per_gene_cost.R` | per-gene smoothing cost | `out/per_gene_cost.csv`, `out/09_per_gene_cost.log` |

Nothing inside the repository was modified.
