# Spec: rank-remapped surrogates for `compareSpatial()`

Addresses `cpp-backend-plan.md` §8: the Viladomat surrogates (smoothed permutation plus Gaussian noise, rescaled to match the variogram) have a near-Gaussian marginal distribution. For a gene detected in few pixels, the permutation distribution of r is heavy-tailed, the surrogates do not reproduce that, and the p-values in the far tail are too small. The maintainer agreed (2026-10-05) to keep `minDetected` as a filter the user can switch off (`minDetected = 0`) and to try rank-remapped surrogates.

## 1. Method

After the engine has built a surrogate s (all N pixels) for a permutation of the source vector x:

1. compute the ranks of s (ties broken by pixel index; ties have probability zero because of the Gaussian noise);
2. replace s by the sorted values of x placed in that rank order: `s_remapped[order(s)] = sort(x)`.

Every surrogate then has exactly the marginal distribution of x, including its zeros, while keeping the spatial arrangement of the Viladomat surrogate (its ranks). This is the amplitude-adjustment step of AAFT surrogates (Theiler et al. 1992); a single step, no iteration.

Properties:
- For x without spatial structure, the smoothed surrogate's ranks are a uniformly random permutation independent of x's values, so the remapped surrogate is a plain permutation of x and the test is the exact permutation test.
- For a smooth x, the remap changes the variogram slightly (the rank transform is monotone, so the spatial pattern is preserved up to a monotone distortion of the amplitudes).
- The null correlation is computed from the remapped surrogate; the remapped surrogate is also what `keep_surrogates` stores.

## 2. Scope

- `compareSpatial()` only, through a new argument `surrogate = c("remap", "gaussian")`. The default is decided from the calibration study below (§4); until then the implementation keeps `"gaussian"` as the default so that existing tests are unchanged.
- The legacy functions are untouched. Legacy mode must stay bit-identical: `bench/validate-published.R` (exported mode) must still report 0 mismatches.
- `minDetected = 0` must be accepted and must skip only constant genes.

## 3. Engine

- A per-task flag `remap` and the task's sorted source values (computed once when the task is defined, from the pool column; NA-free by the pre-checks).
- In the per-permutation step, after `out[i] = out[i] * a1 + e[i] * a0` and before the correlations: fill a thread-local index buffer 0..N-1, sort it by `out` (stable on ties by index), then `out[idx[r]] = sorted_source[r]`.
- Cost: one sort of N doubles per permutation and direction, which is small next to the smoothing (N × m) and the variograms (about 120k pairs).
- Determinism: the step is a fixed sequence of operations per (task, permutation), so results stay independent of threads, batches and chunks.

## 4. Calibration and power study (`bench/calibrate-surrogates.R`, results in `bench/calibration-results.md`)

All runs with the C++ engine; nothing slow.

1. **Null, dense:** all 4,950 `simRanPatternRasts` pairs (or a large random subset), both surrogate modes, default settings. Report P(p ≤ α) for α ∈ {0.1, 0.05, 0.01, 0.001} with binomial standard errors.
2. **Null, sparse:** independent pairs of zero-inflated spatial fields on the AKI coordinates (311 pixels) and on the brain coordinates (2,170 pixels), with detection fractions from about 1% to 100% (reuse the generator of the acceptance review, `scratchpad/acc-stats-api/scripts`, if still available; otherwise: a smooth Gaussian field thinned by a Bernoulli mask, and Poisson counts with low means). Report P(p ≤ α) by detection-fraction bin for both modes, with `minDetected = 0`, and the number of BH false discoveries at the default settings in the reviewer's setting (1,800 independent sparse null genes).
3. **Power:** mixed fields with known correlation ρ ∈ {0.2, 0.4, 0.6} (`stc_mix` recipe of the calibration fixture) at matched settings for both modes; and the real AKI and brain inputs from the cache: the number of significant genes per mode, the overlap, and the overlap with the published results.
4. **Cost:** wall time of both modes on the full AKI and brain inputs (16 threads).
5. **Recommendation:** which mode should be the default, and whether `minDetected` can then default to a smaller value or 0.

## 5. Tests (lean, default suite ≤ ~10 s extra)

- Each remapped surrogate has exactly the multiset of values of x (use the engine's `keep_surrogates` through the internal interface on a small case).
- For an iid x, the remapped null distribution of r matches the plain permutation null (Kolmogorov–Smirnov p > 0.01 on a few hundred permutations), and the exact-permutation p-value of a sparse gene is reproduced to within Monte Carlo error.
- Thread, batch and chunk invariance in remap mode; the global RNG state unchanged.
- `minDetected = 0` skips only constant genes; the argument is documented as such.
- `surrogate = "gaussian"` reproduces the current results exactly.
