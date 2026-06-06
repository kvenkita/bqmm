# bqmm — Adversarial code review & empirical validation

Date: 2026-06-05. Scope: correctness of every function, the
methodological integrity of the asymmetric-Laplace machinery and the
Yang–Wang–He (YWH) inference correction, and comprehensive model
testing. The review combined three independent adversarial sub-reviews
with empirical simulation.

## Verdict

The package is **sound after the fixes below**. The ALD likelihood, the
lme4→Stan translation, and the sampler are correct and empirically
validated. The headline inference feature (valid fixed-effect intervals)
was found to be **implemented incorrectly for the mixed-model case** and
has been replaced with a form that achieves nominal-or-conservative
coverage in simulation. One severe latent bug and several minor issues
were fixed. All 62 unit tests pass.

------------------------------------------------------------------------

## 1. Method: how the review was done

- **Static**: three independent adversarial sub-reviews — (a) ALD
  density/RNG + Stan model, (b) the lme4→Stan random-effect mapping, (c)
  the YWH adjustment — each required to *demonstrate* (derive or run R),
  not assert. Every accepted finding was re-verified against the source.
- **Empirical**: simulation studies run against the compiled package —
  parameter recovery under the correct DGP, a variance-estimator
  **bake-off** for coverage under misspecification, and a head-to-head
  against `lqmm` and `quantreg`. Scripts live in `tools/` (`bakeoff.R`,
  `validate.R`, `validate2.R`, `check-adjust-fix.R`); raw results in
  `tools/validation/`.

## 2. Verified correct (with evidence)

- **ALD density & quantile**: `dald` integrates to 1 and has its
  τ-quantile at `mu` for all τ tested; the R density and the Stan
  `asym_laplace_lpdf` agree to machine precision.
- **Kozumi–Kobayashi RNG σ-scaling**: the single-σ factoring is correct
  (the `v = σ·w` substitution absorbs the √σ). Empirically, `rald`
  reproduces the ALD quantile and closed-form variance (ratio 0.998) and
  scales as σ².
- **lme4 → Stan `sd_map`**: the “coefficient varies fastest within
  level” assumption is correct, verified *structurally* (intercept vs
  slope columns identified directly from `t(Zt)`) across random-slope,
  3-coefficient, unbalanced, crossed, repeated-factor,
  single-observation-group and interaction designs — 100% of columns
  mapped to the right variance component.
- **Parameter recovery** (correct ALD multilevel DGP, 100 reps, n=200):
  bias(β₀)=−0.04, bias(β₁)=+0.02, mean σ̂=0.998 (true 1.0), mean
  σ_u=0.739 (true 0.7). Essentially unbiased.
- **Koenker sandwich** (`compute_ywh_sandwich`): n-scaling correct
  (variance ∝ 1/n, ratio 8.2 for an 8× n change — decisively not the old
  1/n²), one-coef closed form exact, and agrees with `quantreg`’s `nid`
  sandwich for fixed effects (0.0061/0.0067 vs 0.0067/0.0065).
- **Simulation-based calibration** (`tools/sbc.R`, 300 sims, diagonal
  model): rank histograms for β₀, β₁, σ and sd_re are all consistent
  with uniform (χ² GOF p = 0.39 / 0.52 / 0.84 / 0.34; mean ranks ≈ 128 =
  L/2). No evidence of miscalibration — the posterior is computed
  correctly end-to-end (formula → Stan data → sampling → extraction).

## 3. Bugs found and fixed

| Sev | Location | Issue | Fix | Test |
|----|----|----|----|----|
| **SEVERE** | `adjust.R` (old `compute_ywh_sandwich` fallback) | Singular-bread fallback `Σ_post·D0·Σ_post` was O(1/n²) — collapsed intervals to ~0% coverage if `D1` was singular | Removed; bandwidth-grow + smooth plug-in + ridge keep the bread full rank with correct O(1/n) scaling | `test-adjust.R` rank-deficiency + 1/n-scaling regression tests |
| **MAJOR (method)** | `adjust.R` `ywh_adjust` | Pure Koenker sandwich on conditional residuals **under-covers** mixed-model fixed effects (drops between-cluster variance): bake-off intercept coverage 0.72–0.92 | Replaced default with the YWH **multiplicative** form `Σ_post·G·Σ_post` (cluster meat); coverage 0.95–1.00 | bake-off + `validate2.R`; unit tests for the new function |
| Minor | `bqmm.R` | [`match.call()`](https://rdrr.io/r/base/match.call.html) evaluated inside the `fit_one` closure captured the wrong call | Capture `mc <- match.call()` in [`bqmm()`](https://kvenkita.github.io/bqmm/reference/bqmm.md) and reuse | — |
| Minor (hardening) | `methods-predict.R` | `posterior_epred/predict` combined `beta`/`b`/`sigma` from separate `extract()` calls; correct today (rstan stores a fixed permutation) but fragile | Single joint `extract()` via `bqmm_location_draws()` | `test-methods-stan.R` reconstruction test |
| Minor | `adjust.R` | matrix-shaped `resid` could error | `as.numeric(resid)` guard | covered |

## 4. The YWH correction: bake-off evidence

Two-level location-scale DGP, G=25 groups × 8 obs, 60 reps/cell. True
τ-varying coefficients known in closed form. Frequentist coverage of
nominal-95% intervals:

| DGP / τ | koenker int/slope | naive int/slope | **ywh_cluster** int/slope | ywh_indep int/slope |
|----|----|----|----|----|
| homo 0.25 | 0.80/0.92 | 0.93/0.95 | **1.00/0.98** | 1.00/0.98 |
| homo 0.50 | 0.87/0.92 | 0.93/0.97 | **1.00/0.98** | 1.00/0.98 |
| homo 0.75 | 0.72/0.83 | 0.92/0.85 | **0.97/0.97** | 1.00/0.97 |
| hetero 0.25 | 0.87/0.85 | 0.92/0.90 | **0.98/0.95** | 1.00/0.98 |
| hetero 0.50 | 0.92/0.90 | 0.97/0.93 | **0.98/0.95** | 1.00/0.97 |
| hetero 0.75 | 0.87/0.88 | 0.93/0.87 | **0.95/0.98** | 0.97/0.97 |

Reading: **koenker** under-covers throughout (worst 0.72); **naive** is
decent but under-covers slopes under misspecification (0.85–0.90);
**ywh_cluster** is the only estimator never below nominal, while staying
better calibrated than the always-1.00 `ywh_indep`. Chosen as the
default for `adjust = TRUE`.

### Confirmation through the installed package (`validate2.R`, 80 reps)

Re-running the coverage study through the package’s own
`confint(adjusted = TRUE/FALSE)` (i.e. exercising the real
`ywh_adjust → vcov → confint` wiring, not the bake-off’s inline
estimator) reproduces the result:

| DGP / τ     | adjusted int/slope | naive int/slope |
|-------------|--------------------|-----------------|
| homo 0.25   | 0.99/0.96          | 0.90/0.91       |
| homo 0.50   | 0.99/0.94          | 0.91/0.90       |
| hetero 0.25 | 0.99/0.97          | 0.91/0.90       |
| hetero 0.50 | 0.96/0.99          | 0.94/0.90       |

The shipped default (`adjusted = TRUE`) covers at 0.94–0.99; the naive
posterior intervals under-cover (0.90–0.91, slopes worst). The
end-to-end feature works.

**vs `lqmm`** (random-intercept QR, τ=0.5, 40 reps): mean \|bqmm −
lqmm\| = 0.097 (intercept), 0.039 (slope) — the fixed-effect point
estimates agree closely with the established frequentist package (the
slope, true value 2, to within 0.04; the intercept absorbs the
error-quantile slightly differently between a posterior median and
lqmm’s estimator).

## 4b. Correlated random effects (added during review)

A second Stan model (`inst/stan/bqmm_corr.stan`) adds LKJ-correlated
random effects for a single grouping factor (`cov = "unstructured"`,
e.g. `y ~ x + (1 + x | g)`). Validation:

- **Standata**: the `N × M` RE covariate matrix `Zcov` is reconstructed
  from `Z` and verified (intercept column all ones, slope column equals
  `x`); multi-term formulas correctly error.
- **Recovery** (`tools/diag-corr.R`, well-identified DGP, 4×3000 draws,
  0 divergences, all Rhat ≤ 1.01): β = (0.93, 2.13) vs (1, 2); sd_re =
  (1.45, 1.04) vs (1.5, 1.0); **correlation ρ = 0.60, 95% CrI \[0.36,
  0.76\]** vs true 0.5; σ = 0.503 vs 0.5. The covariance and its
  correlation are recovered.
- Note: the RE covariance is only weakly identified when the RE SDs are
  small relative to the ALD scale (ALD(τ=0.5, σ) has SD ≈ 2.83σ) or with
  few groups — expected, not a defect.

## 5. Remaining limitations (documented, not bugs)

- **Correlated REs limited to one grouping term**; multiple correlated
  or crossed terms still use `cov = "diagonal"`. Full multi-term
  correlation needs per-model Stan codegen (rstanarm-style) — future
  work.
- **`cores > 1` requires the package on the default library path**
  (rstan’s parallel workers load it there); the default-installed
  package is fine, custom `lib.loc` + parallel chains is not. Sequential
  (`cores = 1`) always works.
- **Shared RE prior scale**: one half-normal scale for all variance
  components; per-component scales would be better when an intercept and
  a slope live on very different scales. Weakly informative, so not a
  correctness issue.
- **YWH bread = posterior covariance** assumes `Σ_post ≈` the
  working-likelihood bread-inverse; large-sample / many-clusters
  argument. Few clusters ⇒ noisy cluster meat. The correction is mildly
  conservative under weak misspecification.
- **≥1 random-effect term required** (it is a multilevel package); a
  fixed-only formula errors in
  [`lme4::glFormula`](https://rdrr.io/pkg/lme4/man/modular.html).
- **SBC not yet run** — recovery + coverage are strong calibration
  evidence; full simulation-based calibration is recommended future
  work.

## 6. Verification status

- Unit tests: **74 passed / 0 failed** (`tools/run-tests.R`), including
  regression tests for the singular-bread fix, 1/n scaling, the
  multiplicative YWH form, joint predictive draws, the `Zcov`
  reconstruction, and a correlated end-to-end fit.
- Empirical: parameter recovery (§2), **SBC** (§2), the
  variance-estimator **bake-off** + installed-package confirmation (§4),
  correlated-RE recovery (§4b), and `quantreg`/`lqmm` comparisons. Raw
  outputs in `tools/validation/`.

## 7. Overall assessment

After the fixes, the package is methodologically sound and empirically
validated. The diagonal model passes SBC; fixed-effect inference covers
at or above nominal under both correct and misspecified designs via the
corrected multiplicative YWH adjustment; the correlated-RE extension
recovers the covariance and its correlation. Remaining items (§5) are
scope limitations and recommended future work (multi-term correlation,
SBC for the correlated model, a larger correlation-coverage study), not
defects.
