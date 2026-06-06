# bqmm

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/kvenkita/bqmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kvenkita/bqmm/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Bayesian Multilevel Quantile Regression

**bqmm** fits Bayesian mixed-effects (multilevel) quantile regression models in
R using the asymmetric Laplace working likelihood and
[Stan](https://mc-stan.org/). It lets you ask how a predictor relates to *any
quantile* of an outcome — the median, the tails, or a whole grid — while
accounting for clustered or repeated-measures data through random effects, and
it returns full Bayesian uncertainty.

The package fills a genuine gap in the R ecosystem. Existing tools are either
frequentist (`lqmm`, `qrLMM`), Bayesian but single-level (`bayesQR`, `Brq`), or
able to fit multilevel quantile models only awkwardly and with statistically
**invalid** uncertainty (`brms`'s `asym_laplace()`). `bqmm` provides a clean,
quantile-first interface **and** valid fixed-effect inference via the
Yang, Wang & He (2016) correction.

📖 **Full documentation, primer, and articles:**
<https://kvenkita.github.io/bqmm/>

## Installation

```r
# install.packages("remotes")
remotes::install_github("kvenkita/bqmm")
```

`bqmm` compiles Stan models on installation, so a working C++ toolchain is
required (Rtools on Windows, the standard compiler chain on macOS/Linux).

## Quick start

```r
library(bqmm)
data(Orthodont, package = "nlme")

# Conditional median of growth, with a random intercept per child
fit <- bqmm(distance ~ age + (1 | Subject), data = Orthodont, tau = 0.5)

summary(fit)              # fixed effects with valid (adjusted) intervals
VarCorr(fit)             # random-effect standard deviations

# Several quantiles in one call
fit_q <- bqmm(distance ~ age + (1 | Subject), data = Orthodont,
              tau = c(0.1, 0.5, 0.9))
plot(fit_q)              # coefficient-versus-quantile paths
predict(fit_q, noncrossing = "rearrange")   # non-crossing quantiles
```

## Key features

- **Familiar `lme4` formula interface** — `y ~ x + (1 + x | group)`; nested
  *and* crossed random effects work out of the box.
- **One or many quantiles** in a single call (a scalar or vector `tau`).
- **Valid inference** — the Yang, Wang & He (2016) posterior-variance correction
  is applied by default, fixing the well-known invalidity of naive
  asymmetric-Laplace credible intervals. Verified by simulation to cover at or
  above the nominal level.
- **Correlated random effects** — `cov = "unstructured"` adds an LKJ-correlated
  random intercept and slope, with the correlation reported by `VarCorr()`.
- **Non-crossing** — optional post-hoc rearrangement so fitted quantiles never
  cross.
- **Ecosystem citizen** — works with the `posterior` and `bayesplot` stacks via
  `as_draws()`, and ships the usual `lme4`/`rstanarm`-style methods.

## Key functions

| Function | Purpose |
|---|---|
| `bqmm()` | Fit a Bayesian multilevel quantile regression model |
| `bqmm_prior()` | Specify priors (fixed effects, scale, random-effect SDs, LKJ) |
| `ald()` | The asymmetric Laplace family object |
| `summary()`, `fixef()`, `coef()` | Fixed-effect estimates and intervals |
| `ranef()`, `VarCorr()` | Random effects and their (co)variances |
| `vcov(fit, adjusted = TRUE)` | Yang–Wang–He–corrected covariance |
| `predict()`, `fitted()` | Fitted / predicted conditional quantiles |
| `posterior_predict()`, `posterior_epred()` | Posterior predictive draws |
| `as_draws()` | Hand the fit to `posterior` / `bayesplot` |
| `rearrange_quantiles()` | Remove quantile crossing |

## Documentation

- **[Primer](https://kvenkita.github.io/bqmm/articles/bqmm-primer.html)** — a
  thorough, accessible guide for applied researchers: the method, fitting,
  priors, inference, diagnostics, and visualization.
- **[Get started](https://kvenkita.github.io/bqmm/articles/bqmm-intro.html)** —
  a short tour.
- **[Multilevel structure](https://kvenkita.github.io/bqmm/articles/bqmm-multilevel.html)**
  — random intercepts/slopes, correlated and crossed effects.
- **[Valid inference](https://kvenkita.github.io/bqmm/articles/bqmm-inference.html)**
  — why naive intervals fail and how the correction works.
- **[Review & testing report](https://kvenkita.github.io/bqmm/articles/bqmm-review.html)**
  — the validation evidence behind the package.
- **[Function reference](https://kvenkita.github.io/bqmm/reference/)**.

## Citation

If you use `bqmm`, please cite it:

> Venkitasubramanian, K. (2026). *bqmm: Bayesian Multilevel Quantile
> Regression*. R package version 0.1.0.
> <https://github.com/kvenkita/bqmm>

```r
citation("bqmm")
```

Please also cite the underlying methodology where appropriate — Yu & Moyeed
(2001) for the asymmetric Laplace approach and Yang, Wang & He (2016) for the
inference correction.

## Author and license

Created and maintained by **Kailas Venkitasubramanian**. Released under the
[MIT License](LICENSE.md).

## References

- Yu, K. & Moyeed, R. A. (2001). Bayesian quantile regression.
  *Statistics & Probability Letters*, 54(3), 437–447.
- Kozumi, H. & Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
  quantile regression. *J. Stat. Comput. Simul.*, 81(11), 1565–1578.
- Geraci, M. & Bottai, M. (2014). Linear quantile mixed models.
  *Statistics and Computing*, 24(3), 461–479.
- Yang, Y., Wang, H. J. & He, X. (2016). Posterior inference in Bayesian
  quantile regression with asymmetric Laplace likelihood.
  *International Statistical Review*, 84(3), 327–344.
- Chernozhukov, V., Fernández-Val, I. & Galichon, A. (2010). Quantile and
  probability curves without crossing. *Econometrica*, 78(3), 1093–1125.
