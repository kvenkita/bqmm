# Package index

## Model fitting

The main entry point.

- [`bqmm()`](https://kvenkita.github.io/bqmm/reference/bqmm.md) :
  Bayesian multilevel quantile regression

## Priors and family

Prior specification and the asymmetric Laplace family.

- [`bqmm_prior()`](https://kvenkita.github.io/bqmm/reference/bqmm_prior.md)
  : Priors for a Bayesian quantile mixed model
- [`ald()`](https://kvenkita.github.io/bqmm/reference/ald.md) : The
  asymmetric Laplace family for quantile regression

## Inference and uncertainty

The Yang-Wang-He variance correction.

- [`vcov(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/vcov.bqmm.md)
  : Variance-covariance of the fixed effects

## Predictions, methods, and plots

- [`predict(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/predict.bqmm.md)
  : Predictions from a bqmm fit
- [`fitted(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/fitted.bqmm.md)
  : Linear predictor (conditional tau-quantile) at the posterior median
- [`plot(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/plot.bqmm.md)
  : Plot a bqmm fit
- [`as_draws(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/as_draws.bqmm.md)
  : Convert a bqmm fit to a posterior draws object

## Non-crossing

- [`rearrange_quantiles()`](https://kvenkita.github.io/bqmm/reference/rearrange_quantiles.md)
  : Rearrange fitted quantiles to remove crossing

## Package and re-exports

- [`bqmm-package`](https://kvenkita.github.io/bqmm/reference/bqmm-package.md)
  : bqmm: Bayesian Multilevel Quantile Regression
- [`reexports`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`fixef`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`ranef`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`VarCorr`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`as_draws`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`posterior_predict`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`posterior_epred`](https://kvenkita.github.io/bqmm/reference/reexports.md)
  [`log_lik`](https://kvenkita.github.io/bqmm/reference/reexports.md) :
  Objects exported from other packages
