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

## Estimates and summaries

Extract coefficients, random effects, and model summaries.

- [`summary(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/summary.bqmm.md)
  : Summarize a bqmm fit
- [`coef(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/coef.bqmm.md)
  : Extract model coefficients
- [`fixef(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/fixef.bqmm.md)
  : Posterior-median fixed effects
- [`ranef(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/ranef.bqmm.md)
  : Posterior-median random effects
- [`VarCorr(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/VarCorr.bqmm.md)
  : Random-effect standard deviations and correlations
- [`nobs(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/nobs.bqmm.md)
  : Number of observations used in the fit

## Inference and uncertainty

The Yang-Wang-He variance correction and credible intervals.

- [`vcov(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/vcov.bqmm.md)
  : Variance-covariance of the fixed effects
- [`confint(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/confint.bqmm.md)
  : Confidence (credible) intervals for the fixed effects

## Predictions, methods, and plots

- [`predict(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/predict.bqmm.md)
  : Predictions from a bqmm fit
- [`fitted(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/fitted.bqmm.md)
  : Linear predictor (conditional tau-quantile) at the posterior median
- [`residuals(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/residuals.bqmm.md)
  : Residuals from a bqmm fit
- [`posterior_predict(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/posterior_predict.bqmm.md)
  : Draws from the posterior predictive distribution
- [`posterior_epred(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/posterior_epred.bqmm.md)
  : Draws of the expected response (conditional tau-quantile)
- [`log_lik(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/log_lik.bqmm.md)
  : Pointwise log-likelihood draws
- [`plot(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/plot.bqmm.md)
  : Plot a bqmm fit
- [`as_draws(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/as_draws.bqmm.md)
  : Convert a bqmm fit to a posterior draws object
- [`as.matrix(`*`<bqmm>`*`)`](https://kvenkita.github.io/bqmm/reference/as.matrix.bqmm.md)
  : Coerce a bqmm fit to a matrix of posterior draws

## Multiple quantiles

Methods for joint multi-quantile (`bqmm_multi`) fits.

- [`coef(`*`<bqmm_multi>`*`)`](https://kvenkita.github.io/bqmm/reference/coef.bqmm_multi.md)
  : Coefficient-versus-tau matrix for a bqmm_multi fit
- [`summary(`*`<bqmm_multi>`*`)`](https://kvenkita.github.io/bqmm/reference/summary.bqmm_multi.md)
  : Summarize a bqmm_multi fit
- [`plot(`*`<bqmm_multi>`*`)`](https://kvenkita.github.io/bqmm/reference/plot.bqmm_multi.md)
  : Plot coefficient-versus-tau paths for a bqmm_multi fit

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
