# Bayesian multilevel quantile regression

Fits a Bayesian mixed-effects quantile regression model using the
asymmetric Laplace working likelihood and Stan. The interface follows
`lme4`: random effects are written inline in the formula, e.g.
`y ~ x + (1 + x | group)`, and nested or crossed grouping factors are
both supported.

## Usage

``` r
bqmm(
  formula,
  data,
  tau = 0.5,
  family = ald(),
  prior = NULL,
  cov = c("diagonal", "unstructured"),
  adjust = TRUE,
  prior_only = FALSE,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  cores = getOption("mc.cores", 1L),
  seed = NULL,
  control = list(adapt_delta = 0.95),
  ...
)
```

## Arguments

- formula:

  An lme4-style model formula.

- data:

  A data frame containing the variables in `formula`.

- tau:

  Quantile level(s) in (0, 1). Scalar or vector.

- family:

  A `bqmm_family` object; currently only
  [`ald()`](https://kvenkita.github.io/bqmm/reference/ald.md).

- prior:

  A
  [`bqmm_prior()`](https://kvenkita.github.io/bqmm/reference/bqmm_prior.md)
  object, or `NULL` for data-scaled defaults.

- cov:

  Random-effect covariance structure. `"diagonal"` (default) models
  independent random effects and supports any number of nested or
  crossed terms. `"unstructured"` adds an LKJ-correlated covariance but
  currently requires exactly one random-effects term (e.g.
  `y ~ x + (1 + x | g)`).

- adjust:

  Logical; compute the Yang-Wang-He (2016) variance correction so that
  `vcov(fit, adjusted = TRUE)` returns valid fixed-effect uncertainty.
  Default `TRUE`.

- prior_only:

  Logical; sample from the prior predictive distribution.

- chains, iter, warmup, cores, seed:

  Passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

- control:

  A list of sampler control parameters (e.g. `adapt_delta`). Defaults
  raise `adapt_delta` to 0.95 because ALD posteriors are sharp.

- ...:

  Additional arguments forwarded to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

A `bqmm` object (single `tau`) or a `bqmm_multi` object (vector `tau`).

## Details

One or several quantiles may be requested through `tau`. A scalar
returns a single `bqmm` fit; a vector fits each quantile independently
and returns a `bqmm_multi` container.

## Examples

``` r
# \donttest{
# A minimal fit; raise chains/iter for real analyses.
fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
            tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#> Warning: The largest R-hat is 1.07, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Some Rhat > 1.01; chains may not have converged.
#> Warning: Some effective sample sizes < 100; consider more iterations.
summary(fit)
#> Bayesian multilevel quantile regression (tau = 0.5)
#> Formula: distance ~ age + (1 | Subject)
#> Observations: 108
#> 
#> Fixed effects (Yang-Wang-He adjusted intervals):
#>             Estimate Est.Error   Lower   Upper
#> (Intercept)  17.4954    1.1111 15.3176 19.6733
#> age           0.6004    0.0812  0.4412  0.7596
#> 
#> Random-effect SDs (posterior medians):
#> Subject : (Intercept) 
#>                2.3023 
# }
```
