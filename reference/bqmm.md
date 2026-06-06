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
if (FALSE) { # \dontrun{
fit <- bqmm(distance ~ age + (1 | Subject),
            data = nlme::Orthodont, tau = c(0.1, 0.5, 0.9))
summary(fit)
} # }
```
