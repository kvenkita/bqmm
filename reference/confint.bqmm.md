# Confidence (credible) intervals for the fixed effects

Wald-type intervals built from the posterior-median estimates and the
(optionally misspecification-corrected) fixed-effect covariance.

## Usage

``` r
# S3 method for class 'bqmm'
confint(
  object,
  parm,
  level = 0.95,
  adjusted = TRUE,
  method = c("ywh", "ij"),
  cluster = TRUE,
  ...
)
```

## Arguments

- object:

  A `bqmm` fit.

- parm:

  Optional subset of coefficients (names or indices) to return.

- level:

  Interval coverage (default `0.95`).

- adjusted:

  Logical; if `TRUE` (default) use the corrected covariance from
  [`vcov.bqmm()`](https://kvenkita.github.io/bqmm/reference/vcov.bqmm.md),
  otherwise the naive posterior covariance.

- method:

  Correction to use when `adjusted = TRUE`; see
  [`vcov.bqmm()`](https://kvenkita.github.io/bqmm/reference/vcov.bqmm.md).

- cluster:

  Logical; use the cluster-robust form (default `TRUE`).

- ...:

  Unused.

## Value

A matrix with one row per coefficient and lower/upper interval columns.

## Examples

``` r
# \donttest{
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
confint(fit)
#>                   2.5%      97.5%
#> (Intercept) 15.3176323 19.6732556
#> age          0.4412432  0.7596201
# }
```
