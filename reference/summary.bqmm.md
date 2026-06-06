# Summarize a bqmm fit

Produces a fixed-effect coefficient table (estimate, standard error and
interval) together with random-effect standard deviations.

## Usage

``` r
# S3 method for class 'bqmm'
summary(
  object,
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

- level:

  Interval coverage (default `0.95`).

- adjusted:

  Logical; if `TRUE` (default) use the corrected covariance from
  [`vcov.bqmm()`](https://kvenkita.github.io/bqmm/reference/vcov.bqmm.md)
  for the standard errors and intervals.

- method:

  Correction to use when `adjusted = TRUE`; see
  [`vcov.bqmm()`](https://kvenkita.github.io/bqmm/reference/vcov.bqmm.md).

- cluster:

  Logical; use the cluster-robust form (default `TRUE`).

- ...:

  Unused.

## Value

An object of class `summary.bqmm`.

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
