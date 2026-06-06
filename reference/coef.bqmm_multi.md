# Coefficient-versus-tau matrix for a bqmm_multi fit

Coefficient-versus-tau matrix for a bqmm_multi fit

## Usage

``` r
# S3 method for class 'bqmm_multi'
coef(object, ...)
```

## Arguments

- object:

  A `bqmm_multi` fit.

- ...:

  Unused.

## Value

A tau-by-coefficient matrix of posterior-median fixed effects, with one
row per quantile.

## Examples

``` r
# \donttest{
fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
            tau = c(0.25, 0.75), chains = 1, iter = 250,
            refresh = 0, seed = 1)
#> Warning: There were 113 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.19, indicating chains have not mixed.
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
#> Warning: There were 29 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.09, indicating chains have not mixed.
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
coef(fit)
#>      (Intercept)       age
#> 0.25    16.78089 0.6033671
#> 0.75    17.53018 0.6439472
# }
```
