# Summarize a bqmm_multi fit

Summarize a bqmm_multi fit

## Usage

``` r
# S3 method for class 'bqmm_multi'
summary(object, ...)
```

## Arguments

- object:

  A `bqmm_multi` fit.

- ...:

  Passed to the per-quantile
  [`summary()`](https://rdrr.io/r/base/summary.html) method for each
  fit.

## Value

A list of `summary.bqmm` objects, one per quantile.

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
summary(fit)
#> $`0.25`
#> Bayesian multilevel quantile regression (tau = 0.25)
#> Formula: distance ~ age + (1 | Subject)
#> Observations: 108
#> 
#> Fixed effects (Yang-Wang-He adjusted intervals):
#>             Estimate Est.Error   Lower   Upper
#> (Intercept)  16.7809    1.3130 14.2075 19.3543
#> age           0.6034    0.0967  0.4139  0.7928
#> 
#> Random-effect SDs (posterior medians):
#> Subject : (Intercept) 
#>                2.1742 
#> 
#> $`0.75`
#> Bayesian multilevel quantile regression (tau = 0.75)
#> Formula: distance ~ age + (1 | Subject)
#> Observations: 108
#> 
#> Fixed effects (Yang-Wang-He adjusted intervals):
#>             Estimate Est.Error   Lower   Upper
#> (Intercept)  17.5302    1.0126 15.5456 19.5148
#> age           0.6439    0.0744  0.4981  0.7898
#> 
#> Random-effect SDs (posterior medians):
#> Subject : (Intercept) 
#>                2.3617 
#> 
# }
```
