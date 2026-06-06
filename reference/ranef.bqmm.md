# Posterior-median random effects

Posterior-median random effects

## Usage

``` r
# S3 method for class 'bqmm'
ranef(object, ...)
```

## Arguments

- object:

  A `bqmm` fit.

- ...:

  Unused.

## Value

A numeric vector of posterior-median random effects aligned with the
columns of the random-effects design matrix `Z`, or `NULL` if the model
has no random effects.

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
ranef(fit)
#>  [1] -1.07849003 -0.98944964 -0.80020631 -0.58516679 -0.22629122 -0.28227883
#>  [7] -0.07279188 -0.07450183  1.02625922  0.48023175  0.46550229  1.24757566
#> [13]  2.03819092  2.28031929  3.75275728  5.18810027 -5.40758736 -2.57188205
#> [19] -2.67165211 -2.86997744 -1.19461286 -1.00111012 -0.98160755 -0.61369095
#> [25] -0.09422922  0.85262409  2.11130471
# }
```
