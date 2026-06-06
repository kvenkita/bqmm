# Convert a bqmm fit to a posterior draws object

Convert a bqmm fit to a posterior draws object

## Usage

``` r
# S3 method for class 'bqmm'
as_draws(x, ...)
```

## Arguments

- x:

  A `bqmm` fit.

- ...:

  Unused.

## Value

A `draws_array` (from the `posterior` package) with tidy variable names:
`b_<name>` for fixed effects, `sd_<component>` for random-effect SDs,
and `sigma`.

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
as_draws(fit)
#> # A draws_array: 150 iterations, 1 chains, and 275 variables
#> , , variable = b_(Intercept)
#> 
#>          chain
#> iteration  1
#>         1 17
#>         2 18
#>         3 18
#>         4 17
#>         5 17
#> 
#> , , variable = b_age
#> 
#>          chain
#> iteration    1
#>         1 0.64
#>         2 0.61
#>         3 0.62
#>         4 0.60
#>         5 0.61
#> 
#> , , variable = sigma
#> 
#>          chain
#> iteration    1
#>         1 0.54
#>         2 0.43
#>         3 0.41
#>         4 0.42
#>         5 0.44
#> 
#> , , variable = z_std[1]
#> 
#>          chain
#> iteration     1
#>         1 -0.67
#>         2 -0.73
#>         3 -0.94
#>         4 -0.59
#>         5 -0.53
#> 
#> # ... with 145 more iterations, and 271 more variables
# }
```
