# Plot a bqmm fit

Default plot shows fixed-effect posterior intervals. With `bayesplot`
installed, richer MCMC plots are available via
[`as_draws()`](https://mc-stan.org/posterior/reference/draws.html).

## Usage

``` r
# S3 method for class 'bqmm'
plot(x, ...)
```

## Arguments

- x:

  A `bqmm` fit.

- ...:

  Unused.

## Value

Invisibly, `x`.
