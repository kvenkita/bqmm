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
