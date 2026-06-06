# Variance-covariance of the fixed effects

Variance-covariance of the fixed effects

## Usage

``` r
# S3 method for class 'bqmm'
vcov(object, adjusted = TRUE, ...)
```

## Arguments

- object:

  A `bqmm` fit.

- adjusted:

  Logical; if `TRUE` (default) return the Yang-Wang-He corrected
  covariance, otherwise the naive posterior covariance.

- ...:

  Unused.

## Value

A K x K covariance matrix for the fixed effects.
