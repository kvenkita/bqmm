# Predictions from a bqmm fit

Predictions from a bqmm fit

## Usage

``` r
# S3 method for class 'bqmm'
predict(
  object,
  newdata = NULL,
  re.form = NULL,
  noncrossing = c("none", "rearrange"),
  ...
)
```

## Arguments

- object:

  A `bqmm` fit.

- newdata:

  Optional data frame; if omitted, training data are used.

- re.form:

  `NULL` includes random effects (training data only); `NA` gives
  population-level predictions.

- noncrossing:

  One of `"none"` or `"rearrange"`. Rearrangement only has an effect for
  `bqmm_multi` objects (multiple quantiles).

- ...:

  Unused.

## Value

Numeric vector of predicted conditional quantiles.
