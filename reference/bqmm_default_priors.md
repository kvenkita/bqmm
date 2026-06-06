# Resolve data-scaled default priors

Fills in any `NULL` elements of a
[`bqmm_prior()`](https://kvenkita.github.io/bqmm/reference/bqmm_prior.md)
using simple, robust scales derived from the response and design matrix.
The intent is to keep `sigma` identified and the sampler away from
divergences, not to be informative.

## Usage

``` r
bqmm_default_priors(prior, y, K)
```

## Arguments

- prior:

  A `bqmm_prior` (or `NULL` for all defaults).

- y:

  Numeric response vector.

- K:

  Number of fixed-effect columns.

## Value

A fully-specified `bqmm_prior` with numeric hyperparameters.
