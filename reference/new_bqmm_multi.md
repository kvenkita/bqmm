# Construct a bqmm_multi container

Holds a list of independent
[`bqmm()`](https://kvenkita.github.io/bqmm/reference/bqmm.md) fits, one
per quantile, and presents them jointly through S3 methods (e.g.
coefficient-versus-tau paths).

## Usage

``` r
new_bqmm_multi(fits, parsed, formula, call)
```

## Arguments

- fits:

  A list of `bqmm` objects.

- parsed:

  The shared parsed formula.

- formula:

  The model formula.

- call:

  The originating call.

## Value

A `bqmm_multi` object.
