# The asymmetric Laplace family for quantile regression

A lightweight family object describing the asymmetric Laplace
distribution (ALD) working likelihood used by
[`bqmm()`](https://kvenkita.github.io/bqmm/reference/bqmm.md). It
mirrors the role of a
[stats::family](https://rdrr.io/r/stats/family.html) object but is
intentionally minimal in this release.

## Usage

``` r
ald(link = "identity")
```

## Arguments

- link:

  Name of the link for the location (quantile) parameter. Only
  `"identity"` is supported in v0.1.

## Value

An object of class `"bqmm_family"`.

## Examples

``` r
ald()
#> bqmm family: asymmetric Laplace (link = identity )
```
