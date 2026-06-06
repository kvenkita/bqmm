# Numeric core of the Infinitesimal Jackknife variance

Pure-numeric and Stan-free, so it can be unit tested directly.

## Usage

``` r
compute_ij(beta_draws, loglik_draws, groups = NULL)
```

## Arguments

- beta_draws:

  S x K matrix of fixed-effect posterior draws.

- loglik_draws:

  S x n matrix of per-observation log-likelihood draws.

- groups:

  Optional integer cluster index (length n) for the cluster IJ.

## Value

A symmetric K x K covariance matrix.
