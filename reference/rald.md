# Random generation from the asymmetric Laplace distribution

Uses the normal-exponential location-scale mixture representation of the
ALD (Kozumi and Kobayashi, 2011).

## Usage

``` r
rald(n, mu = 0, sigma = 1, tau = 0.5)
```

## Arguments

- n:

  Number of draws.

- mu:

  Location (the `tau`-quantile).

- sigma:

  Positive scale.

- tau:

  Quantile level in (0, 1).

## Value

Numeric vector of length `n`.
