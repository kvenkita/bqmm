# Asymmetric Laplace check (pinball) loss

Asymmetric Laplace check (pinball) loss

## Usage

``` r
rho_tau(u, tau)
```

## Arguments

- u:

  Numeric vector of residuals.

- tau:

  Quantile level in (0, 1).

## Value

Numeric vector of loss values `u * (tau - (u < 0))`.
