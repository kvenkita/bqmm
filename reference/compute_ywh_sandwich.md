# Numeric core of the Yang-Wang-He sandwich

Pure-numeric and Stan-free, so it can be unit tested directly.

## Usage

``` r
compute_ywh_sandwich(
  X,
  resid,
  tau,
  Sigma_post = NULL,
  groups = NULL,
  bandwidth = NULL
)
```

## Arguments

- X:

  Fixed-effect design matrix (N x K).

- resid:

  Residuals `y - mu` at the posterior median fit.

- tau:

  Quantile level.

- Sigma_post:

  Retained for backward compatibility; ignored. (Earlier versions used a
  `Sigma_post`-based fallback for a singular bread, which was
  dimensionally incorrect; the bread is now stabilised directly.)

- groups:

  Optional integer cluster index (length N) for a cluster-robust meat.
  `NULL` uses the independence meat.

- bandwidth:

  Optional Powell bandwidth; default uses Hall-Sheather.

## Value

List with `vcov`, `D0`, `D1`, `bandwidth` (the bandwidth actually used,
which may have been grown to keep the bread full rank).
