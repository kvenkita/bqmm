# Yang-Wang-He multiplicative variance correction

Computes `V = Sigma_post %*% G %*% Sigma_post`, the mixed-model form of
the YWH correction (see
[`ywh_adjust()`](https://kvenkita.github.io/bqmm/reference/ywh_adjust.md)).
The ALD working-likelihood score for observation i is
`s_i = (1/sigma) x_i (tau - 1{resid_i < 0})`, so the meat is
`G = (1/sigma^2) sum_i x_i x_i' psi_i^2` (independence) or
`G = (1/sigma^2) sum_g (sum_{i in g} x_i psi_i)(...)'` (cluster-robust),
with `psi_i = tau - 1{resid_i < 0}`.

## Usage

``` r
compute_ywh_multiplicative(Sigma_post, X, resid, sigma, tau, groups = NULL)
```

## Arguments

- Sigma_post:

  Posterior covariance of the fixed effects (K x K).

- X:

  Fixed-effect design matrix (N x K).

- resid:

  Conditional residuals `y - (X beta + Z b)` at the posterior median
  fit.

- sigma:

  Posterior ALD scale (a positive scalar).

- tau:

  Quantile level.

- groups:

  Optional integer cluster index (length N). `NULL` gives the
  independence meat.

## Value

A symmetric K x K covariance matrix.
