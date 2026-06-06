# Infinitesimal Jackknife (IJ) standard errors

Computes Infinitesimal Jackknife variance estimates for the fixed
effects (Giordano & Broderick, 2023; Ji, Lee & Rabe-Hesketh, 2024) from
a single MCMC run. The IJ influence of observation `i` is
`I_i = n * cov_post(beta, loglik_i)` — the posterior covariance between
the fixed-effect draws and the per-observation (conditional)
log-likelihood draws — and the variance estimate is

## Usage

``` r
ij_vcov(object, cluster = TRUE)
```

## Arguments

- object:

  A fitted `bqmm` object.

- cluster:

  Logical; use the cluster-robust IJ (default `TRUE`, clustering on the
  first grouping factor) or the independence IJ.

## Value

A K x K covariance matrix for the fixed effects.

## Details

V_IJ = (1 / (n (n - 1))) sum_i (I_i - Ibar)(I_i - Ibar)',

with a cluster-robust version that aggregates influences within cluster
`j` as `I_j = (J / n) sum_{i in j} I_i` and replaces `n` by the number
of clusters `J`.

## Caveat for hierarchical models

The influences use the *conditional* per-observation log-likelihood
(given the random effects), whereas the IJ of Ji, Lee & Rabe-Hesketh
(2024) is derived for a *marginal* model. For coefficients identified by
*within*-cluster variation (e.g. slopes) the conditional IJ agrees well
with the Yang-Wang-He sandwich. For coefficients identified by
*between*-cluster variation (the intercept of a random-intercept model)
it **can under-estimate** the variance — up-weighting a cluster mostly
shifts that cluster's random effect rather than the fixed effect — and
it is noisier; how much depends on the random-effect-to-noise ratio. For
valid fixed-effect inference in the hierarchical model fitted by `bqmm`,
prefer the default Yang-Wang-He sandwich
([`ywh_adjust()`](https://kvenkita.github.io/bqmm/reference/ywh_adjust.md));
`method = "ij"` is provided for benchmarking.
