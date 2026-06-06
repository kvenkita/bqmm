# Assemble the Stan data list for the correlated-RE model

Builds the data for `inst/stan/bqmm_corr.stan` (a single grouping factor
with `M` correlated coefficients). Requires the parsed formula to have
exactly one random-effects term.

## Usage

``` r
bqmm_corr_standata(parsed, tau, prior, prior_only = FALSE)
```

## Arguments

- parsed:

  Output of
  [`bqmm_parse_formula()`](https://kvenkita.github.io/bqmm/reference/bqmm_parse_formula.md).

- tau:

  Quantile level in (0, 1).

- prior:

  A fully-resolved
  [`bqmm_prior()`](https://kvenkita.github.io/bqmm/reference/bqmm_prior.md)
  (see
  [`bqmm_default_priors()`](https://kvenkita.github.io/bqmm/reference/bqmm_default_priors.md)).

- prior_only:

  Logical; sample from the prior predictive only.

## Value

A named list for
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
with the correlated-RE model.
