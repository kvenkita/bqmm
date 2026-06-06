# Assemble the Stan data list for one quantile

Turns a parsed formula (from
[`bqmm_parse_formula()`](https://kvenkita.github.io/bqmm/reference/bqmm_parse_formula.md))
and resolved priors into the named list consumed by
`inst/stan/bqmm.stan`.

## Usage

``` r
bqmm_standata(parsed, tau, prior, prior_only = FALSE)
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

A named list suitable for
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
