# Priors for a Bayesian quantile mixed model

Builds the list of prior hyperparameters passed to Stan. Defaults are
weakly informative and scaled to the data (see
[`bqmm_default_priors()`](https://kvenkita.github.io/bqmm/reference/bqmm_default_priors.md));
any element supplied here overrides the default.

## Usage

``` r
bqmm_prior(
  beta_mean = 0,
  beta_sd = NULL,
  sigma_scale = NULL,
  re_scale = NULL,
  lkj = 2
)
```

## Arguments

- beta_mean:

  Numeric scalar or vector: prior mean(s) for the fixed-effect
  coefficients. Recycled to the number of columns of the design matrix.
  Default `0`.

- beta_sd:

  Positive scalar or vector: prior SD(s) for the fixed-effect
  coefficients. `NULL` (default) uses a data-scaled value.

- sigma_scale:

  Positive scalar: half-normal scale for the ALD scale `sigma`. `NULL`
  (default) uses a data-scaled value.

- re_scale:

  Positive scalar: half-normal scale for the random-effect standard
  deviations. `NULL` (default) uses a data-scaled value.

- lkj:

  Positive scalar: LKJ shape parameter for the random-effect correlation
  matrix (used only when `cov = "unstructured"`). `2` favours weak
  correlations; `1` is uniform over correlation matrices.

## Value

An object of class `"bqmm_prior"`.

## Examples

``` r
bqmm_prior(beta_sd = 5)
#> <bqmm_prior>
#>   beta_mean    0
#>   beta_sd      5
#>   sigma_scale  data-scaled default
#>   re_scale     data-scaled default
#>   lkj          2
```
