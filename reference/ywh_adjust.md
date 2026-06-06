# Yang-Wang-He posterior-variance correction (mixed-model form)

The asymmetric Laplace likelihood is a *working* (misspecified)
likelihood, so the naive MCMC posterior covariance of the fixed effects
is not the asymptotic variance of the quantile-regression estimator.
Yang, Wang and He (2016) correct this with a multiplicative sandwich
that re-uses the posterior covariance as the "bread":

## Usage

``` r
ywh_adjust(object, meat = c("cluster", "independence"))
```

## Arguments

- object:

  A fitted `bqmm` object.

- meat:

  Meat estimator: `"cluster"` (default; cluster-robust on the first
  grouping factor) or `"independence"`.

## Value

A list with the adjusted fixed-effect covariance `vcov`, the naive
posterior covariance `vcov_naive`, the meat `G`, and `sigma`.

## Details

V_adj = Sigma_post %*% G %*% Sigma_post,

where `Sigma_post` is the posterior covariance of the fixed effects and
`G` is the meat (variance of the ALD working-likelihood score). For the
mixed model this is the right object to correct: `Sigma_post` already
encodes the multilevel pooling, so the correction retains the
random-effect contribution to fixed-effect uncertainty while fixing the
misspecified ALD scale. Under correct specification `G` approximately
equals `Sigma_post^{-1}` and the correction reduces to `~Sigma_post`.
See
[`compute_ywh_multiplicative()`](https://kvenkita.github.io/bqmm/reference/compute_ywh_multiplicative.md).

The pure Koenker-Powell sandwich (`[compute_ywh_sandwich()]`, valid for
*fixed-effect* quantile regression) was found by simulation to
**under-cover** the fixed effects of a mixed model, because it is
computed on residuals with the random effects removed and therefore
drops the between-cluster variance. The multiplicative form here covers
at or slightly above the nominal level across homoscedastic and
heteroscedastic designs (see `tools/bakeoff.R`). Validity is claimed for
the fixed-effect block only; variance components keep their model-based
posterior.
