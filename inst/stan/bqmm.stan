// bqmm: Bayesian multilevel quantile regression via the asymmetric Laplace
// working likelihood.
//
// This is the v0.1 model. It supports:
//   * arbitrary fixed effects (design matrix X),
//   * arbitrary random effects with nested AND crossed grouping, encoded as a
//     dense design matrix Z plus a mapping from each random-effect column to a
//     variance component (built on the R side from lme4::mkReTrms),
//   * a single quantile `p`, fit with a non-centred parameterisation.
//
// Random effects are modelled as DIAGONAL (independent within a grouping
// factor), which matches lqmm's default `pdDiag` covariance. Correlated random
// effects (LKJ on a Cholesky correlation factor, per term) are the next step
// and require per-model code generation; see SETUP.md / the roadmap.

functions {
#include /include/ald.stanfunctions
}

data {
  int<lower=1> N;                 // number of observations
  int<lower=0> K;                 // number of fixed-effect columns
  matrix[N, K] X;                 // fixed-effect design matrix
  vector[N] y;                    // response

  real<lower=0, upper=1> p;       // target quantile

  int<lower=0> Q;                 // number of random-effect columns (ncol(Z))
  matrix[N, Q] Z;                 // random-effect design matrix (dense)
  int<lower=0> G;                 // number of variance components
  array[Q] int<lower=1, upper=(G > 0 ? G : 1)> sd_map;  // RE column -> component

  // prior hyperparameters (set with sensible scaled defaults on the R side)
  vector[K] prior_beta_mean;
  vector<lower=0>[K] prior_beta_sd;
  real<lower=0> prior_sigma_scale;   // half-normal scale for sigma
  real<lower=0> prior_re_scale;      // half-normal scale for RE sds

  int<lower=0, upper=1> prior_only;  // 1 = sample from the prior predictive
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;
  vector[Q] z_std;                   // standardised random effects
  vector<lower=0>[G > 0 ? G : 0] sd_re;
}

transformed parameters {
  vector[Q] b;                       // random effects on the data scale
  for (q in 1:Q) {
    b[q] = z_std[q] * sd_re[sd_map[q]];
  }
}

model {
  // priors
  beta ~ normal(prior_beta_mean, prior_beta_sd);
  sigma ~ normal(0, prior_sigma_scale);          // half-normal via <lower=0>
  z_std ~ std_normal();
  if (G > 0) {
    sd_re ~ normal(0, prior_re_scale);           // half-normal via <lower=0>
  }

  // likelihood
  if (prior_only == 0) {
    vector[N] mu = X * beta;
    if (Q > 0) {
      mu += Z * b;
    }
    y ~ asym_laplace(mu, sigma, p);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  {
    vector[N] mu = X * beta;
    if (Q > 0) {
      mu += Z * b;
    }
    for (i in 1:N) {
      log_lik[i] = asym_laplace_lpdf_scalar(y[i], mu[i], sigma, p);
      y_rep[i]   = asym_laplace_rng(mu[i], sigma, p);
    }
  }
}
