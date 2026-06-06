// bqmm correlated-random-effects model: a SINGLE grouping factor with M
// correlated random-effect coefficients per level, using an LKJ prior on the
// correlation and a non-centred parameterisation. This complements bqmm.stan
// (which models diagonal/independent random effects and supports any number of
// nested/crossed terms). Selected via cov = "unstructured" when the formula has
// exactly one random-effects term.

functions {
#include /include/ald.stanfunctions
}

data {
  int<lower=1> N;                 // observations
  int<lower=0> K;                 // fixed-effect columns
  matrix[N, K] X;                 // fixed-effect design
  vector[N] y;                    // response
  real<lower=0, upper=1> p;       // quantile

  int<lower=1> M;                 // random-effect coefficients per level
  int<lower=1> L;                 // number of levels (groups)
  array[N] int<lower=1, upper=L> level_id;   // level of each observation
  matrix[N, M] Zcov;              // random-effect covariate rows

  vector[K] prior_beta_mean;
  vector<lower=0>[K] prior_beta_sd;
  real<lower=0> prior_sigma_scale;
  real<lower=0> prior_re_scale;
  real<lower=0> prior_lkj;        // LKJ shape (eta)

  int<lower=0, upper=1> prior_only;
}

parameters {
  vector[K] beta;
  real<lower=0> sigma;
  matrix[M, L] z_std;                 // standardised random effects
  vector<lower=0>[M] sd_re;           // random-effect SDs
  cholesky_factor_corr[M] L_corr;     // correlation Cholesky factor
}

transformed parameters {
  // non-centred: columns are level-specific RE vectors
  matrix[M, L] b = diag_pre_multiply(sd_re, L_corr) * z_std;
}

model {
  beta ~ normal(prior_beta_mean, prior_beta_sd);
  sigma ~ normal(0, prior_sigma_scale);          // half-normal via <lower=0>
  to_vector(z_std) ~ std_normal();
  sd_re ~ normal(0, prior_re_scale);             // half-normal via <lower=0>
  L_corr ~ lkj_corr_cholesky(prior_lkj);

  if (prior_only == 0) {
    vector[N] mu = X * beta;
    for (i in 1:N) {
      mu[i] += Zcov[i] * b[, level_id[i]];       // row_vector * vector = real
    }
    y ~ asym_laplace(mu, sigma, p);
  }
}

generated quantities {
  matrix[M, M] Corr = multiply_lower_tri_self_transpose(L_corr);  // correlation
  vector[N] log_lik;
  vector[N] y_rep;
  {
    vector[N] mu = X * beta;
    for (i in 1:N) {
      mu[i] += Zcov[i] * b[, level_id[i]];
    }
    for (i in 1:N) {
      log_lik[i] = asym_laplace_lpdf_scalar(y[i], mu[i], sigma, p);
      y_rep[i]   = asym_laplace_rng(mu[i], sigma, p);
    }
  }
}
