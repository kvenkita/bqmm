#' Yang-Wang-He posterior-variance correction (mixed-model form)
#'
#' The asymmetric Laplace likelihood is a *working* (misspecified) likelihood,
#' so the naive MCMC posterior covariance of the fixed effects is not the
#' asymptotic variance of the quantile-regression estimator. Yang, Wang and He
#' (2016) correct this with a multiplicative sandwich that re-uses the posterior
#' covariance as the "bread":
#'
#'   V_adj = Sigma_post %*% G %*% Sigma_post,
#'
#' where `Sigma_post` is the posterior covariance of the fixed effects and `G`
#' is the meat (variance of the ALD working-likelihood score). For the mixed
#' model this is the right object to correct: `Sigma_post` already encodes the
#' multilevel pooling, so the correction retains the random-effect contribution
#' to fixed-effect uncertainty while fixing the misspecified ALD scale. Under
#' correct specification `G` approximately equals `Sigma_post^{-1}` and the
#' correction reduces to `~Sigma_post`. See [compute_ywh_multiplicative()].
#'
#' @details
#' The pure Koenker-Powell sandwich (`[compute_ywh_sandwich()]`, valid for
#' *fixed-effect* quantile regression) was found by simulation to **under-cover**
#' the fixed effects of a mixed model, because it is computed on residuals with
#' the random effects removed and therefore drops the between-cluster variance.
#' The multiplicative form here covers at or slightly above the nominal level
#' across homoscedastic and heteroscedastic designs (see `tools/bakeoff.R`).
#' Validity is claimed for the fixed-effect block only; variance components keep
#' their model-based posterior.
#'
#' @param object A fitted `bqmm` object.
#' @param meat Meat estimator: `"cluster"` (default; cluster-robust on the first
#'   grouping factor) or `"independence"`.
#' @return A list with the adjusted fixed-effect covariance `vcov`, the naive
#'   posterior covariance `vcov_naive`, the meat `G`, and `sigma`.
#' @keywords internal
ywh_adjust <- function(object, meat = c("cluster", "independence")) {
  meat <- match.arg(meat)
  draws <- get_fixef_draws(object)         # S x K matrix of fixed-effect draws
  betahat <- apply(draws, 2L, stats::median)
  Sigma_post <- stats::cov(draws)

  X <- object$parsed$X
  resid <- bqmm_residuals_median(object)   # y - (X beta + Z b) at medians
  sigma <- stats::median(rstan::extract(object$stanfit, pars = "sigma",
                                        permuted = TRUE)$sigma)
  groups <- if (identical(meat, "cluster")) ywh_group_index(object) else NULL

  ad <- compute_ywh_multiplicative(Sigma_post, X, resid, sigma,
                                   tau = object$tau, groups = groups)

  list(
    beta       = betahat,
    vcov       = ad$vcov,
    vcov_naive = Sigma_post,
    G          = ad$G,
    sigma      = sigma
  )
}

#' Yang-Wang-He multiplicative variance correction
#'
#' Computes `V = Sigma_post %*% G %*% Sigma_post`, the mixed-model form of the
#' YWH correction (see [ywh_adjust()]). The ALD working-likelihood score for
#' observation i is `s_i = (1/sigma) x_i (tau - 1{resid_i < 0})`, so the meat is
#' `G = (1/sigma^2) sum_i x_i x_i' psi_i^2` (independence) or
#' `G = (1/sigma^2) sum_g (sum_{i in g} x_i psi_i)(...)'` (cluster-robust),
#' with `psi_i = tau - 1{resid_i < 0}`.
#'
#' @param Sigma_post Posterior covariance of the fixed effects (K x K).
#' @param X Fixed-effect design matrix (N x K).
#' @param resid Conditional residuals `y - (X beta + Z b)` at the posterior
#'   median fit.
#' @param sigma Posterior ALD scale (a positive scalar).
#' @param tau Quantile level.
#' @param groups Optional integer cluster index (length N). `NULL` gives the
#'   independence meat.
#' @return A symmetric K x K covariance matrix.
#' @keywords internal
compute_ywh_multiplicative <- function(Sigma_post, X, resid, sigma, tau,
                                       groups = NULL) {
  X <- as.matrix(X)
  resid <- as.numeric(resid)
  psi <- tau - (resid < 0)
  if (is.null(groups)) {
    G <- crossprod(X * psi) / sigma^2
  } else {
    s_g <- rowsum(X * psi, group = groups)
    G <- crossprod(s_g) / sigma^2
  }
  V <- Sigma_post %*% G %*% Sigma_post
  V <- (V + t(V)) / 2
  dimnames(V) <- list(colnames(X), colnames(X))
  list(vcov = V, G = G)
}

#' Numeric core of the Yang-Wang-He sandwich
#'
#' Pure-numeric and Stan-free, so it can be unit tested directly.
#'
#' @param X Fixed-effect design matrix (N x K).
#' @param resid Residuals `y - mu` at the posterior median fit.
#' @param tau Quantile level.
#' @param Sigma_post Retained for backward compatibility; ignored. (Earlier
#'   versions used a `Sigma_post`-based fallback for a singular bread, which was
#'   dimensionally incorrect; the bread is now stabilised directly.)
#' @param groups Optional integer cluster index (length N) for a cluster-robust
#'   meat. `NULL` uses the independence meat.
#' @param bandwidth Optional Powell bandwidth; default uses Hall-Sheather.
#' @return List with `vcov`, `D0`, `D1`, `bandwidth` (the bandwidth actually
#'   used, which may have been grown to keep the bread full rank).
#' @keywords internal
compute_ywh_sandwich <- function(X, resid, tau, Sigma_post = NULL,
                                 groups = NULL, bandwidth = NULL) {
  X <- as.matrix(X)
  resid <- as.numeric(resid)
  n <- nrow(X)
  K <- ncol(X)

  if (is.null(bandwidth)) bandwidth <- hall_sheather_bandwidth(n, tau)
  h <- bandwidth

  # bread D1 = (1/n) sum_i f_i(0) x_i x_i' via a Powell uniform kernel. If too
  # few residuals fall in the band to identify D1 (rank-deficient), grow the
  # bandwidth; as a last resort use a smooth normal-density plug-in, which is
  # full rank whenever X is. (The previous Sigma_post fallback was O(1/n^2) and
  # is removed.)
  bread <- function(h) {
    kdens <- (abs(resid) < h) / (2 * h)
    crossprod(X * sqrt(kdens)) / n
  }
  D1 <- bread(h)
  grow <- 0L
  while (qr(D1)$rank < K && grow < 20L) {
    h <- h * 1.5
    D1 <- bread(h)
    grow <- grow + 1L
  }
  if (qr(D1)$rank < K) {
    s <- stats::sd(resid)
    if (!is.finite(s) || s <= 0) s <- 1
    f0 <- mean(stats::dnorm(resid / s)) / s
    D1 <- f0 * crossprod(X) / n
  }

  # meat D0
  if (is.null(groups)) {
    D0 <- tau * (1 - tau) * crossprod(X) / n
  } else {
    psi <- tau - (resid < 0)               # quantile score (carries tau(1-tau))
    s_g <- rowsum(X * psi, group = groups) # cluster score sums
    D0  <- crossprod(s_g) / n
  }

  D1_inv <- tryCatch(
    solve(D1),
    error = function(e) solve(D1 + diag(1e-7 * mean(diag(D1)), K))
  )
  V <- D1_inv %*% D0 %*% D1_inv / n
  V <- (V + t(V)) / 2                       # enforce exact symmetry
  dimnames(V) <- list(colnames(X), colnames(X))

  list(vcov = V, D0 = D0, D1 = D1, bandwidth = h)
}

#' Hall-Sheather bandwidth for quantile sparsity estimation
#'
#' @param n Sample size.
#' @param tau Quantile level.
#' @param alpha Nominal level for the bandwidth (default 0.05).
#' @return A positive bandwidth.
#' @keywords internal
hall_sheather_bandwidth <- function(n, tau, alpha = 0.05) {
  z  <- stats::qnorm(1 - alpha / 2)
  zt <- stats::qnorm(tau)
  num <- 1.5 * stats::dnorm(zt)^2
  den <- 2 * zt^2 + 1
  n^(-1 / 3) * z^(2 / 3) * (num / den)^(1 / 3)
}

# ---- helpers that touch the fitted object (require a compiled stanfit) -------

#' Extract fixed-effect posterior draws
#' @keywords internal
get_fixef_draws <- function(object) {
  ex <- rstan::extract(object$stanfit, pars = "beta", permuted = TRUE)$beta
  ex <- as.matrix(ex)
  colnames(ex) <- object$parsed$fixed_names
  ex
}

#' Residuals at the posterior-median fit (fixed + random part)
#' @keywords internal
bqmm_residuals_median <- function(object) {
  beta <- apply(get_fixef_draws(object), 2L, stats::median)
  mu <- as.numeric(object$parsed$X %*% beta)
  if (ncol(object$parsed$Z) > 0L) {
    bhat <- bqmm_ranef_vector(object)   # handles diagonal and correlated models
    mu <- mu + as.numeric(object$parsed$Z %*% bhat)
  }
  object$parsed$y - mu
}

#' Cluster index for the cluster-robust meat
#'
#' Uses the first (outermost) grouping factor when one is available.
#' @keywords internal
ywh_group_index <- function(object) {
  flist <- object$parsed$flist
  if (is.null(flist) || length(flist) == 0L) return(NULL)
  as.integer(flist[[1L]])
}
