#' Infinitesimal Jackknife (IJ) standard errors
#'
#' Computes Infinitesimal Jackknife variance estimates for the fixed effects
#' (Giordano & Broderick, 2023; Ji, Lee & Rabe-Hesketh, 2024) from a single
#' MCMC run. The IJ influence of observation `i` is
#' `I_i = n * cov_post(beta, loglik_i)` — the posterior covariance between the
#' fixed-effect draws and the per-observation (conditional) log-likelihood draws
#' — and the variance estimate is
#'
#'   V_IJ = (1 / (n (n - 1))) sum_i (I_i - Ibar)(I_i - Ibar)',
#'
#' with a cluster-robust version that aggregates influences within cluster `j`
#' as `I_j = (J / n) sum_{i in j} I_i` and replaces `n` by the number of clusters
#' `J`.
#'
#' @section Caveat for hierarchical models:
#' The influences use the *conditional* per-observation log-likelihood (given the
#' random effects), whereas the IJ of Ji, Lee & Rabe-Hesketh (2024) is derived
#' for a *marginal* model. For coefficients identified by *within*-cluster
#' variation (e.g. slopes) the conditional IJ agrees well with the Yang-Wang-He
#' sandwich. For coefficients identified by *between*-cluster variation (the
#' intercept of a random-intercept model) it **can under-estimate** the variance
#' — up-weighting a cluster mostly shifts that cluster's random effect rather
#' than the fixed effect — and it is noisier; how much depends on the
#' random-effect-to-noise ratio. For valid fixed-effect inference in the
#' hierarchical model fitted by `bqmm`, prefer the default Yang-Wang-He sandwich
#' ([ywh_adjust()]); `method = "ij"` is provided for benchmarking.
#'
#' @param object A fitted `bqmm` object.
#' @param cluster Logical; use the cluster-robust IJ (default `TRUE`, clustering
#'   on the first grouping factor) or the independence IJ.
#' @return A K x K covariance matrix for the fixed effects.
#' @keywords internal
ij_vcov <- function(object, cluster = TRUE) {
  beta_draws <- get_fixef_draws(object)                       # S x K
  loglik <- rstan::extract(object$stanfit, pars = "log_lik",
                           permuted = TRUE)$log_lik           # S x n
  groups <- if (isTRUE(cluster)) ij_group_index(object) else NULL
  V <- compute_ij(beta_draws, loglik, groups = groups)
  dimnames(V) <- list(object$parsed$fixed_names, object$parsed$fixed_names)
  V
}

#' Numeric core of the Infinitesimal Jackknife variance
#'
#' Pure-numeric and Stan-free, so it can be unit tested directly.
#'
#' @param beta_draws S x K matrix of fixed-effect posterior draws.
#' @param loglik_draws S x n matrix of per-observation log-likelihood draws.
#' @param groups Optional integer cluster index (length n) for the cluster IJ.
#' @return A symmetric K x K covariance matrix.
#' @keywords internal
compute_ij <- function(beta_draws, loglik_draws, groups = NULL) {
  beta_draws <- as.matrix(beta_draws)
  loglik_draws <- as.matrix(loglik_draws)
  S <- nrow(beta_draws)
  n <- ncol(loglik_draws)
  if (nrow(loglik_draws) != S) {
    stop("beta_draws and loglik_draws must have the same number of draws.",
         call. = FALSE)
  }

  # posterior covariance C[k, i] = cov(beta[, k], loglik[, i]); influences I_i.
  bc <- sweep(beta_draws, 2L, colMeans(beta_draws))          # S x K, centred
  lc <- sweep(loglik_draws, 2L, colMeans(loglik_draws))      # S x n, centred
  C <- crossprod(bc, lc) / (S - 1)                           # K x n
  infl <- n * C                                              # K x n: columns I_i

  if (is.null(groups)) {
    Ibar <- rowMeans(infl)
    Ic <- infl - Ibar
    V <- tcrossprod(Ic) / (n * (n - 1))
  } else {
    groups <- as.integer(groups)
    J <- length(unique(groups))
    agg <- t(rowsum(t(infl), group = groups))                # K x J: cluster sums
    infl_cl <- (J / n) * agg
    Ibar_cl <- rowMeans(infl_cl)
    Icc <- infl_cl - Ibar_cl
    V <- tcrossprod(Icc) / (J * (J - 1))
  }
  V <- (V + t(V)) / 2                                         # exact symmetry
  V
}

#' Cluster index for the cluster-robust IJ
#' @keywords internal
ij_group_index <- function(object) {
  flist <- object$parsed$flist
  if (is.null(flist) || length(flist) == 0L) return(NULL)
  as.integer(flist[[1L]])
}
