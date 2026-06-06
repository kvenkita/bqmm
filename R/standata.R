#' Assemble the Stan data list for one quantile
#'
#' Turns a parsed formula (from [bqmm_parse_formula()]) and resolved priors into
#' the named list consumed by `inst/stan/bqmm.stan`.
#'
#' @param parsed Output of [bqmm_parse_formula()].
#' @param tau Quantile level in (0, 1).
#' @param prior A fully-resolved [bqmm_prior()] (see [bqmm_default_priors()]).
#' @param prior_only Logical; sample from the prior predictive only.
#' @return A named list suitable for `rstan::sampling()`.
#' @keywords internal
bqmm_standata <- function(parsed, tau, prior, prior_only = FALSE) {
  stopifnot(tau > 0, tau < 1)

  X <- parsed$X
  Z <- parsed$Z
  N <- nrow(X)
  K <- ncol(X)
  Q <- if (is.null(Z)) 0L else ncol(Z)
  G <- length(parsed$re_components)

  if (Q == 0L) {
    Z <- matrix(0, nrow = N, ncol = 0)
  }

  list(
    N = N,
    K = K,
    X = X,
    y = parsed$y,
    p = tau,

    Q = Q,
    Z = Z,
    G = G,
    sd_map = if (Q > 0L) as.integer(parsed$sd_map) else integer(0),

    prior_beta_mean   = as.numeric(prior$beta_mean),
    prior_beta_sd     = as.numeric(prior$beta_sd),
    prior_sigma_scale = as.numeric(prior$sigma_scale),
    prior_re_scale    = as.numeric(prior$re_scale),

    prior_only = as.integer(isTRUE(prior_only))
  )
}

#' Assemble the Stan data list for the correlated-RE model
#'
#' Builds the data for `inst/stan/bqmm_corr.stan` (a single grouping factor with
#' `M` correlated coefficients). Requires the parsed formula to have exactly one
#' random-effects term.
#'
#' @inheritParams bqmm_standata
#' @return A named list for `rstan::sampling()` with the correlated-RE model.
#' @keywords internal
bqmm_corr_standata <- function(parsed, tau, prior, prior_only = FALSE) {
  stopifnot(tau > 0, tau < 1)
  if (length(parsed$re_terms) != 1L) {
    stop("cov = \"unstructured\" currently supports exactly one random-effects ",
         "term; use cov = \"diagonal\" for multiple or crossed terms.",
         call. = FALSE)
  }
  term <- parsed$re_terms[[1L]]
  M <- length(term$coefs)
  L <- term$n_levels
  level_id <- as.integer(parsed$flist[[1L]])
  if (length(level_id) != nrow(parsed$X)) {
    stop("Internal error: grouping factor length does not match data.",
         call. = FALSE)
  }

  # Reconstruct the N x M random-effect covariate matrix from Z. With lme4's
  # coefficient-fastest-within-level column ordering, the M columns for level l
  # occupy positions ((l-1)*M + 1):(l*M); each observation's nonzero entries in
  # its level block are exactly its RE covariate row.
  Z <- parsed$Z
  Zcov <- matrix(0, nrow = nrow(Z), ncol = M)
  for (i in seq_len(nrow(Z))) {
    cols <- ((level_id[i] - 1L) * M + 1L):(level_id[i] * M)
    Zcov[i, ] <- Z[i, cols]
  }

  X <- parsed$X
  list(
    N = nrow(X), K = ncol(X), X = X, y = parsed$y, p = tau,
    M = M, L = L, level_id = level_id, Zcov = Zcov,
    prior_beta_mean   = as.numeric(prior$beta_mean),
    prior_beta_sd     = as.numeric(prior$beta_sd),
    prior_sigma_scale = as.numeric(prior$sigma_scale),
    prior_re_scale    = as.numeric(prior$re_scale),
    prior_lkj         = as.numeric(prior$lkj),
    prior_only        = as.integer(isTRUE(prior_only))
  )
}
