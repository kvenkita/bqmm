# Prediction / posterior-predictive S3 methods for `bqmm`.

#' Linear predictor (conditional tau-quantile) at the posterior median
#'
#' @param object A `bqmm` fit.
#' @param ... Unused.
#' @return Numeric vector of fitted conditional quantiles.
#' @export
fitted.bqmm <- function(object, ...) {
  beta <- fixef(object)
  mu <- as.numeric(object$parsed$X %*% beta)
  if (ncol(object$parsed$Z) > 0L) {
    bhat <- ranef(object)
    mu <- mu + as.numeric(object$parsed$Z %*% bhat)
  }
  mu
}

#' @export
residuals.bqmm <- function(object, ...) {
  object$parsed$y - fitted(object)
}

#' Predictions from a bqmm fit
#'
#' @param object A `bqmm` fit.
#' @param newdata Optional data frame; if omitted, training data are used.
#' @param re.form `NULL` includes random effects (training data only); `NA`
#'   gives population-level predictions.
#' @param noncrossing One of `"none"` or `"rearrange"`. Rearrangement only has
#'   an effect for `bqmm_multi` objects (multiple quantiles).
#' @param ... Unused.
#' @return Numeric vector of predicted conditional quantiles.
#' @export
predict.bqmm <- function(object, newdata = NULL,
                         re.form = NULL,
                         noncrossing = c("none", "rearrange"), ...) {
  noncrossing <- match.arg(noncrossing)
  if (is.null(newdata)) {
    if (identical(re.form, NA)) {
      return(as.numeric(object$parsed$X %*% fixef(object)))
    }
    return(fitted(object))
  }
  # population-level prediction on new data (random effects set to 0).
  # Drop the (.|.) random-effect bars so grouping variables aren't required.
  fe_form <- lme4::nobars(object$formula)
  tt <- stats::delete.response(stats::terms(fe_form))
  mm <- stats::model.matrix(tt, data = newdata)
  beta <- fixef(object)
  common <- intersect(colnames(mm), names(beta))
  as.numeric(mm[, common, drop = FALSE] %*% beta[common])
}

# Joint draws of the linear predictor (S x N). Extracts beta and b in ONE
# rstan::extract() call so the per-draw correspondence between fixed and random
# effects is preserved (separate permuted=TRUE calls are not guaranteed to share
# an iteration ordering).
bqmm_location_draws <- function(object, extra = character(0)) {
  has_re <- ncol(object$parsed$Z) > 0L
  pars <- c("beta", if (has_re) "b", extra)
  ex <- rstan::extract(object$stanfit, pars = pars, permuted = TRUE)
  beta <- as.matrix(ex$beta)
  loc <- beta %*% t(object$parsed$X)
  if (has_re) {
    bd <- ex$b
    # diagonal model: S x Q. correlated model: S x M x L, flattened to S x Q
    # with column order (level - 1) * M + coef to match Z.
    bQ <- if (length(dim(bd)) == 3L) matrix(bd, nrow = dim(bd)[1L]) else as.matrix(bd)
    loc <- loc + bQ %*% t(object$parsed$Z)
  }
  list(loc = loc, extra = ex[extra])
}

#' @export
posterior_epred.bqmm <- function(object, ...) {
  # the conditional tau-quantile is the ALD location, so epred == location draws
  bqmm_location_draws(object)$loc
}

#' @export
posterior_predict.bqmm <- function(object, ...) {
  d <- bqmm_location_draws(object, extra = "sigma")
  loc <- d$loc
  sigma <- as.numeric(d$extra$sigma)
  tau <- object$tau
  yrep <- matrix(NA_real_, nrow = nrow(loc), ncol = ncol(loc))
  for (s in seq_len(nrow(loc))) {
    yrep[s, ] <- rald(ncol(loc), mu = loc[s, ], sigma = sigma[s], tau = tau)
  }
  yrep
}

#' @export
log_lik.bqmm <- function(object, ...) {
  rstan::extract(object$stanfit, pars = "log_lik", permuted = TRUE)$log_lik
}
