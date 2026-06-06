# Core S3 methods for `bqmm` objects.
# Some methods are fully implemented from the posterior draws; a few are marked
# .NotYetImplemented() and will be completed alongside the compiled model.

#' @export
print.bqmm <- function(x, ...) {
  cat("<bqmm>  Bayesian multilevel quantile regression\n")
  cat("Formula: ", deparse(x$formula), "\n", sep = "")
  cat("Quantile (tau): ", format(x$tau), "\n", sep = "")
  cat("Family: asymmetric Laplace\n")
  cf <- tryCatch(fixef(x), error = function(e) NULL)
  if (!is.null(cf)) {
    cat("\nFixed effects (posterior medians):\n")
    print(round(cf, 4))
  }
  invisible(x)
}

#' Posterior-median fixed effects
#'
#' @param object A `bqmm` fit.
#' @param ... Unused.
#' @return A named numeric vector of posterior-median fixed-effect coefficients.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' fixef(fit)
#' }
#' @export
fixef.bqmm <- function(object, ...) {
  draws <- get_fixef_draws(object)
  stats::setNames(apply(draws, 2L, stats::median), object$parsed$fixed_names)
}

#' Posterior-median random effects as a Q-vector aligned with the columns of Z
#'
#' Handles both the diagonal model (`b` drawn as an S x Q matrix) and the
#' correlated model (`b` drawn as an S x M x L array). For the correlated model
#' the column order `(level - 1) * M + coef` matches the columns of `Z`.
#' @keywords internal
bqmm_ranef_vector <- function(object) {
  if (ncol(object$parsed$Z) == 0L) return(numeric(0))
  b <- rstan::extract(object$stanfit, pars = "b", permuted = TRUE)$b
  if (length(dim(b)) == 3L) {
    as.vector(apply(b, c(2L, 3L), stats::median))     # M x L -> Q (coef fastest)
  } else {
    apply(as.matrix(b), 2L, stats::median)
  }
}

#' Posterior-median random effects
#'
#' @param object A `bqmm` fit.
#' @param ... Unused.
#' @return A numeric vector of posterior-median random effects aligned with the
#'   columns of the random-effects design matrix `Z`, or `NULL` if the model has
#'   no random effects.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' ranef(fit)
#' }
#' @export
ranef.bqmm <- function(object, ...) {
  if (ncol(object$parsed$Z) == 0L) return(NULL)
  bqmm_ranef_vector(object)
}

#' Extract model coefficients
#'
#' Alias for `fixef()`; returns the posterior-median fixed effects.
#'
#' @param object A `bqmm` fit.
#' @param ... Unused.
#' @return A named numeric vector of posterior-median fixed-effect coefficients.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' coef(fit)
#' }
#' @export
coef.bqmm <- function(object, ...) {
  fixef(object)
}

#' Random-effect standard deviations and correlations
#'
#' @param x A `bqmm` fit.
#' @param sigma Ignored; present for compatibility with the generic.
#' @param ... Unused.
#' @return A named numeric vector of posterior-median random-effect standard
#'   deviations (with a posterior-median correlation matrix attached as the
#'   `"correlation"` attribute for unstructured models), or `NULL` if the model
#'   has no random effects.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' VarCorr(fit)
#' }
#' @export
VarCorr.bqmm <- function(x, sigma = 1, ...) {
  if (length(x$parsed$re_components) == 0L) return(NULL)
  sd_draws <- rstan::extract(x$stanfit, pars = "sd_re", permuted = TRUE)$sd_re
  sd_med <- apply(as.matrix(sd_draws), 2L, stats::median)
  out <- stats::setNames(sd_med, x$parsed$re_components)

  # for the correlated model, attach the posterior-median correlation matrix
  if (identical(x$cov, "unstructured")) {
    Corr <- rstan::extract(x$stanfit, pars = "Corr", permuted = TRUE)$Corr
    if (length(dim(Corr)) == 3L) {
      Cmed <- apply(Corr, c(2L, 3L), stats::median)
      coefs <- x$parsed$re_terms[[1L]]$coefs
      dimnames(Cmed) <- list(coefs, coefs)
      attr(out, "correlation") <- Cmed
    }
  }
  out
}

#' Variance-covariance of the fixed effects
#'
#' @param object A `bqmm` fit.
#' @param adjusted Logical; if `TRUE` (default) return a misspecification-
#'   corrected covariance (chosen by `method`), otherwise the naive posterior
#'   covariance.
#' @param method Correction to use when `adjusted = TRUE`: `"ywh"` (default) is
#'   the Yang-Wang-He posterior-covariance sandwich
#'   ([compute_ywh_multiplicative()]); `"ij"` is the Infinitesimal Jackknife
#'   ([ij_vcov()]). Both are cluster-robust by default for a mixed model.
#' @param cluster Logical; use the cluster-robust form (default `TRUE`).
#' @param ... Unused.
#' @return A K x K covariance matrix for the fixed effects.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' vcov(fit)
#' }
#' @export
vcov.bqmm <- function(object, adjusted = TRUE, method = c("ywh", "ij"),
                      cluster = TRUE, ...) {
  if (!isTRUE(adjusted)) {
    return(stats::cov(get_fixef_draws(object)))
  }
  method <- match.arg(method)
  if (method == "ij") {
    return(ij_vcov(object, cluster = cluster))
  }
  # ywh: cache the default (cluster) adjustment; recompute for the independence
  # meat if requested.
  if (isTRUE(cluster)) {
    if (is.null(object$adjustment)) object$adjustment <- ywh_adjust(object)
    return(object$adjustment$vcov)
  }
  ywh_adjust(object, meat = "independence")$vcov
}

#' Confidence (credible) intervals for the fixed effects
#'
#' Wald-type intervals built from the posterior-median estimates and the
#' (optionally misspecification-corrected) fixed-effect covariance.
#'
#' @param object A `bqmm` fit.
#' @param parm Optional subset of coefficients (names or indices) to return.
#' @param level Interval coverage (default `0.95`).
#' @param adjusted Logical; if `TRUE` (default) use the corrected covariance
#'   from [vcov.bqmm()], otherwise the naive posterior covariance.
#' @param method Correction to use when `adjusted = TRUE`; see [vcov.bqmm()].
#' @param cluster Logical; use the cluster-robust form (default `TRUE`).
#' @param ... Unused.
#' @return A matrix with one row per coefficient and lower/upper interval
#'   columns.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' confint(fit)
#' }
#' @export
confint.bqmm <- function(object, parm, level = 0.95, adjusted = TRUE,
                         method = c("ywh", "ij"), cluster = TRUE, ...) {
  method <- match.arg(method)
  beta <- fixef(object)
  V <- vcov(object, adjusted = adjusted, method = method, cluster = cluster)
  se <- sqrt(diag(V))
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(beta - z * se, beta + z * se)
  colnames(ci) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), "%")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' Number of observations used in the fit
#'
#' @param object A `bqmm` fit.
#' @param ... Unused.
#' @return An integer, the number of observations.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' nobs(fit)
#' }
#' @export
nobs.bqmm <- function(object, ...) {
  length(object$parsed$y)
}

#' Summarize a bqmm fit
#'
#' Produces a fixed-effect coefficient table (estimate, standard error and
#' interval) together with random-effect standard deviations.
#'
#' @param object A `bqmm` fit.
#' @param level Interval coverage (default `0.95`).
#' @param adjusted Logical; if `TRUE` (default) use the corrected covariance
#'   from [vcov.bqmm()] for the standard errors and intervals.
#' @param method Correction to use when `adjusted = TRUE`; see [vcov.bqmm()].
#' @param cluster Logical; use the cluster-robust form (default `TRUE`).
#' @param ... Unused.
#' @return An object of class `summary.bqmm`.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' summary(fit)
#' }
#' @export
summary.bqmm <- function(object, level = 0.95, adjusted = TRUE,
                         method = c("ywh", "ij"), cluster = TRUE, ...) {
  method <- match.arg(method)
  beta <- fixef(object)
  V <- vcov(object, adjusted = adjusted, method = method, cluster = cluster)
  se <- sqrt(diag(V))
  z <- stats::qnorm(1 - (1 - level) / 2)
  tab <- cbind(
    Estimate = beta,
    `Est.Error` = se,
    Lower = beta - z * se,
    Upper = beta + z * se
  )
  out <- list(
    formula = object$formula,
    tau = object$tau,
    fixed = tab,
    varcorr = tryCatch(VarCorr.bqmm(object), error = function(e) NULL),
    adjusted = adjusted,
    nobs = nobs(object)
  )
  class(out) <- "summary.bqmm"
  out
}

#' @export
print.summary.bqmm <- function(x, ...) {
  cat("Bayesian multilevel quantile regression (tau = ", format(x$tau), ")\n",
      sep = "")
  cat("Formula: ", deparse(x$formula), "\n", sep = "")
  cat("Observations: ", x$nobs, "\n\n", sep = "")
  cat("Fixed effects",
      if (x$adjusted) " (Yang-Wang-He adjusted intervals):\n" else ":\n",
      sep = "")
  print(round(x$fixed, 4))
  if (!is.null(x$varcorr)) {
    cat("\nRandom-effect SDs (posterior medians):\n")
    print(round(x$varcorr, 4))
  }
  invisible(x)
}
