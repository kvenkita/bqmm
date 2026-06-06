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

#' @export
ranef.bqmm <- function(object, ...) {
  if (ncol(object$parsed$Z) == 0L) return(NULL)
  bqmm_ranef_vector(object)
}

#' @export
coef.bqmm <- function(object, ...) {
  fixef(object)
}

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
#' @param adjusted Logical; if `TRUE` (default) return the Yang-Wang-He
#'   corrected covariance, otherwise the naive posterior covariance.
#' @param ... Unused.
#' @return A K x K covariance matrix for the fixed effects.
#' @export
vcov.bqmm <- function(object, adjusted = TRUE, ...) {
  if (isTRUE(adjusted)) {
    if (is.null(object$adjustment)) object$adjustment <- ywh_adjust(object)
    return(object$adjustment$vcov)
  }
  stats::cov(get_fixef_draws(object))
}

#' @export
confint.bqmm <- function(object, parm, level = 0.95, adjusted = TRUE, ...) {
  beta <- fixef(object)
  V <- vcov(object, adjusted = adjusted)
  se <- sqrt(diag(V))
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(beta - z * se, beta + z * se)
  colnames(ci) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), "%")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
nobs.bqmm <- function(object, ...) {
  length(object$parsed$y)
}

#' @export
summary.bqmm <- function(object, level = 0.95, adjusted = TRUE, ...) {
  beta <- fixef(object)
  V <- vcov(object, adjusted = adjusted)
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
