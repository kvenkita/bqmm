#' Construct a bqmm_multi container
#'
#' Holds a list of independent [bqmm()] fits, one per quantile, and presents
#' them jointly through S3 methods (e.g. coefficient-versus-tau paths).
#'
#' @param fits A list of `bqmm` objects.
#' @param parsed The shared parsed formula.
#' @param formula The model formula.
#' @param call The originating call.
#' @return A `bqmm_multi` object.
#' @keywords internal
new_bqmm_multi <- function(fits, parsed, formula, call) {
  structure(
    list(
      fits    = fits,
      taus    = vapply(fits, function(f) f$tau, numeric(1)),
      parsed  = parsed,
      formula = formula,
      call    = call
    ),
    class = "bqmm_multi"
  )
}

#' @export
print.bqmm_multi <- function(x, ...) {
  cat("<bqmm_multi>  Bayesian multilevel quantile regression\n")
  cat("Formula: ", deparse(x$formula), "\n", sep = "")
  cat("Quantiles: ", paste(format(x$taus), collapse = ", "), "\n", sep = "")
  cat("Fixed effects: ", paste(x$parsed$fixed_names, collapse = ", "), "\n",
      sep = "")
  invisible(x)
}

#' Coefficient-versus-tau matrix for a bqmm_multi fit
#'
#' @param object A `bqmm_multi` fit.
#' @param ... Unused.
#' @return A tau-by-coefficient matrix of posterior-median fixed effects, with
#'   one row per quantile.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = c(0.25, 0.75), chains = 1, iter = 250,
#'             refresh = 0, seed = 1)
#' coef(fit)
#' }
#' @export
coef.bqmm_multi <- function(object, ...) {
  # tau x fixed-coefficient matrix of posterior medians
  mat <- t(vapply(object$fits, function(f) fixef(f), numeric(length(object$parsed$fixed_names))))
  rownames(mat) <- format(object$taus)
  colnames(mat) <- object$parsed$fixed_names
  mat
}

#' Summarize a bqmm_multi fit
#'
#' @param object A `bqmm_multi` fit.
#' @param ... Passed to the per-quantile `summary()` method for each fit.
#' @return A list of `summary.bqmm` objects, one per quantile.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = c(0.25, 0.75), chains = 1, iter = 250,
#'             refresh = 0, seed = 1)
#' summary(fit)
#' }
#' @export
summary.bqmm_multi <- function(object, ...) {
  lapply(object$fits, summary, ...)
}

#' Plot coefficient-versus-tau paths for a bqmm_multi fit
#'
#' @param x A `bqmm_multi` fit.
#' @param ... Unused.
#' @return Invisibly, `x`.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = c(0.25, 0.75), chains = 1, iter = 250,
#'             refresh = 0, seed = 1)
#' plot(fit)
#' }
#' @export
plot.bqmm_multi <- function(x, ...) {
  # coefficient-versus-tau paths (base graphics; bayesplot integration later)
  cf <- coef(x)
  taus <- x$taus
  op <- graphics::par(mfrow = grDevices::n2mfrow(ncol(cf)))
  on.exit(graphics::par(op))
  for (j in seq_len(ncol(cf))) {
    graphics::plot(taus, cf[, j], type = "b",
                   xlab = expression(tau), ylab = colnames(cf)[j],
                   main = colnames(cf)[j])
  }
  invisible(x)
}
