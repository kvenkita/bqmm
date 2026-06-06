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

#' @export
coef.bqmm_multi <- function(object, ...) {
  # tau x fixed-coefficient matrix of posterior medians
  mat <- t(vapply(object$fits, function(f) fixef(f), numeric(length(object$parsed$fixed_names))))
  rownames(mat) <- format(object$taus)
  colnames(mat) <- object$parsed$fixed_names
  mat
}

#' @export
summary.bqmm_multi <- function(object, ...) {
  lapply(object$fits, summary, ...)
}

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
