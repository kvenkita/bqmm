# Plotting methods for `bqmm`.

#' Plot a bqmm fit
#'
#' Default plot shows fixed-effect posterior intervals. With `bayesplot`
#' installed, richer MCMC plots are available via [as_draws()].
#'
#' @param x A `bqmm` fit.
#' @param ... Unused.
#' @return Invisibly, `x`.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' plot(fit)
#' }
#' @export
plot.bqmm <- function(x, ...) {
  beta <- fixef(x)
  V <- vcov(x, adjusted = TRUE)
  se <- sqrt(diag(V))
  idx <- seq_along(beta)
  graphics::plot(beta, idx, yaxt = "n",
                 xlim = range(c(beta - 2 * se, beta + 2 * se)),
                 xlab = "Estimate", ylab = "",
                 main = sprintf("Fixed effects (tau = %s)", format(x$tau)))
  graphics::axis(2, at = idx, labels = names(beta), las = 1)
  graphics::segments(beta - 2 * se, idx, beta + 2 * se, idx)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  invisible(x)
}
