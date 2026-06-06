#' The asymmetric Laplace family for quantile regression
#'
#' A lightweight family object describing the asymmetric Laplace distribution
#' (ALD) working likelihood used by [bqmm()]. It mirrors the role of a
#' [stats::family] object but is intentionally minimal in this release.
#'
#' @param link Name of the link for the location (quantile) parameter.
#'   Only `"identity"` is supported in v0.1.
#'
#' @return An object of class `"bqmm_family"`.
#' @export
#' @examples
#' ald()
ald <- function(link = "identity") {
  link <- match.arg(link, "identity")
  structure(
    list(family = "asymmetric_laplace", link = link),
    class = "bqmm_family"
  )
}

#' @export
print.bqmm_family <- function(x, ...) {
  cat("bqmm family: asymmetric Laplace (link =", x$link, ")\n")
  invisible(x)
}

#' Asymmetric Laplace check (pinball) loss
#'
#' @param u Numeric vector of residuals.
#' @param tau Quantile level in (0, 1).
#' @return Numeric vector of loss values `u * (tau - (u < 0))`.
#' @keywords internal
rho_tau <- function(u, tau) {
  u * (tau - (u < 0))
}

#' Density of the asymmetric Laplace distribution
#'
#' @param x Numeric vector of evaluation points.
#' @param mu Location (the `tau`-quantile).
#' @param sigma Positive scale.
#' @param tau Quantile level in (0, 1).
#' @param log Logical; return the log density.
#' @return Numeric vector of (log) density values.
#' @keywords internal
dald <- function(x, mu = 0, sigma = 1, tau = 0.5, log = FALSE) {
  if (any(sigma <= 0)) stop("`sigma` must be positive.", call. = FALSE)
  if (tau <= 0 || tau >= 1) stop("`tau` must be in (0, 1).", call. = FALSE)
  u <- (x - mu) / sigma
  ld <- log(tau) + log1p(-tau) - log(sigma) - rho_tau(u, tau)
  if (log) ld else exp(ld)
}

#' Random generation from the asymmetric Laplace distribution
#'
#' Uses the normal-exponential location-scale mixture representation of the ALD
#' (Kozumi and Kobayashi, 2011).
#'
#' @param n Number of draws.
#' @inheritParams dald
#' @return Numeric vector of length `n`.
#' @keywords internal
rald <- function(n, mu = 0, sigma = 1, tau = 0.5) {
  if (any(sigma <= 0)) stop("`sigma` must be positive.", call. = FALSE)
  if (tau <= 0 || tau >= 1) stop("`tau` must be in (0, 1).", call. = FALSE)
  theta  <- (1 - 2 * tau) / (tau * (1 - tau))
  tau_sd <- sqrt(2 / (tau * (1 - tau)))
  w <- stats::rexp(n, rate = 1)
  z <- stats::rnorm(n)
  mu + sigma * (theta * w + tau_sd * sqrt(w) * z)
}
