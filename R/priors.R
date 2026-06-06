#' Priors for a Bayesian quantile mixed model
#'
#' Builds the list of prior hyperparameters passed to Stan. Defaults are weakly
#' informative and scaled to the data (see [bqmm_default_priors()]); any element
#' supplied here overrides the default.
#'
#' @param beta_mean Numeric scalar or vector: prior mean(s) for the fixed-effect
#'   coefficients. Recycled to the number of columns of the design matrix.
#'   Default `0`.
#' @param beta_sd Positive scalar or vector: prior SD(s) for the fixed-effect
#'   coefficients. `NULL` (default) uses a data-scaled value.
#' @param sigma_scale Positive scalar: half-normal scale for the ALD scale
#'   `sigma`. `NULL` (default) uses a data-scaled value.
#' @param re_scale Positive scalar: half-normal scale for the random-effect
#'   standard deviations. `NULL` (default) uses a data-scaled value.
#' @param lkj Positive scalar: LKJ shape parameter for the random-effect
#'   correlation matrix (used only when `cov = "unstructured"`). `2` favours
#'   weak correlations; `1` is uniform over correlation matrices.
#'
#' @return An object of class `"bqmm_prior"`.
#' @export
#' @examples
#' bqmm_prior(beta_sd = 5)
bqmm_prior <- function(beta_mean = 0, beta_sd = NULL,
                       sigma_scale = NULL, re_scale = NULL, lkj = 2) {
  structure(
    list(
      beta_mean   = beta_mean,
      beta_sd     = beta_sd,
      sigma_scale = sigma_scale,
      re_scale    = re_scale,
      lkj         = lkj
    ),
    class = "bqmm_prior"
  )
}

#' @export
print.bqmm_prior <- function(x, ...) {
  cat("<bqmm_prior>\n")
  for (nm in names(x)) {
    val <- x[[nm]]
    cat(sprintf("  %-12s %s\n", nm,
                if (is.null(val)) "data-scaled default" else
                  paste(format(val), collapse = ", ")))
  }
  invisible(x)
}

#' Resolve data-scaled default priors
#'
#' Fills in any `NULL` elements of a [bqmm_prior()] using simple, robust scales
#' derived from the response and design matrix. The intent is to keep `sigma`
#' identified and the sampler away from divergences, not to be informative.
#'
#' @param prior A `bqmm_prior` (or `NULL` for all defaults).
#' @param y Numeric response vector.
#' @param K Number of fixed-effect columns.
#' @return A fully-specified `bqmm_prior` with numeric hyperparameters.
#' @keywords internal
bqmm_default_priors <- function(prior, y, K) {
  if (is.null(prior)) prior <- bqmm_prior()
  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1

  if (is.null(prior$beta_sd))     prior$beta_sd     <- 2.5 * y_sd
  if (is.null(prior$sigma_scale)) prior$sigma_scale <- y_sd
  if (is.null(prior$re_scale))    prior$re_scale    <- y_sd
  if (is.null(prior$lkj))         prior$lkj         <- 2

  prior$beta_mean <- rep_len(prior$beta_mean, K)
  prior$beta_sd   <- rep_len(prior$beta_sd, K)
  prior
}
