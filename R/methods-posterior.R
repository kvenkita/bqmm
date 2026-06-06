# Integration with the `posterior` and `bayesplot` ecosystems.

#' Convert a bqmm fit to a posterior draws object
#'
#' @param x A `bqmm` fit.
#' @param ... Unused.
#' @return A `draws_array` (from the `posterior` package) with tidy variable
#'   names: `b_<name>` for fixed effects, `sd_<component>` for random-effect
#'   SDs, and `sigma`.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' as_draws(fit)
#' }
#' @export
as_draws.bqmm <- function(x, ...) {
  arr <- as.array(x$stanfit)        # iterations x chains x parameters
  d <- posterior::as_draws_array(arr)
  posterior::variables(d) <- tidy_param_names(posterior::variables(d), x)
  d
}

#' Coerce a bqmm fit to a matrix of posterior draws
#'
#' @param x A `bqmm` fit.
#' @param ... Unused.
#' @return A matrix of posterior draws (draws in rows, parameters in columns),
#'   using the raw Stan parameter names.
#' @examples
#' \donttest{
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' dim(as.matrix(fit))
#' }
#' @export
as.matrix.bqmm <- function(x, ...) {
  as.matrix(x$stanfit)
}

#' Map raw Stan parameter names to tidy bqmm names
#' @keywords internal
tidy_param_names <- function(nms, object) {
  fixed <- object$parsed$fixed_names
  comps <- object$parsed$re_components
  out <- nms
  for (k in seq_along(fixed)) {
    out[out == sprintf("beta[%d]", k)] <- paste0("b_", fixed[k])
  }
  for (g in seq_along(comps)) {
    out[out == sprintf("sd_re[%d]", g)] <-
      paste0("sd_", gsub("[^[:alnum:]]+", "_", comps[g]))
  }
  out
}
