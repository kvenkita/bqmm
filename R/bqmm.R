#' Bayesian multilevel quantile regression
#'
#' Fits a Bayesian mixed-effects quantile regression model using the asymmetric
#' Laplace working likelihood and Stan. The interface follows `lme4`: random
#' effects are written inline in the formula, e.g. `y ~ x + (1 + x | group)`,
#' and nested or crossed grouping factors are both supported.
#'
#' One or several quantiles may be requested through `tau`. A scalar returns a
#' single `bqmm` fit; a vector fits each quantile independently and returns a
#' `bqmm_multi` container.
#'
#' @param formula An lme4-style model formula.
#' @param data A data frame containing the variables in `formula`.
#' @param tau Quantile level(s) in (0, 1). Scalar or vector.
#' @param family A `bqmm_family` object; currently only [ald()].
#' @param prior A [bqmm_prior()] object, or `NULL` for data-scaled defaults.
#' @param cov Random-effect covariance structure. `"diagonal"` (default)
#'   models independent random effects and supports any number of nested or
#'   crossed terms. `"unstructured"` adds an LKJ-correlated covariance but
#'   currently requires exactly one random-effects term (e.g.
#'   `y ~ x + (1 + x | g)`).
#' @param adjust Logical; compute the Yang-Wang-He (2016) variance correction so
#'   that `vcov(fit, adjusted = TRUE)` returns valid fixed-effect uncertainty.
#'   Default `TRUE`.
#' @param prior_only Logical; sample from the prior predictive distribution.
#' @param chains,iter,warmup,cores,seed Passed to [rstan::sampling()].
#' @param control A list of sampler control parameters (e.g. `adapt_delta`).
#'   Defaults raise `adapt_delta` to 0.95 because ALD posteriors are sharp.
#' @param ... Additional arguments forwarded to [rstan::sampling()].
#'
#' @return A `bqmm` object (single `tau`) or a `bqmm_multi` object (vector
#'   `tau`).
#' @export
#' @examples
#' \donttest{
#' # A minimal fit; raise chains/iter for real analyses.
#' fit <- bqmm(distance ~ age + (1 | Subject), data = nlme::Orthodont,
#'             tau = 0.5, chains = 1, iter = 300, refresh = 0, seed = 1)
#' summary(fit)
#' }
bqmm <- function(formula, data, tau = 0.5,
                 family = ald(), prior = NULL,
                 cov = c("diagonal", "unstructured"),
                 adjust = TRUE, prior_only = FALSE,
                 chains = 4, iter = 2000, warmup = floor(iter / 2),
                 cores = getOption("mc.cores", 1L), seed = NULL,
                 control = list(adapt_delta = 0.95), ...) {

  if (!inherits(family, "bqmm_family")) {
    stop("`family` must be created with ald().", call. = FALSE)
  }
  if (any(tau <= 0 | tau >= 1)) {
    stop("All `tau` values must lie strictly in (0, 1).", call. = FALSE)
  }
  cov <- match.arg(cov)

  mc <- match.call()
  parsed <- bqmm_parse_formula(formula, data)
  resolved_prior <- bqmm_default_priors(prior, parsed$y, ncol(parsed$X))

  model <- if (cov == "unstructured") "bqmm_corr" else "bqmm"

  fit_one <- function(tau_i) {
    standata <- if (cov == "unstructured") {
      bqmm_corr_standata(parsed, tau_i, resolved_prior, prior_only)
    } else {
      bqmm_standata(parsed, tau_i, resolved_prior, prior_only)
    }
    stanfit  <- bqmm_sample(standata, model = model, chains = chains, iter = iter,
                            warmup = warmup, cores = cores, seed = seed,
                            control = control, ...)
    out <- new_bqmm(
      stanfit = stanfit,
      parsed  = parsed,
      tau     = tau_i,
      prior   = resolved_prior,
      family  = family,
      cov     = cov,
      call    = mc,
      formula = formula
    )
    if (isTRUE(adjust)) {
      out$adjustment <- ywh_adjust(out)
    }
    out
  }

  if (length(tau) == 1L) {
    return(fit_one(tau))
  }

  fits <- lapply(sort(tau), fit_one)
  names(fits) <- format(sort(tau))
  new_bqmm_multi(fits, parsed = parsed, formula = formula, call = mc)
}

#' Construct a bqmm object
#' @keywords internal
new_bqmm <- function(stanfit, parsed, tau, prior, family, cov = "diagonal",
                     call, formula) {
  structure(
    list(
      stanfit    = stanfit,
      parsed     = parsed,
      tau        = tau,
      prior      = prior,
      family     = family,
      cov        = cov,
      call       = call,
      formula    = formula,
      adjustment = NULL
    ),
    class = "bqmm"
  )
}
