# Generics defined in other packages and re-exported so that, e.g.,
# `fixef(fit)` and `posterior_predict(fit)` work after `library(bqmm)` without
# attaching lme4/posterior/rstantools. Documented on a single reexports page.

#' @importFrom lme4 fixef
#' @export
lme4::fixef

#' @importFrom lme4 ranef
#' @export
lme4::ranef

#' @importFrom lme4 VarCorr
#' @export
lme4::VarCorr

#' @importFrom posterior as_draws
#' @export
posterior::as_draws

#' @importFrom rstantools posterior_predict
#' @export
rstantools::posterior_predict

#' @importFrom rstantools posterior_epred
#' @export
rstantools::posterior_epred

#' @importFrom rstantools log_lik
#' @export
rstantools::log_lik
