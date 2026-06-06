#' bqmm: Bayesian Multilevel Quantile Regression
#'
#' Fits Bayesian mixed-effects (multilevel) quantile regression models using the
#' asymmetric Laplace working likelihood and Stan, with an `lme4`-style formula
#' interface, one or several quantiles per call, optional LKJ-correlated random
#' effects, post-hoc non-crossing rearrangement, and the Yang, Wang and He
#' (2016) posterior-variance correction for valid fixed-effect inference.
#'
#' @keywords internal
#' @name bqmm-package
#' @aliases bqmm-package
#'
#' @useDynLib bqmm, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats coef confint fitted nobs predict residuals vcov
#'
#' @references
#' Yu, K. and Moyeed, R. A. (2001). Bayesian quantile regression.
#' \emph{Statistics & Probability Letters}, 54(3), 437-447.
#'
#' Geraci, M. and Bottai, M. (2014). Linear quantile mixed models.
#' \emph{Statistics and Computing}, 24(3), 461-479.
#'
#' Yang, Y., Wang, H. J. and He, X. (2016). Posterior inference in Bayesian
#' quantile regression with asymmetric Laplace likelihood.
#' \emph{International Statistical Review}, 84(3), 327-344.
"_PACKAGE"
