#' Rearrange fitted quantiles to remove crossing
#'
#' Post-hoc monotonisation of estimated quantile curves (Chernozhukov,
#' Fernandez-Val and Galichon, 2010): at any covariate point, sorting the fitted
#' values across quantile levels into increasing order yields a valid,
#' non-crossing set of quantiles and never increases estimation error. This is
#' the v0.1 non-crossing strategy; joint constrained estimation is deferred.
#'
#' @param preds A numeric matrix of fitted quantiles with one column per
#'   quantile level, ordered by increasing `tau` (rows = observations).
#' @return A matrix of the same shape with each row sorted increasingly.
#' @export
#' @examples
#' m <- rbind(c(1, 0.5, 2), c(0, 1, 0.8))   # some crossings
#' rearrange_quantiles(m)
rearrange_quantiles <- function(preds) {
  preds <- as.matrix(preds)
  if (ncol(preds) <= 1L) return(preds)
  out <- t(apply(preds, 1L, sort))
  dimnames(out) <- dimnames(preds)
  out
}

#' Check whether quantile predictions are monotone in tau
#'
#' @param preds Matrix of fitted quantiles (columns ordered by increasing tau).
#' @return Logical scalar: `TRUE` if every row is non-decreasing.
#' @keywords internal
is_noncrossing <- function(preds) {
  preds <- as.matrix(preds)
  if (ncol(preds) <= 1L) return(TRUE)
  all(apply(preds, 1L, function(r) all(diff(r) >= 0)))
}
