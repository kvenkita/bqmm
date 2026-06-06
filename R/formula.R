#' Parse an lme4-style mixed-model formula
#'
#' Thin wrapper around [lme4::glFormula()] that extracts everything `bqmm`
#' needs: the fixed-effect design matrix `X`, the random-effect design matrix
#' `Z` (dense), the response `y`, and a mapping from each column of `Z` to a
#' variance component. Reusing lme4's parser means nested *and* crossed random
#' effects are handled for free.
#'
#' @param formula A model formula such as `y ~ x + (1 + x | group)`.
#' @param data A data frame.
#' @param na.action,contrasts Passed through to model-frame construction.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{y}{numeric response vector.}
#'     \item{X}{fixed-effect design matrix (N x K).}
#'     \item{Z}{random-effect design matrix (N x Q), dense.}
#'     \item{sd_map}{integer vector (length Q) mapping each Z column to a
#'       variance component in `1:G`.}
#'     \item{re_components}{character labels for the `G` variance components.}
#'     \item{re_terms}{per-term metadata: grouping factor, coefficient names,
#'       number of levels.}
#'     \item{cnms, flist, Gp}{the raw lme4 random-effect structures.}
#'     \item{fixed_names}{column names of `X`.}
#'   }
#' @keywords internal
bqmm_parse_formula <- function(formula, data, na.action = stats::na.omit,
                               contrasts = NULL) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required to parse the model formula.", call. = FALSE)
  }

  pf <- lme4::glFormula(
    formula  = formula,
    data     = data,
    family   = stats::gaussian(),
    na.action = na.action,
    contrasts = contrasts
  )

  y <- stats::model.response(pf$fr)
  if (is.null(y)) stop("The formula must have a response on the left-hand side.",
                       call. = FALSE)
  y <- as.numeric(y)

  X <- pf$X
  reTrms <- pf$reTrms

  # Z: lme4 stores Zt (Q x N, sparse); we want a dense N x Q matrix.
  Z <- as.matrix(Matrix::t(reTrms$Zt))

  cnms <- reTrms$cnms          # named list: factor -> coefficient names
  flist <- reTrms$flist        # grouping factors
  Gp <- reTrms$Gp              # block boundaries into the columns of Z

  re <- build_sd_map(cnms, Gp)

  list(
    y             = y,
    X             = X,
    Z             = Z,
    sd_map        = re$sd_map,
    re_components = re$components,
    re_terms      = re$terms,
    cnms          = cnms,
    flist         = flist,
    Gp            = Gp,
    fixed_names   = colnames(X)
  )
}

#' Map random-effect columns to variance components
#'
#' lme4 orders the columns of `t(Zt)` block by term (delimited by `Gp`), and
#' within a term block the coefficients vary fastest within each level. This
#' function turns that layout into an integer `sd_map` assigning each column to
#' a `(term, coefficient)` variance component, plus human-readable labels.
#'
#' @param cnms Named list of coefficient names per grouping factor (from
#'   [lme4::mkReTrms()]).
#' @param Gp Integer vector of block boundaries (from `mkReTrms`).
#' @return List with `sd_map`, `components` (labels) and `terms` (metadata).
#' @keywords internal
build_sd_map <- function(cnms, Gp) {
  n_terms <- length(cnms)
  if (n_terms == 0L) {
    return(list(sd_map = integer(0), components = character(0), terms = list()))
  }

  sd_map     <- integer(0)
  components <- character(0)
  terms_meta <- vector("list", n_terms)
  comp_offset <- 0L

  for (i in seq_len(n_terms)) {
    coef_names <- cnms[[i]]
    p_i  <- length(coef_names)
    rows <- (Gp[i] + 1L):Gp[i + 1L]
    n_block <- length(rows)
    n_lev   <- n_block / p_i

    # coefficient index varies fastest within each level
    coef_idx <- ((seq_len(n_block) - 1L) %% p_i) + 1L
    sd_map   <- c(sd_map, comp_offset + coef_idx)

    components <- c(components,
                    paste(names(cnms)[i], coef_names, sep = " : "))
    terms_meta[[i]] <- list(
      group   = names(cnms)[i],
      coefs   = coef_names,
      n_levels = n_lev
    )
    comp_offset <- comp_offset + p_i
  }

  list(sd_map = sd_map, components = components, terms = terms_meta)
}
