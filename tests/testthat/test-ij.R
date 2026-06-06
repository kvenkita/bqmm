# Infinitesimal Jackknife estimator. Pure-numeric core tests always run; the
# end-to-end comparison is gated on the compiled model.

has_stanmodel <- function() {
  exists("stanmodels", where = asNamespace("bqmm"), inherits = FALSE) &&
    is.list(get("stanmodels", envir = asNamespace("bqmm")))
}

test_that("compute_ij returns symmetric PSD matrices of the right size", {
  set.seed(1)
  S <- 1500; K <- 2; n <- 200
  beta <- matrix(rnorm(S * K), S, K)
  a <- matrix(rnorm(n * K, sd = 0.3), K, n)
  loglik <- beta %*% a + matrix(rnorm(S * n, sd = 0.5), S, n)
  V <- bqmm:::compute_ij(beta, loglik)
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values >= -1e-8))
})

test_that("cluster IJ differs from independence IJ", {
  set.seed(2)
  S <- 1500; K <- 2; n <- 240
  beta <- matrix(rnorm(S * K), S, K)
  a <- matrix(rnorm(n * K, sd = 0.3), K, n)
  loglik <- beta %*% a + matrix(rnorm(S * n, sd = 0.5), S, n)
  g <- rep(1:24, each = 10)
  expect_false(isTRUE(all.equal(
    bqmm:::compute_ij(beta, loglik),
    bqmm:::compute_ij(beta, loglik, groups = g))))
})

test_that("compute_ij errors on mismatched draw counts", {
  expect_error(bqmm:::compute_ij(matrix(0, 100, 2), matrix(0, 50, 10)),
               "same number of draws")
})

test_that("vcov method = 'ij' runs and is the same order as YWH", {
  skip_on_cran(); skip_if_not(has_stanmodel(), "no compiled model")
  set.seed(3)
  G <- 25; npg <- 10; n <- G * npg
  g <- factor(rep(seq_len(G), each = npg)); x <- rnorm(n)
  d <- data.frame(y = 1 + 2 * x + rnorm(G)[g] + rnorm(n), x = x, g = g)
  fit <- suppressWarnings(bqmm(y ~ x + (1 | g), d, tau = 0.5,
                              chains = 2, iter = 1000, seed = 3, refresh = 0))
  v_ij  <- vcov(fit, method = "ij")
  v_ywh <- vcov(fit, method = "ywh")
  expect_equal(dim(v_ij), c(2L, 2L))
  expect_true(all(is.finite(v_ij)))
  expect_true(all(eigen(v_ij, only.values = TRUE)$values >= -1e-8))
  # The within-cluster-identified SLOPE agrees with the YWH sandwich (same
  # family). The intercept is identified between clusters and the conditional IJ
  # under-estimates it -- a documented limitation -- so it is excluded here.
  r_slope <- sqrt(v_ij["x", "x"]) / sqrt(v_ywh["x", "x"])
  expect_true(r_slope > 0.5 && r_slope < 2)
})
