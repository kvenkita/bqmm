# The Yang-Wang-He sandwich core is pure-numeric and Stan-free, so we test it
# directly via compute_ywh_sandwich().

test_that("sandwich returns a symmetric PSD matrix of the right size", {
  set.seed(1)
  n <- 200
  X <- cbind(1, rnorm(n))
  resid <- bqmm:::rald(n, 0, 1, 0.5)
  sw <- bqmm:::compute_ywh_sandwich(X, resid, tau = 0.5)
  expect_equal(dim(sw$vcov), c(2L, 2L))
  expect_equal(sw$vcov, t(sw$vcov), tolerance = 1e-8)
  ev <- eigen(sw$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-8))
})

test_that("one-coefficient closed form matches the analytic sandwich", {
  # With X = intercept only, V = tau(1-tau) / (n f0^2), where f0 is the density
  # of the residuals at 0 (uniform-kernel Powell estimate).
  set.seed(2)
  n <- 500
  tau <- 0.5
  X <- matrix(1, n, 1)
  resid <- bqmm:::rald(n, 0, 1, tau)

  sw <- bqmm:::compute_ywh_sandwich(X, resid, tau = tau)
  h <- bqmm:::hall_sheather_bandwidth(n, tau)
  f0 <- mean(abs(resid) < h) / (2 * h)
  analytic <- tau * (1 - tau) / (n * f0^2)

  expect_equal(as.numeric(sw$vcov), analytic, tolerance = 1e-6)
})

test_that("adjustment changes interval width relative to a naive guess", {
  # Under the (misspecified) ALD, the naive posterior SD differs from the
  # sandwich SD; here we only assert they are not identical, i.e. the
  # correction does something.
  set.seed(3)
  n <- 300
  X <- cbind(1, rnorm(n))
  resid <- bqmm:::rald(n, 0, 1, 0.7)
  Sigma_post <- diag(c(0.01, 0.01))
  sw <- bqmm:::compute_ywh_sandwich(X, resid, tau = 0.7,
                                    Sigma_post = Sigma_post)
  expect_false(isTRUE(all.equal(diag(sw$vcov), diag(Sigma_post))))
})

test_that("rank-deficient bread is repaired, not collapsed (regression)", {
  # Previously a singular bread fell back to Sigma_post %*% D0 %*% Sigma_post,
  # which is O(1/n^2) and collapsed intervals. Forcing a near-empty kernel band
  # via a tiny bandwidth must instead grow the bandwidth / use the smooth
  # plug-in and return a finite, full-rank, correctly-scaled covariance.
  set.seed(11)
  n <- 400
  X <- cbind(1, rnorm(n))
  resid <- bqmm:::rald(n, 0, 1, 0.5)

  sw_tiny <- bqmm:::compute_ywh_sandwich(X, resid, tau = 0.5, bandwidth = 1e-8)
  sw_ok   <- bqmm:::compute_ywh_sandwich(X, resid, tau = 0.5)

  expect_true(all(is.finite(sw_tiny$vcov)))
  expect_equal(qr(sw_tiny$vcov)$rank, 2L)
  # same order of magnitude as the well-conditioned sandwich (not ~n smaller)
  ratio <- diag(sw_tiny$vcov) / diag(sw_ok$vcov)
  expect_true(all(ratio > 0.2 & ratio < 5))
})

test_that("sandwich variance scales like 1/n (not 1/n^2)", {
  set.seed(12)
  # intercept-only closed form V = tau(1-tau)/(n f0^2) is stable, so the
  # n-scaling is unambiguous: 8x n => ~8x smaller variance for 1/n (64x for the
  # old 1/n^2 collapse bug).
  v <- function(n) as.numeric(
    bqmm:::compute_ywh_sandwich(matrix(1, n, 1), bqmm:::rald(n, 0, 1, 0.5), 0.5)$vcov)
  v1 <- mean(sapply(1:60, function(i) v(500)))
  v8 <- mean(sapply(1:60, function(i) v(4000)))
  expect_gt(v1 / v8, 5.5)
  expect_lt(v1 / v8, 11)
})

test_that("multiplicative YWH form is symmetric PSD with correct construction", {
  set.seed(20)
  n <- 200
  X <- cbind(1, rnorm(n))
  resid <- bqmm:::rald(n, 0, 1, 0.5)
  Sigma_post <- crossprod(matrix(rnorm(4), 2)) / 50   # a plausible posterior cov
  ad <- bqmm:::compute_ywh_multiplicative(Sigma_post, X, resid, sigma = 1,
                                          tau = 0.5)
  expect_equal(ad$vcov, t(ad$vcov), tolerance = 1e-10)
  ev <- eigen(ad$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-8))
  # matches Sigma_post %*% G %*% Sigma_post exactly
  expect_equal(unname(ad$vcov), unname(Sigma_post %*% ad$G %*% Sigma_post),
               tolerance = 1e-10)
})

test_that("multiplicative meat scales as 1/sigma^2", {
  set.seed(21)
  X <- cbind(1, rnorm(150))
  resid <- bqmm:::rald(150, 0, 1, 0.5)
  Sp <- diag(c(0.02, 0.02))
  v1 <- bqmm:::compute_ywh_multiplicative(Sp, X, resid, sigma = 1, tau = 0.5)$vcov
  v2 <- bqmm:::compute_ywh_multiplicative(Sp, X, resid, sigma = 2, tau = 0.5)$vcov
  # G ~ 1/sigma^2 so V ~ 1/sigma^2; quadrupling implied by sigma 1 -> 2
  expect_equal(as.numeric(v1 / v2), rep(4, length(v1)), tolerance = 1e-8)
})

test_that("cluster-robust meat differs from independence meat", {
  set.seed(4)
  n <- 240
  g <- rep(1:12, each = n / 12)
  X <- cbind(1, rnorm(n))
  # induce within-cluster correlation in residual signs
  base <- bqmm:::rald(12, 0, 1, 0.5)[g]
  resid <- base + bqmm:::rald(n, 0, 0.3, 0.5)

  indep <- bqmm:::compute_ywh_sandwich(X, resid, tau = 0.5)$vcov
  clust <- bqmm:::compute_ywh_sandwich(X, resid, tau = 0.5, groups = g)$vcov
  expect_false(isTRUE(all.equal(indep, clust)))
})
