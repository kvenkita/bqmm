test_that("ALD density integrates to one", {
  for (tau in c(0.25, 0.5, 0.8)) {
    area <- stats::integrate(
      function(x) bqmm:::dald(x, mu = 0, sigma = 1, tau = tau),
      lower = -Inf, upper = Inf
    )$value
    expect_equal(area, 1, tolerance = 1e-5)
  }
})

test_that("the tau-quantile of the ALD location is mu", {
  # P(Y <= mu) should equal tau for ALD(mu, sigma, tau)
  for (tau in c(0.2, 0.5, 0.9)) {
    cdf_at_mu <- stats::integrate(
      function(x) bqmm:::dald(x, mu = 0, sigma = 1, tau = tau),
      lower = -Inf, upper = 0
    )$value
    expect_equal(cdf_at_mu, tau, tolerance = 1e-5)
  }
})

test_that("rald draws have the right empirical quantile", {
  set.seed(42)
  for (tau in c(0.25, 0.75)) {
    y <- bqmm:::rald(2e5, mu = 0, sigma = 1, tau = tau)
    expect_equal(mean(y <= 0), tau, tolerance = 0.01)
  }
})

test_that("check loss is minimised at the location", {
  set.seed(7)
  y <- bqmm:::rald(1e4, mu = 2, sigma = 1, tau = 0.3)
  obj <- function(m) sum(bqmm:::rho_tau(y - m, 0.3))
  est <- stats::optimize(obj, c(-5, 10))$minimum
  expect_equal(est, 2, tolerance = 0.15)
})
