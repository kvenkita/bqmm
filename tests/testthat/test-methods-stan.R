# Tests for posterior/predict S3 methods that need the compiled model.
# Skipped until the Stan model is available (see SETUP.md).

has_stanmodel <- function() {
  exists("stanmodels", where = asNamespace("bqmm"), inherits = FALSE) &&
    is.list(get("stanmodels", envir = asNamespace("bqmm")))
}

make_fit <- function(seed = 1) {
  set.seed(seed)
  g <- factor(rep(1:8, each = 6)); x <- rnorm(48)
  d <- data.frame(y = 1 + 2 * x + rnorm(8)[g] + rnorm(48), x = x, g = g)
  suppressWarnings(bqmm(y ~ x + (1 | g), d, tau = 0.5,
                        chains = 2, iter = 600, warmup = 300, seed = seed, refresh = 0))
}

test_that("posterior_epred/predict have coherent dimensions and finite values", {
  skip_on_cran(); skip_if_not(has_stanmodel(), "no compiled model")
  fit <- make_fit()
  n <- nobs(fit)
  ep <- posterior_epred(fit)
  yp <- posterior_predict(fit)
  expect_equal(ncol(ep), n)
  expect_equal(dim(yp), dim(ep))
  expect_true(all(is.finite(ep)))
  expect_true(all(is.finite(yp)))
})

test_that("posterior_epred equals X*beta + Z*b reconstructed from a joint extract", {
  skip_on_cran(); skip_if_not(has_stanmodel(), "no compiled model")
  fit <- make_fit(2)
  ep <- posterior_epred(fit)
  ex <- rstan::extract(fit$stanfit, pars = c("beta", "b"), permuted = TRUE)
  loc <- as.matrix(ex$beta) %*% t(fit$parsed$X) +
    as.matrix(ex$b) %*% t(fit$parsed$Z)
  expect_equal(unname(ep), unname(loc), tolerance = 1e-10)
})

test_that("as_draws returns a draws object with tidy names", {
  skip_on_cran(); skip_if_not(has_stanmodel(), "no compiled model")
  skip_if_not_installed("posterior")
  fit <- make_fit(3)
  d <- as_draws(fit)
  v <- posterior::variables(d)
  expect_true("b_(Intercept)" %in% v || any(grepl("^b_", v)))
  expect_true(any(grepl("^sd_", v)))
  expect_true("sigma" %in% v)
})

test_that("fitted tracks the central tendency of posterior_epred draws", {
  skip_on_cran(); skip_if_not(has_stanmodel(), "no compiled model")
  fit <- make_fit(4)
  ep <- posterior_epred(fit)
  # fitted uses median(beta)/median(b); the median of sums is not identical to
  # the sum of medians, but the two should be very highly correlated.
  expect_gt(stats::cor(fitted(fit), apply(ep, 2, stats::median)), 0.99)
})
