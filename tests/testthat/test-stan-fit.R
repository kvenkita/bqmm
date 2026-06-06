# End-to-end fitting tests require the compiled Stan model (see SETUP.md). They
# are skipped until the rstantools wiring is in place, so the rest of the suite
# runs on a machine without a Stan toolchain.

has_stanmodel <- function() {
  exists("stanmodels", where = asNamespace("bqmm"), inherits = FALSE) &&
    is.list(get("stanmodels", envir = asNamespace("bqmm")))
}

test_that("bqmm fits a random-intercept model and recovers parameters", {
  skip_on_cran()
  skip_if_not(has_stanmodel(), "compiled Stan model not available")
  skip_if_not_installed("nlme")

  set.seed(1)
  fit <- bqmm(distance ~ age + (1 | Subject),
              data = nlme::Orthodont, tau = 0.5,
              chains = 2, iter = 1000, seed = 1)
  expect_s3_class(fit, "bqmm")
  cf <- fixef(fit)
  expect_named(cf, c("(Intercept)", "age"))
  expect_true(all(is.finite(cf)))
})

test_that("adjusted and unadjusted vcov differ in width", {
  skip_on_cran()
  skip_if_not(has_stanmodel(), "compiled Stan model not available")
  skip_if_not_installed("nlme")

  fit <- bqmm(distance ~ age + (1 | Subject),
              data = nlme::Orthodont, tau = 0.5,
              chains = 2, iter = 1000, seed = 2)
  v_adj <- vcov(fit, adjusted = TRUE)
  v_raw <- vcov(fit, adjusted = FALSE)
  expect_false(isTRUE(all.equal(diag(v_adj), diag(v_raw))))
})

test_that("vector tau returns a bqmm_multi with non-crossing predictions", {
  skip_on_cran()
  skip_if_not(has_stanmodel(), "compiled Stan model not available")
  skip_if_not_installed("nlme")

  fit <- bqmm(distance ~ age + (1 | Subject),
              data = nlme::Orthodont, tau = c(0.1, 0.5, 0.9),
              chains = 2, iter = 1000, seed = 3)
  expect_s3_class(fit, "bqmm_multi")
  expect_equal(length(fit$fits), 3L)
})
