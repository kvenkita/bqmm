# Correlated random-effects model. Pure-R standata tests always run; fit tests
# are gated on the compiled bqmm_corr model.

has_corr_model <- function() {
  exists("stanmodels", where = asNamespace("bqmm"), inherits = FALSE) &&
    "bqmm_corr" %in% names(get("stanmodels", envir = asNamespace("bqmm")))
}

test_that("Zcov is reconstructed correctly from Z for (1 + x | g)", {
  skip_if_not_installed("lme4")
  set.seed(1)
  n <- 60
  d <- data.frame(y = rnorm(n), x = rnorm(n), g = factor(rep(1:6, length.out = n)))
  parsed <- bqmm:::bqmm_parse_formula(y ~ x + (1 + x | g), data = d)
  pr <- bqmm:::bqmm_default_priors(NULL, parsed$y, ncol(parsed$X))
  sd <- bqmm:::bqmm_corr_standata(parsed, tau = 0.5, prior = pr)

  expect_equal(sd$M, 2L)
  expect_equal(sd$L, nlevels(d$g))
  expect_equal(length(sd$level_id), n)
  expect_equal(dim(sd$Zcov), c(n, 2L))
  # column 1 is the intercept (all ones), column 2 is x
  expect_true(all(sd$Zcov[, 1] == 1))
  expect_equal(sd$Zcov[, 2], d$x, tolerance = 1e-12)
})

test_that("unstructured covariance rejects multiple RE terms", {
  skip_if_not_installed("lme4")
  set.seed(2)
  n <- 60
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  g = factor(rep(1:6, length.out = n)),
                  h = factor(rep(1:3, each = n / 3)))
  parsed <- bqmm:::bqmm_parse_formula(y ~ x + (1 | g) + (1 | h), data = d)
  pr <- bqmm:::bqmm_default_priors(NULL, parsed$y, ncol(parsed$X))
  expect_error(bqmm:::bqmm_corr_standata(parsed, 0.5, pr), "one random-effects")
})

test_that("correlated model fits and reports a correlation matrix", {
  skip_on_cran(); skip_if_not(has_corr_model(), "no compiled bqmm_corr model")
  skip_if_not_installed("MASS")
  set.seed(3)
  G <- 25; npg <- 10; n <- G * npg
  g <- factor(rep(seq_len(G), each = npg)); x <- rnorm(n)
  s0 <- 0.8; s1 <- 0.5; rho <- 0.5
  Su <- matrix(c(s0^2, rho*s0*s1, rho*s0*s1, s1^2), 2)
  U <- MASS::mvrnorm(G, c(0, 0), Su)
  y <- 1 + 2 * x + U[g, 1] + U[g, 2] * x + bqmm:::rald(n, 0, 1, 0.5)
  d <- data.frame(y = y, x = x, g = g)

  fit <- suppressWarnings(bqmm(y ~ x + (1 + x | g), d, tau = 0.5,
                              cov = "unstructured", chains = 2, iter = 800,
                              warmup = 400, seed = 3, refresh = 0))
  expect_identical(fit$cov, "unstructured")
  vc <- VarCorr(fit)
  cm <- attr(vc, "correlation")
  expect_equal(dim(cm), c(2L, 2L))
  expect_equal(unname(diag(cm)), c(1, 1), tolerance = 1e-6)
  expect_true(abs(cm[1, 2]) <= 1)                 # valid correlation
  expect_equal(rownames(cm), c("(Intercept)", "x"))
  # fixed effects recovered
  expect_equal(unname(fixef(fit)), c(1, 2), tolerance = 0.4)
  # predictive surface is coherent
  expect_equal(ncol(posterior_epred(fit)), n)
})
