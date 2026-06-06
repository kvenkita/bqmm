# Formula translation is the riskiest reuse point (depends on lme4 internals),
# so these tests check it carefully against hand-built references. They run
# without Stan.

make_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  g <- factor(rep(1:6, length.out = n))
  h <- factor(rep(1:3, each = n / 3))
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = g,
    h = h
  )
}

test_that("random intercept: X, Z, and sd_map have expected shapes", {
  skip_if_not_installed("lme4")
  d <- make_data()
  p <- bqmm_parse_formula(y ~ x + (1 | g), data = d)

  expect_equal(ncol(p$X), 2L)                 # intercept + x
  expect_equal(colnames(p$X), c("(Intercept)", "x"))
  expect_equal(ncol(p$Z), nlevels(d$g))       # one column per group level
  expect_equal(length(p$sd_map), ncol(p$Z))
  expect_equal(length(p$re_components), 1L)   # one variance component
  expect_true(all(p$sd_map == 1L))            # all columns share it
})

test_that("random slope: two variance components, columns split correctly", {
  skip_if_not_installed("lme4")
  d <- make_data()
  p <- bqmm_parse_formula(y ~ x + (1 + x | g), data = d)

  # 2 coefficients x 6 levels = 12 random-effect columns
  expect_equal(ncol(p$Z), 2L * nlevels(d$g))
  expect_equal(length(p$re_components), 2L)
  # each variance component should map to exactly nlevels(g) columns
  expect_equal(as.integer(table(p$sd_map)),
               rep(nlevels(d$g), 2L))
})

test_that("crossed random effects are handled (two grouping factors)", {
  skip_if_not_installed("lme4")
  d <- make_data()
  p <- bqmm_parse_formula(y ~ x + (1 | g) + (1 | h), data = d)

  expect_equal(ncol(p$Z), nlevels(d$g) + nlevels(d$h))
  expect_equal(length(p$re_components), 2L)
  expect_equal(sort(unique(p$sd_map)), c(1L, 2L))
})

test_that("Z matches t(Zt) from lme4 exactly", {
  skip_if_not_installed("lme4")
  d <- make_data()
  ref <- lme4::glFormula(y ~ x + (1 + x | g), data = d,
                         family = stats::gaussian())
  Zref <- as.matrix(Matrix::t(ref$reTrms$Zt))
  p <- bqmm_parse_formula(y ~ x + (1 + x | g), data = d)
  expect_equal(unname(p$Z), unname(Zref))
})

test_that("response is required", {
  skip_if_not_installed("lme4")
  d <- make_data()
  expect_error(bqmm_parse_formula(~ x + (1 | g), data = d))
})
