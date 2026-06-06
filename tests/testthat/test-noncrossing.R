test_that("rearrangement makes every row non-decreasing", {
  m <- rbind(c(1, 0.5, 2),
             c(0, 1, 0.8),
             c(-1, -2, 0))
  out <- rearrange_quantiles(m)
  expect_true(bqmm:::is_noncrossing(out))
})

test_that("rearrangement is the identity when already monotone", {
  m <- rbind(c(0.1, 0.5, 0.9),
             c(-1, 0, 1))
  out <- rearrange_quantiles(m)
  expect_equal(unname(out), unname(m))
})

test_that("rearrangement preserves the multiset of each row", {
  set.seed(10)
  m <- matrix(rnorm(30), nrow = 6)
  out <- rearrange_quantiles(m)
  for (i in seq_len(nrow(m))) {
    expect_equal(sort(m[i, ]), unname(out[i, ]))
  }
})

test_that("single-column input is returned unchanged", {
  m <- matrix(c(3, 1, 2), ncol = 1)
  expect_equal(rearrange_quantiles(m), m)
})
