context("test-binomial")

test_that("fit.binomial works correctly", {

  # only one trial
  expect_equal(fit(binomial(1), c(0,1)), binomial(1, 0.5))

  # many trials
  x <- c(0, rep(100, 99))
  expect_equal(fit(binomial(100), x), binomial(100, 0.99))

  # raises error when data has negative counts
  expect_error(fit(binomial(1), -1))

  # raises error when data has counts greater than binomial's size
  expect_error(fit(binomial(1), 2))
})
