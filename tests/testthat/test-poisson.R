context("test-poisson")

test_that("suff_stat.poisson works correctly", {

  ss <- list(sum = 1, samples = 1)
  expect_equal(suff_stat(poisson(1), 1), ss)

  ss <- list(sum = 5050, samples = 101)
  expect_equal(suff_stat(poisson(1), 0:100), ss)

  expect_error(suff_stat(poisson(1), 0.5))

  expect_error(suff_stat(poisson(1), -1))
})

test_that("fit_mle.poisson works correctly", {

  expect_equal(fit(poisson(1), 2), poisson(2))

  expect_equal(fit(poisson(1), 0:100), poisson(50))

  expect_error(fit(poisson(1), -1))

  expect_error(fit(poisson(1), 0.5))
})
