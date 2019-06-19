context("test-Poisson")

test_that("suff_stat.Poisson works correctly", {

  ss <- list(sum = 1, samples = 1)
  expect_equal(suff_stat(Poisson(1), 1), ss)

  ss <- list(sum = 5050, samples = 101)
  expect_equal(suff_stat(Poisson(1), 0:100), ss)

  expect_error(suff_stat(Poisson(1), 0.5))

  expect_error(suff_stat(Poisson(1), -1))
})

test_that("fit_mle.Poisson works correctly", {

  expect_equal(fit_mle(Poisson(1), 2), Poisson(2))

  expect_equal(fit_mle(Poisson(1), 0:100), Poisson(50))

  expect_error(fit_mle(Poisson(1), -1))

  expect_error(fit_mle(Poisson(1), 0.5))
})
