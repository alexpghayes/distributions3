context("test-bernoulli")

test_that("fit.Bernoulli works correctly", {

  expect_equal(fit(Bernoulli(), c(0,1)), Bernoulli(0.5))

})

test_that("suff_stats.Bernoulli works correctly", {

  ss <- list(successes = 3, failures = 2)
  expect_equal(suff_stat(Bernoulli(), c(1,1,1,0,0)), ss)

  expect_error(suff_stat(Bernoulli(), 2))

  expect_error(suff_stat(Bernoulli(), -1))
})


test_that("likelihood.Bernoulli and log_likelihood.Bernoulli work correctly", {

  b <- Bernoulli(0.1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(b, 1), 0.1)
  expect_equal(likelihood(b, x), 0.009)

  expect_equal(log_likelihood(b, 1), log(0.1))
  expect_equal(log_likelihood(b, x), log(0.009))

})
