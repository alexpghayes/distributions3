context("test-bernoulli")

test_that("fit.bernoulli works correctly", {

  expect_equal(fit(bernoulli(), c(0,1)), bernoulli(0.5))

})

test_that("suff_stats.bernoulli works correctly", {

  ss <- list(successes = 3, failures = 2)
  expect_equal(suff_stat(bernoulli(), c(1,1,1,0,0)), ss)

  expect_error(suff_stat(bernoulli(), 2))

  expect_error(suff_stat(bernoulli(), -1))
})
