context("test-exponential")

test_that("exponential() works correctly", {

})

test_that("suff_stat.exponential works correctly", {

  ss <- list(sum = 1, samples = 1)
  expect_equal(suff_stat(exponential(), 1), ss)

  ss <- list(sum = 25, samples = 5)
  expect_equal(suff_stat(exponential(), rep(5,5)), ss)

  expect_error(suff_stat(exponential, -1))

})

test_that("fit_mle.exponential works correctly", {

  expect_equal(fit(exponential(), 1), exponential(1))

  expect_error(fit(exponential(), -1))

  expect_true(is.numeric(fit(exponential(), rexp(100))$rate))
})
