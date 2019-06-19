context("test-exponential")

test_that("suff_stat.Exponential works correctly", {

  ss <- list(sum = 1, samples = 1)
  expect_equal(suff_stat(Exponential(), 1), ss)

  ss <- list(sum = 25, samples = 5)
  expect_equal(suff_stat(Exponential(), rep(5,5)), ss)

  expect_error(suff_stat(Exponential, -1))

})

test_that("fit_mle.Exponential works correctly", {

  expect_equal(fit_mle(Exponential(), 1), Exponential(1))

  expect_error(fit_mle(Exponential(), -1))

  expect_true(is.numeric(fit_mle(Exponential(), rexp(100))$rate))
})
