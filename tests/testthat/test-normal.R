context("test-normal")

test_that("suff_stat.Normal works correctly", {

  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(Normal(), c(0, 0)), ss)

  expect_error(suff_stat(Normal(), "abc"))
})

test_that("fit_mle.Normal works correctly", {

  expect_equal(fit_mle(Normal(), c(0, 0)), Normal(0, 0))

})
