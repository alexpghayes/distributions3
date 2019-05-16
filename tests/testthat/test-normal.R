context("test-normal")

test_that("suff_stat.normal works correctly", {

  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(normal(), c(0,0)), ss)

  expect_error(suff_stat(normal(), "abc"))
})

test_that("fit_mle.normal works correctly", {

  expect_equal(fit(normal(), c(0,0)), normal(0, 0))

})
