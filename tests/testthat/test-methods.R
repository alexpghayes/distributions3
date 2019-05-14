context("test-methods")

test_that("pmf works", {
  N <- normal()
  B <- bernoulli()

  expect_equal(pmf(N, 0), pdf(N, 0))
  expect_equal(pmf(B, 0), pdf(B, 0))
})
