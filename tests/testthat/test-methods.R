context("test-methods")

test_that("pmf() works", {
  N <- Normal()
  B <- Bernoulli()

  expect_equal(pmf(N, 0), pdf(N, 0))
  expect_equal(pmf(B, 0), pdf(B, 0))
})

test_that("fit_mle() works", {
  N <- Normal()
  bern <- Bernoulli()
  B <- Binomial(10)

  expect_equal(Normal(1, 0), fit_mle(N, c(1, 1)))
  expect_equal(Bernoulli(), fit_mle(bern, c(0, 1)))
  expect_equal(Binomial(10, 0.5), fit_mle(B, c(5, 5)))
})
