context("test-methods")

test_that("pmf works", {
  N <- normal()
  B <- Bernoulli()

  expect_equal(pmf(N, 0), pdf(N, 0))
  expect_equal(pmf(B, 0), pdf(B, 0))
})

test_that("fit works", {
  N <- normal()
  bern <- Bernoulli()
  B <- Binomial(10)

  expect_equal(normal(1, 0), fit(N, c(1, 1)))
  expect_equal(Bernoulli(), fit(bern, c(0, 1)))
  expect_equal(Binomial(10, 0.5), fit(B, c(5,5)))

})
