context("test-Uniform")

test_that("print.Uniform works", {
  expect_output(print(Uniform(1, 1)), regexp = "Uniform distribution")
})

test_that("likelihood.Uniform and log_likelihood.Uniform work correctly", {
  u <- Uniform()
  x <- c(1, 1, 0)

  expect_equal(likelihood(u, 1), dunif(1))
  expect_equal(likelihood(u, x), dunif(1) * dunif(1) * dunif(0))

  expect_equal(log_likelihood(u, 1), log(dunif(1)))
  expect_equal(log_likelihood(u, x), log(dunif(1) * dunif(1) * dunif(0)))
})

test_that("random.Uniform work correctly", {
  u <- Uniform()

  expect_length(random(u), 1)
  expect_length(random(u, 100), 100)
  expect_length(random(u, 0), 0)
  expect_error(random(u, -2))
})

test_that("pdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(pdf(u, 0), dunif(0, 0, 1))
  expect_equal(pdf(u, 1), dunif(1, 0, 1))

  expect_length(pdf(u, seq_len(0)), 0)
  expect_length(pdf(u, seq_len(1)), 1)
  expect_length(pdf(u, seq_len(10)), 10)
})

test_that("log_pdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(log_pdf(u, 0), log(dunif(0, 0, 1)))
  expect_equal(log_pdf(u, 1), log(dunif(1, 0, 1)))

  expect_length(log_pdf(u, seq_len(0)), 0)
  expect_length(log_pdf(u, seq_len(1)), 1)
  expect_length(log_pdf(u, seq_len(10)), 10)
})

test_that("cdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(cdf(u, 0), punif(0, 0, 1))
  expect_equal(cdf(u, 1), punif(1, 0, 1))


  expect_length(cdf(u, seq_len(0)), 0)
  expect_length(cdf(u, seq_len(1)), 1)
  expect_length(cdf(u, seq_len(10)), 10)
})

test_that("quantile.Uniform work correctly", {
  u <- Uniform()

  expect_equal(quantile(u, 0), qunif(0, 0, 1))
  expect_equal(quantile(u, 1), qunif(1, 0, 1))


  expect_length(quantile(u, seq_len(0)), 0)
  expect_length(quantile(u, c(0, 1)), 2)
})

test_that("{moments}.Uniform work correctly", {
  u <- Uniform()

  expect_equal(mean(u), 0.5)
  expect_equal(variance(u), 1/12)
  expect_equal(skewness(u), 0)
  expect_equal(kurtosis(u), -6/5)
})
