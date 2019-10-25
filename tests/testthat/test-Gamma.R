context("test-Gamma")

test_that("print.Gamma works", {
  expect_output(print(Gamma(1, 1)), regexp = "Gamma distribution")
})

test_that("likelihood.Gamma and log_likelihood.Gamma work correctly", {
  cau <- Gamma(1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dgamma(1, 1, 1))
  expect_equal(likelihood(cau, x), dgamma(1, 1, 1) * dgamma(1, 1, 1) * dgamma(0, 1, 1))

  expect_equal(log_likelihood(cau, 1), log(dgamma(1, 1, 1)))
  expect_equal(log_likelihood(cau, x), log(dgamma(1, 1, 1) * dgamma(1, 1, 1) * dgamma(0, 1, 1)))
})

test_that("random.Gamma work correctly", {
  cau <- Gamma(1, 1)

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(pdf(cau, 0), dgamma(0, 0, 1))
  expect_equal(pdf(cau, 1), dgamma(1, 0, 1))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(log_pdf(cau, 0), log(dgamma(0, 0, 1)))
  expect_equal(log_pdf(cau, 1), log(dgamma(1, 0, 1)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(cdf(cau, 0), pgamma(0, 0, 1))
  expect_equal(cdf(cau, 1), pgamma(1, 0, 1))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(quantile(cau, 0), qgamma(0, 0, 1))
  expect_equal(quantile(cau, 1), qgamma(1, 0, 1))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})
