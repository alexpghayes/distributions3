context("test-LogNormal")

test_that("print.LogNormal works", {
  expect_output(print(LogNormal()), regexp = "Lognormal distribution")
})

test_that("likelihood.LogNormal and log_likelihood.LogNormal work correctly", {
  cau <- LogNormal()
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dlnorm(1))
  expect_equal(likelihood(cau, x), dlnorm(1) * dlnorm(1) * dlnorm(0))

  expect_equal(log_likelihood(cau, 1), log(dlnorm(1)))
  expect_equal(log_likelihood(cau, x), log(dlnorm(1) * dlnorm(1) * dlnorm(0)))
})

test_that("random.LogNormal work correctly", {
  cau <- LogNormal()

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.LogNormal work correctly", {
  cau <- LogNormal()

  expect_equal(pdf(cau, 0), dlnorm(0, 0, 1))
  expect_equal(pdf(cau, 1), dlnorm(1, 0, 1))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.LogNormal work correctly", {
  cau <- LogNormal()

  expect_equal(log_pdf(cau, 0), log(dlnorm(0, 0, 1)))
  expect_equal(log_pdf(cau, 1), log(dlnorm(1, 0, 1)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.LogNormal work correctly", {
  cau <- LogNormal()

  expect_equal(cdf(cau, 0), plnorm(0, 0, 1))
  expect_equal(cdf(cau, 1), plnorm(1, 0, 1))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.LogNormal work correctly", {
  cau <- LogNormal()

  expect_equal(quantile(cau, 0), qlnorm(0, 0, 1))
  expect_equal(quantile(cau, 1), qlnorm(1, 0, 1))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})
