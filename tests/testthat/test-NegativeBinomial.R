context("test-NegativeBinomial")

test_that("print.NegativeBinomial works", {
  expect_output(print(NegativeBinomial(1, 1)), regexp = "Negative Binomial distribution")
})

test_that("likelihood.NegativeBinomial and log_likelihood.NegativeBinomial work correctly", {
  cau <- NegativeBinomial(1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dnbinom(1, 1, 1))
  expect_equal(likelihood(cau, x), dnbinom(1, 1, 1) * dnbinom(1, 1, 1) * dnbinom(0, 1, 1))

  expect_equal(log_likelihood(cau, 1), log(dnbinom(1, 1, 1)))
  expect_equal(log_likelihood(cau, x), log(dnbinom(1, 1, 1) * dnbinom(1, 1, 1) * dnbinom(0, 1, 1)))
})

test_that("random.NegativeBinomial work correctly", {
  cau <- NegativeBinomial(1, 1)

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.NegativeBinomial work correctly", {
  cau <- NegativeBinomial(1, 1)

  expect_equal(pdf(cau, 0), dnbinom(0, 0, 1))
  expect_equal(pdf(cau, 1), dnbinom(1, 0, 1))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.NegativeBinomial work correctly", {
  cau <- NegativeBinomial(1, 1)

  expect_equal(log_pdf(cau, 0), log(dnbinom(0, 0, 1)))
  expect_equal(log_pdf(cau, 1), log(dnbinom(1, 0, 1)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.NegativeBinomial work correctly", {
  cau <- NegativeBinomial(1, 1)

  expect_equal(cdf(cau, 0), pnbinom(0, 0, 1))
  expect_equal(cdf(cau, 1), pnbinom(1, 0, 1))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.NegativeBinomial work correctly", {
  cau <- NegativeBinomial(1, 1)

  expect_equal(quantile(cau, 0), qnbinom(0, 0, 1))
  expect_equal(quantile(cau, 1), qnbinom(1, 0, 1))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})
