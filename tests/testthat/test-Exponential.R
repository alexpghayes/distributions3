context("test-Exponential")

test_that("fit_mle.Exponential works correctly", {
  expect_equal(fit_mle(Exponential(), 1), Exponential(1))

  expect_error(fit_mle(Exponential(), -1))

  expect_true(is.numeric(fit_mle(Exponential(), rexp(100))$rate))
})

test_that("print.Beta works", {
  expect_output(print(Exponential()), regexp = "Exponential distribution")
})

test_that("random.Exponential work correctly", {
  e <- Exponential()

  expect_length(random(e), 1)
  expect_length(random(e, 100), 100)
  expect_length(random(e, 0), 0)
  expect_error(random(e, -2))
})

test_that("pdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(pdf(e, 0), 1)
  expect_equal(pdf(e, 1), 1 / exp(1))
  expect_equal(pdf(e, -12), 0)

  expect_length(pdf(e, seq_len(0)), 0)
  expect_length(pdf(e, seq_len(1)), 1)
  expect_length(pdf(e, seq_len(10)), 10)
})

test_that("pdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(log_pdf(e, 0), log(1))
  expect_equal(log_pdf(e, 1), log(1 / exp(1)))
  expect_equal(log_pdf(e, -12), log(0))

  expect_length(log_pdf(e, seq_len(0)), 0)
  expect_length(log_pdf(e, seq_len(1)), 1)
  expect_length(log_pdf(e, seq_len(10)), 10)
})

test_that("cdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(cdf(e, 0), 0)
  expect_equal(cdf(e, 1), 1 - 1 / exp(1))


  expect_length(cdf(e, seq_len(0)), 0)
  expect_length(cdf(e, seq_len(1)), 1)
  expect_length(cdf(e, seq_len(10)), 10)
})

test_that("quantile.Exponential work correctly", {
  e <- Exponential()

  expect_equal(quantile(e, 0), 0)
  expect_equal(quantile(e, 1), Inf)


  expect_length(quantile(e, seq_len(0)), 0)
  expect_length(quantile(e, c(0, 1)), 2)
})

test_that("{moments}.Exponential work correctly", {
  e <- Exponential()

  expect_equal(mean(e), 1)
  expect_equal(variance(e), 1)
  expect_equal(skewness(e), 2)
  expect_equal(kurtosis(e), 6)
})
