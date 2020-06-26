context("test-Beta")

test_that("print.Beta works", {
  expect_output(print(Beta()), regexp = "Beta distribution")
})

test_that("random.Beta work correctly", {
  b <- Beta(1, 5)

  expect_length(random(b), 1)
  expect_length(random(b, 100), 100)
  expect_length(random(b, 0), 0)
  expect_error(random(b, -2))
})

test_that("pdf.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(pdf(b, 0), 2)
  expect_equal(pdf(b, 0.5), 1)
  expect_equal(pdf(b, -12), 0)

  expect_length(pdf(b, seq_len(0)), 0)
  expect_length(pdf(b, seq_len(1)), 1)
  expect_length(pdf(b, seq_len(10)), 10)
})

test_that("log_pdf.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(log_pdf(b, 0), log(2))
  expect_equal(log_pdf(b, 0.5), log(1))
  expect_equal(log_pdf(b, -12), log(0))

  expect_length(log_pdf(b, seq_len(0)), 0)
  expect_length(log_pdf(b, seq_len(1)), 1)
  expect_length(log_pdf(b, seq_len(10)), 10)
})

test_that("cdf.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(cdf(b, 0), 0)
  expect_equal(cdf(b, 0.5), 0.75)
  expect_equal(cdf(b, 1), 1)


  expect_length(cdf(b, seq_len(0)), 0)
  expect_length(cdf(b, seq_len(1)), 1)
  expect_length(cdf(b, seq_len(10)), 10)
})

test_that("quantile.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(quantile(b, 0), 0)
  expect_equal(quantile(b, 0.75), 0.5)
  expect_equal(quantile(b, 1), 1)


  expect_length(quantile(b, seq_len(0)), 0)
  expect_length(quantile(b, c(0, 1)), 2)
})

test_that("{moments}.Beta work correctly", {
  b <- Beta(1, 1)

  expect_equal(mean(b), 0.5)
  expect_equal(variance(b), 1/12)
  expect_equal(skewness(b), 0)
  expect_equal(kurtosis(b), -6/5)
})
