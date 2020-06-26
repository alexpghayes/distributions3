context("test-Cauchy")

test_that("print.Cauchy works", {
  expect_output(print(Cauchy(1, 1)), regexp = "Cauchy distribution")
})

test_that("likelihood.Cauchy and log_likelihood.Cauchy work correctly", {
  cau <- Cauchy()
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dcauchy(1))
  expect_equal(likelihood(cau, x), dcauchy(1) * dcauchy(1) * dcauchy(0))

  expect_equal(log_likelihood(cau, 1), log(dcauchy(1)))
  expect_equal(log_likelihood(cau, x), log(dcauchy(1) * dcauchy(1) * dcauchy(0)))
})

test_that("random.Cauchy work correctly", {
  cau <- Cauchy()

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.Cauchy work correctly", {
  cau <- Cauchy()

  expect_equal(pdf(cau, 0), dcauchy(0, 0, 1))
  expect_equal(pdf(cau, 1), dcauchy(1, 0, 1))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.Cauchy work correctly", {
  cau <- Cauchy()

  expect_equal(log_pdf(cau, 0), log(dcauchy(0, 0, 1)))
  expect_equal(log_pdf(cau, 1), log(dcauchy(1, 0, 1)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.Cauchy work correctly", {
  cau <- Cauchy()

  expect_equal(cdf(cau, 0), pcauchy(0, 0, 1))
  expect_equal(cdf(cau, 1), pcauchy(1, 0, 1))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.Cauchy work correctly", {
  cau <- Cauchy()

  expect_equal(quantile(cau, 0), qcauchy(0, 0, 1))
  expect_equal(quantile(cau, 1), qcauchy(1, 0, 1))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})

test_that("{moments}.Cauchy work correctly", {
  cau <- Cauchy()

  expect_true(is.nan(mean(cau)))
  expect_true(is.nan(variance(cau)))
  expect_true(is.nan(skewness(cau)))
  expect_true(is.nan(kurtosis(cau)))
})
