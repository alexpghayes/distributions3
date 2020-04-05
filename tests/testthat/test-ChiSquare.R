test_that("print.ChiSquare works", {
  expect_output(print(ChiSquare(df = 1)), regexp = "Chi Square distribution")
})

test_that("random.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_length(random(cs), 1)
  expect_length(random(cs, 100), 100)
  expect_length(random(cs, 0), 0)
  expect_error(random(cs, -2))
})

test_that("pdf.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(pdf(cs, 0), Inf)
  expect_equal(pdf(cs, 1), dchisq(1, 1))
  expect_equal(pdf(cs, -12), 0)

  expect_length(pdf(cs, seq_len(0)), 0)
  expect_length(pdf(cs, seq_len(1)), 1)
  expect_length(pdf(cs, seq_len(10)), 10)
})

test_that("log_pdf.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(log_pdf(cs, 0), Inf)
  expect_equal(log_pdf(cs, 1), log(dchisq(1, 1)))
  expect_equal(log_pdf(cs, -12), log(0))

  expect_length(log_pdf(cs, seq_len(0)), 0)
  expect_length(log_pdf(cs, seq_len(1)), 1)
  expect_length(log_pdf(cs, seq_len(10)), 10)
})

test_that("cdf.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(cdf(cs, 0), 0)
  expect_equal(cdf(cs, 1), pchisq(1, 1))


  expect_length(cdf(cs, seq_len(0)), 0)
  expect_length(cdf(cs, seq_len(1)), 1)
  expect_length(cdf(cs, seq_len(10)), 10)
})

test_that("quantile.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(quantile(cs, 0), 0)
  expect_equal(quantile(cs, 0.5), qchisq(0.5, 1))


  expect_length(quantile(cs, seq_len(0)), 0)
  expect_length(quantile(cs, c(0, 1)), 2)
})

test_that("{moments}.Normal work correctly", {
  df <- 5
  cs <- ChiSquare(5)

  expect_equal(mean(cs), df)
  expect_equal(variance(cs), 2 * df)
  expect_equal(skewness(cs), sqrt(8 / df))
  expect_equal(kurtosis(cs), 12 / df)
})
