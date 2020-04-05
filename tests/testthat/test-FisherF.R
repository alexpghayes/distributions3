context("test-FisherF")

test_that("print.FisherF works", {
  expect_output(print(FisherF(1, 1)), regexp = "Fisher's F distribution")
})

test_that("likelihood.FisherF and log_likelihood.FisherF work correctly", {
  s <- FisherF(1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(s, 1), df(1, 1, 1, 0))
  expect_equal(likelihood(s, x), df(1, 1, 1, 0) * df(1, 1, 1, 0) * df(0, 1, 1, 0))

  expect_equal(log_likelihood(s, 1), log(df(1, 1, 1, 0)))
  expect_equal(log_likelihood(s, x), log(df(1, 1, 1, 0) * df(1, 1, 1, 0) * df(0, 1, 1, 0)))
})

test_that("random.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_length(random(s), 1)
  expect_length(random(s, 100), 100)
  expect_length(random(s, 0), 0)
  expect_error(random(s, -2))
})

test_that("pdf.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(pdf(s, 0), Inf)
  expect_equal(pdf(s, 1), df(1, 1, 1, 0))
  expect_equal(pdf(s, -12), 0)

  expect_length(pdf(s, seq_len(0)), 0)
  expect_length(pdf(s, seq_len(1)), 1)
  expect_length(pdf(s, seq_len(10)), 10)
})

test_that("log_pdf.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(log_pdf(s, 0), Inf)
  expect_equal(log_pdf(s, 1), log(df(1, 1, 1, 0)))
  expect_equal(log_pdf(s, -12), -Inf)

  expect_length(log_pdf(s, seq_len(0)), 0)
  expect_length(log_pdf(s, seq_len(1)), 1)
  expect_length(log_pdf(s, seq_len(10)), 10)
})

test_that("cdf.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(cdf(s, 0), pf(0, 1, 1, 0))
  expect_equal(cdf(s, 1), pf(1, 1, 1, 0))


  expect_length(cdf(s, seq_len(0)), 0)
  expect_length(cdf(s, seq_len(1)), 1)
  expect_length(cdf(s, seq_len(10)), 10)
})

test_that("quantile.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(quantile(s, 0), qf(0, 1, 1, 0))
  expect_equal(quantile(s, 1), qf(1, 1, 1, 0))


  expect_length(quantile(s, seq_len(0)), 0)
  expect_length(quantile(s, c(0, 1)), 2)
})
