context("test-Weibull")

test_that("print.Weibull works", {
  expect_output(print(Weibull(1, 1)), regexp = "Weibull distribution")
})

test_that("likelihood.Weibull and log_likelihood.Weibull work correctly", {
  w <- Weibull(1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(w, 1), 1 / exp(1))
  expect_equal(likelihood(w, x), (1 / exp(1))^2)

  expect_equal(log_likelihood(w, 1), log(1 / exp(1)))
  expect_equal(log_likelihood(w, x), log((1 / exp(1))^2))
})

test_that("random.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_length(random(w), 1)
  expect_length(random(w, 100), 100)
  expect_length(random(w, 0), 0)
  expect_error(random(w, -2))
})

test_that("pdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(pdf(w, 0), 1)
  expect_equal(pdf(w, 1), 1 / exp(1))
  expect_equal(pdf(w, -12), 0)

  expect_length(pdf(w, seq_len(0)), 0)
  expect_length(pdf(w, seq_len(1)), 1)
  expect_length(pdf(w, seq_len(10)), 10)
})

test_that("log_pdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(log_pdf(w, 0), 0)
  expect_equal(log_pdf(w, 1), log(1 / exp(1)))
  expect_equal(log_pdf(w, -12), -Inf)

  expect_length(log_pdf(w, seq_len(0)), 0)
  expect_length(log_pdf(w, seq_len(1)), 1)
  expect_length(log_pdf(w, seq_len(10)), 10)
})

test_that("cdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(cdf(w, 0), 0)
  expect_equal(cdf(w, 1), 1 - 1 / exp(1))


  expect_length(cdf(w, seq_len(0)), 0)
  expect_length(cdf(w, seq_len(1)), 1)
  expect_length(cdf(w, seq_len(10)), 10)
})

test_that("quantile.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(quantile(w, 0), 0)
  expect_equal(quantile(w, 1), Inf)


  expect_length(quantile(w, seq_len(0)), 0)
  expect_length(quantile(w, c(0, 1)), 2)
})
