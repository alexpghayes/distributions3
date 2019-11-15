context("test-Laplace")

test_that("print.Laplace works", {
  expect_output(print(Laplace()), regexp = "Laplace distribution")
})

theta <- 1
phi <- 2
L <- Laplace(theta, phi)

test_that("random.Laplace works correctly", {
  expect_length(random(L), 1)
  expect_length(random(L, 100), 100)
  expect_length(random(L, 0), 0)
  expect_error(random(L, -2))
})

test_that("pdf.Laplace works correctly", {
  expect_equal(pdf(L, 1), 1 / (2 * phi))
  expect_equal(pdf(L, phi + theta), exp(-1) / (2 * phi))
  expect_equal(pdf(L, -Inf), 0)
  expect_equal(pdf(L, Inf), 0)

  expect_length(pdf(L, seq_len(0)), 0)
  expect_length(pdf(L, seq_len(1)), 1)
  expect_length(pdf(L, seq_len(10)), 10)
})

test_that("log_pdf.Laplace works correctly", {
  expect_equal(log_pdf(L, 1), -log(2 * phi))
  expect_equal(log_pdf(L, phi + theta), -1 - log(2 * phi))
  expect_equal(log_pdf(L, -Inf), -Inf)
  expect_equal(log_pdf(L, Inf), -Inf)

  expect_length(log_pdf(L, seq_len(0)), 0)
  expect_length(log_pdf(L, seq_len(1)), 1)
  expect_length(log_pdf(L, seq_len(10)), 10)
})

test_that("cdf.Laplace works correctly", {
  expect_equal(cdf(L, -Inf), 0)
  expect_equal(cdf(L, theta), 0.5)
  expect_equal(cdf(L, Inf), 1)

  expect_length(cdf(L, seq_len(0)), 0)
  expect_length(cdf(L, seq_len(1)), 1)
  expect_length(cdf(L, seq_len(10)), 10)
})

test_that("quantile.Laplace works correctly", {
  expect_equal(quantile(L, 0), -Inf)
  expect_equal(quantile(L, 0.5), theta)
  expect_equal(quantile(L, 1), Inf)

  expect_length(quantile(L, seq_len(0)), 0)
  expect_length(quantile(L, c(0, 1)), 2)
})

test_that("cdf.Laplace and quantile.Laplace are consistent", {
  pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)
  expect_equal(cdf(L, quantile(L, pvec)), pvec)
})
