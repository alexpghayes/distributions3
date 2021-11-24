context("test-Erlang")
e <- Erlang(lambda = 0.5, k = 3)

test_that("print.Erlang works", {
  expect_output(print(e), regexp = "Erlang distribution")
})

test_that("random.Erlang works correctly", {
  expect_length(random(e), 1)
  expect_length(random(e, 100), 100)
  expect_length(random(e, 0), 0)
  expect_error(random(e, -2))
})

test_that("pdf.Erlang works correctly", {
  expect_length(pdf(e, seq_len(0)), 0)
  expect_length(pdf(e, seq_len(1)), 1)
  expect_length(pdf(e, seq_len(10)), 10)
  expect_error(pdf(e, -42))
})

test_that("log_pdf.Erlang works correctly", {
  expect_length(log_pdf(e, seq_len(0)), 0)
  expect_length(log_pdf(e, seq_len(1)), 1)
  expect_length(log_pdf(e, seq_len(10)), 10)
  expect_error(log_pdf(e, -42))
})

test_that("cdf.Erlang works correctly", {
  expect_equal(cdf(e, 0), 0)
  expect_length(cdf(e, seq_len(0)), 0)
  expect_length(cdf(e, seq_len(1)), 1)
  expect_length(cdf(e, seq_len(10)), 10)
  expect_error(cdf(e, -42))
})

test_that("quantile.Erlang works correctly", {
  expect_equal(quantile(e, 0), 0)
  expect_length(quantile(e, seq_len(1)), 1)
  expect_length(quantile(e, seq(0.1, 0.9, by = 0.1)), 9)
  expect_error(quantile(e, p = -42))
  expect_error(quantile(e, p = 42))
  expect_error(quantile(e, p = 0.5, interval = 42))
  expect_error(quantile(e, p = 0.5, interval = "Hi"))
  expect_error(quantile(e, p = 0.5, tol = "Hi"))
})

test_that("support.Erlang works correctly", {
  expect_equal(support(e), c(0, Inf))
})
