context("test-Poisson")

test_that("print.Poisson works", {
  expect_output(print(Poisson(1)), regexp = "Poisson distribution")
})

test_that("suff_stat.Poisson works correctly", {
  ss <- list(sum = 1, samples = 1)
  expect_equal(suff_stat(Poisson(1), 1), ss)

  ss <- list(sum = 5050, samples = 101)
  expect_equal(suff_stat(Poisson(1), 0:100), ss)

  expect_error(suff_stat(Poisson(1), 0.5))

  expect_error(suff_stat(Poisson(1), -1))
})

test_that("fit_mle.Poisson works correctly", {
  expect_equal(fit_mle(Poisson(1), 2), Poisson(2))

  expect_equal(fit_mle(Poisson(1), 0:100), Poisson(50))

  expect_error(fit_mle(Poisson(1), -1))

  expect_error(fit_mle(Poisson(1), 0.5))
})

test_that("random.Poisson work correctly", {
  p <- Poisson(1)

  expect_length(random(p), 1)
  expect_length(random(p, 100), 100)
  expect_length(random(p, 0), 0)
  expect_error(random(p, -2))
})

test_that("pdf.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(pdf(p, 0), 1 / exp(1))
  expect_equal(pdf(p, 1), 1 / exp(1))
  expect_equal(pdf(p, -12), 0)

  expect_warning(pdf(p, 0.5))

  expect_length(pdf(p, seq_len(0)), 0)
  expect_length(pdf(p, seq_len(1)), 1)
  expect_length(pdf(p, seq_len(10)), 10)
})

test_that("log_pdf.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(log_pdf(p, 0), -1)
  expect_equal(log_pdf(p, 1), -1)
  expect_equal(log_pdf(p, -12), -Inf)

  expect_warning(log_pdf(p, 0.5))

  expect_length(log_pdf(p, seq_len(0)), 0)
  expect_length(log_pdf(p, seq_len(1)), 1)
  expect_length(log_pdf(p, seq_len(10)), 10)
})

test_that("cdf.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(cdf(p, 0), 1 / exp(1))
  expect_equal(cdf(p, 1), 1 / exp(1) * 2)


  expect_length(cdf(p, seq_len(0)), 0)
  expect_length(cdf(p, seq_len(1)), 1)
  expect_length(cdf(p, seq_len(10)), 10)
})

test_that("quantile.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(quantile(p, 0), 0)
  expect_equal(quantile(p, 0.5), 1)


  expect_length(quantile(p, seq_len(0)), 0)
  expect_length(quantile(p, c(0, 1)), 2)
})
