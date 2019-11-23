context("test-Pareto")

test_that("print.Pareto works", {
  expect_output(print(Pareto()), regexp = "Pareto distribution")
})

k <- 1
a <- 2
P <- Pareto(k, a)

set.seed(27)

test_that("suff_stat.Pareto works correctly", {
  n <- 10
  x <- random(P, n)
  ss <- list(minimum = min(x), mean_log = mean(log(x)), samples = n)
  expect_equal(suff_stat(Pareto(), x), ss)

  expect_error(suff_stat(Pareto(), "abc"))
})

test_that("fit_mle.Pareto works correctly", {
  n <- 10
  x <- random(P, n)
  khat <- min(x)
  ahat <- n / sum(log(x / khat))
  expect_equal(fit_mle(Pareto(), x), Pareto(khat, ahat))
})

test_that("random.Pareto work correctly", {
  expect_length(random(P), 1)
  expect_length(random(P, 100), 100)
  expect_length(random(P, 0), 0)
  expect_error(random(P, -2))
})

test_that("pdf.Pareto work correctly", {
  expect_equal(pdf(P, k), a / k)
  expect_equal(pdf(P, k * 2),  a / k / 2 ^ (a + 1))
  expect_equal(pdf(P, k - 0.1), 0)
  expect_equal(pdf(P, -12), 0)

  expect_length(pdf(P, seq_len(0)), 0)
  expect_length(pdf(P, seq_len(1)), 1)
  expect_length(pdf(P, seq_len(10)), 10)
})

test_that("log_pdf.Pareto work correctly", {
  expect_equal(log_pdf(P, k), log(a) - log(k))
  expect_equal(log_pdf(P, k * 2), log(a) - log(k) - (a + 1) * log(2))
  expect_equal(log_pdf(P, k - 0.1), log(0))
  expect_equal(log_pdf(P, -12), log(0))

  expect_length(log_pdf(P, seq_len(0)), 0)
  expect_length(log_pdf(P, seq_len(1)), 1)
  expect_length(log_pdf(P, seq_len(10)), 10)
})

test_that("cdf.Pareto work correctly", {
  expect_equal(cdf(P, k), 0)
  expect_equal(cdf(P, 2 * k), 1 - 0.5 ^ a)
  expect_equal(cdf(P, Inf), 1)

  expect_length(cdf(P, seq_len(0)), 0)
  expect_length(cdf(P, seq_len(1)), 1)
  expect_length(cdf(P, seq_len(10)), 10)
})

test_that("quantile.Pareto work correctly", {
  expect_equal(quantile(P, 0), k)
  expect_equal(quantile(P, 0.5), k * 0.5 ^ (-1 / a))
  expect_equal(quantile(P, 1), Inf)

  expect_length(quantile(P, seq_len(0)), 0)
  expect_length(quantile(P, c(0, 1)), 2)

  expect_error(quantile(L, -0.1))
  expect_error(quantile(L, -1.1))
})

test_that("cdf.Pareto and quantile.Pareto are consistent", {
  pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)
  expect_equal(cdf(P, quantile(P, pvec)), pvec)
})
