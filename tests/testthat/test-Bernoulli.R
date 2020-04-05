context("test-Bernoulli")

test_that("print.Bernoulli works", {
  expect_output(print(Bernoulli()), regexp = "Bernoulli distribution")
})

test_that("fit_mle.Bernoulli works correctly", {
  expect_equal(fit_mle(Bernoulli(), c(0, 1)), Bernoulli(0.5))
})

test_that("suff_stats.Bernoulli works correctly", {
  ss <- list(successes = 3, failures = 2)

  expect_equal(suff_stat(Bernoulli(), c(1, 1, 1, 0, 0)), ss)

  expect_error(suff_stat(Bernoulli(), 2))

  expect_error(suff_stat(Bernoulli(), -1))
})


test_that("likelihood.Bernoulli and log_likelihood.Bernoulli work correctly", {
  b <- Bernoulli(0.1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(b, 1), 0.1)
  expect_equal(likelihood(b, x), 0.009)

  expect_equal(log_likelihood(b, 1), log(0.1))
  expect_equal(log_likelihood(b, x), log(0.009))
})

test_that("random.Bernoulli work correctly", {
  b <- Bernoulli()

  expect_length(random(b), 1)
  expect_length(random(b, 100), 100)
  expect_length(random(b, 0), 0)
  expect_error(random(b, -2))
})

test_that("pdf.Bernoulli work correctly", {
  b <- Bernoulli(0.1)

  expect_equal(pdf(b, 0), 0.9)
  expect_equal(pdf(b, 1), 0.1)
  expect_equal(pdf(b, -12), 0)

  expect_warning(pdf(b, 0.5))

  expect_length(pdf(b, seq_len(0)), 0)
  expect_length(pdf(b, seq_len(1)), 1)
  expect_length(pdf(b, seq_len(10)), 10)
})

test_that("cdf.Bernoulli work correctly", {
  b <- Bernoulli(0.1)

  expect_equal(cdf(b, 0), 0.9)
  expect_equal(cdf(b, 1), 1)


  expect_length(cdf(b, seq_len(0)), 0)
  expect_length(cdf(b, seq_len(1)), 1)
  expect_length(cdf(b, seq_len(10)), 10)
})

test_that("quantile.Bernoulli work correctly", {
  b <- Bernoulli(0.1)

  expect_equal(quantile(b, 0), 0)
  expect_equal(quantile(b, 1), 1)


  expect_length(quantile(b, seq_len(0)), 0)
  expect_length(quantile(b, c(0, 1)), 2)
})

test_that("{moments}.Bernoulli work correctly", {
  n <- Bernoulli()

  expect_equal(mean(n), 0.5)
  expect_equal(variance(n), 0.25)
  expect_equal(skewness(n), 0)
  expect_equal(kurtosis(n), -2)
})
