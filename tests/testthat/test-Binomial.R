context("test-Binomial")

test_that("print.Binomial works", {
  expect_output(print(Binomial(1)), regexp = "Binomial distribution")
})

test_that("fit_mle.Binomial works correctly", {

  # only one trial
  expect_equal(fit_mle(Binomial(1), c(0, 1)), Binomial(1, 0.5))

  # many trials
  x <- c(0, rep(100, 99))
  expect_equal(fit_mle(Binomial(100), x), Binomial(100, 0.99))

  # raises error when data has negative counts
  expect_error(fit_mle(Binomial(1), -1))

  # raises error when data has counts greater than Binomial's size
  expect_error(fit_mle(Binomial(1), 2))
})

test_that("suff_stat.Binomial works correctly", {
  ss_1 <- list(successes = 1, experiments = 2, trials = 1)
  expect_equal(suff_stat(Binomial(1), c(0, 1)), ss_1)

  ss_2 <- list(successes = 9, experiments = 3, trials = 5)
  expect_equal(suff_stat(Binomial(5), c(3, 3, 3)), ss_2)
})

test_that("likelihood.Binomial and log_likelihood.Binomial work correctly", {
  b <- Binomial(size = 10, p = 0.1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(b, 1), dbinom(1, 10, 0.1))
  expect_equal(likelihood(b, x), prod(dbinom(x, 10, 0.1)))

  expect_equal(log_likelihood(b, 1), dbinom(1, 10, 0.1, log = TRUE))
  expect_equal(log_likelihood(b, x), sum(dbinom(x, 10, 0.1, log = TRUE)))
})

test_that("random.Binomial work correctly", {
  b <- Binomial(size = 10, p = 0.1)

  expect_length(random(b), 1)
  expect_length(random(b, 100), 100)
  expect_length(random(b, 0), 0)
  expect_error(random(b, -2))
})


test_that("pdf.Bernoulli work correctly", {
  b <- Binomial(size = 2, p = 0.1)

  expect_equal(pdf(b, 0), 0.9^2)
  expect_equal(pdf(b, 2), 0.1^2)
  expect_equal(pdf(b, -12), 0)

  expect_warning(pdf(b, 0.5))

  expect_length(pdf(b, seq_len(0)), 0)
  expect_length(pdf(b, seq_len(1)), 1)
  expect_length(pdf(b, seq_len(10)), 10)
})

test_that("cdf.Bernoulli work correctly", {
  b <- Binomial(size = 2, p = 0.1)

  expect_equal(cdf(b, 0), 0.9^2)
  expect_equal(cdf(b, 2), 1)


  expect_length(cdf(b, seq_len(0)), 0)
  expect_length(cdf(b, seq_len(1)), 1)
  expect_length(cdf(b, seq_len(10)), 10)
})

test_that("quantile.Bernoulli work correctly", {
  b <- Binomial(size = 2, p = 0.1)

  expect_equal(quantile(b, 0), 0)
  expect_equal(quantile(b, 1), 2)


  expect_length(quantile(b, seq_len(0)), 0)
  expect_length(quantile(b, c(0, 1)), 2)
})
