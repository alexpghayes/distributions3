context("test-binomial")

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
