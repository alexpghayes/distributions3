context("test-binomial")

test_that("fit.binomial works correctly", {

  # only one trial
  expect_equal(fit(binomial(1), c(0,1)), binomial(1, 0.5))

  # many trials
  x <- c(0, rep(100, 99))
  expect_equal(fit(binomial(100), x), binomial(100, 0.99))

  # raises error when data has negative counts
  expect_error(fit(binomial(1), -1))

  # raises error when data has counts greater than binomial's size
  expect_error(fit(binomial(1), 2))
})

test_that("suff_stat.binomial works correctly", {

  ss_1 <- list(successes = 1, experiments = 2, trials = 1)
  expect_equal(suff_stat(binomial(1), c(0,1)), ss_1)

  ss_2 <- list(successes = 9, experiments = 3, trials = 5)
  expect_equal(suff_stat(binomial(5), c(3,3,3)), ss_2)
})

test_that("likelihood.binomial and log_likelihood.binomial work correctly", {

  b <- binomial(size = 10, p = 0.1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(b, 1), dbinom(1, 10, 0.1))
  expect_equal(likelihood(b, x), prod(dbinom(x, 10, 0.1)))

  expect_equal(log_likelihood(b, 1), dbinom(1, 10, 0.1, log = TRUE))
  expect_equal(log_likelihood(b, x), sum(dbinom(x, 10, 0.1, log = TRUE)))

})
