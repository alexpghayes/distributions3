context("test-Normal")

test_that("print.Normal works", {
  expect_output(print(Normal()), regexp = "Normal distribution")
})


test_that("suff_stat.Normal works correctly", {
  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(Normal(), c(0, 0)), ss)

  expect_error(suff_stat(Normal(), "abc"))
})

test_that("fit_mle.Normal works correctly", {
  expect_equal(fit_mle(Normal(), c(0, 0)), Normal(0, 0))
})

test_that("random.Normal work correctly", {
  n <- Normal()

  expect_length(random(n), 1)
  expect_length(random(n, 100), 100)
  expect_length(random(n, 0), 0)
  expect_error(random(n, -2))
})


test_that("pdf.Normal work correctly", {
  n <- Normal()

  expect_equal(pdf(n, 0), dnorm(0, 0, 1))
  expect_equal(pdf(n, 1), dnorm(1, 0, 1))

  expect_length(pdf(n, seq_len(0)), 0)
  expect_length(pdf(n, seq_len(1)), 1)
  expect_length(pdf(n, seq_len(10)), 10)
})

test_that("log_pdf.Normal work correctly", {
  n <- Normal()

  expect_equal(log_pdf(n, 0), log(dnorm(0, 0, 1)))
  expect_equal(log_pdf(n, 1), log(dnorm(1, 0, 1)))

  expect_length(log_pdf(n, seq_len(0)), 0)
  expect_length(log_pdf(n, seq_len(1)), 1)
  expect_length(log_pdf(n, seq_len(10)), 10)
})


test_that("cdf.Normal work correctly", {
  n <- Normal()

  expect_equal(cdf(n, 0), 0.5)

  expect_length(cdf(n, seq_len(0)), 0)
  expect_length(cdf(n, seq_len(1)), 1)
  expect_length(cdf(n, seq_len(10)), 10)
})

test_that("quantile.Normal work correctly", {
  n <- Normal()

  expect_equal(quantile(n, 0), -Inf)
  expect_equal(quantile(n, 0.5), 0)
  expect_equal(quantile(n, 1), Inf)


  expect_length(quantile(n, seq_len(0)), 0)
  expect_length(quantile(n, c(0, 1)), 2)
})

test_that("{moments}.Normal work correctly", {
  n <- Normal()

  expect_equal(mean(n), 0)
  expect_equal(variance(n), 1)
  expect_equal(skewness(n), 0)
  expect_equal(kurtosis(n), 0)
})

test_that("Normal ops work correctly", {
  x <- Normal(1,1)
  y <- Normal(2,sqrt(3))

  s <- x + y

  expect_equal(s$mu, 3)
  expect_equal(s$sigma, 2)

  d <- x - y

  expect_equal(d$mu, -1)
  expect_equal(d$sigma, 2)

})

