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

test_that("vectorization of a Normal distribution work correctly", {
  d <- Normal(c(0, 10), c(1, 1))
  d1 <- d[1]
  d2 <- d[2]

  expect_equal(mean(d), c(mean(d1), mean(d2)))
  expect_equal(variance(d), c(variance(d1), variance(d2)))
  expect_equal(skewness(d), c(skewness(d1), skewness(d2)))
  expect_equal(kurtosis(d), c(kurtosis(d1), kurtosis(d2)))

  set.seed(123); r1 <- random(d)
  set.seed(123); r2 <- c(random(d1), random(d2))
  expect_equal(r1, r2)

  expect_equal(pdf(d, 0), c(pdf(d1, 0), pdf(d2, 0)))
  expect_equal(log_pdf(d, 0), c(log_pdf(d1, 0), log_pdf(d2, 0)))
  expect_equal(cdf(d, 0.5), c(cdf(d1, 0.5), cdf(d2, 0.5)))

  expect_equal(quantile(d, 0.5), c(quantile(d1, 0.5), quantile(d2, 0.5)))
  expect_equal(quantile(d, c(0.5, 0.5)), c(quantile(d1, 0.5), quantile(d2, 0.5)))
  expect_equal(
    quantile(d, c(0.1, 0.5, 0.9)),
    matrix(
      c(quantile(d1, c(0.1, 0.5, 0.9)), quantile(d2, c(0.1, 0.5, 0.9))),
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("q_0.1", "q_0.5", "q_0.9"))
    )
  )

  expect_equal(
    support(d),
    matrix(
      c(support(d1), support(d2)),
      nrow = 2,
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("min", "max"))
    )
  )
})
