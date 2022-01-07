context("test-Exponential")

test_that("fit_mle.Exponential works correctly", {
  expect_equal(fit_mle(Exponential(), 1), Exponential(1))

  expect_error(fit_mle(Exponential(), -1))

  expect_true(is.numeric(fit_mle(Exponential(), rexp(100))$rate))
})

test_that("print.Beta works", {
  expect_output(print(Exponential()), regexp = "Exponential distribution")
})

test_that("random.Exponential work correctly", {
  e <- Exponential()

  expect_length(random(e), 1)
  expect_length(random(e, 100), 100)
  expect_length(random(e, 0), 0)
  expect_error(random(e, -2))
})

test_that("pdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(pdf(e, 0), 1)
  expect_equal(pdf(e, 1), 1 / exp(1))
  expect_equal(pdf(e, -12), 0)

  expect_length(pdf(e, seq_len(0)), 0)
  expect_length(pdf(e, seq_len(1)), 1)
  expect_length(pdf(e, seq_len(10)), 10)
})

test_that("pdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(log_pdf(e, 0), log(1))
  expect_equal(log_pdf(e, 1), log(1 / exp(1)))
  expect_equal(log_pdf(e, -12), log(0))

  expect_length(log_pdf(e, seq_len(0)), 0)
  expect_length(log_pdf(e, seq_len(1)), 1)
  expect_length(log_pdf(e, seq_len(10)), 10)
})

test_that("cdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(cdf(e, 0), 0)
  expect_equal(cdf(e, 1), 1 - 1 / exp(1))


  expect_length(cdf(e, seq_len(0)), 0)
  expect_length(cdf(e, seq_len(1)), 1)
  expect_length(cdf(e, seq_len(10)), 10)
})

test_that("quantile.Exponential work correctly", {
  e <- Exponential()

  expect_equal(quantile(e, 0), 0)
  expect_equal(quantile(e, 1), Inf)


  expect_length(quantile(e, seq_len(0)), 0)
  expect_length(quantile(e, c(0, 1)), 2)
})

test_that("{moments}.Exponential work correctly", {
  e <- Exponential()

  expect_equal(mean(e), 1)
  expect_equal(variance(e), 1)
  expect_equal(skewness(e), 2)
  expect_equal(kurtosis(e), 6)
})

test_that("vectorization of a Exponential distribution work correctly", {
  d <- Exponential(c(1, 2))
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

