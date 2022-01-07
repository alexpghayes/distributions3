context("test-FisherF")

test_that("print.FisherF works", {
  expect_output(print(FisherF(1, 1)), regexp = "Fisher's F distribution")
})

test_that("likelihood.FisherF and log_likelihood.FisherF work correctly", {
  s <- FisherF(1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(s, 1), df(1, 1, 1, 0))
  expect_equal(likelihood(s, x), df(1, 1, 1, 0) * df(1, 1, 1, 0) * df(0, 1, 1, 0))

  expect_equal(log_likelihood(s, 1), log(df(1, 1, 1, 0)))
  expect_equal(log_likelihood(s, x), log(df(1, 1, 1, 0) * df(1, 1, 1, 0) * df(0, 1, 1, 0)))
})

test_that("random.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_length(random(s), 1)
  expect_length(random(s, 100), 100)
  expect_length(random(s, 0), 0)
  expect_error(random(s, -2))
})

test_that("pdf.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(pdf(s, 0), Inf)
  expect_equal(pdf(s, 1), df(1, 1, 1, 0))
  expect_equal(pdf(s, -12), 0)

  expect_length(pdf(s, seq_len(0)), 0)
  expect_length(pdf(s, seq_len(1)), 1)
  expect_length(pdf(s, seq_len(10)), 10)
})

test_that("log_pdf.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(log_pdf(s, 0), Inf)
  expect_equal(log_pdf(s, 1), log(df(1, 1, 1, 0)))
  expect_equal(log_pdf(s, -12), -Inf)

  expect_length(log_pdf(s, seq_len(0)), 0)
  expect_length(log_pdf(s, seq_len(1)), 1)
  expect_length(log_pdf(s, seq_len(10)), 10)
})

test_that("cdf.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(cdf(s, 0), pf(0, 1, 1, 0))
  expect_equal(cdf(s, 1), pf(1, 1, 1, 0))


  expect_length(cdf(s, seq_len(0)), 0)
  expect_length(cdf(s, seq_len(1)), 1)
  expect_length(cdf(s, seq_len(10)), 10)
})

test_that("quantile.FisherF work correctly", {
  s <- FisherF(1, 1)

  expect_equal(quantile(s, 0), qf(0, 1, 1, 0))
  expect_equal(quantile(s, 1), qf(1, 1, 1, 0))


  expect_length(quantile(s, seq_len(0)), 0)
  expect_length(quantile(s, c(0, 1)), 2)
})

test_that("vectorization of a FisherF distribution work correctly", {
  d <- FisherF(c(5, 10), c(5, 10))
  d1 <- d[1]
  d2 <- d[2]

  expect_equal(mean(d), c(mean(d1), mean(d2)))
  expect_equal(variance(d), c(variance(d1), variance(d2)))
  expect_equal(skewness(d), c(skewness(d1), skewness(d2)))
  expect_equal(kurtosis(d), c(kurtosis(d1), kurtosis(d2)))

  set.seed(123); r1 <- random(d)
  set.seed(123); r2 <- rf(2, df1 = c(5, 10), df2 = c(5, 10), ncp = 0)
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

