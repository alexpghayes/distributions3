context("test-Beta")

test_that("print.Beta works", {
  expect_output(print(Beta()), regexp = "Beta distribution")
})

test_that("random.Beta work correctly", {
  b <- Beta(1, 5)

  expect_length(random(b), 1)
  expect_length(random(b, 100), 100)
  expect_length(random(b, 0), 0)
  expect_error(random(b, -2))
})

test_that("pdf.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(pdf(b, 0), 2)
  expect_equal(pdf(b, 0.5), 1)
  expect_equal(pdf(b, -12), 0)

  expect_length(pdf(b, seq_len(0)), 0)
  expect_length(pdf(b, seq_len(1)), 1)
  expect_length(pdf(b, seq_len(10)), 10)
})

test_that("log_pdf.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(log_pdf(b, 0), log(2))
  expect_equal(log_pdf(b, 0.5), log(1))
  expect_equal(log_pdf(b, -12), log(0))

  expect_length(log_pdf(b, seq_len(0)), 0)
  expect_length(log_pdf(b, seq_len(1)), 1)
  expect_length(log_pdf(b, seq_len(10)), 10)
})

test_that("cdf.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(cdf(b, 0), 0)
  expect_equal(cdf(b, 0.5), 0.75)
  expect_equal(cdf(b, 1), 1)


  expect_length(cdf(b, seq_len(0)), 0)
  expect_length(cdf(b, seq_len(1)), 1)
  expect_length(cdf(b, seq_len(10)), 10)
})

test_that("quantile.Beta work correctly", {
  b <- Beta(1, 2)

  expect_equal(quantile(b, 0), 0)
  expect_equal(quantile(b, 0.75), 0.5)
  expect_equal(quantile(b, 1), 1)


  expect_length(quantile(b, seq_len(0)), 0)
  expect_length(quantile(b, c(0, 1)), 2)
})

test_that("{moments}.Beta work correctly", {
  b <- Beta(1, 1)

  expect_equal(mean(b), 0.5)
  expect_equal(variance(b), 1/12)
  expect_equal(skewness(b), 0)
  expect_equal(kurtosis(b), -6/5)
})

test_that("vectorization of a Beta distribution work correctly", {
  d <- Beta(c(1, 1), c(2, 5))
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

