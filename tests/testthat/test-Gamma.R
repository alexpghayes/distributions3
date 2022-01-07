context("test-Gamma")

test_that("print.Gamma works", {
  expect_output(print(Gamma(1, 1)), regexp = "Gamma distribution")
})

test_that("likelihood.Gamma and log_likelihood.Gamma work correctly", {
  cau <- Gamma(1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dgamma(1, 1, 1))
  expect_equal(likelihood(cau, x), dgamma(1, 1, 1) * dgamma(1, 1, 1) * dgamma(0, 1, 1))

  expect_equal(log_likelihood(cau, 1), log(dgamma(1, 1, 1)))
  expect_equal(log_likelihood(cau, x), log(dgamma(1, 1, 1) * dgamma(1, 1, 1) * dgamma(0, 1, 1)))
})

test_that("random.Gamma work correctly", {
  cau <- Gamma(1, 1)

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(pdf(cau, 0), dgamma(0, 0, 1))
  expect_equal(pdf(cau, 1), dgamma(1, 0, 1))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(log_pdf(cau, 0), log(dgamma(0, 0, 1)))
  expect_equal(log_pdf(cau, 1), log(dgamma(1, 0, 1)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(cdf(cau, 0), pgamma(0, 0, 1))
  expect_equal(cdf(cau, 1), pgamma(1, 0, 1))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.Gamma work correctly", {
  cau <- Gamma(0, 1)

  expect_equal(quantile(cau, 0), qgamma(0, 0, 1))
  expect_equal(quantile(cau, 1), qgamma(1, 0, 1))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})

test_that("vectorization of a Gamma distribution work correctly", {
  d <- Gamma(c(0, 10), c(1, 1))
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
