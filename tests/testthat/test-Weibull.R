context("test-Weibull")

test_that("print.Weibull works", {
  expect_output(print(Weibull(1, 1)), regexp = "Weibull distribution")
})

test_that("likelihood.Weibull and log_likelihood.Weibull work correctly", {
  w <- Weibull(1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(w, 1), 1 / exp(1))
  expect_equal(likelihood(w, x), (1 / exp(1))^2)

  expect_equal(log_likelihood(w, 1), log(1 / exp(1)))
  expect_equal(log_likelihood(w, x), log((1 / exp(1))^2))
})

test_that("random.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_length(random(w), 1)
  expect_length(random(w, 100), 100)
  expect_length(random(w, 0), 0)
  expect_error(random(w, -2))
})

test_that("pdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(pdf(w, 0), 1)
  expect_equal(pdf(w, 1), 1 / exp(1))
  expect_equal(pdf(w, -12), 0)

  expect_length(pdf(w, seq_len(0)), 0)
  expect_length(pdf(w, seq_len(1)), 1)
  expect_length(pdf(w, seq_len(10)), 10)
})

test_that("log_pdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(log_pdf(w, 0), 0)
  expect_equal(log_pdf(w, 1), log(1 / exp(1)))
  expect_equal(log_pdf(w, -12), -Inf)

  expect_length(log_pdf(w, seq_len(0)), 0)
  expect_length(log_pdf(w, seq_len(1)), 1)
  expect_length(log_pdf(w, seq_len(10)), 10)
})

test_that("cdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(cdf(w, 0), 0)
  expect_equal(cdf(w, 1), 1 - 1 / exp(1))


  expect_length(cdf(w, seq_len(0)), 0)
  expect_length(cdf(w, seq_len(1)), 1)
  expect_length(cdf(w, seq_len(10)), 10)
})

test_that("quantile.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(quantile(w, 0), 0)
  expect_equal(quantile(w, 1), Inf)


  expect_length(quantile(w, seq_len(0)), 0)
  expect_length(quantile(w, c(0, 1)), 2)
})

test_that("vectorization of a Weibull distribution work correctly", {
  d <- Weibull(c(1, 0.3), c(1, 2))
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

