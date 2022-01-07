context("test-HyperGeometric")

test_that("HyperGeometric works as intended when k > n + m", {
  expect_error(HyperGeometric(1,1,3))
})

test_that("print.HyperGeometric works", {
  expect_output(print(HyperGeometric(1, 1, 1)), regexp = "HyperGeometric distribution")
})


test_that("likelihood.HyperGeometric and log_likelihood.HyperGeometric work correctly", {
  h <- HyperGeometric(1, 1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(h, 1), 0.5)
  expect_equal(likelihood(h, x), 0.5^3)

  expect_equal(log_likelihood(h, 1), log(0.5))
  expect_equal(log_likelihood(h, x), log(0.5) * 3)
})

test_that("random.HyperGeometric work correctly", {
  h <- HyperGeometric(1, 1, 1)

  expect_length(random(h), 1)
  expect_length(random(h, 100), 100)
  expect_length(random(h, 0), 0)
  expect_error(random(h, -2))
})

test_that("pdf.HyperGeometric work correctly", {
  h <- HyperGeometric(1, 1, 1)

  expect_equal(pdf(h, 0), 0.5)
  expect_equal(pdf(h, 1), 0.5)
  expect_equal(pdf(h, -12), 0)

  expect_warning(pdf(h, 0.5))

  expect_length(pdf(h, seq_len(0)), 0)
  expect_length(pdf(h, seq_len(1)), 1)
  expect_length(pdf(h, seq_len(10)), 10)
})

test_that("cdf.HyperGeometric work correctly", {
  h <- HyperGeometric(1, 1, 1)

  expect_equal(cdf(h, 0), 0.5)
  expect_equal(cdf(h, 1), 1)


  expect_length(cdf(h, seq_len(0)), 0)
  expect_length(cdf(h, seq_len(1)), 1)
  expect_length(cdf(h, seq_len(10)), 10)
})

test_that("quantile.HyperGeometric work correctly", {
  h <- HyperGeometric(1, 1, 1)

  expect_equal(quantile(h, 0), 0)
  expect_equal(quantile(h, 1), 1)


  expect_length(quantile(h, seq_len(0)), 0)
  expect_length(quantile(h, c(0, 1)), 2)
})

test_that("vectorization of a HyperGeometric distribution work correctly", {
  d <- HyperGeometric(1, c(1, 3), 1)
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

