context("test-Uniform")

test_that("print.Uniform works", {
  expect_output(print(Uniform(1, 1)), regexp = "Uniform distribution")
})

test_that("likelihood.Uniform and log_likelihood.Uniform work correctly", {
  u <- Uniform()
  x <- c(1, 1, 0)

  expect_equal(likelihood(u, 1), dunif(1))
  expect_equal(likelihood(u, x), dunif(1) * dunif(1) * dunif(0))

  expect_equal(log_likelihood(u, 1), log(dunif(1)))
  expect_equal(log_likelihood(u, x), log(dunif(1) * dunif(1) * dunif(0)))
})

test_that("random.Uniform work correctly", {
  u <- Uniform()

  expect_length(random(u), 1)
  expect_length(random(u, 100), 100)
  expect_length(random(u, 0), 0)
  expect_error(random(u, -2))
})

test_that("pdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(pdf(u, 0), dunif(0, 0, 1))
  expect_equal(pdf(u, 1), dunif(1, 0, 1))

  expect_length(pdf(u, seq_len(0)), 0)
  expect_length(pdf(u, seq_len(1)), 1)
  expect_length(pdf(u, seq_len(10)), 10)
})

test_that("log_pdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(log_pdf(u, 0), log(dunif(0, 0, 1)))
  expect_equal(log_pdf(u, 1), log(dunif(1, 0, 1)))

  expect_length(log_pdf(u, seq_len(0)), 0)
  expect_length(log_pdf(u, seq_len(1)), 1)
  expect_length(log_pdf(u, seq_len(10)), 10)
})

test_that("cdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(cdf(u, 0), punif(0, 0, 1))
  expect_equal(cdf(u, 1), punif(1, 0, 1))


  expect_length(cdf(u, seq_len(0)), 0)
  expect_length(cdf(u, seq_len(1)), 1)
  expect_length(cdf(u, seq_len(10)), 10)
})

test_that("quantile.Uniform work correctly", {
  u <- Uniform()

  expect_equal(quantile(u, 0), qunif(0, 0, 1))
  expect_equal(quantile(u, 1), qunif(1, 0, 1))


  expect_length(quantile(u, seq_len(0)), 0)
  expect_length(quantile(u, c(0, 1)), 2)
})

test_that("{moments}.Uniform work correctly", {
  u <- Uniform()

  expect_equal(mean(u), 0.5)
  expect_equal(variance(u), 1/12)
  expect_equal(skewness(u), 0)
  expect_equal(kurtosis(u), -6/5)
})

test_that("vectorization of a Uniform distribution work correctly", {
  d <- Uniform(c(0, 10), c(1, 20))
  d1 <- d[1]
  d2 <- d[2]

  expect_equal(mean(d), c(mean(d1), mean(d2)))
  expect_equal(variance(d), c(variance(d1), variance(d2)))
  expect_equal(skewness(d), c(skewness(d1), skewness(d2)))
  expect_equal(kurtosis(d), c(kurtosis(d1), kurtosis(d2)))

  set.seed(123); r1 <- random(d)
  set.seed(123); r2 <- c(random(d1), random(d2))
  expect_equal(r1, r2)

  #set.seed(123); r3 <- random(d, 10)
  #set.seed(123); r4 <- c(random(d1, 10), random(d2, 10))
  #expect_equal(unname(r3[1, ]), r4[1:10])
  #expect_equal(unname(r3[2, ]), r4[11:20])

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

