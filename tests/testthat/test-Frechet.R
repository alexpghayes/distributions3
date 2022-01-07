context("test-Frechet")

test_that("print.Frechet works", {
  expect_output(print(Frechet()), regexp = "Frechet distribution")
})

## Example distributions

# GEV parameter values

mu <- 0
sigma <- 1
xi1 <- 0.1

# Equivalent Frechet parameter values

alpha <- 1 / xi1
s <- alpha * sigma
m <- mu - s

g1 <- Frechet(m, s, alpha)
low <- m

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x1 <- c(low, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.Frechet works correctly", {
  expect_length(random(g1), 1)
  expect_length(random(g1, 100), 100)
  expect_length(random(g1, 0), 0)
  expect_error(random(g1, -2))
})

test_that("pdf.Frechet works correctly", {
  p <- pvec[2:4]
  expect_equal(pdf(g1, x1), c(0, 0, exp(-1), 0, NA))
  expect_equal(pdf(g1, quantile(g1, p)), (-log(p)) ^ (1 + xi1) * p)
  expect_length(pdf(g1, seq_len(0)), 0)
  expect_length(pdf(g1, seq_len(1)), 1)
  expect_length(pdf(g1, seq_len(10)), 10)
})

test_that("log_pdf.Frechet works correctly", {
  expect_equal(log_pdf(g1, x1), c(-Inf, -Inf, -1, -Inf, NA))
  expect_length(log_pdf(g1, seq_len(0)), 0)
  expect_length(log_pdf(g1, seq_len(1)), 1)
  expect_length(log_pdf(g1, seq_len(10)), 10)
})

test_that("cdf.Frechet works correctly", {
  expect_equal(cdf(g1, x1), c(0, 0, exp(-1), 1, NA))
  expect_length(cdf(g1, seq_len(0)), 0)
  expect_length(cdf(g1, seq_len(1)), 1)
  expect_length(cdf(g1, seq_len(10)), 10)
})

test_that("quantile.Frechet works correctly", {
  q1 <- ((-log(pvec[2:4])) ^ (-xi1) - 1) / xi1
  expect_equal(quantile(g1, pvec), c(low, q1, Inf, NA))
  expect_length(quantile(g1, seq_len(0)), 0)
  expect_length(quantile(g1, c(0, 1)), 2)
  expect_length(quantile(g1, seq_len(10) / 10), 10)
})

test_that("cdf.Frechet and quantile.Frechet are consistent", {
  expect_equal(cdf(g1, quantile(g1, pvec)), pvec)
})

test_that("vectorization of a Frechet distribution work correctly", {
  d <- Frechet(c(-10, -5), c(10, 5), c(10, 5))
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
})

