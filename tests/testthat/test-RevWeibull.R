context("test-RevWeibull")

test_that("print.RevWeibull works", {
  expect_output(print(RevWeibull()), regexp = "RevWeibull distribution")
})

## Example distributions

# GEV parameter values

mu <- 0
sigma <- 1
xi3 <- -0.1

# Equivalent reversed Weibull parameter values

alpha <- -1 / xi3
s <- alpha * sigma
m <- mu + s

g3 <- RevWeibull(m, s, alpha)
up <- m

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.RevWeibull works correctly", {
  expect_length(random(g3), 1)
  expect_length(random(g3, 100), 100)
  expect_length(random(g3, 0), 0)
  expect_error(random(g3, -2))
})

test_that("pdf.RevWeibull works correctly", {
  p <- pvec[2:4]
  expect_equal(pdf(g3, x3), c(0, 0, exp(-1), 0, NA))
  expect_equal(pdf(g3, quantile(g3, p)), (-log(p)) ^ (1 + xi3) * p)
  expect_length(pdf(g3, seq_len(0)), 0)
  expect_length(pdf(g3, seq_len(1)), 1)
  expect_length(pdf(g3, seq_len(10)), 10)
})

test_that("log_pdf.RevWeibull works correctly", {
  expect_equal(log_pdf(g3, x3), c(-Inf, -Inf, -1, -Inf, NA))
  expect_length(log_pdf(g3, seq_len(0)), 0)
  expect_length(log_pdf(g3, seq_len(1)), 1)
  expect_length(log_pdf(g3, seq_len(10)), 10)
})

test_that("cdf.RevWeibull works correctly", {
  expect_equal(cdf(g3, x3), c(1, 0, exp(-1), 1, NA))
  expect_length(cdf(g3, seq_len(0)), 0)
  expect_length(cdf(g3, seq_len(1)), 1)
  expect_length(cdf(g3, seq_len(10)), 10)
})

test_that("quantile.RevWeibull works correctly", {
  q3 <- ((-log(pvec[2:4])) ^ (-xi3) - 1) / xi3
  expect_equal(quantile(g3, pvec), c(-Inf, q3, up, NA))
  expect_length(quantile(g3, seq_len(0)), 0)
  expect_length(quantile(g3, c(0, 1)), 2)
  expect_length(quantile(g3, seq_len(10) / 10), 10)
})

test_that("cdf.RevWeibull and quantile.RevWeibull are consistent", {
  expect_equal(cdf(g3, quantile(g3, pvec)), pvec)
})
