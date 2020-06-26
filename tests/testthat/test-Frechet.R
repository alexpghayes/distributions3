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
