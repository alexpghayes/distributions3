context("test-GP")

test_that("print.GP works", {
  expect_output(print(GP()), regexp = "GP distribution")
})

## Example distributions

# Positive shape, finite lower end point
xi1 <- 0.1
g1 <- GP(0, 1, xi1)

# Zero shape
xi2 <- 0
g2 <- GP(0, 1, xi2)

# Negative shape, finite upper end point
xi3 <- -1e-7
g3 <- GP(0, 1, xi3)
up <- -1 / xi3

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(0, Inf, NA)
x1 <- xvec
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.GP works correctly", {
  expect_length(random(g1), 1)
  expect_length(random(g1, 100), 100)
  expect_length(random(g1, 0), 0)
  expect_error(random(g1, -2))

  expect_length(random(g2), 1)
  expect_length(random(g2, 100), 100)
  expect_length(random(g2, 0), 0)
  expect_error(random(g2, -2))

  expect_length(random(g3), 1)
  expect_length(random(g3, 100), 100)
  expect_length(random(g3, 0), 0)
  expect_error(random(g3, -2))
})

test_that("pdf.GP works correctly", {
  p <- pvec[2:4]
  expect_equal(pdf(g1, x1), c(1, 0, NA))
  expect_equal(pdf(g1, quantile(g1, p)), (1 - p) ^ (1 + xi1))
  expect_length(pdf(g1, seq_len(0)), 0)
  expect_length(pdf(g1, seq_len(1)), 1)
  expect_length(pdf(g1, seq_len(10)), 10)

  expect_equal(pdf(g2, x2), c(1, 0, NA))
  expect_equal(pdf(g2, quantile(g2, p)), (1 - p) ^ (1 + xi2))
  expect_length(pdf(g2, seq_len(0)), 0)
  expect_length(pdf(g2, seq_len(1)), 1)
  expect_length(pdf(g2, seq_len(10)), 10)

  expect_equal(pdf(g3, x3), c(0, 1, 0, NA))
  expect_equal(pdf(g3, quantile(g3, p)), (1 - p) ^ (1 + xi3))
  expect_length(pdf(g3, seq_len(0)), 0)
  expect_length(pdf(g3, seq_len(1)), 1)
  expect_length(pdf(g3, seq_len(10)), 10)
})

test_that("log_pdf.GP works correctly", {
  expect_equal(log_pdf(g1, x1), c(0, -Inf, NA))
  expect_length(log_pdf(g1, seq_len(0)), 0)
  expect_length(log_pdf(g1, seq_len(1)), 1)
  expect_length(log_pdf(g1, seq_len(10)), 10)

  expect_equal(log_pdf(g2, x2), c(0, -Inf, NA))
  expect_length(log_pdf(g2, seq_len(0)), 0)
  expect_length(log_pdf(g2, seq_len(1)), 1)
  expect_length(log_pdf(g2, seq_len(10)), 10)

  expect_equal(log_pdf(g3, x3), c(-Inf, 0, -Inf, NA))
  expect_length(log_pdf(g3, seq_len(0)), 0)
  expect_length(log_pdf(g3, seq_len(1)), 1)
  expect_length(log_pdf(g3, seq_len(10)), 10)
})

test_that("cdf.GP works correctly", {
  expect_equal(cdf(g1, x1), c(0, 1, NA))
  expect_length(cdf(g1, seq_len(0)), 0)
  expect_length(cdf(g1, seq_len(1)), 1)
  expect_length(cdf(g1, seq_len(10)), 10)

  expect_equal(cdf(g2, x2), c(0, 1, NA))
  expect_length(cdf(g2, seq_len(0)), 0)
  expect_length(cdf(g2, seq_len(1)), 1)
  expect_length(cdf(g2, seq_len(10)), 10)

  expect_equal(cdf(g3, x3), c(1, 0, 1, NA))
  expect_length(cdf(g3, seq_len(0)), 0)
  expect_length(cdf(g3, seq_len(1)), 1)
  expect_length(cdf(g3, seq_len(10)), 10)
})

test_that("quantile.GP works correctly", {
  q1 <- ((1 - pvec[2:4]) ^ (-xi1) - 1) / xi1
  expect_equal(quantile(g1, pvec), c(0, q1, Inf, NA))
  expect_length(quantile(g1, seq_len(0)), 0)
  expect_length(quantile(g1, c(0, 1)), 2)
  expect_length(quantile(g1, seq_len(10) / 10), 10)

  q2 <- -log(1 - pvec[2:4])
  expect_equal(quantile(g2, pvec), c(0, q2, Inf, NA))
  expect_length(quantile(g2, seq_len(0)), 0)
  expect_length(quantile(g2, c(0, 1)), 2)
  expect_length(quantile(g2, seq_len(10) / 10), 10)

  q3 <- ((1 - pvec[2:4]) ^ (-xi3) - 1) / xi3
  expect_equal(quantile(g3, pvec), c(0, q3, up, NA))
  expect_length(quantile(g3, seq_len(0)), 0)
  expect_length(quantile(g3, c(0, 1)), 2)
  expect_length(quantile(g3, seq_len(10) / 10), 10)
})

test_that("cdf.GP and quantile.GP are consistent", {
  expect_equal(cdf(g1, quantile(g1, pvec)), pvec)
  expect_equal(cdf(g2, quantile(g2, pvec)), pvec)
  expect_equal(cdf(g3, quantile(g3, pvec)), pvec)
})
