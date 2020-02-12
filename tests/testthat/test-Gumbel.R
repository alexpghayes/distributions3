context("test-Gumbel")

test_that("print.Gumbel works", {
  expect_output(print(Gumbel()), regexp = "Gumbel distribution")
})

## Example distribution (from test-GeneralisedExtremeValue.R)

g2 <- Gumbel(0, 1)
xi2 <- 0

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x2 <- xvec
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.Gumbel works correctly", {
  expect_length(random(g2), 1)
  expect_length(random(g2, 100), 100)
  expect_length(random(g2, 0), 0)
  expect_error(random(g2, -2))
})

test_that("pdf.Gumbel works correctly", {
  p <- pvec[2:4]
  expect_equal(pdf(g2, x2), c(0, exp(-1), 0, NA))
  expect_equal(pdf(g2, quantile(g2, p)), (-log(p)) ^ (1 + xi2) * p)
  expect_length(pdf(g2, seq_len(0)), 0)
  expect_length(pdf(g2, seq_len(1)), 1)
  expect_length(pdf(g2, seq_len(10)), 10)
})

test_that("log_pdf.Gumbel works correctly", {
  expect_equal(log_pdf(g2, x2), c(-Inf, -1, -Inf, NA))
  expect_length(log_pdf(g2, seq_len(0)), 0)
  expect_length(log_pdf(g2, seq_len(1)), 1)
  expect_length(log_pdf(g2, seq_len(10)), 10)
})

test_that("cdf.Gumbel works correctly", {
  expect_equal(cdf(g2, x2), c(0, exp(-1), 1, NA))
  expect_length(cdf(g2, seq_len(0)), 0)
  expect_length(cdf(g2, seq_len(1)), 1)
  expect_length(cdf(g2, seq_len(10)), 10)
})

test_that("quantile.Gumbel works correctly", {
  q2 <- -log(-log(pvec[2:4]))
  expect_equal(quantile(g2, pvec), c(-Inf, q2, Inf, NA))
  expect_length(quantile(g2, seq_len(0)), 0)
  expect_length(quantile(g2, c(0, 1)), 2)
  expect_length(quantile(g2, seq_len(10) / 10), 10)
})

test_that("cdf.Gumbel and quantile.Gumbel are consistent", {
  expect_equal(cdf(g2, quantile(g2, pvec)), pvec)
})
