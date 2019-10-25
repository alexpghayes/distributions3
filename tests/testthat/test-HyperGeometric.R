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
