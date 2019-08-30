context("test-Geometric")

test_that("print.Geometric works", {
  expect_output(print(Geometric()), regexp = "Geometric distribution")
})

test_that("likelihood.Geometric and log_likelihood.Geometric work correctly", {
  cau <- Geometric()
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dgeom(1, 0.5))
  expect_equal(likelihood(cau, x), dgeom(1, 0.5) * dgeom(1, 0.5) * dgeom(0, 0.5))

  expect_equal(log_likelihood(cau, 1), log(dgeom(1, 0.5)))
  expect_equal(log_likelihood(cau, x), log(dgeom(1, 0.5) * dgeom(1, 0.5) * dgeom(0, 0.5)))
})

test_that("random.Geometric work correctly", {
  cau <- Geometric()

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.Geometric work correctly", {
  cau <- Geometric()

  expect_equal(pdf(cau, 0), dgeom(0, 0.5))
  expect_equal(pdf(cau, 1), dgeom(1, 0.5))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.Geometric work correctly", {
  cau <- Geometric()

  expect_equal(log_pdf(cau, 0), log(dgeom(0, 0.5)))
  expect_equal(log_pdf(cau, 1), log(dgeom(1, 0.5)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.Geometric work correctly", {
  cau <- Geometric()

  expect_equal(cdf(cau, 0), pgeom(0, 0.5))
  expect_equal(cdf(cau, 1), pgeom(1, 0.5))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.Geometric work correctly", {
  cau <- Geometric()

  expect_equal(quantile(cau, 0), qgeom(0, 0.5))
  expect_equal(quantile(cau, 1), qgeom(1, 0.5))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})
