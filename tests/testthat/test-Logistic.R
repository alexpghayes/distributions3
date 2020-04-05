context("test-Logistic")

test_that("print.Logistic works", {
  expect_output(print(Logistic(1, 1)), regexp = "Logistic distribution")
})

test_that("likelihood.Logistic and log_likelihood.Logistic work correctly", {
  cau <- Logistic()
  x <- c(1, 1, 0)

  expect_equal(likelihood(cau, 1), dlogis(1))
  expect_equal(likelihood(cau, x), dlogis(1) * dlogis(1) * dlogis(0))

  expect_equal(log_likelihood(cau, 1), log(dlogis(1)))
  expect_equal(log_likelihood(cau, x), log(dlogis(1) * dlogis(1) * dlogis(0)))
})

test_that("random.Logistic work correctly", {
  cau <- Logistic()

  expect_length(random(cau), 1)
  expect_length(random(cau, 100), 100)
  expect_length(random(cau, 0), 0)
  expect_error(random(cau, -2))
})

test_that("pdf.Logistic work correctly", {
  cau <- Logistic()

  expect_equal(pdf(cau, 0), dlogis(0, 0, 1))
  expect_equal(pdf(cau, 1), dlogis(1, 0, 1))

  expect_length(pdf(cau, seq_len(0)), 0)
  expect_length(pdf(cau, seq_len(1)), 1)
  expect_length(pdf(cau, seq_len(10)), 10)
})

test_that("log_pdf.Logistic work correctly", {
  cau <- Logistic()

  expect_equal(log_pdf(cau, 0), log(dlogis(0, 0, 1)))
  expect_equal(log_pdf(cau, 1), log(dlogis(1, 0, 1)))

  expect_length(log_pdf(cau, seq_len(0)), 0)
  expect_length(log_pdf(cau, seq_len(1)), 1)
  expect_length(log_pdf(cau, seq_len(10)), 10)
})

test_that("cdf.Logistic work correctly", {
  cau <- Logistic()

  expect_equal(cdf(cau, 0), plogis(0, 0, 1))
  expect_equal(cdf(cau, 1), plogis(1, 0, 1))


  expect_length(cdf(cau, seq_len(0)), 0)
  expect_length(cdf(cau, seq_len(1)), 1)
  expect_length(cdf(cau, seq_len(10)), 10)
})

test_that("quantile.Logistic work correctly", {
  cau <- Logistic()

  expect_equal(quantile(cau, 0), qlogis(0, 0, 1))
  expect_equal(quantile(cau, 1), qlogis(1, 0, 1))


  expect_length(quantile(cau, seq_len(0)), 0)
  expect_length(quantile(cau, c(0, 1)), 2)
})
