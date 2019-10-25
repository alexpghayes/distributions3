context("test-StudentsT")

test_that("print.StudentsT works", {
  expect_output(print(StudentsT(1)), regexp = "Student's T distribution")
})

test_that("likelihood.StudentsT and log_likelihood.StudentsT work correctly", {
  s <- StudentsT(1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(s, 1), dt(1, 1))
  expect_equal(likelihood(s, x), dt(1, 1) * dt(1, 1) * dt(0, 1))

  expect_equal(log_likelihood(s, 1), log(dt(1, 1)))
  expect_equal(log_likelihood(s, x), log(dt(1, 1) * dt(1, 1) * dt(0, 1)))
})

test_that("random.StudentsT work correctly", {
  s <- StudentsT(1)

  expect_length(random(s), 1)
  expect_length(random(s, 100), 100)
  expect_length(random(s, 0), 0)
  expect_error(random(s, -2))
})

test_that("pdf.StudentsT work correctly", {
  s <- StudentsT(1)

  expect_equal(pdf(s, 0), dt(0, 1))
  expect_equal(pdf(s, 1), dt(1, 1))
  expect_equal(pdf(s, -12), dt(-12, 1))

  expect_length(pdf(s, seq_len(0)), 0)
  expect_length(pdf(s, seq_len(1)), 1)
  expect_length(pdf(s, seq_len(10)), 10)
})

test_that("log_pdf.StudentsT work correctly", {
  s <- StudentsT(1)

  expect_equal(log_pdf(s, 0), log(dt(0, 1)))
  expect_equal(log_pdf(s, 1), log(dt(1, 1)))
  expect_equal(log_pdf(s, -12), log(dt(-12, 1)))

  expect_length(log_pdf(s, seq_len(0)), 0)
  expect_length(log_pdf(s, seq_len(1)), 1)
  expect_length(log_pdf(s, seq_len(10)), 10)
})

test_that("cdf.StudentsT work correctly", {
  s <- StudentsT(1)

  expect_equal(cdf(s, 0), pt(0, 1))
  expect_equal(cdf(s, 1), pt(1, 1))


  expect_length(cdf(s, seq_len(0)), 0)
  expect_length(cdf(s, seq_len(1)), 1)
  expect_length(cdf(s, seq_len(10)), 10)
})

test_that("quantile.StudentsT work correctly", {
  s <- StudentsT(1)

  expect_equal(quantile(s, 0), qt(0, 1))
  expect_equal(quantile(s, 1), qt(1, 1))


  expect_length(quantile(s, seq_len(0)), 0)
  expect_length(quantile(s, c(0, 1)), 2)
})
