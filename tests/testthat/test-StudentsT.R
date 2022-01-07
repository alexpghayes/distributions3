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

test_that("vectorization of a StudentsT distribution work correctly", {
  d <- StudentsT(c(1, 3))
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
