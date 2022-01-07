context("test-Tukey")

test_that("print.Tukey works", {
  expect_output(print(Tukey(1, 2, 2)), regexp = "Tukey distribution")
})

test_that("vectorization of a Normal distribution work correctly", {
  d <- Tukey(c(3L, 4L), 16L, 2L)
  d1 <- d[1]
  d2 <- d[2]

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
