context("test-Tukey")

test_that("print.Tukey works", {
  expect_output(print(Tukey(1, 2, 2)), regexp = "Tukey distribution")
})
