context("test-utils")

test_that("is_distribution works", {

  expect_true(is_distribution(normal()))

  expect_false(is_distribution(123))
})
