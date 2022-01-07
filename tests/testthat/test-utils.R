context("test-utils")

test_that("is_distribution() works", {
  expect_true(is_distribution(Normal()))
  expect_false(is_distribution(123))
})


test_that("{methods}.dstribution work correctly", {
  n <- Normal(c(0, 10), c(1, 1))

  expect_null(dim.distribution(n))
  expect_equal(length.distribution(n), 2)
  expect_equal(n[1], Normal(0, 1))
  expect_length(format.distribution(n), 2L)
  expect_length(print.distribution(n), 2L)
  expect_null(names.distribution(n))

  expect_silent(n <- `names<-.distribution`(n, c("a", "b")))
  expect_named(n, c("a", "b"))
  expect_silent(names(n) <- NULL)
  expect_equal(as.matrix.distribution(n), as.matrix(data.frame(mu = c(0, 10), sigma = c(1, 1)))) 
  df <- data.frame(n = 1:2)
  df$n <- n
  expect_equal(as.data.frame.distribution(n), df)
  expect_equal(as.list.distribution(n), list(mu = c(0, 10), sigma = c(1, 1)))
  expect_equal(n, c.distribution(n[1], n[2]))
  expect_equal(
    capture_output(print(summary(n))),
    capture_output({
    cat("Normal distribution: \n")
    print(summary.data.frame(as.data.frame(as.matrix(n))))
    })
  )
})

