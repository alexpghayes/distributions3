context("test-utils")

test_that("is_distribution() works", {
  expect_true(is_distribution(Normal()))
  expect_false(is_distribution(123))
})


test_that("{methods}.dstribution work", {
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

test_that("apply_dpqr() applied to 'random' works", {
  N <- Normal(c(0, 10), c(1, 1))

  expect_length(random(N[1], 1), 1)
  expect_length(random(N[1], 2), 2)
  expect_true(is.data.frame(random(N[1], 1, drop = FALSE)))
  expect_true(is.data.frame(random(N[1], 2, drop = FALSE)))
  expect_equal(dim(random(N[1], 1, drop = FALSE)), c(1L, 1L))
  expect_equal(dim(random(N[1], 2, drop = FALSE)), c(1L, 2L))
  expect_error(random(N[1], c(1, 2)))
  expect_error(random(N[1], cbind(1, 2)))
  expect_error(random(N[1], rbind(1, 2)))

  expect_length(random(N, 1), length(N))
  expect_equal(dim(random(N, 2)), c(2L, 2L))
  expect_true(is.matrix(random(N, 2)))
  expect_false(is.data.frame(random(N, 2)))
  expect_true(is.data.frame(random(N, 1, drop = FALSE)))
  expect_true(is.data.frame(random(N, 2, drop = FALSE)))
  expect_equal(dim(random(N, 1, drop = FALSE)), c(2L, 1L))
  expect_equal(dim(random(N, 2, drop = FALSE)), c(2L, 2L))
  expect_error(random(N, c(1, 2)))
  expect_error(random(N, cbind(1, 2)))
  expect_error(random(N, rbind(1, 2)))
})

test_that("apply_dpqr() applied to 'pdf' works", {
  N <- Normal(c(0, 100), c(1, 1))

  expect_length(pdf(N[1], 0.5), 1)
  expect_length(pdf(N[1], c(0.2, 0.5)), 2)
  expect_true(is.data.frame(pdf(N[1], 0.5, drop = FALSE)))
  expect_true(is.data.frame(pdf(N[1], c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(pdf(N[1], 0.5, drop = FALSE)), c(1L, 1L))
  expect_equal(dim(pdf(N[1], c(0.2, 0.5), drop = FALSE)), c(1L, 2L))
  expect_equal(pdf(N[1], cbind(0.2, 0.5)), pdf(N[1], c(0.2, 0.5)))
  expect_equal(pdf(N[1], cbind(0.2, 0.5), drop = FALSE), pdf(N[1], c(0.2, 0.5), drop = FALSE))
  expect_equal(pdf(N[1], rbind(0.2, 0.5)), matrix(c(pdf(N[1], 0.2), pdf(N[1], 0.5))))
  expect_equal(
    pdf(N[1], rbind(0.2, 0.5), drop = FALSE),
    data.frame(density = c(pdf(N[1], 0.2), pdf(N[1], 0.5)))
  )

  expect_length(pdf(N, 0.5), length(N))
  expect_equal(dim(pdf(N, c(0.2, 0.5))), c(2L, 2L))
  expect_true(is.matrix(pdf(N, c(0.2, 0.5))))
  expect_false(is.data.frame(pdf(N, c(0.2, 0.5))))
  expect_true(is.data.frame(pdf(N, 0.2, drop = FALSE)))
  expect_true(is.data.frame(pdf(N, c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(pdf(N, 0.5, drop = FALSE)), c(2L, 1L))
  expect_equal(dim(pdf(N, c(0.3, 0.5), drop = FALSE)), c(2L, 2L))
  expect_equal(pdf(N, cbind(0.2, 0.5)), pdf(N, c(0.2, 0.5)))
  expect_equal(pdf(N, cbind(0.2, 0.5), drop = FALSE), pdf(N, c(0.2, 0.5), drop = FALSE))
  expect_equal(pdf(N, rbind(0.2, 0.5)), matrix(c(pdf(N[1], 0.2), pdf(N[2], 0.5))))
  expect_equal(pdf(N, rbind(0.2, 0.5)), matrix(c(pdf(N[1], 0.2), pdf(N[2], 0.5))))
  expect_equal(
    pdf(N, rbind(0.2, 0.5), drop = FALSE),
    data.frame(density = c(pdf(N[1], 0.2), pdf(N[2], 0.5)))
  )
})

test_that("apply_dpqr() applied to 'log_pdf' works", {
  N <- Normal(c(0, 100), c(1, 1))

  expect_length(log_pdf(N[1], 0.5), 1)
  expect_length(log_pdf(N[1], c(0.2, 0.5)), 2)
  expect_true(is.data.frame(log_pdf(N[1], 0.5, drop = FALSE)))
  expect_true(is.data.frame(log_pdf(N[1], c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(log_pdf(N[1], 0.5, drop = FALSE)), c(1L, 1L))
  expect_equal(dim(log_pdf(N[1], c(0.2, 0.5), drop = FALSE)), c(1L, 2L))
  expect_equal(log_pdf(N[1], cbind(0.2, 0.5)), log_pdf(N[1], c(0.2, 0.5)))
  expect_equal(log_pdf(N[1], cbind(0.2, 0.5), drop = FALSE), log_pdf(N[1], c(0.2, 0.5), drop = FALSE))
  expect_equal(log_pdf(N[1], rbind(0.2, 0.5)), matrix(c(log_pdf(N[1], 0.2), log_pdf(N[1], 0.5))))
  expect_equal(
    log_pdf(N[1], rbind(0.2, 0.5), drop = FALSE),
    data.frame("logLik" = c(log_pdf(N[1], 0.2), log_pdf(N[1], 0.5)))
  )

  expect_length(log_pdf(N, 0.5), length(N))
  expect_equal(dim(log_pdf(N, c(0.2, 0.5))), c(2L, 2L))
  expect_true(is.matrix(log_pdf(N, c(0.2, 0.5))))
  expect_false(is.data.frame(log_pdf(N, c(0.2, 0.5))))
  expect_true(is.data.frame(log_pdf(N, 0.2, drop = FALSE)))
  expect_true(is.data.frame(log_pdf(N, c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(log_pdf(N, 0.5, drop = FALSE)), c(2L, 1L))
  expect_equal(dim(log_pdf(N, c(0.3, 0.5), drop = FALSE)), c(2L, 2L))
  expect_equal(log_pdf(N, cbind(0.2, 0.5)), log_pdf(N, c(0.2, 0.5)))
  expect_equal(log_pdf(N, cbind(0.2, 0.5), drop = FALSE), log_pdf(N, c(0.2, 0.5), drop = FALSE))
  expect_equal(log_pdf(N, rbind(0.2, 0.5)), matrix(c(log_pdf(N[1], 0.2), log_pdf(N[2], 0.5))))
  expect_equal(log_pdf(N, rbind(0.2, 0.5)), matrix(c(log_pdf(N[1], 0.2), log_pdf(N[2], 0.5))))
  expect_equal(
    log_pdf(N, rbind(0.2, 0.5), drop = FALSE),
    data.frame("logLik" = c(log_pdf(N[1], 0.2), log_pdf(N[2], 0.5)))
  )
})

test_that("apply_dpqr() applied to 'cdf' works", {
  N <- Normal(c(0, 100), c(1, 1))

  expect_length(cdf(N[1], 0.5), 1)
  expect_length(cdf(N[1], c(0.2, 0.5)), 2)
  expect_true(is.data.frame(cdf(N[1], 0.5, drop = FALSE)))
  expect_true(is.data.frame(cdf(N[1], c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(cdf(N[1], 0.5, drop = FALSE)), c(1L, 1L))
  expect_equal(dim(cdf(N[1], c(0.2, 0.5), drop = FALSE)), c(1L, 2L))
  expect_equal(cdf(N[1], cbind(0.2, 0.5)), cdf(N[1], c(0.2, 0.5)))
  expect_equal(cdf(N[1], cbind(0.2, 0.5), drop = FALSE), cdf(N[1], c(0.2, 0.5), drop = FALSE))
  expect_equal(cdf(N[1], rbind(0.2, 0.5)), matrix(c(cdf(N[1], 0.2), cdf(N[1], 0.5))))
  expect_equal(
    cdf(N[1], rbind(0.2, 0.5), drop = FALSE),
    data.frame(probability = c(cdf(N[1], 0.2), cdf(N[1], 0.5)))
  )

  expect_length(cdf(N, 0.5), length(N))
  expect_equal(dim(cdf(N, c(0.2, 0.5))), c(2L, 2L))
  expect_true(is.matrix(cdf(N, c(0.2, 0.5))))
  expect_false(is.data.frame(cdf(N, c(0.2, 0.5))))
  expect_true(is.data.frame(cdf(N, 0.2, drop = FALSE)))
  expect_true(is.data.frame(cdf(N, c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(cdf(N, 0.5, drop = FALSE)), c(2L, 1L))
  expect_equal(dim(cdf(N, c(0.3, 0.5), drop = FALSE)), c(2L, 2L))
  expect_equal(cdf(N, cbind(0.2, 0.5)), cdf(N, c(0.2, 0.5)))
  expect_equal(cdf(N, cbind(0.2, 0.5), drop = FALSE), cdf(N, c(0.2, 0.5), drop = FALSE))
  expect_equal(cdf(N, rbind(0.2, 0.5)), matrix(c(cdf(N[1], 0.2), cdf(N[2], 0.5))))
  expect_equal(cdf(N, rbind(0.2, 0.5)), matrix(c(cdf(N[1], 0.2), cdf(N[2], 0.5))))
  expect_equal(
    cdf(N, rbind(0.2, 0.5), drop = FALSE),
    data.frame(probability = c(cdf(N[1], 0.2), cdf(N[2], 0.5)))
  )
})

test_that("apply_dpqr() applied to 'quantile' works", {
  N <- Normal(c(0, 100), c(1, 1))

  expect_length(quantile(N[1], 0.5), 1)
  expect_length(quantile(N[1], c(0.2, 0.5)), 2)
  expect_true(is.data.frame(quantile(N[1], 0.5, drop = FALSE)))
  expect_true(is.data.frame(quantile(N[1], c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(quantile(N[1], 0.5, drop = FALSE)), c(1L, 1L))
  expect_equal(dim(quantile(N[1], c(0.2, 0.5), drop = FALSE)), c(1L, 2L))
  expect_equal(quantile(N[1], cbind(0.2, 0.5)), quantile(N[1], c(0.2, 0.5)))
  expect_equal(quantile(N[1], cbind(0.2, 0.5), drop = FALSE), quantile(N[1], c(0.2, 0.5), drop = FALSE))
  expect_equal(quantile(N[1], rbind(0.2, 0.5)), matrix(c(quantile(N[1], 0.2), quantile(N[1], 0.5))))
  expect_equal(
    quantile(N[1], rbind(0.2, 0.5), drop = FALSE),
    data.frame(quantile = c(quantile(N[1], 0.2), quantile(N[1], 0.5)))
  )

  expect_length(quantile(N, 0.5), length(N))
  expect_equal(dim(quantile(N, c(0.2, 0.5))), c(2L, 2L))
  expect_true(is.matrix(quantile(N, c(0.2, 0.5))))
  expect_false(is.data.frame(quantile(N, c(0.2, 0.5))))
  expect_true(is.data.frame(quantile(N, 0.2, drop = FALSE)))
  expect_true(is.data.frame(quantile(N, c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(quantile(N, 0.5, drop = FALSE)), c(2L, 1L))
  expect_equal(dim(quantile(N, c(0.3, 0.5), drop = FALSE)), c(2L, 2L))
  expect_equal(quantile(N, cbind(0.2, 0.5)), quantile(N, c(0.2, 0.5)))
  expect_equal(quantile(N, cbind(0.2, 0.5), drop = FALSE), quantile(N, c(0.2, 0.5), drop = FALSE))
  expect_equal(quantile(N, rbind(0.2, 0.5)), matrix(c(quantile(N[1], 0.2), quantile(N[2], 0.5))))
  expect_equal(quantile(N, rbind(0.2, 0.5)), matrix(c(quantile(N[1], 0.2), quantile(N[2], 0.5))))
  expect_equal(
    quantile(N, rbind(0.2, 0.5), drop = FALSE),
    data.frame(quantile = c(quantile(N[1], 0.2), quantile(N[2], 0.5)))
  )
})
