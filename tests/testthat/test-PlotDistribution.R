context("test-PlotDistribution")

test_that("Error thrown for incorrect input to plot.distribution", {
  expect_error(plot.distribution(1))
})

test_that("Binomial plot for scalar parameters", {
  B <- Binomial(20, 0.7)
  x <- plot(B)
  expect_equal(x$size, 20)
})

test_that("Binomial plot for vector parameters, all = FALSE", {
  B2 <- Binomial(20, c(0.1, 0.5, 0.9))
  x <- plot(B2, cdf = TRUE)
  expect_equal(x$size, c(20, 20, 20))
})

test_that("Binomial plot for vector parameters, all = TRUE", {
  B2 <- Binomial(20, c(0.1, 0.5, 0.9))
  x <- plot(B2, cdf = TRUE, all = TRUE)
  expect_equal(x$size, c(20, 20, 20))
})

test_that("Gamma plot for vector parameters, all = FALSE", {
  G <- Gamma(c(1, 3), 1:2)
  x <- plot(G)
  expect_equal(x$rate, G$rate)
  expect_equal(x$shape, G$shape)
})

test_that("Gamma plot for vector parameters, all = TRUE", {
  G <- Gamma(c(1, 3), 1:2)
  x <- plot(G, all = TRUE)
  expect_equal(x$rate, c(1, 1, 2, 2))
  expect_equal(x$shape, c(1, 3, 1, 3))
})
