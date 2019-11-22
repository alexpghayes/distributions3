context("test-PlotDistribution")

test_that("Error thrown for incorrect input to plot.distribution", {
  expect_error(plot.distribution(1))
})

test_that("Error thrown for unsupported distribution", {
  x <- 1
  class(x) <- c("unsupported_name", "distribution")
  expect_error(plot(x))
})

test_that("Binomial plot for scalar parameters", {
  B <- Binomial(20, 0.7)
  expect_equal(plot(B)$size, 20)
})

test_that("Binomial pmf plot for vector parameters, all = FALSE", {
  B2 <- Binomial(20, c(0.1, 0.5, 0.9))
  expect_equal(plot(B2)$size, c(20, 20, 20))
})

test_that("Binomial cdf plot for vector parameters, all = TRUE", {
  B <- Binomial(20, c(0.1, 0.5, 0.9))
  expect_equal(plot(B, cdf = TRUE, all = TRUE)$size, c(20, 20, 20))
})

test_that("Gamma plot for vector parameters, all = FALSE", {
  G <- Gamma(c(1, 3), 1:2)
  x <- plot(G, xlim = c(-1, 7))
  expect_equal(x$rate, G$rate)
  expect_equal(x$shape, G$shape)
})

test_that("Gamma plot for vector parameters, all = TRUE", {
  G <- Gamma(c(1, 3), 1:2)
  x <- plot(G, cdf = TRUE, all = TRUE)
  expect_equal(x$rate, c(1, 1, 2, 2))
  expect_equal(x$shape, c(1, 3, 1, 3))
})
