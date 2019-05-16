#' Create a normal distribution
#'
#' The normal distribution is ubiquituous in statistics, partially because
#' of the central limit theorem, which states that sums of i.i.d. random
#' variables eventually become normal. Linear transformations of normal
#' random variables result in new random variables that are also normal. If
#' you are taking an intro stats course, you'll likely use the normal
#' distribution for Z-tests and in simple linear regression. Under
#' regularity conditions, maximum likelihood estimators are
#' asymptotically normal. The normal distribution is also called the
#' gaussian distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks,
#'   which is also the mean of the distribution. Can be any real number.
#'   Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks,
#'   which is also the **standard deviation** of the distribution. Can be any
#'   positive number. Defaults to `1`. If you would like a normal
#'   distribution with **variance** \eqn{\sigma^2}, be sure to take the
#'   square root, as this is a common source of errors.
#'
#' @return A `normal` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a normal random variable with mean
#'   `mu` = \eqn{\mu} and standard deviation `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{\mathbb{R}}{R}, the set of all real numbers
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{\sigma^2}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2}
#'   }{
#'     f(x) = 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2))
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   The cumulative distribution function has the form
#'
#'   \deqn{
#'     F(t) = \int_{-\infty}^t \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2} dx
#'   }{
#'     F(t) = integral_{-\infty}^t 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2)) dx
#'   }
#'
#'   but this integral does not have a closed form solution and must be
#'   approximated numerically. The c.d.f. of a standard normal is sometimes
#'   called the "error function". The notation \eqn{\Phi(t)} also stands
#'   for the c.d.f. of a standard normal evaluated at \eqn{t}. Z-tables
#'   list the value of \eqn{\Phi(t)} for various \eqn{t}.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     \mathbb{E}(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2}
#'   }{
#'     E(e^(tX)) = e^(\mu t + \sigma^2 t^2 / 2)
#'   }
#'
#' @examples
#'
#' n <- normal(5, 2)
#' n
#'
#' random(n, 10)
#' pdf(n, 2)
#' cdf(n, 4)
#' quantile(n, 0.7)
#'
#' ### example: calculating p-values for two-sided Z-test
#'
#' # here the null hypothesis is H_0: mu = 3
#' # and we assume sigma = 2
#'
#' # exactly the same as: Z <- normal(0, 1)
#' Z <- normal()
#'
#' # data to test
#' x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' nx <- length(x)
#'
#' # calculate the z-statistic
#' z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
#' z_stat
#'
#' # calculate the two-sided p-value
#' 1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#'
#' # exactly equivalent to the above
#' 2 * cdf(Z, -abs(z_stat))
#'
#' # p-value for one-sided test
#' # H_0: mu <= 3   vs   H_A: mu > 3
#' 1 - cdf(Z, z_stat)
#'
#' # p-value for one-sided test
#' # H_0: mu >= 3   vs   H_A: mu < 3
#' cdf(Z, z_stat)
#'
#' ### example: calculating a 88 percent Z CI for a mean
#'
#' # same `x` as before, still assume `sigma = 2`
#'
#' # lower-bound
#' mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # upper-bound
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # equivalent to
#' mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # also equivalent to
#' mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ### generating random samples and plugging in ks.test()
#'
#' set.seed(27)
#'
#' # generate a random sample
#' normal_sample <- random(normal(3, 7), 26)
#'
#' # test if sample is normal(3, 7)
#' ks.test(normal_sample, pnorm, mean = 3, sd = 7)
#'
#' # test if sample is gamma(8, 3) using base R pgamma()
#' ks.test(normal_sample, pgamma, shape = 8, rate = 3)
#'
#'
#' ###
#'
#' # note that the cdf() and quantile() functions are inverses
#' cdf(n, quantile(n, 0.7))
#' quantile(n, cdf(n, 7))
#'
normal <- function(mu = 0, sigma = 1) {
  d <- list(mu = mu, sigma = sigma)
  class(d) <- "normal"
  d
}

#' @export
print.normal <- function(x, ...) {
  cat(glue("normal distribution (mu = {x$mu}, sigma = {x$sigma})"))
}

#' Draw a random sample from a normal distribution
#'
#' Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit normal examples
#'
#' @param d A `normal` object created by a call to [normal()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family normal distribution
#'
#' @return A numeric vector of length `n`.
#' @export
#'
#'
random.normal <- function(d, n = 1L, ...) {
  rnorm(n = n, mean = d$mu, sd = d$sigma)
}

#' Evaluate the probability mass function of a normal distribution
#'
#' Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.normal <- function(d, x, ...) {
  dnorm(x = x, mean = d$mu, sd = d$sigma)
}

#' Evaluate the cumulative distribution function of a normal distribution
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.normal <- function(d, x, ...) {
  pnorm(q = x, mean = d$mu, sd = d$sigma)
}

#' Determine quantiles of a normal distribution
#'
#' Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family normal distribution
#'
quantile.normal <- function(d, p, ...) {
  qnorm(p = p, mean = d$mu, sd = d$sigma)
}

#' Fit a normal distribution to data
#'
#' @param d A `normal` object created by a call to [normal()].
#' @param x A vector to fit the distribution to.
#'
#' @family normal distribution
#'
#' @return A `normal` object.
#' @export
fit_mle.normal <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  normal(ss$mu, ss$sigma)
}


#' Compute the sufficient statistics for a normal distribution from data
#'
#' @inheritParams normal
#' @export
suff_stat.normal <- function(d, x, ...) {
  valid_x <- is.numeric(x)
  if(!valid_x) stop("`x` must be a numeric vector")
  list(mu = mean(x), sigma = sd(x), samples = length(x))
}
