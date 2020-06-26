#' Create a Normal distribution
#'
#' The Normal distribution is ubiquituous in statistics, partially because
#' of the central limit theorem, which states that sums of i.i.d. random
#' variables eventually become Normal. Linear transformations of Normal
#' random variables result in new random variables that are also Normal. If
#' you are taking an intro stats course, you'll likely use the Normal
#' distribution for Z-tests and in simple linear regression. Under
#' regularity conditions, maximum likelihood estimators are
#' asymptotically Normal. The Normal distribution is also called the
#' gaussian distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks,
#'   which is also the mean of the distribution. Can be any real number.
#'   Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks,
#'   which is also the **standard deviation** of the distribution. Can be any
#'   positive number. Defaults to `1`. If you would like a Normal
#'   distribution with **variance** \eqn{\sigma^2}, be sure to take the
#'   square root, as this is a common source of errors.
#'
#' @return A `Normal` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Normal random variable with mean
#'   `mu` = \eqn{\mu} and standard deviation `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
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
#'   approximated numerically. The c.d.f. of a standard Normal is sometimes
#'   called the "error function". The notation \eqn{\Phi(t)} also stands
#'   for the c.d.f. of a standard Normal evaluated at \eqn{t}. Z-tables
#'   list the value of \eqn{\Phi(t)} for various \eqn{t}.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2}
#'   }{
#'     E(e^(tX)) = e^(\mu t + \sigma^2 t^2 / 2)
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Normal(5, 2)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
#' ### example: calculating p-values for two-sided Z-test
#'
#' # here the null hypothesis is H_0: mu = 3
#' # and we assume sigma = 2
#'
#' # exactly the same as: Z <- Normal(0, 1)
#' Z <- Normal()
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
#' ns <- random(Normal(3, 7), 26)
#'
#' # test if sample is Normal(3, 7)
#' ks.test(ns, pnorm, mean = 3, sd = 7)
#'
#' # test if sample is gamma(8, 3) using base R pgamma()
#' ks.test(ns, pgamma, shape = 8, rate = 3)
#'
#' ### MISC
#'
#' # note that the cdf() and quantile() functions are inverses
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
#'
Normal <- function(mu = 0, sigma = 1) {
  d <- list(mu = mu, sigma = sigma)
  class(d) <- c("Normal", "distribution")
  d
}

#' @export
print.Normal <- function(x, ...) {
  cat(glue("Normal distribution (mu = {x$mu}, sigma = {x$sigma})"), "\n")
}

#' @export
mean.Normal <- function(d, ...) d$mu

#' @export
variance.Normal <- function(d, ...) d$sigma ^ 2

#' @export
skewness.Normal <- function(d, ...) 0

#' @export
kurtosis.Normal <- function(d, ...) 0

#' Draw a random sample from a Normal distribution
#'
#' Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Normal examples
#'
#' @param d A `Normal` object created by a call to [Normal()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
#'
random.Normal <- function(d, n = 1L, ...) {
  rnorm(n = n, mean = d$mu, sd = d$sigma)
}

#' Evaluate the probability mass function of a Normal distribution
#'
#' Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Normal examples
#' @inheritParams random.Normal
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Normal <- function(d, x, ...) {
  dnorm(x = x, mean = d$mu, sd = d$sigma)
}

#' @rdname pdf.Normal
#' @export
#'
log_pdf.Normal <- function(d, x, ...) {
  dnorm(x = x, mean = d$mu, sd = d$sigma, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Normal distribution
#'
#' @inherit Normal examples
#' @inheritParams random.Normal
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Normal <- function(d, x, ...) {
  pnorm(q = x, mean = d$mu, sd = d$sigma)
}

#' Determine quantiles of a Normal distribution
#'
#' Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Normal examples
#' @inheritParams random.Normal
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family Normal distribution
#'
quantile.Normal <- function(d, p, ...) {
  qnorm(p = p, mean = d$mu, sd = d$sigma)
}

#' Fit a Normal distribution to data
#'
#' @param d A `Normal` object created by a call to [Normal()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Normal distribution
#'
#' @return A `Normal` object.
#' @export
fit_mle.Normal <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Normal(ss$mu, ss$sigma)
}


#' Compute the sufficient statistics for a Normal distribution from data
#'
#' @inheritParams fit_mle.Normal
#'
#' @return A named list of the sufficient statistics of the normal
#'   distribution:
#'
#'   - `mu`: The sample mean of the data.
#'   - `sigma`: The sample standard deviation of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
suff_stat.Normal <- function(d, x, ...) {
  valid_x <- is.numeric(x)
  if (!valid_x) stop("`x` must be a numeric vector")
  list(mu = mean(x), sigma = sd(x), samples = length(x))
}

#' Return the support of the Normal distribution
#'
#' @param d An `Normal` object created by a call to [Normal()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Normal <- function(d){
  if(!is_distribution(d)){
    message("d has to be a disitrubtion")
    stop()
  }
  return(c(-Inf, Inf))
}
