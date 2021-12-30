#' Create a Gamma distribution
#'
#' Several important distributions are special cases of the Gamma
#' distribution. When the shape parameter is `1`, the Gamma is an
#' exponential distribution with parameter \eqn{1/\beta}. When the
#' \eqn{shape = n/2} and \eqn{rate = 1/2}, the Gamma is a equivalent to
#' a chi squared distribution with n degrees of freedom. Moreover, if
#' we have \eqn{X_1} is \eqn{Gamma(\alpha_1, \beta)} and
#' \eqn{X_2} is \eqn{Gamma(\alpha_2, \beta)}, a function of these two variables
#' of the form \eqn{\frac{X_1}{X_1 + X_2}} \eqn{Beta(\alpha_1, \alpha_2)}.
#' This last property frequently appears in another distributions, and it
#' has extensively been used in multivariate methods. More about the Gamma
#' distribution will be added soon.
#'
#' @param shape The shape parameter. Can be any positive number.
#' @param rate The rate parameter. Can be any positive number. Defaults
#'   to `1`.
#'
#' @return A `Gamma` object.
#' @export

#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Gamma random variable
#'   with parameters
#'   `shape` = \eqn{\alpha} and
#'   `rate` = \eqn{\beta}.
#'
#'   **Support**: \eqn{x \in (0, \infty)}
#'
#'   **Mean**: \eqn{\frac{\alpha}{\beta}}
#'
#'   **Variance**: \eqn{\frac{\alpha}{\beta^2}}
#'
#'   **Probability density function (p.m.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\beta^{\alpha}}{\Gamma(\alpha)} x^{\alpha - 1} e^{-\beta x}
#'   }{
#'     f(x) = \frac{\beta^{\alpha}}{\Gamma(\alpha)} x^{\alpha - 1} e^{-\beta x}
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\Gamma(\alpha, \beta x)}{\Gamma{\alpha}}
#'   }{
#'     f(x) = \frac{\Gamma(\alpha, \beta x)}{\Gamma{\alpha}}
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \Big(\frac{\beta}{ \beta - t}\Big)^{\alpha}, \thinspace t < \beta
#'   }{
#'     E(e^(tX)) = \Big(\frac{\beta}{ \beta - t}\Big)^{\alpha}, \thinspace t < \beta
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Gamma(5, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Gamma <- function(shape, rate = 1) {
  d <- list(shape = shape, rate = rate)
  class(d) <- c("Gamma", "distribution")
  d
}

#' @export
print.Gamma <- function(x, ...) {
  cat(glue("Gamma distribution (shape = {x$shape}, rate = {x$rate})"), "\n")
}

#' @export
mean.Gamma <- function(x, ...) {
  ellipsis::check_dots_used()
  x$shape / x$rate
}

#' @export
variance.Gamma <- function(x, ...) x$shape / x$rate^2

#' @export
skewness.Gamma <- function(x, ...) 2 / sqrt(x$shape)

#' @export
kurtosis.Gamma <- function(x, ...) 6 / x$shape

#' Draw a random sample from a Gamma distribution
#'
#' @inherit Gamma examples
#'
#' @param x A `Gamma` object created by a call to [Gamma()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Gamma <- function(x, n = 1L, ...) {
  rgamma(n = n, shape = x$shape, rate = x$rate)
}

#' Evaluate the probability mass function of a Gamma distribution
#'
#' @inherit Gamma examples
#'
#' @param d A `Gamma` object created by a call to [Gamma()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Gamma <- function(d, x, ...) {
  dgamma(x = x, shape = d$shape, rate = d$rate)
}

#' @rdname pdf.Gamma
#' @export
#'
log_pdf.Gamma <- function(d, x, ...) {
  dgamma(x = x, shape = d$shape, rate = d$rate, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Gamma distribution
#'
#' @inherit Gamma examples
#'
#' @param d A `Gamma` object created by a call to [Gamma()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Gamma <- function(d, x, ...) {
  pgamma(q = x, shape = d$shape, rate = d$rate)
}

#' Determine quantiles of a Gamma distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Gamma examples
#' @inheritParams random.Gamma
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.Gamma <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qgamma(p = probs, shape = x$shape, rate = x$rate)
}

#' Fit a Gamma distribution to data
#'
#' @param d A `Gamma` object created by a call to [Gamma()].
#' @param x A vector to fit the Gamma distribution to.
#' @param ... Unused.
#'
#' @return a `Gamma` object
#' @export
fit_mle.Gamma <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  stop("`fit_mle` is not implemented for the Gamma distribution yet")
}

#' Compute the sufficient statistics for a Gamma distribution from data
#'
#'   - `sum`: The sum of the data.
#'   - `log_sum`: The log of the sum of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @inherit fit_mle.Gamma
#' @export
suff_stat.Gamma <- function(d, x, ...) {
  if (any(x < 0)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), log_sum = sum(log(x)), samples = length(x))
}

#' Return the support of the Gamma distribution
#'
#' @param d An `Gamma` object created by a call to [Gamma()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Gamma <- function(d) c(0, Inf)
