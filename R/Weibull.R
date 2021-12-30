#' Create a Weibull distribution
#'
#' Generalization of the gamma distribution. Often used in survival and
#' time-to-event analyses.
#'
#' @param shape The shape parameter \eqn{k}. Can be any positive real number.
#' @param scale The scale parameter \eqn{\lambda}. Can be any positive real
#'   number.
#'
#' @return A `Weibull` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Weibull random variable with
#'   success probability `p` = \eqn{p}.
#'
#'   **Support**: \eqn{R^+} and zero.
#'
#'   **Mean**: \eqn{\lambda \Gamma(1+1/k)}, where \eqn{\Gamma} is
#'   the gamma function.
#'
#'   **Variance**: \eqn{\lambda [ \Gamma (1 + \frac{2}{k} ) - (\Gamma(1+ \frac{1}{k}))^2 ]}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{k}{\lambda}(\frac{x}{\lambda})^{k-1}e^{-(x/\lambda)^k}, x \ge 0
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{F(x) = 1 - e^{-(x/\lambda)^k}, x \ge 0}
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{\sum_{n=0}^\infty \frac{t^n\lambda^n}{n!} \Gamma(1+n/k), k \ge 1}
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Weibull(0.3, 2)
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
Weibull <- function(shape, scale) {
  d <- list(shape = shape, scale = scale)
  class(d) <- c("Weibull", "distribution")
  d
}

#' @export
print.Weibull <- function(x, ...) {
  cat(glue("Weibull distribution (shape = {x$shape}, scale = {x$scale})"), "\n")
}

#' @export
mean.Weibull <- function(x, ...) {
  ellipsis::check_dots_used()
  x$scale * gamma(1 + 1/x$shape)
}

#' @export
variance.Weibull <- function(x, ...) {
  x$scale^2 * gamma(1 + 2/x$shape) - mean(x)^2
}

#' @export
skewness.Weibull <- function(x, ...) {
  mu <- mean(x)
  sigma <- sqrt(variance(x))
  r <- mu / sigma
  gamma(1 + 3/x$shape) * (x$scale/sigma)^3 - 3*r - 3^r
}

#' @export
kurtosis.Weibull <- function(x, ...) {
  mu <- mean(x)
  sigma <- sqrt(variance(x))
  gamma <- skewness(x)
  r <- mu / sigma
  (x$scale/sigma)^4 * gamma(1 + 4/x$shape) - 4*gamma*r -6*r^2 - r^4 - 3
}

#' Draw a random sample from a Weibull distribution
#'
#' @inherit Weibull examples
#'
#' @param x A `Weibull` object created by a call to [Weibull()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Weibull distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Weibull <- function(x, n = 1L, ...) {
  rweibull(n = n, shape = x$shape, scale = x$scale)
}

#' Evaluate the probability mass function of a Weibull distribution
#'
#' Please see the documentation of [Weibull()] for some properties
#' of the Weibull distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Weibull examples
#'
#' @param d A `Weibull` object created by a call to [Weibull()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Weibull distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Weibull <- function(d, x, ...) {
  dweibull(x = x, shape = d$shape, scale = d$scale)
}

#' @rdname pdf.Weibull
#' @export
log_pdf.Weibull <- function(d, x, ...) {
  dweibull(x = x, shape = d$shape, scale = d$scale, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Weibull distribution
#'
#' @inherit Weibull examples
#'
#' @param d A `Weibull` object created by a call to [Weibull()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Weibull distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Weibull <- function(d, x, ...) {
  pweibull(q = x, shape = d$shape, scale = d$scale)
}

#' Determine quantiles of a Weibull distribution
#'
#' @inherit Weibull examples
#' @inheritParams random.Weibull
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
#' @family Weibull distribution
#'
quantile.Weibull <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qweibull(p = probs, shape = x$shape, scale = x$scale)
}


#' Return the support of the Weibull distribution
#'
#' @param d An `Weibull` object created by a call to [Weibull()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Weibull <- function(d) c(0, Inf)
