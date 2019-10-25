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
#'   <https://alexpghayes.github.io/distributions3>, where the math
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
  cat(glue("Weibull distribution (shape = {x$shape}, scale = {x$scale})\n"))
}

#' Draw a random sample from a Weibull distribution
#'
#' @inherit Weibull examples
#'
#' @param d A `Weibull` object created by a call to [Weibull()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Weibull distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Weibull <- function(d, n = 1L, ...) {
  rweibull(n = n, shape = d$shape, scale = d$scale)
}

#' Evaluate the probability mass function of a Weibull distribution
#'
#' Please see the documentation of [Weibull()] for some properties
#' of the Weibull distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Weibull examples
#' @inheritParams random.Weibull
#'
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
#' @inheritParams random.Weibull
#'
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
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family Weibull distribution
#'
quantile.Weibull <- function(d, p, ...) {
  qweibull(p = p, shape = d$shape, scale = d$scale)
}
