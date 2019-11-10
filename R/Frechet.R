#' Create a Frechet distribution
#'
#' The Frechet distribution is a special case of the `\link{GEV}` distribution,
#' obtained when the GEV shape parameter \eqn{\xi} is positive.
#' It may be referred to as a type II extreme value distribution.
#'
#' @param location The location (minimum) parameter \eqn{m}.
#'   `location` can be any real number.  Defaults to `0`.
#' @param scale The scale parameter \eqn{s}.
#'   `scale` can be any positive number.  Defaults to `1`.
#' @param shape The scale parameter \eqn{\alpha}.
#'   `shape` can be any positive number.  Defaults to `1`.
#'
#' @return A `Frechet` object.
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
#'   In the following, let \eqn{X} be a Frechet random variable with location
#'   parameter  `location` = \eqn{m}, scale parameter `scale` =
#'   \eqn{s}, and shape parameter `shape` = \eqn{\alpha}.
#'   A Frechet(\eqn{m, s, \alpha}) distribution is equivalent to a
#'   `\link{GEV}`(\eqn{m + s, s / \alpha, 1 / \alpha}) distribution.
#'
#'   **Support**: \eqn{(m, \infty)}.
#'
#'   **Mean**: \eqn{m + s\Gamma(1 - 1/\alpha)}, for \eqn{\alpha > 1}; undefined
#'   otherwise.
#'
#'   **Median**: \eqn{m + (\ln 2)^{-1/\alpha}}{m + (\ln 2)^(-1/\alpha)}.
#'
#'   **Variance**:
#'   \eqn{s^2 [\Gamma(1 - 2 / \alpha) - \Gamma(1 - 1 / \alpha)^2]}
#'   for \eqn{\alpha > 2}; undefined otherwise.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{f(x) = \alpha s ^ {-1} [(x - m) / s] ^ {-(1 + \alpha)}%
#'         \exp\{-[(x - m) / s] ^ {-\alpha} \}}{%
#'        f(x) = (\alpha / s) [(x - m) / s] ^ [-(1 + \alpha)]%
#'         exp{-[(x - m) / s] ^ (-\alpha)}}
#'   for \eqn{x > m}.  The p.d.f. is 0 outside the support.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{F(x) = \exp\{-[(x - m) / s] ^ {-\alpha} \}}{%
#'        F(x) = exp{-[(x - m) / s] ^ (-\alpha)}}
#'   for \eqn{x > m}.  The c.d.f. is 0 below the support and 1 above the
#'   support.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Frechet(0, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 0.7)
#' log_pdf(X, 0.7)
#'
#' cdf(X, 0.7)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
Frechet <- function(location = 0, scale = 1, shape = 1) {
  d <- list(location = location, scale = scale, shape = shape)
  class(d) <- c("Frechet", "distribution")
  d
}

#' @export
print.Frechet <- function(x, ...) {
  cat(glue("Frechet distribution (location = {x$location},
           scale = {x$scale}, shape = {x$shape})\n"))
}

#' Draw a random sample from a Frechet distribution
#'
#' @inherit Frechet examples
#'
#' @param d A `Frechet` object created by a call to [Frechet()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Frechet <- function(d, n = 1L, ...) {
  # Convert to the GEV parameterisation
  loc <- d$location + d$scale
  scale <- d$scale / d$shape
  shape <- 1 / d$shape
  revdbayes::rgev(n = n, loc = loc, scale = scale, shape = shape)
}

#' Evaluate the probability mass function of a Frechet distribution
#'
#' @inherit Frechet examples
#' @inheritParams random.Frechet
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Frechet <- function(d, x, ...) {
  # Convert to the GEV parameterisation
  loc <- d$location + d$scale
  scale <- d$scale / d$shape
  shape <- 1 / d$shape
  revdbayes::dgev(x = x, loc = loc, scale = scale, shape = shape)
}

#' @rdname pdf.Frechet
#' @export
#'
log_pdf.Frechet <- function(d, x, ...) {
  # Convert to the GEV parameterisation
  loc <- d$location + d$scale
  scale <- d$scale / d$shape
  shape <- 1 / d$shape
  revdbayes::dgev(x = x, loc = loc, scale = scale, shape = shape, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Frechet distribution
#'
#' @inherit Frechet examples
#' @inheritParams random.Frechet
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Frechet <- function(d, x, ...) {
  # Convert to the GEV parameterisation
  loc <- d$location + d$scale
  scale <- d$scale / d$shape
  shape <- 1 / d$shape
  revdbayes::pgev(q = x, loc = loc, scale = scale, shape = shape)
}

#' Determine quantiles of a Frechet distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Frechet examples
#' @inheritParams random.Frechet
#'
#' @param p A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Frechet <- function(d, p, ...) {
  # Convert to the GEV parameterisation
  loc <- d$location + d$scale
  scale <- d$scale / d$shape
  shape <- 1 / d$shape
  revdbayes::qgev(p = p, loc = loc, scale = scale, shape = shape)
}
