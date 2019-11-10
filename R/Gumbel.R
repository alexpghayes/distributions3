#' Create a Gumbel distribution
#'
#' The Gumbel distribution is a special case of the `\link{GEV}` distribution.
#' It may be referred to as a type I extreme value distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   `mu` can be any real number.  Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   `sigma` can be any positive number.  Defaults to `1`.
#'
#' @return A `Gumbel` object.
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
#'   In the following, let \eqn{X} be a Gumbel random variable with location
#'   parameter  `mu` = \eqn{\mu}, scale parameter `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers.
#'
#'   **Mean**: \eqn{\mu + \sigma\gamma}, where \eqn{\gamma} is Euler's
#'   constant, approximately equal to 0.57722; undefined otherwise.
#'
#'   **Median**: \eqn{\mu - \sigma\ln(\ln 2)}{\mu - \sigma ln(ln 2)}.
#'
#'   **Variance**: \eqn{\sigma^2 \pi^2 / 6}.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{f(x) = \sigma ^ {-1} \exp[-(x - \mu) / \sigma]%
#'         \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'        f(x) = (1 / \sigma) exp[-(x - \mu) / \sigma]%
#'         exp{-exp[-(x - \mu) / \sigma]}}
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   In the \eqn{\xi = 0} (Gumbel) special case
#'   \deqn{F(x) = \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'         F(x) = exp{ - exp[-(x - \mu) / \sigma]} }
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Gumbel(1, 2)
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
Gumbel <- function(mu = 0, sigma = 1) {
  d <- list(mu = mu, sigma = sigma)
  class(d) <- c("Gumbel", "distribution")
  d
}

#' @export
print.Gumbel <- function(x, ...) {
  cat(glue("Gumbel distribution (mu = {x$mu}, sigma = {x$sigma})\n"))
}

#' Draw a random sample from a Gumbel distribution
#'
#' @inherit Gumbel examples
#'
#' @param d A `Gumbel` object created by a call to [Gumbel()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Gumbel <- function(d, n = 1L, ...) {
  revdbayes::rgev(n = n, loc = d$mu, scale = d$sigma, shape = 0)
}

#' Evaluate the probability mass function of a Gumbel distribution
#'
#' @inherit Gumbel examples
#' @inheritParams random.Gumbel
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Gumbel <- function(d, x, ...) {
  revdbayes::dgev(x = x, loc = d$mu, scale = d$sigma, shape = 0)
}

#' @rdname pdf.Gumbel
#' @export
#'
log_pdf.Gumbel <- function(d, x, ...) {
  revdbayes::dgev(x = x, loc = d$mu, scale = d$sigma, shape = 0, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Gumbel distribution
#'
#' @inherit Gumbel examples
#' @inheritParams random.Gumbel
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Gumbel <- function(d, x, ...) {
  revdbayes::pgev(q = x, loc = d$mu, scale = d$sigma, shape = 0)
}

#' Determine quantiles of a Gumbel distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Gumbel examples
#' @inheritParams random.Gumbel
#'
#' @param p A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Gumbel <- function(d, p, ...) {
  revdbayes::qgev(p = p, loc = d$mu, scale = d$sigma, shape = 0)
}
