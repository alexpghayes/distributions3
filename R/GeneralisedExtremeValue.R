#' Create a Generalised Extreme Value (GEV) distribution
#'
#' The GEV distribution arises from the Extremal Types Theorem, which is rather
#' like the Central Limit Theorem (see \code{\link{Normal}}) but it relates to
#' the *maximum* of \eqn{n} i.i.d. random variables rather than to the sum.
#' If, after a suitable linear rescaling, the distribution of this maximum
#' tends to a non-degenerate limit as \eqn{n} tends to infinity then this limit
#' must be a GEV distribution. The requirement that the variables are independent
#' can be relaxed substantially. Therefore, the GEV distribution is often used
#' to model the maximum of a large number of random variables.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   `mu` can be any real number.  Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   `sigma` can be any positive number.  Defaults to `1`.
#' @param xi The shape parameter, written \eqn{\xi} in textbooks.
#'   `xi` can be any real number.  Defaults to `0`, which corresponds to a
#'   Gumbel distribution.
#'
#' @return A `GEV` object.
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
#'   In the following, let \eqn{X} be a GEV random variable with location
#'   parameter  `mu` = \eqn{\mu}, scale parameter `sigma` = \eqn{\sigma} and
#'   shape parameter `xi` = \eqn{\xi}.
#'
#'   **Support**:
#'   \eqn{(-\infty, \mu - \sigma / \xi]} for \eqn{\xi < 0};
#'   \eqn{[\mu - \sigma / \xi, \infty)} for \eqn{\xi > 0};
#'   and \eqn{R}, the set of all real numbers, for \eqn{\xi = 0}.
#'
#'   **Mean**: \eqn{\mu + \sigma[\Gamma(1 - \xi) - 1]/\xi} for
#'   \eqn{\xi < 1, \xi \neq 0}{\xi < 1, \xi != 0};
#'   \eqn{\mu + \sigma\gamma} for \eqn{\xi = 0}, where \eqn{\gamma}
#'   is Euler's constant, approximately equal to 0.57722; undefined otherwise.
#'
#'   **Median**: \eqn{\mu + \sigma[(\ln 2) ^ {-\xi} - 1]/\xi}{%
#'   \mu + \sigma[(ln 2)^(-\xi) - 1] / \xi} for \eqn{\xi \neq 0}{\xi != 0};
#'   \eqn{\mu - \sigma\ln(\ln 2)}{\mu - \sigma ln(ln 2)} for \eqn{\xi = 0}.
#'
#'   **Variance**:
#'   \eqn{\sigma^2 [\Gamma(1 - 2 \xi) - \Gamma(1 - \xi)^2] / \xi^2}
#'   for \eqn{\xi < 1 / 2, \xi \neq 0}{\xi < 1 / 2, \xi != 0};
#'   \eqn{\sigma^2 \pi^2 / 6} for \eqn{\xi = 0}; undefined otherwise.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{f(x) = \sigma ^ {-1} [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}%
#'          \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        f(x) = (1 / \sigma) [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}%
#'          exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)}}
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The p.d.f. is 0 outside the
#'   support.
#'
#'   In the \eqn{\xi = 0} (Gumbel) special case
#'   \deqn{f(x) = \sigma ^ {-1} \exp[-(x - \mu) / \sigma]%
#'         \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'        f(x) = (1 / \sigma) exp[-(x - \mu) / \sigma]%
#'         exp{-exp[-(x - \mu) / \sigma]}}
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{F(x) = \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        F(x) = exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The c.d.f. is 0 below the
#'   support and 1 above the support.
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
#' X <- GEV(1, 2, 0.1)
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
GEV <- function(mu = 0, sigma = 1, xi = 0) {
  d <- list(mu = mu, sigma = sigma, xi = xi)
  class(d) <- c("GEV", "distribution")
  d
}

#' @export
print.GEV <- function(x, ...) {
  cat(glue("GEV distribution (mu = {x$mu}, sigma = {x$sigma}, xi = {x$xi})\n"))
}

#' Draw a random sample from a GEV distribution
#'
#' @inherit GEV examples
#'
#' @param d A `GEV` object created by a call to [GEV()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector containing values in ??? of length `n`.
#' @export
#'
random.GEV <- function(d, n = 1L, ...) {
  revdbayes::rgev(n = n, loc = d$mu, scale = d$sigma, shape = d$xi)
}

#' Evaluate the probability mass function of a GEV distribution
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.GEV <- function(d, x, ...) {
  revdbayes::dgev(x = x, loc = d$mu, scale = d$sigma, shape = d$xi)
}

#' @rdname pdf.GEV
#' @export
#'
log_pdf.GEV <- function(d, x, ...) {
  revdbayes::dgev(x = x, loc = d$mu, scale = d$sigma, shape = d$xi, log = TRUE)
}

#' Evaluate the cumulative distribution function of a GEV distribution
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GEV <- function(d, x, ...) {
  revdbayes::pgev(q = x, loc = d$mu, scale = d$sigma, shape = d$xi)
}

#' Determine quantiles of a GEV distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.GEV <- function(d, p, ...) {
  revdbayes::qgev(p = p, loc = d$mu, scale = d$sigma, shape = d$xi)
}