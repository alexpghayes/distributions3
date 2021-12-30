#' Create a Generalised Pareto (GP) distribution
#'
#' The GP distribution has a link to the `\link{GEV}` distribution.
#' Suppose that the maximum of \eqn{n} i.i.d. random variables has
#' approximately a GEV distribution. For a sufficiently large threshold
#' \eqn{u}, the conditional distribution of the amount (the threshold
#' excess) by which a variable exceeds \eqn{u} given that it exceeds \eqn{u}
#' has approximately a GP distribution.  Therefore, the GP distribution is
#' often used to model the threshold excesses of a high threshold \eqn{u}.
#' The requirement that the variables are independent can be relaxed
#' substantially, but then exceedances of \eqn{u} may cluster.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   `mu` can be any real number.  Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   `sigma` can be any positive number.  Defaults to `1`.
#' @param xi The shape parameter, written \eqn{\xi} in textbooks.
#'   `xi` can be any real number.  Defaults to `0`, which corresponds to a
#'   Gumbel distribution.
#'
#' @return A `GP` object.
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
#'   In the following, let \eqn{X} be a GP random variable with location
#'   parameter  `mu` = \eqn{\mu}, scale parameter `sigma` = \eqn{\sigma} and
#'   shape parameter `xi` = \eqn{\xi}.
#'
#'   **Support**:
#'   \eqn{[\mu, \mu - \sigma / \xi]} for \eqn{\xi < 0};
#'   \eqn{[\mu, \infty)} for \eqn{\xi \geq 0}{\xi >= 0}.
#'
#'   **Mean**: \eqn{\mu + \sigma/(1 - \xi)} for
#'   \eqn{\xi < 1}; undefined otherwise.
#'
#'   **Median**: \eqn{\mu + \sigma[2 ^ \xi - 1]/\xi}{%
#'   \mu + \sigma[2^\xi - 1] / \xi} for \eqn{\xi \neq 0}{\xi != 0};
#'   \eqn{\mu + \sigma\ln 2}{\mu + \sigma ln2} for \eqn{\xi = 0}.
#'
#'   **Variance**:
#'   \eqn{\sigma^2 / (1 - \xi)^2 (1 - 2\xi)}
#'   for \eqn{\xi < 1 / 2}; undefined otherwise.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{f(x) = \sigma^{-1} [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}}{%
#'        f(x) = (1 / \sigma) [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}}
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The p.d.f. is 0 outside the
#'   support.
#'
#'   In the \eqn{\xi = 0} special case
#'   \deqn{f(x) = \sigma ^ {-1} \exp[-(x - \mu) / \sigma]}{%
#'        f(x) = (1 / \sigma) exp[-(x - \mu) / \sigma]}
#'   for \eqn{x} in [\eqn{\mu, \infty}).  The p.d.f. is 0 outside the support.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{F(x) = 1 - \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        F(x) = 1 - exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The c.d.f. is 0 below the
#'   support and 1 above the support.
#'
#'   In the \eqn{\xi = 0} special case
#'   \deqn{F(x) = 1 - \exp[-(x - \mu) / \sigma] \}}{%
#'         F(x) = 1 - exp[-(x - \mu) / \sigma] }
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GP(0, 2, 0.1)
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
GP <- function(mu = 0, sigma = 1, xi = 0) {
  if (sigma <= 0) {
    stop("sigma must be positive")
  }
  d <- list(mu = mu, sigma = sigma, xi = xi)
  class(d) <- c("GP", "distribution")
  d
}

#' @export
print.GP <- function(x, ...) {
  cat(glue("GP distribution (mu = {x$mu}, sigma = {x$sigma}, xi = {x$xi})\n"))
}

#' @export
mean.GP <- function(x, ...) {
  ellipsis::check_dots_used()
  mu <- x$mu
  sigma <- x$sigma
  xi <- x$xi

  if (xi < 1) mu + sigma / (1 - xi)
  else Inf
}

#' @export
variance.GP <- function(x, ...) {
  sigma <- x$sigma
  xi <- x$xi

  if (xi < 1/2) sigma^2 / ((1 - xi)^2 - (1 - 2*xi))
  else Inf
}

#' @export
skewness.GP <- function(x, ...) {
  xi <- x$xi

  if (xi < 1/3) 2*(1 + xi) * sqrt(1 - 2*xi) / (1 - 3*xi)
  else Inf
}

#' @export
kurtosis.GP <- function(x, ...) {
  xi <- x$xi

  if (xi < 1/4) {
    k1 <- (1 - 2*xi) * (2*xi^2 + xi + 3)
    k2 <- (1 - 3*xi) * (1 - 4*xi)
    3 * k1 / k2 - 3
  } else {
    Inf
  }
}

#' Draw a random sample from a GP distribution
#'
#' @inherit GP examples
#'
#' @param x A `GP` object created by a call to [GP()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.GP <- function(x, n = 1L, ...) {
  revdbayes::rgp(n = n, loc = x$mu, scale = x$sigma, shape = x$xi)
}

#' Evaluate the probability mass function of a GP distribution
#'
#' @inherit GP examples
#'
#' @param d A `GP` object created by a call to [GP()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.GP <- function(d, x, ...) {
  revdbayes::dgp(x = x, loc = d$mu, scale = d$sigma, shape = d$xi)
}

#' @rdname pdf.GP
#' @export
#'
log_pdf.GP <- function(d, x, ...) {
  revdbayes::dgp(x = x, loc = d$mu, scale = d$sigma, shape = d$xi, log = TRUE)
}

#' Evaluate the cumulative distribution function of a GP distribution
#'
#' @inherit GP examples
#'
#' @param d A `GP` object created by a call to [GP()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.GP <- function(d, x, ...) {
  revdbayes::pgp(q = x, loc = d$mu, scale = d$sigma, shape = d$xi)
}

#' Determine quantiles of a GP distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GP examples
#' @inheritParams random.GP
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.GP <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  revdbayes::qgp(p = probs, loc = x$mu, scale = x$sigma, shape = x$xi)
}
