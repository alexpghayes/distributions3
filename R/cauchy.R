#' Create a cauchy distribution
#'
#' Note that the cauchy distribution is the student's t distribution
#' with one degree of freedom. The cauchy distribution does not have
#' a well defined mean or variance. Cauchy distributions often appear
#' as priors in Bayesian contexts due to their heavy tails.
#'
#' @param location The location parameter. Can be any real number. Defaults
#'   to `0`.
#' @param scale The scale parameter. Must be greater than zero (?). Defaults
#'   to `1`.
#'
#' @return A `cauchy` object.
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
#'   In the following, let \eqn{X} be a Cauchy variable with mean
#'   `location =` \eqn{x_0} and `scale` = \eqn{\gamma}.
#'
#'   **Support**: \eqn{\mathbb{R}}{R}, the set of all real numbers
#'
#'   **Mean**: Undefined.
#'
#'   **Variance**: Undefined.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\pi \gamma \left[1 + \left(\frac{x - x_0}{\gamma} \right)^2 \right]}
#'   }{
#'     f(x) = 1 / (\pi \gamma (1 + ((x - x_0) / \gamma)^2)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(t) = \frac{1}{\pi} \arctan \left( \frac{t - x_0}{\gamma} \right) +
#'       \frac{1}{2}
#'   }{
#'     F(t) = arctan((t - x_0) / \gamma) / \pi + 1/2
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   Does not exist.
#'
#' @examples
#'
#' c <- cauchy(10, 0.2)
#' c
#'
#' random(c, 10)
#' pdf(c, 2)
#' log_pdf(c, 2)
#' cdf(c, 2)
#' quantile(c, 0.7)
#'
#' cdf(c, quantile(c, 0.7))
#' quantile(c, cdf(c, 7))
#'
cauchy <- function(location = 0, scale = 1) {
  d <- list(location = location, scale = scale)
  class(d) <- "cauchy"
  d
}

#' @export
print.cauchy <- function(x, ...) {
  cat(glue("Cauchy distribution (location = {x$location}, scale = {x$scale})"))
}

#' Draw a random sample from a cauchy distribution
#'
#' @inherit cauchy examples
#'
#' @param d A `cauchy` object created by a call to [cauchy()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.cauchy <- function(d, n = 1L, ...) {
  rcauchy(n = n, location = d$location, scale = d$scale)
}

#' Evaluate the probability mass function of a cauchy distribution
#'
#' @inherit cauchy examples
#' @inheritParams random.cauchy
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.cauchy <- function(d, x, ...) {
  dcauchy(x = x, location = d$location, scale = d$scale)
}

#' @rdname pdf.cauchy
#' @export
#'
log_pdf.cauchy <- function(d, x, ...) {
  dcauchy(x = x, location = d$location, scale = d$scale, log = TRUE)
}

#' Evaluate the cumulative distribution function of a cauchy distribution
#'
#' @inherit cauchy examples
#' @inheritParams random.cauchy
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.cauchy <- function(d, x, ...) {
  pcauchy(q = x, location = d$location, scale = d$scale)
}

#' Determine quantiles of a cauchy distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit cauchy examples
#' @inheritParams random.cauchy
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.cauchy <- function(d, p, ...) {
  qcauchy(p = p, location = d$location, scale = d$scale)
}
