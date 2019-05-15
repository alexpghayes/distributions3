#' Create a geometric distribution
#'
#' The geometric distribution can be thought of as a generalization
#' of the [bernoulli()] distribution where ask: "if I keep flipping a
#' coin with probability `p` of heads, what is the probability I need
#' \eqn{k} flips before I get my first heads?" The geometric
#' distribution is a special case of negative binomial distribution
#' with parameters TODO.
#'
#' @param p The success probability for the distribution. `p` can be
#'   any value in [0, 1], and defaults to `0.5`.
#'
#' @return A `geometric` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a geometric random variable with
#'   success probability `p` = \eqn{p}.
#'
#'   TODO: multiple parameterizations BLEH
#'
#'   **Support**: TODO
#'
#'   **Mean**: TODO
#'
#'   **Variance**: TODO
#'
#'   **Probability density function (p.d.f)**:
#'
#'   TODO
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   TODO
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   TODO
#'
#' @examples
#'
#' g <- geometric(0.3)
#' g
#'
#' random(g, 10)
#' pdf(g, 2)
#' cdf(g, 4)
#' quantile(g, 0.7)
#'
geometric <- function(p = 0.5) {
  d <- list(p = p)
  class(d) <- "geometric"
  d
}

#' @export
print.geometric <- function(x, ...) {
  cat(glue("Geometric distribution (p = {x$p})"))
}

#' Draw a random sample from a geometric distribution
#'
#' Please see the documentation of [geometric()] for some properties
#' of the geometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit geometric examples
#'
#' @param d A `geometric` object created by a call to [geometric()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family geometric distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.geometric <- function(d, n = 1L, ...) {
  rgeom(n = n, prob = d$p)
}

#' Evaluate the probability mass function of a geometric distribution
#'
#' Please see the documentation of [geometric()] for some properties
#' of the geometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit geometric examples
#' @inheritParams random.geometric
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family geometric distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.geometric <- function(d, x, ...) {
  dgeom(x = x, prob = d$p)
}

#' Evaluate the cumulative distribution function of a geometric distribution
#'
#' @inherit geometric examples
#' @inheritParams random.geometric
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family geometric distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.geometric <- function(d, x, ...) {
  pgeom(q = x, prob = d$p)
}

#' Determine quantiles of a geometric distribution
#'
#' @inherit geometric examples
#' @inheritParams random.geometric
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family geometric distribution
#'
quantile.geometric <- function(d, p, ...) {
  qgeom(p = p, prob = d$p)
}
