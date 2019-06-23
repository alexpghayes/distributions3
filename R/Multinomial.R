#' Create a Multinomial distribution
#'
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the Multinomial distribution reduces to the
#'   categorical distribution (also called the discrete uniform).
#'   Often called `n` in textbooks.
#' @param p A vector of success probabilities for each trial. `p` can
#'   take on any positive value, and the vector is normalized internally.
#'
#' @return A `Multinomial` object.
#' @export
#'
#' @family discrete distributions
#' @family multivariate distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Multinomial random variable with
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
#' X <- Multinomial(size = 5, p = c(0.3, 0.4, 0.2, 0.1))
#' X
#'
#' random(X, 10)
#'
#' # pdf(X, 2)
#' # log_pdf(X, 2)
#'
#'
Multinomial <- function(size, p) {
  d <- list(size = size, p = p)
  class(d) <- c("Multinomial", "multivariate", "distribution")
  d
}

#' @export
print.Multinomial <- function(x, ...) {
  cat(glue("Multinomial distribution (size = {x$size}, p = {x$p})"))
}

#' Draw a random sample from a Multinomial distribution
#'
#' @inherit Multinomial examples
#'
#' @param d A `Multinomial` object created by a call to [Multinomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Multinomial distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Multinomial <- function(d, n = 1L, ...) {
  rmultinom(n = n, size = d$size, prob = d$size)
}

#' Evaluate the probability mass function of a Multinomial distribution
#'
#' Please see the documentation of [Multinomial()] for some properties
#' of the Multinomial distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Multinomial examples
#' @inheritParams random.Multinomial
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Multinomial distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Multinomial <- function(d, x, ...) {
  dmultinom(x = x, size = d$size, prob = d$size)
}

#' @rdname pdf.Multinomial
#' @export
log_pdf.Multinomial <- function(d, x, ...) {
  dmultinom(x = x, size = d$size, prob = d$size, log = TRUE)
}
