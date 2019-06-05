#' Create a negative_binomial distribution
#'
#' TODO
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the binomial distribution reduces to the
#'   bernoulli distribution. Oftened called `n` in textbooks.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `negative_binomial` object.
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
#'   In the following, let \eqn{X} be a negative_binomial random variable with
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
#' nb <- negative_binomial(10, 0.3)
#' nb
#'
#' random(nb, 10)
#' pdf(nb, 2)
#' log_pdf(nb, 2)
#' cdf(nb, 4)
#' quantile(nb, 0.7)
#'
negative_binomial <- function(size, p = 0.5) {
  d <- list(size = size, p = p)
  class(d) <- c("negative_binomial", "distribution")
  d
}

#' @export
print.negative_binomial <- function(x, ...) {
  cat(glue("Negative Binomial distribution (size = {x$size}, p = {x$p})"))
}

#' Draw a random sample from a negative binomial distribution
#'
#' @inherit negative_binomial examples
#'
#' @param d A `negative_binomial` object created by a call to
#'   [negative_binomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family negative_binomial distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.negative_binomial <- function(d, n = 1L, ...) {
  rnbinom(n = n, size = d$size, prob = d$p)
}

#' Evaluate the probability mass function of a negative_binomial distribution
#'
#' @inherit negative_binomial examples
#' @inheritParams random.negative_binomial
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family negative_binomial distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.negative_binomial <- function(d, x, ...) {
  dnbinom(x = x, size = d$size, prob = d$p)
}

#' @rdname pdf.negative_binomial
#' @export
#'
log_pdf.negative_binomial <- function(d, x, ...) {
  dnbinom(x = x, size = d$size, prob = d$p, log = TRUE)
}

#' Evaluate the cumulative distribution function of a negative binomial distribution
#'
#' @inherit negative_binomial examples
#' @inheritParams random.negative_binomial
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family negative_binomial distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.negative_binomial <- function(d, x, ...) {
  pnbinom(q = x, size = d$size, prob = d$p)
}

#' Determine quantiles of a negative_binomial distribution
#'
#' @inherit negative_binomial examples
#' @inheritParams random.negative_binomial
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family negative_binomial distribution
#'
quantile.negative_binomial <- function(d, p, ...) {
  qnbinom(p = p, size = d$size, prob = d$p)
}
