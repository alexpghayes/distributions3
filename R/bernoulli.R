#' Create a bernoulli distribution
#'
#' Bernoulli distributions are used to represent events like coin flips
#' when there is single trial that is either successful or unsuccessful.
#' The bernoulli distribution is a special case of the [binomial()]
#' distribution with `n = 1`.
#'
#' @param p The success probability for the distribution. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `bernoulli` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the
#'   math will render nicely.
#'
#'   In the following, let \eqn{X} be a Bernoulli random variable with parameter
#'   `p` = \eqn{p}. Some textbooks also define \eqn{q = 1 - p}, or use
#'   \eqn{\pi} instead of \eqn{p}.
#'
#'   **Support**: \eqn{0, 1}
#'
#'   **Mean**: \eqn{p}
#'
#'   **Variance**: \eqn{p \cdot (1 - p) = p \cdot q}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{P(X = x) = p^x (1 - p)^{1-x} = p^x q^{1-x}}
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   TODO ugh formatting
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{\mathbb{E}(e^{tX}) = (1 - p) + p e^t}
#'
#' @examples
#'
#' b <- bernoulli(0.7)
#' b
#'
#' random(b, 10)
#' pdf(b, 1)
#' cdf(b, 0)
#' quantile(b, 0.7)
#'
#' # TODO: looks like I don't quite have quantiles right since the
#' # inverses don't hold like I'd want
#'
#' cdf(b, quantile(b, 0.7))
#' quantile(b, cdf(b, 0.7))
#'
bernoulli <- function(p = 0.5) {

  # TODO: check that 0 <= p <= 1

  d <- list(p = p)
  class(d) <- "bernoulli"
  d
}

#' @export
print.bernoulli <- function(d) {
  cat(glue("Bernoulli distribution (p = {d$p})"))
}

#' Draw a random sample from a bernoulli distribution
#'
#' @inherit bernoulli examples
#'
#' @param d A `bernoulli` object created by a call to [bernoulli()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return An integer vector of zeros and ones of length `n`.
#' @export
#'
random.bernoulli <- function(d, n = 1L, ...) {
  rbinom(n = n, size = 1, prob = d$p)
}

#' Evaluate the probability mass function of a bernoulli distribution
#'
#' @inherit bernoulli examples
#' @inheritParams random.bernoulli
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.bernoulli <- function(d, x, ...) {
  # TODO: call as.integer()? on doubles
  # TODO: out of support: return zero, return zero and warn, or stop?
  if (x == 0) 1 - d$p else d$p
}

#' Evaluate the cumulative distribution function of a bernoulli distribution
#'
#' @inherit bernoulli examples
#' @inheritParams random.bernoulli
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.bernoulli <- function(d, x, ...) {
  pbinom(q = x, size = 1, prob = d$p)
}

#' Determine quantiles of a bernoulli distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit bernoulli examples
#' @inheritParams random.bernoulli
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.bernoulli <- function(d, p, ...) {
  qbinom(p = p, size = 1, prob = d$p)
}
