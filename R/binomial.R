#' Create a binomial distribution
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the binomial distribution reduces to the
#'   bernoulli distribution. Oftened called `n` in textbooks.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'   greater than zero. Defaults to `1`.
#'
#' @return A `binomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @examples
#'
#' b <- binomial(10, 0.2)
#' b
#'
#' random(b, 10)
#' pdf(b, 2L)
#' cdf(b, 4L)
#' quantile(b, 7)
#'
#' cdf(b, quantile(b, 7))
#' quantile(b, cdf(b, 7))
#'
binomial <- function(size, p) {

  stopifnot(is.integer(size))

  d <- list(size = size, p = p)
  class(d) <- "binomial"
  d
}

print.binomial <- function(d) {
  cat(glue("Binomial distribution (size = {d$size}, p = {d$p})"))
}

#' Draw a random sample from a binomial distribution
#'
#' @inherit binomial examples
#'
#' @param d A `binomial` object created by a call to [binomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#'
#' @return An integer vector containing values between `0` and `d$size`
#'   of length `n`.
#' @export
#'
random.binomial <- function(d, n = 1L) {
  rbinom(n = n, size = d$size, prob = d$p)
}

#' Evaluate the probability mass function of a binomial distribution
#'
#' @inherit binomial examples
#' @inheritParams random.binomial
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.binomial <- function(d, x) {
  dbinom(x = x, size = d$size, prob = d$p)
}

#' Evaluate the cumulative distribution function of a binomial distribution
#'
#' @inherit binomial examples
#' @inheritParams random.binomial
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.beta <- function(d, x) {
  pbinom(x = x, size = d$size, prob = d$p)
}

#' Determine quantiles of a binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit binomial examples
#' @inheritParams random.binomial
#'
#' @param p A vector of probabilites.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.binomial <- function(d, p) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qbinom(x = x, size = d$size, prob = d$p)
}
