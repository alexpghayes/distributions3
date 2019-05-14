# things to sort out with the generics
#  - can i get stats::generics() to use ellipsis::check_dots_used()?
#  - pdf() conflict with grDevices::pdf()

#' Draw a random sample from a probability distribution
#'
#' @param d A probability distribution object such as those created by
#'   a call to [bernoulli()], [beta()], or [binomial()].
#' @param n The number of samples to draw. Should be a positive
#'   integer. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @export
random <- function(d, n = 1L, ...) {
  ellipsis::check_dots_used()
  UseMethod("random")
}

#' Evaluate the probability density of a probability distribution
#'
#' For discrete distributions, the probabilty mass function. `pmf()`
#' is an alias.
#'
#' @inheritParams random
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("pdf")
}

#' @rdname pdf
#' @export
pmf <- function(d, x, ...) {
  pdf(d, x, ...)
}

#' Evaluate the probability density of a probability distribution
#'
#' For discrete distributions, the probabilty mass function.
#'
#' @inheritParams random
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("cdf")
}

#' Find the quantile of a probability distribution
#'
#' TODO: Note that this current masks the [stats::quantile()] generic
#' to allow for consistent argument names and warnings when arguments
#' disappear into `...`.
#'
#' @inheritParams random
#'
#' @param p A vector of probabilites.
#'
#' @return A vector of quantiles, one for each element of `p`.
#'
#' @export
quantile <- function(d, p, ...) {
  ellipsis::check_dots_used()
  UseMethod("quantile")
}


#' Fit a distribution to data
#'
#' Approximates an empirical distribution with a theoretical one
#'
#' @param d A probability distribution object such as those created by
#'   a call to [bernoulli()], [beta()], or [binomial()].
#' @param x A vector of data to estimate the parameters of the
#'   distribution with
#'
#' @return an object the same class as \code{d} with updated parameters
#' @export
fit <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("fit")
}
