# things to sort out with the generics
#  - can i get stats::generics() to use ellipsis::check_dots_used()?
#  - pdf() conflict with grDevices::pdf()

#' Draw a random sample from a probability distribution
#'
#' @param x A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#' @param n The number of samples to draw. Should be a positive
#'   integer. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return Random samples drawn from the distriubtion `x`.
#'
#' @examples
#'
#' X <- Normal()
#'
#' random(X, 10)
#' @export
random <- function(x, n = 1L, drop = TRUE, ...) {
  ellipsis::check_dots_used()
  UseMethod("random")
}

#' Evaluate the probability density of a probability distribution
#'
#' For discrete distributions, the probability mass function. `pmf()`
#' is an alias.
#'
#' @inheritParams random
#'
#' @param d A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return Probabilities corresponding to the vector `x`.
#'
#' @examples
#'
#' X <- Normal()
#'
#' pdf(X, c(1, 2, 3, 4, 5))
#' pmf(X, c(1, 2, 3, 4, 5))
#'
#' log_pdf(X, c(1, 2, 3, 4, 5))
#' @export
pdf <- function(d, x, drop = TRUE, ...) {
  ellipsis::check_dots_used()
  UseMethod("pdf")
}

#' @rdname pdf
#' @export
log_pdf <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("log_pdf")
}

#' @rdname pdf
#' @export
pmf <- function(d, x, ...) {
  pdf(d, x, ...)
}

#' Evaluate the probability density of a probability distribution
#'
#' For discrete distributions, the probability mass function.
#'
#' @inheritParams random
#'
#' @param d A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return Probabilities corresponding to the vector `x`.
#'
#' @examples
#'
#' X <- Normal()
#'
#' cdf(X, c(1, 2, 3, 4, 5))
#' @export
cdf <- function(d, x, drop = TRUE, ...) {
  ellipsis::check_dots_used()
  UseMethod("cdf")
}

#' Compute the moments of a probability distribution
#'
#' The functions \code{variance}, \code{skewness}, and \code{kurtosis} are new
#' generic functions for computing moments of probability distributions such as
#' those provided in this package. Additionally, the probability distributions
#' from \pkg{distributions3} all have methods for the \code{\link[base]{mean}}
#' generic. Moreover, quantiles can be computed with methods for
#' \code{\link[stats]{quantile}}. For examples illustrating the usage with
#' probability distribution objects, see the manual pages of the respective
#' distributions, e.g., \code{\link{Normal}} or \code{\link{Binomial}} etc.
#'
#' @param x An object. The package provides methods for probability
#' distribution objects, e.g., those created by [Normal()] or [Beta()] etc.
#' @param ... Further arguments passed to or from other methods. Unevaluated
#' arguments will generate a warning to catch mispellings or other possible
#' errors.
#'
#' @return A numeric scalar
#' @seealso \code{\link[base]{mean}}, \code{\link[stats]{quantile}},
#' \code{\link{cdf}}, \code{\link{random}}
#' @export
#'
variance <- function(x, ...) {
  ellipsis::check_dots_used()
  UseMethod("variance")
}

#' @rdname variance
#' @export
skewness <- function(x, ...) {
  ellipsis::check_dots_used()
  UseMethod("skewness")
}

#' @rdname variance
#' @export
kurtosis <- function(x, ...) {
  ellipsis::check_dots_used()
  UseMethod("kurtosis")
}


#' Compute the likelihood of a probability distribution given data
#'
#' @param d A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#' @param x A vector of data to compute the likelihood.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return the likelihood
#'
#' @examples
#'
#' X <- Normal()
#'
#' likelihood(X, c(-1, 0, 0, 0, 3))
#' @export
likelihood <- function(d, x, ...) {
  exp(log_likelihood(d, x, ...))
}

#' Compute the log-likelihood of a probability distribution given data
#'
#' @inheritParams likelihood
#'
#' @return the log-likelihood
#'
#' @examples
#'
#' X <- Normal()
#'
#' log_likelihood(X, c(-1, 0, 0, 0, 3))
#' @export
log_likelihood <- function(d, x, ...) {
  sum(log_pdf(d, x, ...))
}

#' Fit a distribution to data
#'
#' Approximates an empirical distribution with a theoretical one
#'
#' @inheritParams likelihood
#'
#' @return A distribution (the same kind as `d`) where the parameters
#'   are the MLE estimates based on `x`.
#'
#' @examples
#'
#' X <- Normal()
#'
#' fit_mle(X, c(-1, 0, 0, 0, 3))
#' @export
fit_mle <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("fit_mle")
}

#' Compute the sufficient statistics of a distribution from data
#'
#' @inheritParams fit_mle
#'
#' @return a named list of sufficient statistics
#' @export
suff_stat <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("suff_stat")
}

#' Return the support of a distribution
#'
#' @param d A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#' @return A vector with two elements indicating the range of the support.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @export
support <- function(d, drop = TRUE) {
  if (!is_distribution(d)) {
    stop("d must be a supported distribution object")
  }

  UseMethod("support")
}
