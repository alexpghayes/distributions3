# things to sort out with the generics
#  - can i get stats::generics() to use ellipsis::check_dots_used()?
#  - pdf() conflict with grDevices::pdf()

#' Draw a random sample from a probability distribution
#'
#' @param d A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#' @param n The number of samples to draw. Should be a positive
#'   integer. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @examples
#'
#' X <- Normal()
#'
#' random(X, 10)
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
#'
#' @examples
#'
#' X <- Normal()
#'
#' pdf(X, c(1, 2, 3, 4, 5))
#' pmf(X, c(1, 2, 3, 4, 5))
#'
#' log_pdf(X, c(1, 2, 3, 4, 5))
#'
#' @export
pdf <- function(d, x, ...) {
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
#' For discrete distributions, the probabilty mass function.
#'
#' @inheritParams random
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#'
#' @return A vector of probabilities, one for each element of `x`.
#'
#' @examples
#'
#' X <- Normal()
#'
#' cdf(X, c(1, 2, 3, 4, 5))
#'
#' @export
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
#' @examples
#'
#' X <- Normal()
#'
#' cdf(X, c(0.2, 0.4, 0.6, 0.8))
#'
#' @export
quantile <- function(d, p, ...) {
  ellipsis::check_dots_used()
  UseMethod("quantile")
}

#' Compute the moments of a probability distribution
#'
#' @param d A probability distribution object such as those created by
#'   a call to [Bernoulli()], [Beta()], or [Binomial()].
#'
#' @return A numeric scalar
#' @export
#'
variance <- function(d, ...) {
  ellipsis::check_dots_used()
  UseMethod("variance")
}

#' @rdname variance
#' @export
skewness <- function(d, ...) {
  ellipsis::check_dots_used()
  UseMethod("skewness")
}

#' @rdname variance
kurtosis <- function(d, ...) {
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
#'
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
#'
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
#'
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
suff_stat <- function(d, x, ...) {
  ellipsis::check_dots_used()
  UseMethod("suff_stat")
}
