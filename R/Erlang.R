#' Create an Erlang distribution
#'
#' The Erlang distribution is a two-parameter family of continuous probability
#' distributions with support \eqn{x \in [0,\infty)}.
#' The two parameters are a positive integer shape parameter \eqn{k} and a
#' positive real rate parameter \eqn{\lambda}.
#' The Erlang distribution with shape parameter \eqn{k = 1} simplifies to the
#' exponential distribution, and it is a special case of the gamma distribution.
#' It corresponds to a sum of \eqn{k} independent exponential variables with mean
#' \eqn{1 / \lambda} each.
#'
#' @param k The shape parameter. Can be any positive integer number.
#' @param lambda The rate parameter. Can be any positive number.
#' @return An `Erlang` object.
#' @export
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Erlang(5, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Erlang <- function(k, lambda) {
  d <- list(k = k, lambda = lambda)
  class(d) <- c("Erlang", "distribution")
  d
}

#' @export
print.Erlang <- function(x, ...) {
  cat(glue("Erlang distribution (k = {x$k}, lambda = {x$lambda})\n"))
}

#' Draw a random sample from an Erlang distribution
#'
#' @inherit Erlang examples
#'
#' @param d An `Erlang` object created by a call to [Erlang()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Erlang <- function(d, n = 1L, ...) {
  replicate(n = n, expr = {
    -1 / d$lambda * log(prod(runif(n = d$k)))
  })
}

#' Evaluate the probability mass function of an Erlang distribution
#'
#' @inherit Erlang examples
#' @inheritParams random.Erlang
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Erlang <- function(d, x, ...) {
  if (any(x < 0)) stop("'x' must be non-negative", call. = FALSE)
  (d$lambda^d$k) * (x^(d$k - 1)) * exp(-d$lambda * x) / factorial(d$k - 1)
}

#' @rdname pdf.Erlang
#' @export
#'
log_pdf.Erlang <- function(d, x, ...) {
  if (any(x < 0)) stop("'x' must be non-negative", call. = FALSE)
  d$k * log(d$lambda) + (d$k - 1) * log(x) - d$lambda * x - log(factorial(d$k - 1))
}

#' Evaluate the cumulative distribution function of an Erlang distribution
#'
#' @inherit Erlang examples
#' @inheritParams random.Erlang
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Erlang <- function(d, x, ...) {
  if (any(x < 0)) stop("'x' must be non-negative", call. = FALSE)
  internal <- Vectorize(FUN = function(d, x, ...) {
    summation <- vector(mode = "numeric", length = d$k)
    n <- 0:(d$k - 1)
    summation <- 1 / factorial(n) * exp(-d$lambda * x) * (d$lambda * x)^n
    return(1 - sum(summation))
  }, vectorize.args = "x")
  internal(d = d, x = x, ...)
}

#' Determine quantiles of an Erlang distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Erlang examples
#' @inheritParams random.Erlang
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#' @param interval Interval being used to search for the quantile using numerical root finding. Defaults to (0, 1e6)
#' @param tol Tolerance of the root finding algorithm. Defaults to `.Machine$double.eps`
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Erlang <- function(d, p, ..., interval = c(0, 1e6), tol = .Machine$double.eps) {
  if (any(p < 0) | any(p > 1)) stop("'p' must be between 0 and 1.", call. = TRUE)
  p[p == 1] <- (1 - .Machine$double.eps^0.25)
  internal <- Vectorize(FUN = function(d, p, ..., interval, tol) {
    qf <- function(x) cdf(d = d, x = x) - p
    root <- stats::uniroot(qf, interval = interval, tol = tol, check.conv = TRUE)
    return(root$root)
  }, vectorize.args = "p")
  internal(d = d, p = p, ..., interval = interval, tol = tol)
}

#' Return the support of the Erlang distribution
#'
#' @param d An `Erlang` object created by a call to [Erlang()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Erlang <- function(d) c(0, Inf)
