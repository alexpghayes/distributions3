#' Create a exponential distribution
#'
#' @param rate The rate parameter. Can be any positive number. Defaults
#'   to `1`.
#'
#' @return A `exponential` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' e <- exponential(5)
#' e
#'
#' random(e, 10)
#' pdf(e, 2)
#' log_pdf(e, 2)
#' cdf(e, 4)
#' quantile(e, 0.7)
#'
#' cdf(e, quantile(e, 0.7))
#' quantile(e, cdf(e, 7))
#'
exponential <- function(rate = 1) {
  d <- list(rate = rate)
  class(d) <- "exponential"
  d
}

#' @export
print.exponential <- function(x, ...) {
  cat(glue("Exponential distribution (rate = {x$rate})"))
}

#' Draw a random sample from a exponential distribution
#'
#' @inherit exponential examples
#'
#' @param d A `exponential` object created by a call to [exponential()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.exponential <- function(d, n = 1L, ...) {
  rexp(n = n, rate = d$rate)
}

#' Evaluate the probability mass function of a exponential distribution
#'
#' @inherit exponential examples
#' @inheritParams random.exponential
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.exponential <- function(d, x, ...) {
  dexp(x = x, rate = d$rate)
}

#' @rdname pdf.exponential
#' @export
#'
log_pdf.exponential <- function(d, x, ...) {
  dexp(x = x, rate = d$rate, log = TRUE)
}

#' Evaluate the cumulative distribution function of a exponential distribution
#'
#' @inherit exponential examples
#' @inheritParams random.exponential
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.exponential <- function(d, x, ...) {
  pexp(q = x, rate = d$rate)
}

#' Determine quantiles of a exponential distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit exponential examples
#' @inheritParams random.exponential
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.exponential <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qexp(p = p, rate = d$rate)
}

#' Fit an exponential distribution to data
#'
#' @param d An `exponential` object created by a call to [exponential()].
#' @param x A vector of data.
#'
#' @family exponential distribution
#'
#' @return An `exponential` object.
#' @export
fit_mle.exponential <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  exponential(ss$sum / ss$samples)
}


#' Compute the sufficient statistics of an exponential distribution from data
#'
#' @inheritParams fit_mle.exponential
#'
#' @return A named list of the sufficient statistics of exponential distribution
#'   \describe{
#'     \item{\code{sum}}{The sum of the data}
#'     \item{\code{samples}}{The number of samples in the data}
#'   }
#'
#' @export
suff_stat.exponential <- function(d, x, ...) {
  valid_x <- (x > 0)
  if(any(!valid_x)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), samples = length(x))
}
