#' Create a Exponential distribution
#'
#' @param rate The rate parameter. Can be any positive number. Defaults
#'   to `1`.
#'
#' @return A `Exponential` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' X <- Exponential(5)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
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
#'
Exponential <- function(rate = 1) {
  d <- list(rate = rate)
  class(d) <- c("Exponential", "distribution")
  d
}

#' @export
print.Exponential <- function(x, ...) {
  cat(glue("Exponential distribution (rate = {x$rate})"))
}

#' @export
mean.Exponential <- function(d, ...) d$rate ^ -1

#' @export
variance.Exponential <- function(d, ...) d$rate ^ 2

#' @export
skewness.Exponential <- function(d, ...) 2

#' @export
kurtosis.Exponential <- function(d, ...) 6

#' Draw a random sample from a Exponential distribution
#'
#' @inherit Exponential examples
#'
#' @param d A `Exponential` object created by a call to [Exponential()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Exponential <- function(d, n = 1L, ...) {
  rexp(n = n, rate = d$rate)
}

#' Evaluate the probability mass function of a Exponential distribution
#'
#' @inherit Exponential examples
#' @inheritParams random.Exponential
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Exponential <- function(d, x, ...) {
  dexp(x = x, rate = d$rate)
}

#' @rdname pdf.Exponential
#' @export
#'
log_pdf.Exponential <- function(d, x, ...) {
  dexp(x = x, rate = d$rate, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Exponential distribution
#'
#' @inherit Exponential examples
#' @inheritParams random.Exponential
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Exponential <- function(d, x, ...) {
  pexp(q = x, rate = d$rate)
}

#' Determine quantiles of a Exponential distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Exponential examples
#' @inheritParams random.Exponential
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Exponential <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qexp(p = p, rate = d$rate)
}

#' Fit an Exponential distribution to data
#'
#' @param d An `Exponential` object created by a call to [Exponential()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Exponential distribution
#'
#' @return An `Exponential` object.
#' @export
fit_mle.Exponential <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Exponential(ss$sum / ss$samples)
}


#' Compute the sufficient statistics of an Exponential distribution from data
#'
#' @inheritParams fit_mle.Exponential
#'
#' @return A named list of the sufficient statistics of the exponential
#'   distribution:
#'
#'   - `sum`: The sum of the observations.
#'   - `samples`: The number of observations.
#'
#' @export
suff_stat.Exponential <- function(d, x, ...) {
  valid_x <- (x > 0)
  if(any(!valid_x)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), samples = length(x))
}
