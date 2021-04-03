#' Create an Exponential distribution
#'
#' Exponential distributions are frequently used for modeling the amount
#' of time that passes until a specific event occurs. For example, exponential
#' distributions could be used to model the time between two earthquakes,
#' the amount of delay between internet packets, or the amount of time a piece
#' of machinery can run before needing repair.
#'
#'
#' @param rate The rate parameter, written \eqn{\lambda} in textbooks.
#'   Can be any positive number. Defaults to `1`.
#'
#' @return An `Exponential` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be an Exponential random variable with
#'   rate parameter `rate` = \eqn{\lambda}.
#'
#'   **Support**: {x in [0, \eqn{\infty})}
#'
#'   **Mean**: 1 / \eqn{\lambda}
#'
#'   **Variance**: 1 / \eqn{\lambda^2}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \lambdae^{-\lambdax}
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(x) = 1 - e^{-\lambdax}
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     \frac{\lambda}{\lambda - t}, for t < \lambda
#'   }
#'
#' @examples
#'
#' set.seed(27)
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
Exponential <- function(rate = 1) {
  d <- list(rate = rate)
  class(d) <- c("Exponential", "distribution")
  d
}

#' @export
print.Exponential <- function(x, ...) {
  cat(glue("Exponential distribution (rate = {x$rate})"), "\n")
}

#' @export
mean.Exponential <- function(d, ...) d$rate ^ -1

#' @export
variance.Exponential <- function(d, ...) d$rate ^ 2

#' @export
skewness.Exponential <- function(d, ...) 2

#' @export
kurtosis.Exponential <- function(d, ...) 6

#' Draw a random sample from an Exponential distribution
#'
#' @inherit Exponential examples
#'
#' @param d An `Exponential` object created by a call to [Exponential()].
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

#' Evaluate the probability density function of an Exponential distribution
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

#' Evaluate the cumulative distribution function of an Exponential distribution
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

#' Determine quantiles of an Exponential distribution
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
  if (any(!valid_x)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), samples = length(x))
}

#' Return the support of the Exponential distribution
#'
#' @param d An `Exponential` object created by a call to [Exponential()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Exponential <- function(d){
  if(!is_distribution(d)){
    message("d has to be a disitrubtion")
    stop()
  }
  return(c(0, Inf))
}
