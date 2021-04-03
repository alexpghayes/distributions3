#' Create a Poisson distribution
#'
#' Poisson distributions are frequently used to model counts.
#'
#' @param lambda The shape parameter, which is also the mean and the
#'   variance of the distribution. Can be any positive number.
#'
#' @return A `Poisson` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Poisson random variable with parameter
#'   `lambda` = \eqn{\lambda}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}{{0, 1, 2, 3, ...}}
#'
#'   **Mean**: \eqn{\lambda}
#'
#'   **Variance**: \eqn{\lambda}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}
#'   }{
#'     P(X = k) = \lambda^k e^(-\lambda) / k!
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) = e^{-\lambda}
#'     \sum_{i = 0}^{\lfloor k \rfloor} \frac{\lambda^i}{i!}
#'   }{
#'     P(X \le k) = e^(-\lambda)
#'     \sum_{i = 0}^k \lambda^i / i!
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\lambda (e^t - 1)}
#'   }{
#'     E(e^(tX)) = e^(\lambda (e^t - 1))
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Poisson(2)
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
#'
Poisson <- function(lambda) {
  d <- list(lambda = lambda)
  class(d) <- c("Poisson", "distribution")
  d
}

#' @export
print.Poisson <- function(x, ...) {
  cat(glue("Poisson distribution (lambda = {x$lambda})"), "\n")
}

#' @export
mean.Poisson <- function(d, ...) d$lambda

#' @export
variance.Poisson <- function(d, ...) d$lambda

#' @export
skewness.Poisson <- function(d, ...) 1 / sqrt(d$lambda)

#' @export
kurtosis.Poisson <- function(d, ...) 1 / d$lambda

#' Draw a random sample from a Poisson distribution
#'
#' @inherit Poisson examples
#'
#' @param d A `Poisson` object created by a call to [Poisson()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Poisson <- function(d, n = 1L, ...) {
  rpois(n = n, lambda = d$lambda)
}

#' Evaluate the probability mass function of a Poisson distribution
#'
#' @inherit Poisson examples
#' @inheritParams random.Poisson
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Poisson <- function(d, x, ...) {
  dpois(x = x, lambda = d$lambda)
}

#' @rdname pdf.Poisson
#' @export
#'
log_pdf.Poisson <- function(d, x, ...) {
  dpois(x = x, lambda = d$lambda, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Poisson distribution
#'
#' @inherit Poisson examples
#' @inheritParams random.Poisson
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Poisson <- function(d, x, ...) {
  ppois(q = x, lambda = d$lambda)
}

#' Determine quantiles of a Poisson distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Poisson examples
#' @inheritParams random.Poisson
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Poisson <- function(d, p, ...) {
  qpois(p = p, lambda = d$lambda)
}

#' Fit an Poisson distribution to data
#'
#' @param d An `Poisson` object created by a call to [Poisson()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Poisson distribution
#'
#' @return An `Poisson` object.
#' @export
fit_mle.Poisson <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Poisson(ss$sum / ss$samples)
}


#' Compute the sufficient statistics of an Poisson distribution from data
#'
#' @inheritParams fit_mle.Poisson
#'
#' @return A named list of the sufficient statistics of the Poisson
#'   distribution:
#'
#'   - `sum`: The sum of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
suff_stat.Poisson <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x %% 1 == 0)
  if (any(!valid_x)) stop("`x` must only contain positive integers")
  list(sum = sum(x), samples = length(x))
}


#' Return the support of the Poisson distribution
#'
#' @param d An `Poisson` object created by a call to [Poisson()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Poisson <- function(d){
  if(!is_distribution(d)){
    message("d has to be a disitrubtion")
    stop()
  }
  return(c(0, Inf))
}
