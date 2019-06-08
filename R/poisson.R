#' Create a poisson distribution
#'
#' Poisson distributions are frequently used to model counts.
#'
#' @param lambda The shape parameter, which is also the mean and the
#'   variance of the distribution. Can be any positive number.
#'
#' @return A `poisson` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Poisson random variable with parameter
#'   `lamdba` = \eqn{\lambda}.
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
#'     \mathbb{E}(e^{tX}) = e^{\lambda (e^t - 1)}
#'   }{
#'     E(e^(tX)) = e^(\lambda (e^t - 1))
#'   }
#'
#' @examples
#'
#' p <- poisson(2)
#' p
#'
#' random(p, 10)
#' pdf(p, 2)
#' log_pdf(p, 2)
#' cdf(p, 4)
#' quantile(p, 0.7)
#'
#' cdf(p, quantile(p, 0.7))
#' quantile(p, cdf(p, 7))
#'
poisson <- function(lambda) {
  d <- list(lambda = lambda)
  class(d) <- c("poisson", "distribution")
  d
}

#' @export
print.poisson <- function(x, ...) {
  cat(glue("Poisson distribution (lambda = {x$lambda})"))
}

#' Draw a random sample from a poisson distribution
#'
#' @inherit poisson examples
#'
#' @param d A `poisson` object created by a call to [poisson()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.poisson <- function(d, n = 1L, ...) {
  rpois(n = n, lambda = d$lambda)
}

#' Evaluate the probability mass function of a poisson distribution
#'
#' @inherit poisson examples
#' @inheritParams random.poisson
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.poisson <- function(d, x, ...) {
  dpois(x = x, lambda = d$lambda)
}

#' @rdname pdf.poisson
#' @export
#'
log_pdf.poisson <- function(d, x, ...) {
  dpois(x = x, lambda = d$lambda, log = TRUE)
}

#' Evaluate the cumulative distribution function of a poisson distribution
#'
#' @inherit poisson examples
#' @inheritParams random.poisson
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.poisson <- function(d, x, ...) {
  ppois(q = x, lambda = d$lambda)
}

#' Determine quantiles of a poisson distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit poisson examples
#' @inheritParams random.poisson
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.poisson <- function(d, p, ...) {
  qpois(p = p, lambda = d$lambda)
}

#' Fit an poisson distribution to data
#'
#' @param d An `poisson` object created by a call to [poisson()].
#' @param x A vector of data.
#'
#' @family poisson distribution
#'
#' @return An `poisson` object.
#' @export
fit_mle.poisson <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  poisson(ss$sum / ss$samples)
}


#' Compute the sufficient statistics of an poisson distribution from data
#'
#' @inheritParams fit_mle.poisson
#'
#' @return A named list of the sufficient statistics of the poisson distribution
#'   \describe{
#'     \item{\code{sum}}{The sum of the data}
#'     \item{\code{samples}}{The number of samples in the data}
#'   }
#'
#' @export
suff_stat.poisson <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x %% 1 == 0)
  if(any(!valid_x)) stop("`x` must only contain positive integers")
  list(sum = sum(x), samples = length(x))
}

