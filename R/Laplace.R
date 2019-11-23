#' Create a Laplace distribution
#'
#' The Laplace distribution is related to the exponential distribution (see
#' `\link{Exponential}`).  If \eqn{X} has a Laplace distribution with location
#' parameter \eqn{\theta} and scale parameter \eqn{\phi} then
#' \eqn{| X - \theta |} has an exponential distribution with rate
#' parameter \eqn{1 / \phi}.  If \eqn{\theta} is zero then a Laplace
#' distribution may be viewed as a composite of two identical exponential
#' distributions placed back-to-back.  Therefore, it may be referred to as a
#' *double exponential*, *two-tailed exponential* or *bilateral exponential*
#' distribution.
#'
#' @param theta The location parameter \eqn{\theta}.
#'   `theta` can be any real number.  Defaults to `0`.
#' @param phi The scale parameter \eqn{\phi}.
#'  `phi` can be any value strictly greater than zero.  Defaults to `1`.
#'
#' @return A `Laplace` object.
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
#'   In the following, let \eqn{X} be a Laplace random variable with location
#'   parameter `theta` = \eqn{\theta} and scale parameter `phi` = \eqn{\phi}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean, median, mode**: \eqn{\theta}
#'
#'   **Variance**: \eqn{2 \phi^2}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{2} \phi^{-1} \exp(- |x - \theta)| / \phi)
#'   }{
#'     f(x) = (1 / 2 \phi) exp(- |x - \theta)| / \phi)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(x) = \frac{1}{2} + \frac{1}{2} {\rm ~sgn}(x - \theta)%
#'      [1 - \exp(-|x - \theta| / \phi)]
#'   }{
#'     F(x) = 1/2 + (1/2) sgn(x - \theta) [1 - exp(-|x - \theta| / \phi)]
#'   }
#'   where sgn() is the sign function: sgn(\eqn{x}) is equal to \eqn{-1} if
#'   \eqn{x} is negative, \eqn{0} if \eqn{x} is zero and \eqn{1} if \eqn{x} is
#'   positive.
#'
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \exp(\theta t) / (1 - \phi^2 t^2)
#'   }{
#'     E(e^(tX)) = exp(\theta t) / (1 - \phi^2 t^2)
#'   }
#'   for \eqn{|t| < 1 / \phi}.
#' @examples
#'
#' set.seed(27)
#'
#' X <- Laplace(1, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, c(0.3, 1.7))
#' log_pdf(X, c(0.3, 1.7))
#'
#' cdf(X, c(0.3, 1, 1.7))
#' quantile(X, c(0.25, 0.5, 0.75))
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
#'
#' x <- random(X, 100)
#' fit_mle(X, x)
Laplace <- function(theta = 0, phi = 1) {
  d <- list(theta = theta, phi = phi)
  class(d) <- c("Laplace", "distribution")
  d
}

#' @export
print.Laplace <- function(x, ...) {
  cat(glue("Laplace distribution (theta = {x$theta}, phi = {x$phi})\n"))
}

#' Draw a random sample from a Laplace distribution
#'
#' @inherit Laplace examples
#'
#' @param d A `Laplace` object created by a call to [Laplace()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Laplace <- function(d, n = 1L, ...) {
  # Mimic the vectorised behaviour of R's base random generation functions
  max_len <- ifelse(length(n) > 1L, length(n), n)
  d$theta <- rep_len(d$theta, max_len)
  d$phi <- rep_len(d$phi, max_len)
  quantile(d = d, p = runif(n))
}

#' Evaluate the probability mass function of a Laplace distribution
#'
#' @inherit Laplace examples
#' @inheritParams random.Laplace
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Laplace <- function(d, x, ...) {
  # Mimic the vectorised behaviour of R's base random generation functions
  if (length(x) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(x), length(d$theta), length(d$phi))
  x <- rep_len(x, max_len)
  d$theta <- rep_len(d$theta, max_len)
  d$phi <- rep_len(d$phi, max_len)
  dexp(x = abs(x - d$theta), rate = 1 / d$phi) / 2
}

#' @rdname pdf.Laplace
#' @export
#'
log_pdf.Laplace <- function(d, x, ...) {
  if (length(x) == 0) {
    return(numeric(0))
  }
  # Mimic the vectorised behaviour of R's base random generation functions
  max_len <- max(length(x), length(d$theta), length(d$phi))
  x <- rep_len(x, max_len)
  d$theta <- rep_len(d$theta, max_len)
  d$phi <- rep_len(d$phi, max_len)
  dexp(x = abs(x - d$theta), rate = 1 / d$phi, log = TRUE) - log(2)
}

#' Evaluate the cumulative distribution function of a Laplace distribution
#'
#' @inherit Laplace examples
#' @inheritParams random.Laplace
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Laplace <- function(d, x, ...) {
  if (length(x) == 0) {
    return(numeric(0))
  }
  # Mimic the vectorised behaviour of R's base random generation functions
  max_len <- max(length(x), length(d$theta), length(d$phi))
  x <- rep_len(x, max_len)
  d$theta <- rep_len(d$theta, max_len)
  d$phi <- rep_len(d$phi, max_len)
  0.5 + 0.5 * sign(x - d$theta) * pexp(q = abs(x - d$theta), rate = 1 / d$phi)
}

#' Determine quantiles of a Laplace distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Laplace examples
#' @inheritParams random.Laplace
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Laplace <- function(d, p, ...) {
  if (length(p) == 0) {
    return(numeric(0))
  }
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("invalid p: p must be in [0, 1].")
  }
  # Mimic the vectorised behaviour of R's base random generation functions
  max_len <- max(length(p), length(d$theta), length(d$phi))
  p <- rep_len(p, max_len)
  d$theta <- rep_len(d$theta, max_len)
  d$phi <- rep_len(d$phi, max_len)
  d$theta + sign(p - 0.5) * qexp(2 * abs(p - 0.5), rate = 1 / d$phi)
}

#' Fit a Laplace distribution to data
#'
#' @param d A `Laplace` object created by a call to [Laplace()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Laplace distribution
#'
#' @return A `Laplace` object.
#' @export
fit_mle.Laplace <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Laplace(ss$theta, ss$phi)
}

#' Compute the sufficient statistics for a Laplace distribution from data
#'
#' @inheritParams fit_mle.Laplace
#'
#' @return A named list of the sufficient statistics of the Laplace
#'   distribution:
#'
#'   - `theta`: The sample median.
#'   - `sigma`: The sample mean absolute deviation from the sample median.
#'   - `samples`: The number of samples in the data.
#'
#' @export
suff_stat.Laplace <- function(d, x, ...) {
  valid_x <- is.numeric(x)
  if (!valid_x) stop("`x` must be a numeric vector")
  thetahat <- median(x)
  phihat <- mean(abs(x - median(x)))
  list(theta = thetahat, phi = phihat, samples = length(x))
}
