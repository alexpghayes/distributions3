#' Create a chi square distribution
#'
#' Chi square distributions show up often in frequentist settings
#' as the sampling distribution of test statistics, especially
#' in maximum likelihood estimation settings. The chi square
#' distribution with a single degree of freedom is equivalent
#' to a squared standard normal distribution. Summing chi square
#' distributions results in a new chi square distribution with
#' added degrees of freedom. A chi square distribution is a special
#' case of the gamma distribution with `shape = 2` (double check this).
#' TODO: add relationship to F and T distributions.
#'
#' @param df Degrees of freedom. Must be positive.
#'
#' @return A `ChiSquare` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a \eqn{\chi^2} random variable with
#'   `df` = \eqn{k}.
#'
#'   **Support**: \eqn{R^+}, the set of positive real numbers
#'
#'   **Mean**: \eqn{k}
#'
#'   **Variance**: \eqn{2k}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2}
#'   }{
#'     f(x) = 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2))
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   The cumulative distribution function has the form
#'
#'   \deqn{
#'     F(t) = \int_{-\infty}^t \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2} dx
#'   }{
#'     F(t) = integral_{-\infty}^t 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2)) dx
#'   }
#'
#'   but this integral does not have a closed form solution and must be
#'   approximated numerically. The c.d.f. of a standard normal is sometimes
#'   called the "error function". The notation \eqn{\Phi(t)} also stands
#'   for the c.d.f. of a standard normal evaluated at \eqn{t}. Z-tables
#'   list the value of \eqn{\Phi(t)} for various \eqn{t}.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2}
#'   }{
#'     E(e^(tX)) = e^(\mu t + \sigma^2 t^2 / 2)
#'   }
#'
#' @examples
#'
#' X <- ChiSquare(5)
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
ChiSquare <- function(df) {
  d <- list(df = df)
  class(d) <- c("ChiSquare", "distribution")
  d
}

#' @export
print.ChiSquare <- function(x, ...) {
  cat(glue("Chi Square distribution (df = {x$df})"))
}

#' @export
mean.ChiSquare <- function(d, ...) d$df

#' @export
variance.ChiSquare <- function(d, ...) d$df * 2

#' @export
skewness.ChiSquare <- function(d, ...) sqrt(8 / d$df)

#' @export
kurtosis.ChiSquare <- function(d, ...) 12 / d$df

#' Draw a random sample from a chi square distribution
#'
#' @inherit ChiSquare examples
#'
#' @param d A `ChiSquare` object created by a call to [ChiSquare()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.ChiSquare <- function(d, n = 1L, ...) {
  rchisq(n = n, df = d$df)
}

#' Evaluate the probability mass function of a chi square distribution
#'
#' @inherit ChiSquare examples
#' @inheritParams random.ChiSquare
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.ChiSquare <- function(d, x, ...) {
  dchisq(x = x, df = d$df)
}

#' @rdname pdf.ChiSquare
#' @export
#'
log_pdf.ChiSquare <- function(d, x, ...) {
  dchisq(x = x, df = d$df, log = TRUE)
}

#' Evaluate the cumulative distribution function of a chi square distribution
#'
#' @inherit ChiSquare examples
#' @inheritParams random.ChiSquare
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.ChiSquare <- function(d, x, ...) {
  pchisq(q = x, df = d$df)
}

#' Determine quantiles of a chi square distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit ChiSquare examples
#' @inheritParams random.ChiSquare
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.ChiSquare <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qchisq(p = p, df = d$df)
}
