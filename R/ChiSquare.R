#' Create a Chi-Square distribution
#'
#' Chi-square distributions show up often in frequentist settings
#' as the sampling distribution of test statistics, especially
#' in maximum likelihood estimation settings.
#'
#' @param df Degrees of freedom. Must be positive.
#'
#' @return A `ChiSquare` object.
#' @export
#'
#' @family continuous distributions
#'
#' @section Transformations:
#'
#' A squared standard [Normal()] distribution is equivalent to a
#' \eqn{\chi^2_1} distribution with one degree of freedom. The
#' \eqn{\chi^2} distribution is a special case of the [Gamma()]
#' distribution with shape (TODO: check this) parameter equal
#' to a half. Sums of \eqn{\chi^2} distributions
#' are also distributed as \eqn{\chi^2} distributions, where the
#' degrees of freedom of the contributing distributions get summed.
#' The ratio of two \eqn{\chi^2} distributions is a [FisherF()]
#' distribution. The ratio of a [Normal()] and the square root
#' of a scaled [ChiSquare()] is a [StudentsT()] distribution.
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3/>, where the math
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
#' set.seed(27)
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
ChiSquare <- function(df) {
  d <- list(df = df)
  class(d) <- c("ChiSquare", "distribution")
  d
}

#' @export
print.ChiSquare <- function(x, ...) {
  cat(glue("Chi Square distribution (df = {x$df})"), "\n")
}

#' @export
mean.ChiSquare <- function(x, ...) {
  ellipsis::check_dots_used()
  x$df
}


#' @export
variance.ChiSquare <- function(x, ...) x$df * 2

#' @export
skewness.ChiSquare <- function(x, ...) sqrt(8 / x$df)

#' @export
kurtosis.ChiSquare <- function(x, ...) 12 / x$df

#' Draw a random sample from a chi square distribution
#'
#' @inherit ChiSquare examples
#'
#' @param x A `ChiSquare` object created by a call to [ChiSquare()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.ChiSquare <- function(x, n = 1L, ...) {
  rchisq(n = n, df = x$df)
}

#' Evaluate the probability mass function of a chi square distribution
#'
#' @inherit ChiSquare examples
#'
#' @param d A `ChiSquare` object created by a call to [ChiSquare()].
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
#'
#' @param d A `ChiSquare` object created by a call to [ChiSquare()].
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
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.ChiSquare <- function(x, probs, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  ellipsis::check_dots_used()
  qchisq(p = probs, df = x$df)
}


#' Return the support of the ChiSquare distribution
#'
#' @param d An `ChiSquare` object created by a call to [ChiSquare()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.ChiSquare <- function(d){
  c(0, Inf)
}

