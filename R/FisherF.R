#' Create an F distribution
#'
#' @param df1 Numerator degrees of freedom. Can be any positive number.
#' @param df2 Denominator degrees of freedom. Can be any positive number.
#' @param lambda Non-centrality parameter. Can be any positive number.
#'   Defaults to `0`.
#'
#' @return A `FisherF` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   TODO
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- FisherF(5, 10, 0.2)
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
FisherF <- function(df1, df2, lambda = 0) {
  d <- list(df1 = df1, df2 = df2, lambda = lambda)
  class(d) <- c("FisherF", "distribution")
  d
}

#' @export
print.FisherF <- function(x, ...) {
  cat(glue("Fisher's F distribution (df1 = {x$df1}, df2 = {x$df2}, lambda = {x$lambda})"), "\n")
}


#' @export
mean.FisherF <- function(x, ...) {
  ellipsis::check_dots_used()
  # The k-th moment of an F(df1, df2) distribution exists and
  # is finite only when 2k < d2
  d1 <- x$df1
  d2 <- x$df2
  if (d2 > 2) d2 / (d2 - 2) else NaN
}

#' @export
variance.FisherF <- function(x, ...) {
  d1 <- x$df1
  d2 <- x$df2
  if (d2 > 4) {
    (2 * d2^2 * (d1 + d2 - 2)) / (d1 * (d2 - 2)^2 * (d2 - 4))
  } else {
    NaN
  }
}

#' @export
skewness.FisherF <- function(x, ...) {
  d1 <- x$df1
  d2 <- x$df2
  if (d2 > 6) {
    a <- (2 * d1 + d2 - 2) * sqrt(8 * (d2 - 4))
    b <- (d2 - 6) * sqrt(d1 * (d1 + d2 - 2))
    a / b
  } else {
    NaN
  }
}

#' @export
kurtosis.FisherF <- function(x, ...) {
  d1 <- x$df1
  d2 <- x$df2
  if (d2 > 8) {
    a <- d1 * (5 * d2 - 22) * (d1 + d2 - 2) + (d2 - 4) * (d2 - 2)^2
    b <- d1 * (d2 - 6) * (d2 - 8) * (d1 + d2 - 2)
    12 * a / b
  } else {
    NaN
  }
}

#' Draw a random sample from an F distribution
#'
#' @inherit FisherF examples
#'
#' @param x A `FisherF` object created by a call to [FisherF()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.FisherF <- function(x, n = 1L, ...) {
  rf(n = n, df1 = x$df1, df2 = x$df2, ncp = x$lambda)
}

#' Evaluate the probability mass function of an F distribution
#'
#' @inherit FisherF examples
#'
#' @param d A `FisherF` object created by a call to [FisherF()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.FisherF <- function(d, x, ...) {
  df(x = x, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}

#' @rdname pdf.FisherF
#' @export
#'
log_pdf.FisherF <- function(d, x, ...) {
  df(x = x, df1 = d$df1, df2 = d$df2, ncp = d$lambda, log = TRUE)
}

#' Evaluate the cumulative distribution function of an F distribution
#'
#' @inherit FisherF examples
#'
#' @param d A `FisherF` object created by a call to [FisherF()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.FisherF <- function(d, x, ...) {
  pf(q = x, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}

#' Determine quantiles of an F distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit FisherF examples
#' @inheritParams random.FisherF
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.FisherF <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qf(p = probs, df1 = x$df1, df2 = x$df2, ncp = x$lambda)
}

#' Return the support of the FisherF distribution
#'
#' @param d An `FisherF` object created by a call to [FisherF()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.FisherF <- function(d){
  return(c(0, Inf))
}
