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
#' @examples
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
#'
FisherF <- function(df1, df2, lambda = 0) {
  d <- list(df1 = df1, df2 = df2, lambda = lambda)
  class(d) <- c("FisherF", "distribution")
  d
}

#' @export
print.FisherF <- function(x, ...) {
  cat(glue("Fisher's F distribution (df1 = {x$df1}, df2 = {x$df2}, lambda = {x$lambda})"))
}

#' Draw a random sample from an F distribution
#'
#' @inherit FisherF examples
#'
#' @param d A `FisherF` object created by a call to [FisherF()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.FisherF <- function(d, n = 1L, ...) {
  rf(n = n, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}

#' Evaluate the probability mass function of an F distribution
#'
#' @inherit FisherF examples
#' @inheritParams random.FisherF
#'
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
#' @inheritParams random.FisherF
#'
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
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.FisherF <- function(d, p, ...) {
  qf(p = p, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}
