#' Create an F distribution
#'
#' @param df1 Numerator degrees of freedom. Can be any positive number.
#' @param df2 Denominator degrees of freedom. Can be any positive number.
#' @param lambda Non-centrality parameter. Can be any positive number.
#'   Defaults to `0`.
#'
#' @return An `f` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' f <- fisher_f(5, 10, 0.2)
#' f
#'
#' random(f, 10)
#' pdf(f, 2)
#' log_pdf(f, 2)
#' cdf(f, 4)
#' quantile(f, 0.7)
#'
#' cdf(f, quantile(f, 0.7))
#' quantile(f, cdf(f, 7))
#'
fisher_f <- function(df1, df2, lambda = 0) {
  d <- list(df1 = df1, df2 = df2, lambda = lambda)
  class(d) <- "fisher_f"
  d
}

#' @export
print.fisher_f <- function(x, ...) {
  cat(glue("Fisher's F distribution (rate = {x$rate})"))
}

#' Draw a random sample from an F distribution
#'
#' @inherit fisher_f examples
#'
#' @param d A `fisher_f` object created by a call to [fisher_f()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.fisher_f <- function(d, n = 1L, ...) {
  rf(n = n, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}

#' Evaluate the probability mass function of an F distribution
#'
#' @inherit fisher_f examples
#' @inheritParams random.fisher_f
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.fisher_f <- function(d, x, ...) {
  df(x = x, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}

#' @rdname pdf.fisher_f
#' @export
#'
log_pdf.fisher_f <- function(d, x, ...) {
  df(x = x, df1 = d$df1, df2 = d$df2, ncp = d$lambda, log = TRUE)
}

#' Evaluate the cumulative distribution function of an F distribution
#'
#' @inherit fisher_f examples
#' @inheritParams random.fisher_f
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.fisher_f <- function(d, x, ...) {
  pf(q = x, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}

#' Determine quantiles of an F distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit fisher_f examples
#' @inheritParams random.fisher_f
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.fisher_f <- function(d, p, ...) {
  qf(p = p, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
}
