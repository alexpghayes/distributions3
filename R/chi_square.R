#' Create a chi square distribution
#'
#' @param df Degrees of freedom. Must be positive.
#'
#' @return A `chi_square` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' c <- chi_square(5)
#' c
#'
#' random(c, 10)
#' pdf(c, 2)
#' cdf(c, 4)
#' quantile(c, 0.7)
#'
#' cdf(c, quantile(c, 0.7))
#' quantile(c, cdf(c, 7))
#'
chi_square <- function(df) {
  d <- list(df = df)
  class(d) <- "chi_square"
  d
}

#' @export
print.chi_square <- function(d) {
  cat(glue("chi square distribution (df = {d$df})"))
}

#' Draw a random sample from a chi square distribution
#'
#' @inherit chi_square examples
#'
#' @param d A `chi_square` object created by a call to [chi_square()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.chi_square <- function(d, n = 1L, ...) {
  rchisq(n = n, df = d$df)
}

#' Evaluate the probability mass function of a chi square distribution
#'
#' @inherit chi_square examples
#' @inheritParams random.chi_square
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.chi_square <- function(d, x, ...) {
  dchisq(x = x, df = d$df)
}

#' Evaluate the cumulative distribution function of a chi square distribution
#'
#' @inherit chi_square examples
#' @inheritParams random.chi_square
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.chi_square <- function(d, x, ...) {
  pchisq(q = x, df = d$df)
}

#' Determine quantiles of a chi square distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit chi_square examples
#' @inheritParams random.chi_square
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.chi_square <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qchisq(p = p, df = d$df)
}
