#' Create a cauchy distribution
#'
#' Note that the cauchy distribution is the student's t distribution
#' with one degree of freedom. The cauchy distribution does not have
#' a well defined mean or variance.
#'
#' @param location The location parameter. Can be any real number. Defaults
#'   to `0`.
#' @param scale The scale parameter. Must be greater than zero (?). Defaults
#'   to `1`.
#'
#' @return A `cauchy` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' c <- cauchy(10, 0.2)
#' c
#'
#' random(c, 10)
#' pdf(c, 2)
#' cdf(c, 2)
#' quantile(c, 0.7)
#'
#' cdf(c, quantile(c, 0.7))
#' quantile(c, cdf(c, 7))
#'
cauchy <- function(location = 0, scale = 1) {
  d <- list(location = location, scale = scale)
  class(d) <- "cauchy"
  d
}

#' @export
print.cauchy <- function(d, ...) {
  cat(glue("Cauchy distribution (location = {d$location}, scale = {d$scale})"))
}

#' Draw a random sample from a cauchy distribution
#'
#' @inherit cauchy examples
#'
#' @param d A `cauchy` object created by a call to [cauchy()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.cauchy <- function(d, n = 1L, ...) {
  rcauchy(n = n, location = d$location, scale = d$scale)
}

#' Evaluate the probability mass function of a cauchy distribution
#'
#' @inherit cauchy examples
#' @inheritParams random.cauchy
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.cauchy <- function(d, x, ...) {
  dcauchy(x = x, location = d$location, scale = d$scale)
}

#' Evaluate the cumulative distribution function of a cauchy distribution
#'
#' @inherit cauchy examples
#' @inheritParams random.cauchy
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.cauchy <- function(d, x, ...) {
  pcauchy(q = x, location = d$location, scale = d$scale)
}

#' Determine quantiles of a cauchy distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit cauchy examples
#' @inheritParams random.cauchy
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.cauchy <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qcauchy(p = p, location = d$location, scale = d$scale)
}
