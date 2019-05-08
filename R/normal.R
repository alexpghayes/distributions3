#' Create a normal distribution
#'
#' @param mu The location parameter, which is also the mean of the
#'   distribution. Can be any real number. Defaults to `0`.
#' @param sigma The scale parameter, which is also the standard
#'   deviation of the distribution. Can be any positive number.
#'   Defaults to `1`.
#'
#' @return A `normal` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' n <- normal(5, 2)
#' n
#'
#' random(n, 10)
#' pdf(n, 2)
#' cdf(n, 4)
#' quantile(n, 0.7)
#'
#' cdf(n, quantile(n, 0.7))
#' quantile(n, cdf(n, 7))
#'
normal <- function(mu = 0, sigma = 1) {
  d <- list(mu = mu, sigma = sigma)
  class(d) <- "normal"
  d
}

#' @export
print.normal <- function(d) {
  cat(glue("normal distribution (mu = {d$mu}, sigma = {d$sigma})"))
}

#' Draw a random sample from a normal distribution
#'
#' @inherit normal examples
#'
#' @param d A `normal` object created by a call to [normal()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.normal <- function(d, n = 1L, ...) {
  rnorm(n = n, mean = d$mu, sd = d$sigma)
}

#' Evaluate the probability mass function of a normal distribution
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.normal <- function(d, x, ...) {
  dnorm(x = x, mean = d$mu, sd = d$sigma)
}

#' Evaluate the cumulative distribution function of a normal distribution
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.normal <- function(d, x, ...) {
  pnorm(q = x, mean = d$mu, sd = d$sigma)
}

#' Determine quantiles of a normal distribution
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @examples
#'
#' # TODO: ggplot example illustrating that you always get the lower
#' # tail probability
#'
quantile.normal <- function(d, p, ...) {
  qnorm(p = p, mean = d$mu, sd = d$sigma)
}
