#' Create a gamma distribution
#'
#' Several important distributions are special cases of the gamma
#' distribution. When the shape parameter is `1`, the gamma is an
#' gamma. When the (scale?) parametre is `2` the gamma is a
#' chi square.
#'
#' @param shape The shape parameter. Can be any positive number.
#' @param rate The rate parameter. Can be any positive number. Defaults
#'   to `1`.
#'
#' @return A `gamma` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' g <- gamma(5, 2)
#' g
#'
#' random(g, 10)
#' pdf(g, 2)
#' log_pdf(g, 2)
#' cdf(g, 4)
#' quantile(g, 0.7)
#'
#' cdf(g, quantile(g, 0.7))
#' quantile(g, cdf(g, 7))
#'
gamma <- function(shape, rate = 1) {
  d <- list(shape = shape, rate = rate)
  class(d) <- "gamma"
  d
}

#' @export
print.gamma <- function(x, ...) {
  cat(glue("Gamma distribution (shape = {x$shape}, rate = {x$rate})"))
}

#' Draw a random sample from a gamma distribution
#'
#' @inherit gamma examples
#'
#' @param d A `gamma` object created by a call to [gamma()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.gamma <- function(d, n = 1L, ...) {
  rgamma(n = n, shape = d$shape, rate = d$rate)
}

#' Evaluate the probability mass function of a gamma distribution
#'
#' @inherit gamma examples
#' @inheritParams random.gamma
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.gamma <- function(d, x, ...) {
  dgamma(x = x, shape = d$shape, rate = d$rate)
}

#' @rdname pdf.gamma
#' @export
#'
log_pdf.gamma <- function(d, x, ...) {
  dgamma(x = x, shape = d$shape, rate = d$rate, log = TRUE)
}

#' Evaluate the cumulative distribution function of a gamma distribution
#'
#' @inherit gamma examples
#' @inheritParams random.gamma
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.gamma <- function(d, x, ...) {
  pgamma(q = x, shape = d$shape, rate = d$rate)
}

#' Determine quantiles of a gamma distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit gamma examples
#' @inheritParams random.gamma
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.gamma <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qgamma(p = p, shape = d$shape, rate = d$rate)
}

#' Compute the likelihood of a distribution given data
#'
#' @export
likelihood.gamma <- function(d, x, ...) {

}


#' Fit a gamma distribution to data
#'
#' @param d A `gamma` object created by a call to [gamma()].
#' @param x A vector to fit the gamma distribution to.
#'
#' @return a `gamma` object
#' @export
fit_mle.gamma <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  gamma()
}

#' Compute the sufficient statistics for a bernoulli distribution from data
#'
#' @inherit gamma
#' @export
suff_stat.gamma <- function(d, x, ...) {
  if(any(x < 0)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), log_sum = sum(log(x)), samples = length(x))
}
