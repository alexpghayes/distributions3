#' Create a Beta distribution
#'
#'
#' @param alpha The alpha parameter. `alpha` can be any value strictly
#'   greater than zero. Defaults to `1`.
#' @param beta The beta parameter. `beta` can be any value strictly
#'   greater than zero. Defaults to `1`.
#'
#' @return A `beta` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Beta(1, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 0.7)
#' log_pdf(X, 0.7)
#'
#' cdf(X, 0.7)
#' quantile(X, 0.7)
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
Beta <- function(alpha = 1, beta = 1) {
  d <- list(alpha = alpha, beta = beta)
  class(d) <- c("Beta", "distribution")
  d
}

#' @export
print.Beta <- function(x, ...) {
  cat(glue("Beta distribution (alpha = {x$alpha}, beta = {x$beta})"), "\n")
}

#' @export
mean.Beta <- function(x, ...) {
  ellipsis::check_dots_used()
  x$alpha / (x$alpha + x$beta)
}

#' @export
variance.Beta <- function(x, ...) {
  a <- x$alpha
  b <- x$beta
  (a * b) /  ((a + b)^2 * (a + b + 1))
}

#' @export
skewness.Beta <- function(x, ...) {
  a <- x$alpha
  b <- x$beta
  2 * (b - a) * sqrt(a + b + 1) / (a + b + 2) * sqrt(a * b)
}

#' @export
kurtosis.Beta <- function(x, ...) {
  a <- x$alpha
  b <- x$beta
  num <- 6 * ((a - b)^2 * (a + b + 1) - (a * b) * (a + b + 2))
  denom <- a * b * (a + b + 2) * (a + b + 3)
  num / denom
}

#' Draw a random sample from a Beta distribution
#'
#' @inherit Beta examples
#'
#' @param x A `Beta` object created by a call to [Beta()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector containing values in `[0, 1]` of length `n`.
#' @export
#'
random.Beta <- function(x, n = 1L, ...) {
  rbeta(n = n, shape1 = x$alpha, shape2 = x$beta)
}

#' Evaluate the probability mass function of a Beta distribution
#'
#' @inherit Beta examples
#'
#' @param d A `Beta` object created by a call to [Beta()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Beta <- function(d, x, ...) {
  dbeta(x = x, shape1 = d$alpha, shape2 = d$beta)
}

#' @rdname pdf.Beta
#' @export
#'
log_pdf.Beta <- function(d, x, ...) {
  dbeta(x = x, shape1 = d$alpha, shape2 = d$beta, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Beta distribution
#'
#' @inherit Beta examples
#'
#' @param d A `Beta` object created by a call to [Beta()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Beta <- function(d, x, ...) {
  pbeta(q = x, shape1 = d$alpha, shape2 = d$beta)
}

#' Determine quantiles of a Beta distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Beta examples
#' @inheritParams random.Beta
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.Beta <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qbeta(p = probs, shape1 = x$alpha, shape2 = x$beta)
}


#' Return the support of the Beta distribution
#'
#' @param d An `Beta` object created by a call to [Beta()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Beta <- function(d){
  return(c(0, 1))
}
