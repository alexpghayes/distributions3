#' Create a Continuous Uniform distribution
#'
#' A distribution with constant density on an interval. The
#' continuous analogue to the [Categorical()] distribution.
#'
#' @param a The a parameter. `a` can be any value in the set of real
#'   numbers. Defaults to `0`.
#' @param b The a parameter. `b` can be any value in the set of real
#'   numbers. It should be strictly bigger than `a`, but if is not, the
#'   order of the parameters is inverted. Defaults to `1`.
#'
#' @return A `Uniform` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Uniform(1, 2)
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
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
#'
Uniform <- function(a = 0, b = 1) {
  d <- list(a = a, b = b)
  class(d) <- c("Uniform", "distribution")
  d
}

#' @export
print.Uniform <- function(x, ...) {
  if (x$a > x$b) names(x) <- c("b", "a")
  cat(glue("Continuous Uniform distribution (a = {x$a}, b = {x$b})"), "\n")
}

#' @export
mean.Uniform <- function(x, ...) {
  ellipsis::check_dots_used()
  (x$a + x$b) / 2
}

#' @export
variance.Uniform <- function(x, ...) (1 / 12) * (x$b - x$a) ^ 2

#' @export
skewness.Uniform <- function(x, ...) 0

#' @export
kurtosis.Uniform <- function(x, ...) -6/5

#' Draw a random sample from a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param x A `Uniform` object created by a call to [Uniform()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector containing values in `[a, b]` of length `n`.
#' @export
#'
random.Uniform <- function(x, n = 1L, ...) {
  runif(n = n, min = min(x$a, x$b), max = max(x$a, x$b))
}

#' Evaluate the probability mass function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Uniform <- function(d, x, ...) {
  dunif(x = x, min = min(d$a, d$b), max = max(d$a, d$b))
}

#' @rdname pdf.Uniform
#' @export
#'
log_pdf.Uniform <- function(d, x, ...) {
  dunif(x = x, min = min(d$a, d$b), max = max(d$a, d$b), log = TRUE)
}

#' Evaluate the cumulative distribution function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Uniform <- function(d, x, ...) {
  punif(q = x, min = min(d$a, d$b), max = max(d$a, d$b))
}

#' Determine quantiles of a continuous Uniform  distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Uniform examples
#' @inheritParams random.Uniform
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.Uniform <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qunif(p = probs, min = min(x$a, x$b), max = max(x$a, x$b))
}


#' Return the support of the Uniform distribution
#'
#' @param d An `Uniform` object created by a call to [Uniform()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Uniform <- function(d) c(d$a, d$b)
