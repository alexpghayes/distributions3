#' Create a continuous uniform distribution
#'
#'
#' @param a The a parameter. `a` can be any value in the set of real
#'   numbers. Defaults to `0`.
#' @param b The a parameter. `b` can be any value in the set of real
#'   numbers. It should be strictly bigger than `a`, but if is not, the
#'   order of the parameters is inverted. Defaults to `1`.
#'
#' @return A `uniform` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' b <- uniform(1, 2)
#' b
#'
#' random(b, 10)
#' pdf(b, 0.7)
#' log_pdf(b, 0.7)
#' cdf(b, 0.7)
#' quantile(b, 0.7)
#'
#' cdf(b, quantile(b, 0.7))
#' quantile(b, cdf(b, 0.7))
#'
uniform <- function(a = 0, b = 1) {
  d <- list(a = a, b = b)
  class(d) <- "uniform"
  d
}

#' @export
print.uniform <- function(x, ...) {
  if(x$a >  x$b) names(x) <- c("b", "a")
  cat(glue("Continuous uniform distribution (a = {x$a}, b = {x$b})"))
}

#' Draw a random sample from a continuous uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `uniform` object created by a call to [uniform()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector containing values in `[a, b]` of length `n`.
#' @export
#'
random.uniform <- function(d, n = 1L, ...) {
  runif(n = n, min = min(d$a, d$b), max = max(d$a, d$b))
}

#' Evaluate the probability mass function of a continuous uniform distribution
#'
#' @inherit uniform examples
#' @inheritParams random.uniform
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.uniform <- function(d, x, ...) {
  dunif(x = x, min = min(d$a, d$b), max = max(d$a, d$b))
}

#' @rdname pdf.uniform
#' @export
#'
log_pdf.uniform <- function(d, x, ...) {
  dunif(x = x, min = min(d$a, d$b), max = max(d$a, d$b), log = TRUE)
}

#' Evaluate the cumulative distribution function of a continuous uniform distribution
#'
#' @inherit uniform examples
#' @inheritParams random.uniform
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.uniform <- function(d, x, ...) {
  punif(q = x, min = min(d$a, d$b), max = max(d$a, d$b))
}

#' Determine quantiles of a continuous uniform  distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit uniform examples
#' @inheritParams random.uniform
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.uniform <- function(d, p, ...) {
  qunif(p = p, min = min(d$a, d$b), max = max(d$a, d$b))
}
