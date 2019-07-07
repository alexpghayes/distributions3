#' Create a continuous uniform distribution
#'
#' TODO: A categorical distribution that can be used for discrete uniforms.
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
  class(d) <- "Uniform"
  d
}

#' @export
print.Uniform <- function(x, ...) {
  if(x$a >  x$b) names(x) <- c("b", "a")
  cat(glue("Continuous Uniform distribution (a = {x$a}, b = {x$b})"))
}

#' @export
mean.Uniform <- function(d, ...) (d$a + d$b) / 2

#' @export
variance.Uniform <- function(d, ...) (1 / 12) * (d$b - d$a) ^ 2

#' @export
skewness.Uniform <- function(d, ...) 0

#' @export
kurtosis.Uniform <- function(d, ...) -6/5

#' Draw a random sample from a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector containing values in `[a, b]` of length `n`.
#' @export
#'
random.Uniform <- function(d, n = 1L, ...) {
  runif(n = n, min = min(d$a, d$b), max = max(d$a, d$b))
}

#' Evaluate the probability mass function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#' @inheritParams random.Uniform
#'
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
#' @inheritParams random.Uniform
#'
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
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Uniform <- function(d, p, ...) {
  qunif(p = p, min = min(d$a, d$b), max = max(d$a, d$b))
}
