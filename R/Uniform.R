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
Uniform <- function(a = 0, b = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(a) == length(b) | length(a) == 1 | length(b) == 1
  )
  d <- data.frame(a = a, b = b)
  class(d) <- c("Uniform", "distribution")
  d
}

#' @export
mean.Uniform <- function(x, ...) {
  ellipsis::check_dots_used()
  (x$a + x$b) / 2
}

#' @export
variance.Uniform <- function(x, ...) (1 / 12) * (x$b - x$a)^2

#' @export
skewness.Uniform <- function(x, ...) rep.int(0, length(x))

#' @export
kurtosis.Uniform <- function(x, ...) rep(-6 / 5, length(x))

#' Draw a random sample from a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param x A `Uniform` object created by a call to [Uniform()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Values in `[a, b]`. In case of a single distribution object, a numeric
#' vector of length `n` (if `drop = TRUE`, default) or a `data.frame`
#' with `n` columns. In case of a vectorized distribution
#' object, either a matrix (if `drop = TRUE`, default) or a `data.frame`
#' with `n` columns.
#' @export
#'
random.Uniform <- function(x, n = 1L, drop = TRUE, ...) {
  FUN <- function(at, d) {
    runif(
      n = length(d),
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max)
    )
  }
  apply_dpqr(d = x, FUN = FUN, at = rep.int(1, n), type_prefix = "r", drop = drop)
}

#' Evaluate the probability mass function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{dunif}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, a numeric
#'   vector of probabilities of length `x` (if `drop = TRUE`, default)
#'   or a `data.frame` with `n` columns. In case of a vectorized distribution
#'   object, either a matrix (if `drop = TRUE`, default) or a `data.frame`
#'   with `n` columns, containing all possible combinations.
#' @export
#'
pdf.Uniform <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) {
    dunif(
      x = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      ...
    )
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type_prefix = "d", drop = drop)
}

#' @rdname pdf.Uniform
#' @export
#'
log_pdf.Uniform <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) {
    dunif(
      x = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      log = TRUE
    )
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type_prefix = "l", drop = drop)
}

#' Evaluate the cumulative distribution function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{punif}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, a numeric
#'   vector of cumulative probabilities of length `x` (if `drop = TRUE`, default)
#'   or a `data.frame` with `n` columns. In case of a vectorized distribution
#'   object, either a matrix (if `drop = TRUE`, default) or a `data.frame`
#'   with `n` columns, containing all possible combinations.
#' @export
#'
cdf.Uniform <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) {
    punif(
      q = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      ...
    )
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type_prefix = "p", drop = drop)
}

#' Determine quantiles of a continuous Uniform  distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Uniform examples
#' @inheritParams random.Uniform
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{qunif}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, a numeric
#'   vector of quantiles of length `probs` (if `drop = TRUE`, default)
#'   or a `data.frame` with `n` columns. In case of a vectorized distribution
#'   object, either a matrix (if `drop = TRUE`, default) or a `data.frame`
#'   with `n` columns, containing all possible combinations.
#' @export
#'
quantile.Uniform <- function(x, probs, drop = TRUE, ...) {
  ellipsis::check_dots_used()
  FUN <- function(at, d) {
    qunif(
      p = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      ...
    )
  }
  apply_dpqr(d = x, FUN = FUN, at = probs, type_prefix = "q", drop = drop)
}


#' Return the support of the Uniform distribution
#'
#' @param d An `Uniform` object created by a call to [Uniform()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Uniform <- function(d, drop = TRUE) {
  stopifnot("d must be a supported distribution object" = is_distribution(d))
  stopifnot(is.logical(drop))

  min <- d$a
  max <- d$b

  make_support(min, max, drop = drop)
}
