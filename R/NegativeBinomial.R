#' Create a Negative Binomial distribution
#'
#' A generalization of the geometric distribution. It is the number
#' of successes in a sequence of i.i.d. Bernoulli trials before
#' a specified number (\eqn{r}) of failures occurs.
#'
#'
#' @param size The number of failures (an integer greater than \eqn{0})
#'   until the experiment is stopped. Denoted \eqn{r} below.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `NegativeBinomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Negative Binomial random variable with
#'   success probability `p` = \eqn{p}.
#'
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}
#'
#'   **Mean**: \eqn{\frac{p r}{1-p}}
#'
#'   **Variance**: \eqn{\frac{pr}{(1-p)^2}}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'      f(k) = {k + r - 1 \choose k} \cdot (1-p)^r p^k
#'   }{
#'      f(k) = (k+r-1)!/(k!(r-1)!) (1-p)^r p^k
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   Omitted for now.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'      \left(\frac{1-p}{1-pe^t}\right)^r, t < -\log p
#'   }{
#'      \frac{(1-p)^r}{(1-pe^t)^r}, t < -\log p
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- NegativeBinomial(10, 0.3)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
NegativeBinomial <- function(size, p = 0.5) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(size) == length(p) | length(size) == 1 | length(p) == 1
  )
  d <- data.frame(size = size, p = p)
  class(d) <- c("NegativeBinomial", "distribution")
  d
}

#' @export
mean.NegativeBinomial <- function(x, ...) {
  ellipsis::check_dots_used()
  rval <- x$p * x$size / (1 - x$p)
  setNames(rval, names(x))
}

#' @export
variance.NegativeBinomial <- function(x, ...) {
  rval <- (x$p * x$size) / (1 - x$p)^2
  setNames(rval, names(x))
}

#' @export
skewness.NegativeBinomial <- function(x, ...) {
  rval <- (1 + x$p) / sqrt(x$p * x$size)
  setNames(rval, names(x))
}

#' @export
kurtosis.NegativeBinomial <- function(x, ...) {
  rval <- 6 / x$size + (1 - x$p)^2 / x$size * x$p
  setNames(rval, names(x))
}

#' Draw a random sample from a negative binomial distribution
#'
#' @inherit NegativeBinomial examples
#'
#' @param x A `NegativeBinomial` object created by a call to
#'   [NegativeBinomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family NegativeBinomial distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.NegativeBinomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rnbinom(n = length(d), size = d$size, prob = d$p)
  apply_dpqr(d = x, FUN = FUN, at = matrix(1, ncol = n), type = "random", drop = drop)
}

#' Evaluate the probability mass function of a NegativeBinomial distribution
#'
#' @inherit NegativeBinomial examples
#'
#' @param d A `NegativeBinomial` object created by a call to
#'   [NegativeBinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{dnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family NegativeBinomial distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.NegativeBinomial <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) dnbinom(x = at, size = d$size, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop)
}

#' @rdname pdf.NegativeBinomial
#' @export
#'
log_pdf.NegativeBinomial <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) dnbinom(x = at, size = d$size, prob = d$p, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop)
}

#' Evaluate the cumulative distribution function of a negative binomial distribution
#'
#' @inherit NegativeBinomial examples
#'
#' @param d A `NegativeBinomial` object created by a call to
#'   [NegativeBinomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{pnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family NegativeBinomial distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.NegativeBinomial <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) pnbinom(q = at, size = d$size, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop)
}

#' Determine quantiles of a NegativeBinomial distribution
#'
#' @inherit NegativeBinomial examples
#' @inheritParams random.NegativeBinomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{qnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#' @export
#'
#' @family NegativeBinomial distribution
#'
quantile.NegativeBinomial <- function(x, probs, drop = TRUE, ...) {
  ellipsis::check_dots_used()
  FUN <- function(at, d) qnbinom(p = at, size = x$size, prob = x$p, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop)
}


#' Return the support of the NegativeBinomial distribution
#'
#' @param d An `NegativeBinomial` object created by a call to [NegativeBinomial()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.NegativeBinomial <- function(d, drop = TRUE) {
  stopifnot("d must be a supported distribution object" = is_distribution(d))
  stopifnot(is.logical(drop))

  min <- rep(0, length(d))
  max <- rep(Inf, length(d))

  make_support(min, max, d, drop = drop)
}
