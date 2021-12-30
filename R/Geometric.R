#' Create a Geometric distribution
#'
#' The Geometric distribution can be thought of as a generalization
#' of the [Bernoulli()] distribution where we ask: "if I keep flipping a
#' coin with probability `p` of heads, what is the probability I need
#' \eqn{k} flips before I get my first heads?" The Geometric
#' distribution is a special case of Negative Binomial distribution.
#'
#' @param p The success probability for the distribution. `p` can be
#'   any value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `Geometric` object.
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
#'   In the following, let \eqn{X} be a Geometric random variable with
#'   success probability `p` = \eqn{p}. Note that there are multiple
#'   parameterizations of the Geometric distribution.
#'
#'   **Support**: 0 < p < 1, \eqn{x = 0, 1, \dots}
#'
#'   **Mean**: \eqn{\frac{1-p}{p}}
#'
#'   **Variance**: \eqn{\frac{1-p}{p^2}}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = x) = p(1-p)^x,
#'    }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le x) = 1 - (1-p)^{x+1}
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \frac{pe^t}{1 - (1-p)e^t}
#'   }{
#'     E(e^{tX}) = \frac{pe^t}{1 - (1-p)e^t}
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Geometric(0.3)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
Geometric <- function(p = 0.5) {
  d <- list(p = p)
  class(d) <- c("Geometric", "distribution")
  d
}

#' @export
print.Geometric <- function(x, ...) {
  cat(glue("Geometric distribution (p = {x$p})"), "\n")
}

#' @export
mean.Geometric <- function(x, ...) {
  ellipsis::check_dots_used()
  1 / x$p
}
#' @export
variance.Geometric <- function(x, ...) (1 - x$p) / x$p^2

#' @export
skewness.Geometric <- function(x, ...) (2 - x$p) / sqrt(1 - x$p)

#' @export
kurtosis.Geometric <- function(x, ...) 6 + (x$p^2 / (1 - x$p))

#' Draw a random sample from a Geometric distribution
#'
#' Please see the documentation of [Geometric()] for some properties
#' of the Geometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Geometric examples
#'
#' @param x A `Geometric` object created by a call to [Geometric()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Geometric distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Geometric <- function(x, n = 1L, ...) {
  rgeom(n = n, prob = x$p)
}

#' Evaluate the probability mass function of a Geometric distribution
#'
#' Please see the documentation of [Geometric()] for some properties
#' of the Geometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Geometric examples
#'
#' @param d A `Geometric` object created by a call to [Geometric()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Geometric distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Geometric <- function(d, x, ...) {
  dgeom(x = x, prob = d$p)
}

#' @rdname pdf.Geometric
#' @export
#'
log_pdf.Geometric <- function(d, x, ...) {
  dgeom(x = x, prob = d$p, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Geometric distribution
#'
#' @inherit Geometric examples
#'
#' @param d A `Geometric` object created by a call to [Geometric()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Geometric distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Geometric <- function(d, x, ...) {
  pgeom(q = x, prob = d$p)
}

#' Determine quantiles of a Geometric distribution
#'
#' @inherit Geometric examples
#' @inheritParams random.Geometric
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
#' @family Geometric distribution
#'
quantile.Geometric <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qgeom(p = probs, prob = x$p)
}

#' Fit a Geometric distribution to data
#'
#' @param d A `Geometric` object.
#' @param x A vector of zeroes and ones.
#' @param ... Unused.
#'
#' @return a `Geometric` object
#' @export
#'
fit_mle.Geometric <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Geometric(1 / (ss$trials / ss$experiments + 1))
}

#' Compute the sufficient statistics for the Geometric distribution from data
#'
#' @inheritParams fit_mle.Geometric
#'
#' @return A named list of the sufficient statistics of the Geometric
#'   distribution:
#'
#'   - `trials`: The total number of trials ran until the first success.
#'   - `experiments`: The number of experiments run.
#'
#' @export
#'
suff_stat.Geometric <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x %% 1 == 0)
  if (any(!valid_x)) {
    stop("`x` must be a vector of positive discrete numbers")
  }
  list(trials = sum(x), experiments = length(x))
}

#' Return the support of the Geometric distribution
#'
#' @param d An `Geometric` object created by a call to [Geometric()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Geometric <- function(d){
  c(0, Inf)
}
