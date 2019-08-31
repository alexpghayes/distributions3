#' Create a HyperGeometric distribution
#'
#' To understand the HyperGeometric distribution, consider a set of
#' \eqn{r} objects, of which \eqn{m} are of the type I and
#' \eqn{n} are of the type II. A sample with size \eqn{k} (\eqn{k<r})
#'  with no replacement is randomly chosen. The number of observed
#'  type I elements observed in this sample is set to be our random
#'  variable \eqn{X}. For example, consider that in a set of 20
#'  car parts, there are 4 that are defective (type I).
#'  If we take a sample of size 5 from those car parts, the
#'  probability of finding 2 that are defective will be given by
#'  the HyperGeometric distribution (needs double checking).
#'
#'
#' @param m The number of type I elements available.
#' @param n The number of type II elements available.
#' @param k The size of the sample taken.
#'
#' @return A `HyperGeometric` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a HyperGeometric random variable with
#'   success probability `p` = \eqn{p = m/(m+n)}.
#'
#'   **Support**: \eqn{x \in { \{\max{(0, k-(n-m)}, \dots, \min{(k,m)}}\}}
#'
#'   **Mean**: \eqn{\frac{km}{n+m} = kp}
#'
#'   **Variance**: \eqn{\frac{km(n)(n+m-k)}{(n+m)^2 (n+m-1)} =
#'   kp(1-p)(1 - \frac{k-1}{m+n-1})}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
#'   }{
#'     P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) \approx \Phi\Big(\frac{x - kp}{\sqrt{kp(1-p)}}\Big)
#'  }
#'   **Moment generating function (m.g.f)**:
#'
#'   Not useful.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- HyperGeometric(4, 5, 8)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
HyperGeometric <- function(m, n, k) {
  d <- list(m = m, n = n, k = k)
  class(d) <- c("HyperGeometric", "distribution")
  d
}

#' @export
print.HyperGeometric <- function(x, ...) {
  cat(glue("HyperGeometric distribution (m = {x$m}, n = {x$n}, k = {x$k})\n"))
}

#' Draw a random sample from a HyperGeometric distribution
#'
#' Please see the documentation of [HyperGeometric()] for some properties
#' of the HyperGeometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit HyperGeometric examples
#'
#' @param d A `HyperGeometric` object created by a call to [HyperGeometric()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family HyperGeometric distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.HyperGeometric <- function(d, n = 1L, ...) {
  rhyper(nn = n, m = d$m, n = d$n, k = d$k)
}

#' Evaluate the probability mass function of a HyperGeometric distribution
#'
#' Please see the documentation of [HyperGeometric()] for some properties
#' of the HyperGeometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit HyperGeometric examples
#' @inheritParams random.HyperGeometric
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family HyperGeometric distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.HyperGeometric <- function(d, x, ...) {
  dhyper(x = x, m = d$m, n = d$n, k = d$k)
}

#' @rdname pdf.HyperGeometric
#' @export
log_pdf.HyperGeometric <- function(d, x, ...) {
  dhyper(x = x, m = d$m, n = d$n, k = d$k, log = TRUE)
}

#' Evaluate the cumulative distribution function of a HyperGeometric distribution
#'
#' @inherit HyperGeometric examples
#' @inheritParams random.HyperGeometric
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family HyperGeometric distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.HyperGeometric <- function(d, x, ...) {
  phyper(q = x, m = d$m, n = d$n, k = d$k)
}

#' Determine quantiles of a HyperGeometric distribution
#'
#' @inherit HyperGeometric examples
#' @inheritParams random.HyperGeometric
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family HyperGeometric distribution
#'
quantile.HyperGeometric <- function(d, p, ...) {
  qhyper(p = p, m = d$m, n = d$n, k = d$k)
}
