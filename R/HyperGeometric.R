#' Create a HyperGeometric distribution
#'
#' TODO
#'
#' @param m TODO
#' @param n TODO
#' @param k TODO
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
#'   success probability `p` = \eqn{p}.
#'
#'   **Support**: TODO
#'
#'   **Mean**: TODO
#'
#'   **Variance**: TODO
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   TODO
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   TODO
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   TODO
#'
#' @examples
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
#'
HyperGeometric <- function(m, n, k) {
  d <- list(m = m, n = n, k = k)
  class(d) <- c("HyperGeometric", "distribution")
  d
}

#' @export
print.HyperGeometric <- function(x, ...) {
  cat(glue("HyperGeometric distribution (m = {x$m}, n = {x$n}, k = {x$k})"))
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

