#' Create a Logistic distribution
#'
#' TODO
#'
#' be sure to include connection to logistic regression
#'
#' @param location TODO
#' @param scale TODO
#'
#' @return A `Logistic` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Logistic random variable with
#'   success probability `p` = \eqn{p}.
#'
#'   TODO: multiple parameterizations BLEH
#'
#'   **Support**: TODO
#'
#'   **Mean**: TODO
#'
#'   **Variance**: TODO
#'
#'   **Probability density function (p.d.f)**:
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
#' X <- Logistic(2, 4)
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
Logistic <- function(location = 0, scale = 1) {
  d <- list(location = location, scale = scale)
  class(d) <- c("Logistic", "distribution")
  d
}

#' @export
print.Logistic <- function(x, ...) {
  cat(
    glue("Logistic distribution (location = {x$location}, scale = {x$scale})")
  )
}

#' Draw a random sample from a Logistic distribution
#'
#' @inherit Logistic examples
#'
#' @param d A `Logistic` object created by a call to [Logistic()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Logistic <- function(d, n = 1L, ...) {
  rlogis(n = n, location = d$location, scale = d$scale)
}

#' Evaluate the probability mass function of a Logistic distribution
#'
#' Please see the documentation of [Logistic()] for some properties
#' of the Logistic distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Logistic examples
#' @inheritParams random.Logistic
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Logistic <- function(d, x, ...) {
  dlogis(x = x, location = d$location, scale = d$scale)
}

#' @rdname pdf.Logistic
#' @export
log_pdf.Logistic <- function(d, x, ...) {
  dlogis(x = x, location = d$location, scale = d$scale, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Logistic distribution
#'
#' @inherit Logistic examples
#' @inheritParams random.Logistic
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Logistic <- function(d, x, ...) {
  plogis(q = x, location = d$location, scale = d$scale)
}

#' Determine quantiles of a Logistic distribution
#'
#' @inherit Logistic examples
#' @inheritParams random.Logistic
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family Logistic distribution
#'
quantile.Logistic <- function(d, p, ...) {
  qlogis(p = p, location = d$location, scale = d$scale)
}
