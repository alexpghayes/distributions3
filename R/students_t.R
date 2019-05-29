#' Create a Student's T distribution
#'
#' TODO
#'
#' @param df Degrees of freedom. Can be any positive number.
#'
#' @return A `students_t` object.
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
#'   In the following, let \eqn{X} be a students_t random variable with mean
#'   `mu` = \eqn{\mu} and standard deviation `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{\mathbb{R}}{R}, the set of all real numbers
#'
#'   **Mean**: Undefined unless \eqn{\nu \ge 2}, in which case the mean is
#'     zero.
#'
#'   **Variance**:
#'
#'   \deqn{
#'     \frac{\nu}{\nu - 2}
#'   }
#'
#'   Undefined if \eqn{\nu < 1}, infinite when \eqn{1 < \n \le 2}.
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
#' @examples
#'
#' # TODO
#'
students_t <- function(df) {
  d <- list(df = df)
  class(d) <- "students_t"
  d
}

#' @export
print.students_t <- function(x, ...) {
  cat(glue("Student's T distribution (df = {x$df})"))
}

#' Draw a random sample from a students_t distribution
#'
#' Please see the documentation of [students_t()] for some properties
#' of the T distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit students_t examples
#'
#' @param d A `students_t` object created by a call to [students_t()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family students_t distribution
#'
#' @return A numeric vector of length `n`.
#' @export
#'
#'
random.students_t <- function(d, n = 1L, ...) {
  rt(n = n, df = d$df)
}

#' Evaluate the probability mass function of a students_t distribution
#'
#' Please see the documentation of [students_t()] for some properties
#' of the students_t distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit students_t examples
#' @inheritParams random.students_t
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family students_t distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.students_t <- function(d, x, ...) {
  dt(x = x, df = d$df)
}

#' @rdname pdf.students_t
#' @export
#'
log_pdf.students_t <- function(d, x, ...) {
  dt(x = x, df = d$df, log = TRUE)
}

#' Evaluate the cumulative distribution function of a students_t distribution
#'
#' @inherit students_t examples
#' @inheritParams random.students_t
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family students_t distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.students_t <- function(d, x, ...) {
  pt(q = x, df = d$df)
}

#' Determine quantiles of a students_t distribution
#'
#' Please see the documentation of [students_t()] for some properties
#' of the students_t distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [students_t()] for some properties
#' of the students_t distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit students_t examples
#' @inheritParams random.students_t
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family students_t distribution
#'
quantile.students_t <- function(d, p, ...) {
  qt(p = p, df = d$df)
}
