#' Create a Student's T distribution
#'
#' TODO
#'
#' @param df Degrees of freedom. Can be any positive number.
#'
#' @return A `StudentsT` object.
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
#'   In the following, let \eqn{X} be a StudentsT random variable with
#'   `df` = \eqn{\nu}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
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
#'   Undefined if \eqn{\nu < 1}, infinite when \eqn{1 < \nu \le 2}.
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
StudentsT <- function(df) {
  d <- list(df = df)
  class(d) <- c("StudentsT", "distribution")
  d
}

#' @export
print.StudentsT <- function(x, ...) {
  cat(glue("Student's T distribution (df = {x$df})"))
}

#' Draw a random sample from a StudentsT distribution
#'
#' Please see the documentation of [StudentsT()] for some properties
#' of the T distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit StudentsT examples
#'
#' @param d A `StudentsT` object created by a call to [StudentsT()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family StudentsT distribution
#'
#' @return A numeric vector of length `n`.
#' @export
#'
#'
random.StudentsT <- function(d, n = 1L, ...) {
  rt(n = n, df = d$df)
}

#' Evaluate the probability mass function of a StudentsT distribution
#'
#' Please see the documentation of [StudentsT()] for some properties
#' of the StudentsT distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit StudentsT examples
#' @inheritParams random.StudentsT
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family StudentsT distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.StudentsT <- function(d, x, ...) {
  dt(x = x, df = d$df)
}

#' @rdname pdf.StudentsT
#' @export
#'
log_pdf.StudentsT <- function(d, x, ...) {
  dt(x = x, df = d$df, log = TRUE)
}

#' Evaluate the cumulative distribution function of a StudentsT distribution
#'
#' @inherit StudentsT examples
#' @inheritParams random.StudentsT
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family StudentsT distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.StudentsT <- function(d, x, ...) {
  pt(q = x, df = d$df)
}

#' Determine quantiles of a StudentsT distribution
#'
#' Please see the documentation of [StudentsT()] for some properties
#' of the StudentsT distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of
#' [StudentsT()] for some properties
#' of the StudentsT distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit StudentsT examples
#' @inheritParams random.StudentsT
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family StudentsT distribution
#'
quantile.StudentsT <- function(d, p, ...) {
  qt(p = p, df = d$df)
}
