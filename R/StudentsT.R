#' Create a Student's T distribution
#'
#' The Student's T distribution is closely related to the [Normal()]
#' distribution, but has heavier tails. As \eqn{\nu} increases to \eqn{\infty},
#' the Student's T converges to a Normal. The T distribution appears
#' repeatedly throughout classic frequentist hypothesis testing when
#' comparing group means.
#'
#' @param df Degrees of freedom. Can be any positive number. Often
#'   called \eqn{\nu} in textbooks.
#'
#' @return A `StudentsT` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Students T random variable with
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
#'   }{
#'     \nu / (\nu - 2)
#'   }
#'
#'   Undefined if \eqn{\nu < 1}, infinite when \eqn{1 < \nu \le 2}.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\Gamma(\frac{\nu + 1}{2})}{\sqrt{\nu \pi} \Gamma(\frac{\nu}{2})} (1 + \frac{x^2}{\nu} )^{- \frac{\nu + 1}{2}}
#'   }{
#'     f(x) = \Gamma((\nu + 1) / 2) / (\sqrt(\nu \pi) \Gamma(\nu / 2)) (1 + x^2 / \nu)^(- (\nu + 1) / 2)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   Nasty, omitted.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   Undefined.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- StudentsT(3)
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
#' ### example: calculating p-values for two-sided T-test
#'
#' # here the null hypothesis is H_0: mu = 3
#'
#' # data to test
#' x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' nx <- length(x)
#'
#' # calculate the T-statistic
#' t_stat <- (mean(x) - 3) / (sd(x) / sqrt(nx))
#' t_stat
#'
#' # null distribution of statistic depends on sample size!
#' T <- StudentsT(df = nx - 1)
#'
#' # calculate the two-sided p-value
#' 1 - cdf(T, abs(t_stat)) + cdf(T, -abs(t_stat))
#'
#' # exactly equivalent to the above
#' 2 * cdf(T, -abs(t_stat))
#'
#' # p-value for one-sided test
#' # H_0: mu <= 3   vs   H_A: mu > 3
#' 1 - cdf(T, t_stat)
#'
#' # p-value for one-sided test
#' # H_0: mu >= 3   vs   H_A: mu < 3
#' cdf(T, t_stat)
#'
#' ### example: calculating a 88 percent T CI for a mean
#'
#' # lower-bound
#' mean(x) - quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#'
#' # upper-bound
#' mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#'
#' # equivalent to
#' mean(x) + c(-1, 1) * quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#'
#' # also equivalent to
#' mean(x) + quantile(T, 0.12 / 2) * sd(x) / sqrt(nx)
#' mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
StudentsT <- function(df) {
  d <- list(df = df)
  class(d) <- c("StudentsT", "distribution")
  d
}

#' @export
print.StudentsT <- function(x, ...) {
  cat(glue("Student's T distribution (df = {x$df})"), "\n")
}

#' @export
mean.StudentsT <- function(d, ...) if (d$df > 1) 0 else NaN

#' @export
variance.StudentsT <- function(d, ...) {
  if (d$df > 2) d$df / (d$df - 2)
  else if (d$df > 1) Inf
  else NaN
}

#' @export
skewness.StudentsT <- function(d, ...) if(d$df > 3) 0 else NaN

#' @export
kurtosis.StudentsT <- function(d, ...) {
  if (d$df > 4) 6 / (d$df - 4)
  else if (d$df > 2) Inf
  else NaN
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


#' Return the support of the StudentsT distribution
#'
#' @param d An `StudentsT` object created by a call to [StudentsT()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.StudentsT <- function(d){
  return(c(-Inf, Inf))
}
