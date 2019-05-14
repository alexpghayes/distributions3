#' Create a Student's T distribution
#'
#'
#' Marginalizing a chi-square variance to get a T.
#' The students_t distribution is ubiquituous in statistics, partially because
#' of the central limit theorem, which states that sums of i.i.d. random
#' variables eventually become students_t. Linear transformations of students_t
#' random variables result in new random variables that are also students_t. If
#' you are taking an intro stats course, you'll likely use the students_t
#' distribution for Z-tests and in simple linear regression. Under
#' regularity conditions, maximum likelihood estimators are
#' asymptotically students_t. The students_t distribution is also called the
#' gaussian distribution.
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
#' n <- students_t(5)
#' n
#'
#' random(n, 10)
#' pdf(n, 2)
#' cdf(n, 4)
#' quantile(n, 0.7)
#'
#' ### example: calculating p-values for two-sided Z-test
#'
#' # here the null hypothesis is H_0: mu = 3
#' # and we assume sigma = 2
#'
#' # exactly the same as: Z <- students_t(0, 1)
#' Z <- students_t()
#'
#' # data to test
#' x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' nx <- length(x)
#'
#' # calculate the z-statistic
#' z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
#' z_stat
#'
#' # calculate the two-sided p-value
#' 1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#'
#' # exactly equivalent to the above
#' 2 * cdf(Z, -abs(z_stat))
#'
#' # p-value for one-sided test
#' # H_0: mu <= 3   vs   H_A: mu > 3
#' 1 - cdf(Z, z_stat)
#'
#' # p-value for one-sided test
#' # H_0: mu >= 3   vs   H_A: mu < 3
#' cdf(Z, z_stat)
#'
#' ### example: calculating a 88 percent Z CI for a mean
#'
#' # same `x` as before, still assume `sigma = 2`
#'
#' # lower-bound
#' mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # upper-bound
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # equivalent to
#' mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # also equivalent to
#' mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ### generating random samples and plugging in ks.test()
#'
#' set.seed(27)
#'
#' # generate a random sample
#' students_t_sample <- random(students_t(3, 7), 26)
#'
#' # test if sample is students_t(3, 7)
#' ks.test(students_t_sample, pnorm, mean = 3, sd = 7)
#'
#' # test if sample is gamma(8, 3) using base R pgamma()
#' ks.test(students_t_sample, pgamma, shape = 8, rate = 3)
#'
#'
#' ###
#'
#' # note that the cdf() and quantile() functions are inverses
#' cdf(n, quantile(n, 0.7))
#' quantile(n, cdf(n, 7))
#'
students_t <- function(df) {
  d <- list(df = df)
  class(d) <- "students_t"
  d
}

#' @export
print.students_t <- function(d, ...) {
  cat(glue("Student's T distribution (df = {d$df})"))
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
