#' Create a Tukey distribution
#'
#' Tukey's studentized range distribution, used mostly for Tukey's
#' honestly significant differences test in ANOVA
#'
#' @param nmeans TODO
#' @param df TODO
#' @param nranges TODO
#'
#' @return A `Tukey` object.
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
#'   In the following, let \eqn{X} be a Tukey random variable with
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
#' X <- Tukey(4L, 16L, 2L)
#' X
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
Tukey <- function(nmeans, df, nranges) {
  d <- list(nmeans = nmeans, df = df, nranges = nranges)
  class(d) <- c("Tukey", "distribution")
  d
}

#' @export
print.Tukey <- function(x, ...) {
  cat(
    glue("Tukey distribution (nmeans = {x$nmeans}",
         "df = {x$df}, nranges = {x$nranges})"))
}

#' Evaluate the cumulative distribution function of a Tukey distribution
#'
#' @inherit Tukey examples
#'
#' @param d A `Tukey` distribution created by a call to [Tukey()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#'
#' @family Tukey distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Tukey <- function(d, x, ...) {
  ptukey(q = x, nmeans = d$nmeans, df = d$nmeans, nranges = d$nranges)
}

#' Determine quantiles of a Tukey distribution
#'
#' @inherit Tukey examples
#' @inheritParams cdf.Tukey
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family Tukey distribution
#'
quantile.Tukey <- function(d, p, ...) {
  qtukey(p = p, nmeans = d$nmeans, df = d$nmeans, nranges = d$nranges)
}
