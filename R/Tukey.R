#' Create a Tukey distribution
#'
#' Tukey's studentized range distribution, used for Tukey's
#' honestly significant differences test in ANOVA.
#'
#' @param nmeans Sample size for each range.
#' @param df Degrees of freedom.
#' @param nranges Number of groups being compared.
#'
#' @return A `Tukey` object.
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
#'   **Support**: \eqn{R^+}, the set of positive real numbers.
#'
#'   Other properties of Tukey's Studentized Range Distribution
#'   are omitted, largely because the distribution is not fun
#'   to work with.
#'
#' @examples
#'
#' set.seed(27)
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
    glue(
      "Tukey distribution (nmeans = {x$nmeans},",
      "df = {x$df}, nranges = {x$nranges})\n"
    ),
    "\n"
  )
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


#' Return the support of the Tukey distribution
#'
#' @param d An `Tukey` object created by a call to [Tukey()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Tukey <- function(d){
  return(c(0, Inf))
}
