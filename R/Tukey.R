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
#'   <https://alexpghayes.github.io/distributions3/>, where the math
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
Tukey <- function(nmeans, df, nranges) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(nmeans) == length(df) & length(nmeans) == length(nranges) |
        sum(c(length(nmeans) == 1, length(df) == 1, length(nranges) == 1)) >= 2 |
        length(nmeans) == length(df) & length(nranges) == 1 |
        length(nmeans) == length(nranges) & length(df) == 1 |
        length(df) == length(nranges) & length(nmeans) == 1
  )

  d <- data.frame(nmeans = nmeans, df = df, nranges = nranges)
  class(d) <- c("Tukey", "distribution")
  d
}

#' Evaluate the cumulative distribution function of a Tukey distribution
#'
#' @inherit Tukey examples
#'
#' @param d A `Tukey` distribution created by a call to [Tukey()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{ptukey}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Tukey distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Tukey <- function(d, x, drop = TRUE, ...) {
  FUN <- function(at, d) ptukey(q = at, nmeans = d$nmeans, df = d$df, nranges = d$nranges, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop)
}

#' Determine quantiles of a Tukey distribution
#'
#' @inherit Tukey examples
#' @inheritParams cdf.Tukey
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments to be passed to \code{\link[stats]{qtukey}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#' @export
#'
#' @family Tukey distribution
#'
quantile.Tukey <- function(x, probs, drop = TRUE, ...) {
  ellipsis::check_dots_used()
  FUN <- function(at, d) qtukey(p = at, nmeans = x$nmeans, df = x$df, nranges = x$nranges, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop)
}


#' Return the support of the Tukey distribution
#'
#' @param d An `Tukey` object created by a call to [Tukey()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Tukey <- function(d, drop = TRUE) {
  stopifnot("d must be a supported distribution object" = is_distribution(d))
  stopifnot(is.logical(drop))

  min <- rep(0, length(d))
  max <- rep(Inf, length(d))

  make_support(min, max, d, drop = drop)
}
