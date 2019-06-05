#' Create a beta distribution
#'
#'
#' @param alpha The alpha parameter. `alpha` can be any value strictly
#'   greater than zero. Defaults to `1`.
#' @param beta The beta parameter. `beta` can be any value strictly
#'   greater than zero. Defaults to `1`.
#'
#' @return A `beta` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' b <- beta(1, 2)
#' b
#'
#' random(b, 10)
#' pdf(b, 0.7)
#' log_pdf(b, 0.7)
#' cdf(b, 0.7)
#' quantile(b, 0.7)
#'
#' cdf(b, quantile(b, 0.7))
#' quantile(b, cdf(b, 0.7))
#'
beta <- function(alpha = 1, beta = 1) {
  d <- list(alpha = alpha, beta = beta)
  class(d) <- c("beta", "distribution")
  d
}

#' @export
print.beta <- function(x, ...) {
  cat(glue("Beta distribution (alpha = {x$alpha}, beta = {x$beta})"))
}

#' Draw a random sample from a beta distribution
#'
#' @inherit beta examples
#'
#' @param d A `beta` object created by a call to [beta()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector containing values in `[0, 1]` of length `n`.
#' @export
#'
random.beta <- function(d, n = 1L, ...) {
  rbeta(n = n, shape1 = d$alpha, shape2 = d$beta)
}

#' Evaluate the probability mass function of a beta distribution
#'
#' @inherit beta examples
#' @inheritParams random.beta
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.beta <- function(d, x, ...) {
  dbeta(x = x, shape1 = d$alpha, shape2 = d$beta)
}

#' @rdname pdf.beta
#' @export
#'
log_pdf.beta <- function(d, x, ...) {
  dbeta(x = x, shape1 = d$alpha, shape2 = d$beta, log = TRUE)
}

#' Evaluate the cumulative distribution function of a beta distribution
#'
#' @inherit beta examples
#' @inheritParams random.beta
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.beta <- function(d, x, ...) {
  pbeta(q = x, shape1 = d$alpha, shape2 = d$beta)
}

#' Determine quantiles of a beta distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit beta examples
#' @inheritParams random.beta
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.beta <- function(d, p, ...) {
  qbeta(p = p, shape1 = d$alpha, shape2 = d$beta)
}
