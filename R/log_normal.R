#' Create a log_normal distribution
#'
#' @param log_mu TODO explain what the hell this is
#' @param log_sigma TODO explain what the hell this is
#'
#' @return A `log_normal` object.
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
#'   In the following, let \eqn{X} be a log_normal random variable with
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
#' ln <- log_normal(0.3, 2)
#' ln
#'
#' random(ln, 10)
#' pdf(ln, 2)
#' log_pdf(ln, 2)
#' cdf(ln, 4)
#' quantile(ln, 0.7)
#'
log_normal <- function(log_mu = 0, log_sigma = 1) {
  d <- list(log_mu = log_mu, log_sigma = log_sigma)
  class(d) <- c("log_normal", "distribution")
  d
}

#' @export
print.log_normal <- function(x, ...) {
  cat(glue("Lognormal distribution (log_mu = {x$log_mu}, log_sigma = {x$log_sigma})"))
}

#' Draw a random sample from a log_normal distribution
#'
#' @inherit log_normal examples
#'
#' @param d A `log_normal` object created by a call to [log_normal()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family log_normal distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.log_normal <- function(d, n = 1L, ...) {
  rlnorm(n = n, meanlog = d$log_mu, sdlog = d$log_sigma)
}

#' Evaluate the probability mass function of a log_normal distribution
#'
#' Please see the documentation of [log_normal()] for some properties
#' of the log_normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit log_normal examples
#' @inheritParams random.log_normal
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family log_normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.log_normal <- function(d, x, ...) {
  dlnorm(x = x, meanlog = d$log_mu, sdlog = d$log_sigma)
}

#' @rdname pdf.log_normal
#' @export
#'
log_pdf.log_normal <- function(d, x, ...) {
  dlnorm(x = x, meanlog = d$log_mu, sdlog = d$log_sigma, log = TRUE)
}

#' Evaluate the cumulative distribution function of a log_normal distribution
#'
#' @inherit log_normal examples
#' @inheritParams random.log_normal
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family log_normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.log_normal <- function(d, x, ...) {
  plnorm(q = x, meanlog = d$log_mu, sdlog = d$log_sigma)
}

#' Determine quantiles of a log_normal distribution
#'
#' @inherit log_normal examples
#' @inheritParams random.log_normal
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family log_normal distribution
#'
quantile.log_normal <- function(d, p, ...) {
  qlnorm(p = p, meanlog = d$log_mu, sdlog = d$log_sigma)
}

#' Fit a Log-normal distribution to data
#'
#' @param d A `log_normal` object created by a call to [log_normal()].
#' @param x A vector of data.
#'
#' @family log_normal distribution
#'
#' @return A `log_normal` object.
#' @export
#'
fit_mle.log_normal <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  log_normal(ss$mu, ss$sigma)
}

#' Compute the sufficient statistics for a Log-normal distribution from data
#'
#' @inheritParams fit_mle.log_normal
#'
#' @return A named list of the sufficient statistics of the normal distribution:
#'
#'   - `mu`: The sample mean of the log of the data.
#'   - `sigma`: The sample standard deviation of the log of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
#'
suff_stat.log_normal <- function(d, x, ...) {
  valid_x <- x > 0
  if(any(!valid_x)) stop("`x` must be a vector of positive real numbers")
  log_x <- log(x)
  list(mu = mean(log_x), sigma = sd(log_x), samples = length(x))
}
