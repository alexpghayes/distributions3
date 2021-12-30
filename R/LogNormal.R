#' Create a LogNormal distribution
#'
#' A random variable created by exponentiating a [Normal()]
#' distribution. Taking the log of LogNormal data returns in
#' [Normal()] data.
#'
#' @param log_mu The location parameter, written \eqn{\mu} in textbooks.
#'   Can be any real number. Defaults to `0`.
#' @param log_sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   Can be any positive real number. Defaults to `1`.
#'
#' @return A `LogNormal` object.
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
#'   In the following, let \eqn{X} be a LogNormal random variable with
#'   success probability `p` = \eqn{p}.
#'
#'   **Support**: \eqn{R^+}
#'
#'   **Mean**: \eqn{\exp(\mu + \sigma^2/2)}
#'
#'   **Variance**: \eqn{[\exp(\sigma^2)-1]\exp(2\mu+\sigma^2)}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp \left(-\frac{(\log x - \mu)^2}{2 \sigma^2} \right)
#'   }{
#'     f(x) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp (-\frac{(\log x - \mu)^2}{2 \sigma^2})
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{F(x) = \frac{1}{2} + \frac{1}{2\sqrt{pi}}\int_{-x}^x e^{-t^2} dt}
#'
#'   **Moment generating function (m.g.f)**:
#'   Undefined.
#'
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- LogNormal(0.3, 2)
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
LogNormal <- function(log_mu = 0, log_sigma = 1) {
  d <- list(log_mu = log_mu, log_sigma = log_sigma)
  class(d) <- c("LogNormal", "distribution")
  d
}

#' @export
print.LogNormal <- function(x, ...) {
  cat(glue("Lognormal distribution (log_mu = {x$log_mu}, log_sigma = {x$log_sigma})"), "\n")
}

#' @export
mean.LogNormal <- function(x, ...) {
  ellipsis::check_dots_used()
  mu <- x$log_mu
  sigma <- x$log_sigma
  exp(mu + sigma^2 / 2)
}

#' @export
variance.LogNormal <- function(x, ...) {
  mu <- x$log_mu
  sigma <- x$log_sigma
  (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)
}

#' @export
skewness.LogNormal <- function(x, ...) {
  mu <- x$log_mu
  sigma <- x$log_sigma
  (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
}

#' @export
kurtosis.LogNormal <- function(x, ...) {
  mu <- x$log_mu
  sigma <- x$log_sigma
  exp(4 * sigma^2) + 2 * exp(3 * sigma^2) + 3 * exp(2 * sigma^2) - 6
}

#' Draw a random sample from a LogNormal distribution
#'
#' @inherit LogNormal examples
#'
#' @param x A `LogNormal` object created by a call to [LogNormal()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family LogNormal distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.LogNormal <- function(x, n = 1L, ...) {
  rlnorm(n = n, meanlog = x$log_mu, sdlog = x$log_sigma)
}

#' Evaluate the probability mass function of a LogNormal distribution
#'
#' Please see the documentation of [LogNormal()] for some properties
#' of the LogNormal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit LogNormal examples
#'
#' @param d A `LogNormal` object created by a call to [LogNormal()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family LogNormal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.LogNormal <- function(d, x, ...) {
  dlnorm(x = x, meanlog = d$log_mu, sdlog = d$log_sigma)
}

#' @rdname pdf.LogNormal
#' @export
log_pdf.LogNormal <- function(d, x, ...) {
  dlnorm(x = x, meanlog = d$log_mu, sdlog = d$log_sigma, log = TRUE)
}

#' Evaluate the cumulative distribution function of a LogNormal distribution
#'
#' @inherit LogNormal examples
#'
#' @param d A `LogNormal` object created by a call to [LogNormal()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family LogNormal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.LogNormal <- function(d, x, ...) {
  plnorm(q = x, meanlog = d$log_mu, sdlog = d$log_sigma)
}

#' Determine quantiles of a LogNormal distribution
#'
#' @inherit LogNormal examples
#' @inheritParams random.LogNormal
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
#' @family LogNormal distribution
#'
quantile.LogNormal <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qlnorm(p = probs, meanlog = x$log_mu, sdlog = x$log_sigma)
}

#' Fit a Log Normal distribution to data
#'
#' @param d A `LogNormal` object created by a call to [LogNormal()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family LogNormal distribution
#'
#' @return A `LogNormal` object.
#' @export
#'
fit_mle.LogNormal <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  LogNormal(ss$mu, ss$sigma)
}

#' Compute the sufficient statistics for a Log-normal distribution from data
#'
#' @inheritParams fit_mle.LogNormal
#'
#' @return A named list of the sufficient statistics of the normal distribution:
#'
#'   - `mu`: The sample mean of the log of the data.
#'   - `sigma`: The sample standard deviation of the log of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
#'
suff_stat.LogNormal <- function(d, x, ...) {
  valid_x <- x > 0
  if (any(!valid_x)) stop("`x` must be a vector of positive real numbers")
  log_x <- log(x)
  list(mu = mean(log_x), sigma = sd(log_x), samples = length(x))
}

#' Return the support of the LogNormal distribution
#'
#' @param d An `LogNormal` object created by a call to [LogNormal()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.LogNormal <- function(d){
  if(!is_distribution(d)){
    message("d has to be a disitrubtion")
    stop()
  }
  return(c(0, Inf))
}
