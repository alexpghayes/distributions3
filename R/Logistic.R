#' Create a Logistic distribution
#'
#' A continuous distribution on the real line. For binary outcomes
#' the model given by \eqn{P(Y = 1 | X) = F(X \beta)} where
#' \eqn{F} is the Logistic [cdf()] is called *logistic regression*.
#'
#' @param location The location parameter for the distribution. For Logistic
#'   distributions, the location parameter is the mean, median and also mode.
#'   Defaults to zero.
#'
#' @param scale The scale parameter for the distribution. Defaults to one.
#'
#' @return A `Logistic` object.
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
#'   In the following, let \eqn{X} be a Logistic random variable with
#'   `location` = \eqn{\mu} and `scale` = \eqn{s}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{s^2 \pi^2 / 3}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{e^{-(\frac{x - \mu}{s})}}{s [1 + \exp(-(\frac{x - \mu}{s})) ]^2}
#'   }{
#'     f(x) = e^(-(t - \mu) / s) / (s (1 + e^(-(t - \mu) / s))^2)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(t) = \frac{1}{1 + e^{-(\frac{t - \mu}{s})}}
#'   }{
#'     F(t) = 1 / (1 +  e^(-(t - \mu) / s))
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\mu t} \beta(1 - st, 1 + st)
#'   }{
#'     E(e^(tX)) = = e^(\mu t) \beta(1 - st, 1 + st)
#'   }
#'
#'   where \eqn{\beta(x, y)} is the Beta function.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Logistic(2, 4)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
Logistic <- function(location = 0, scale = 1) {
  d <- list(location = location, scale = scale)
  class(d) <- c("Logistic", "distribution")
  d
}

#' @export
print.Logistic <- function(x, ...) {
  cat(
    glue("Logistic distribution (location = {x$location}, scale = {x$scale})", "\n")
  )
}

#' @export
mean.Logistic <- function(x, ...) {
  ellipsis::check_dots_used()
  x$location
}

#' @export
variance.Logistic <- function(x, ...) x$scale^2 * pi^2 / 3

#' @export
skewness.Logistic <- function(x, ...) 0

#' @export
kurtosis.Logistic <- function(x, ...) 6 / 5

#' Draw a random sample from a Logistic distribution
#'
#' @inherit Logistic examples
#'
#' @param x A `Logistic` object created by a call to [Logistic()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Logistic <- function(x, n = 1L, ...) {
  rlogis(n = n, location = x$location, scale = x$scale)
}

#' Evaluate the probability mass function of a Logistic distribution
#'
#' Please see the documentation of [Logistic()] for some properties
#' of the Logistic distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Logistic examples
#'
#' @param d A `Logistic` object created by a call to [Logistic()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Logistic <- function(d, x, ...) {
  dlogis(x = x, location = d$location, scale = d$scale)
}

#' @rdname pdf.Logistic
#' @export
log_pdf.Logistic <- function(d, x, ...) {
  dlogis(x = x, location = d$location, scale = d$scale, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Logistic distribution
#'
#' @inherit Logistic examples
#'
#' @param d A `Logistic` object created by a call to [Logistic()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Logistic <- function(d, x, ...) {
  plogis(q = x, location = d$location, scale = d$scale)
}

#' Determine quantiles of a Logistic distribution
#'
#' @inherit Logistic examples
#' @inheritParams random.Logistic
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
#' @family Logistic distribution
#'
quantile.Logistic <- function(x, probs, ...) {
  ellipsis::check_dots_used()
  qlogis(p = probs, location = x$location, scale = x$scale)
}

#' Return the support of the Logistic distribution
#'
#' @param d An `Logistic` object created by a call to [Logistic()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Logistic <- function(d){
  if(!is_distribution(d)){
    message("d has to be a disitrubtion")
    stop()
  }
  return(c(-Inf, Inf))
}
