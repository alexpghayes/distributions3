#' Create a Pareto distribution
#'
#' The (Type I) Pareto distribution is a power-law distribution originating
#' from Pareto's Law, which was used to describe the distribution of wealth in
#' a society.  It is also known as *Pareto distribution of the first kind*.
#'
#' @param k The scale parameter. `k` can be any positive number.
#'   Defaults to `1`.
#' @param a The shape parameter. `a` can be any positive number.
#'   Defaults to `1`.
#'
#' @return A `Pareto` object.
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
#'   In the following, let \eqn{X} be a Pareto random variable with scale
#'   parameter \eqn{k} and scale parameter \eqn{a}.
#'
#'   **Support**: \eqn{[k, \infty)}.
#'
#'   **Mean**: \eqn{a k / (a - 1)}, for \eqn{a > 1}.
#'     Undefined otherwise.
#'
#'   **Median**: \eqn{2 ^ {1 / a} k}{2 ^ (1 / a) k}.
#'
#'   **Mode**: \eqn{k}.
#'
#'   **Variance**: \eqn{a k ^ 2 / (a - 1)^2 (a - 2)}, for \eqn{a > 2}.
#'     Undefined otherwise.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = a k^a / x^{a + 1}
#'   }{
#'     f(x) = a k^a / x^(a + 1)
#'   }
#'   for \eqn{x \geq k}{x >= k}.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(x) = 1 - (k / x)^a
#'   }{
#'     F(x) = 1 - (k / x)^a
#'   }
#'   for \eqn{x \geq k}{x >= k}.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = a(-kt)^a \Gamma(-a, -kt)
#'   }{
#'     E(e^(tX)) = a(-kt)^a \Gamma(-a, -kt)
#'   }
#'   for \eqn{t < 0}, where \eqn{\Gamma(a, x)} is the (upper) incomplete gamma
#'   function. For \eqn{t = 0}, \eqn{E(e^{tX}) = 1}.  The m.g.f. is undefined
#'   for \eqn{t > 0}.
#' @examples
#'
#' set.seed(27)
#'
#' X <- Pareto(1, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 1.7)
#' log_pdf(X, 1.7)
#'
#' cdf(X, 1.7)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 1.7))
#'
#' x <- random(X, 100)
#' fit_mle(X, x)
Pareto <- function(k = 1, a = 1) {
  d <- list(k = k, a = a)
  class(d) <- c("Pareto", "distribution")
  d
}

#' @export
print.Pareto <- function(x, ...) {
  cat(glue("Pareto distribution (k = {x$k}, a = {x$a})\n"))
}

#' Draw a random sample from a Pareto distribution
#'
#' @inherit Pareto examples
#'
#' @param d A `Pareto` object created by a call to [Pareto()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
random.Pareto <- function(d, n = 1L, ...) {
  # Mimic the vectorised behaviour of R's base random generation functions
  max_len <- ifelse(length(n) > 1L, length(n), n)
  d$k <- rep_len(d$k, max_len)
  d$a <- rep_len(d$a, max_len)
  quantile(d = d, p = runif(n))
}

#' Evaluate the probability mass function of a Pareto distribution
#'
#' @inherit Pareto examples
#' @inheritParams random.Pareto
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Pareto <- function(d, x, ...) {
  # Mimic the vectorised behaviour of R's base random generation functions
  if (length(x) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(x), length(d$k), length(d$a))
  x <- rep_len(x, max_len)
  d$k <- rep_len(d$k, max_len)
  d$a <- rep_len(d$a, max_len)
  lnd <- ifelse(x < d$k, -Inf, log(d$a) + d$a * log(d$k) -  (d$a + 1) * log(x))
  exp(lnd)
}

#' @rdname pdf.Pareto
#' @export
#'
log_pdf.Pareto <- function(d, x, ...) {
  # Mimic the vectorised behaviour of R's base random generation functions
  if (length(x) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(x), length(d$k), length(d$a))
  x <- rep_len(x, max_len)
  d$k <- rep_len(d$k, max_len)
  d$a <- rep_len(d$a, max_len)
  ifelse(x < d$k, -Inf, log(d$a) + d$a * log(d$k) -  (d$a + 1) * log(x))
}

#' Evaluate the cumulative distribution function of a Pareto distribution
#'
#' @inherit Pareto examples
#' @inheritParams random.Pareto
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Pareto <- function(d, x, ...) {
  # Mimic the vectorised behaviour of R's base random generation functions
  if (length(x) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(x), length(d$k), length(d$a))
  x <- rep_len(x, max_len)
  d$k <- rep_len(d$k, max_len)
  d$a <- rep_len(d$a, max_len)
  ln_survivor <- ifelse(x < d$k, 0, d$a * (log(d$k) - log(x)))
  1 - exp(ln_survivor)
}

#' Determine quantiles of a Pareto distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Pareto examples
#' @inheritParams random.Pareto
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Pareto <- function(d, p, ...) {
  if (length(p) == 0) {
    return(numeric(0))
  }
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("invalid p: p must be in [0, 1].")
  }
  # Mimic the vectorised behaviour of R's base random generation functions
  max_len <- max(length(p), length(d$k), length(d$a))
  p <- rep_len(p, max_len)
  d$k <- rep_len(d$k, max_len)
  d$a <- rep_len(d$a, max_len)
  lnq <- log(d$k) - log(1 - p) / d$a
  exp(lnq)
}

#' Fit a Pareto distribution to data
#'
#' @param d A `Pareto` object created by a call to [Pareto()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Pareto distribution
#'
#' @return A `Pareto` object.
#' @export
fit_mle.Pareto <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  khat <- ss$minimum
  ahat <- 1 / (ss$mean_log - log(khat))
  Pareto(khat, ahat)
}

#' Compute the sufficient statistics for a Pareto distribution from data
#'
#' @inheritParams fit_mle.Pareto
#'
#' @return A named list of the sufficient statistics of the Pareto
#'   distribution:
#'
#'   - `minimum`: The sample minimum.
#'   - `mean_log`: The sample mean of the natural logs of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
suff_stat.Pareto <- function(d, x, ...) {
  valid_x <- is.numeric(x)
  if (!valid_x) stop("`x` must be a numeric vector")
  list(minimum = min(x), mean_log = mean(log(x)), samples = length(x))
}
