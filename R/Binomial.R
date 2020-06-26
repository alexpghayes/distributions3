#' Create a Binomial distribution
#'
#' Binomial distributions are used to represent situations can that can
#' be thought as the result of \eqn{n} Bernoulli experiments (here the
#' \eqn{n} is defined as the `size` of the experiment). The classical
#' example is \eqn{n} independent coin flips, where each coin flip has
#' probability `p` of success. In this case, the individual probability of
#' flipping heads or tails is given by the  Bernoulli(p) distribution,
#' and the probability of having \eqn{x} equal results (\eqn{x} heads,
#' for example), in \eqn{n} trials is given by the Binomial(n, p) distribution.
#' The equation of the Binomial distribution is directly derived from
#' the equation of the Bernoulli distribution.
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the Binomial distribution reduces to the
#'   bernoulli distribution. Often called `n` in textbooks.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `Binomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   The Binomial distribution comes up when you are interested in the portion
#'   of people who do a thing. The Binomial distribution
#'   also comes up in the sign test, sometimes called the Binomial test
#'   (see [stats::binom.test()]), where you may need the Binomial C.D.F. to
#'   compute p-values.
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions3>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Binomial random variable with parameter
#'   `size` = \eqn{n} and `p` = \eqn{p}. Some textbooks define \eqn{q = 1 - p},
#'   or called \eqn{\pi} instead of \eqn{p}.
#'
#'   **Support**: \eqn{\{0, 1, 2, ..., n\}}{{0, 1, 2, ..., n}}
#'
#'   **Mean**: \eqn{np}
#'
#'   **Variance**: \eqn{np \cdot (1 - p) = np \cdot q}{np (1 - p)}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = k) = {n \choose k} p^k (1 - p)^{n-k}
#'   }{
#'     P(X = k) = choose(n, k) p^k (1 - p)^(n - k)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) = \sum_{i=0}^{\lfloor k \rfloor} {n \choose i} p^i (1 - p)^{n-i}
#'   }{
#'     P(X \le k) = \sum_{i=0}^k choose(n, i) p^i (1 - p)^(n-i)
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = (1 - p + p e^t)^n
#'   }{
#'     E(e^(tX)) = (1 - p + p e^t)^n
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Binomial(10, 0.2)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' random(X, 10)
#'
#' pdf(X, 2L)
#' log_pdf(X, 2L)
#'
#' cdf(X, 4L)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Binomial <- function(size, p = 0.5) {
  d <- list(size = size, p = p)
  class(d) <- c("Binomial", "distribution")
  d
}

#' @export
print.Binomial <- function(x, ...) {
  cat(glue("Binomial distribution (size = {x$size}, p = {x$p})"), "\n")
}

#' @export
mean.Binomial <- function(d, ...) d$size * d$p

#' @export
variance.Binomial <- function(d, ...) d$size * d$p * (1 - d$p)

#' @export
skewness.Binomial <- function(d, ...) {
  n <- d$size
  p <- d$p
  q <- 1 - d$p
  (1 - (2 * p)) / sqrt(n * p * q)
}

#' @export
kurtosis.Binomial <- function(d, ...) {
  n <- d$size
  p <- d$p
  q <- 1 - d$p
  (1 - (6 * p * q)) / (n * p * q)
}

#' Draw a random sample from a Binomial distribution
#'
#' @inherit Binomial examples
#'
#' @param d A `Binomial` object created by a call to [Binomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return An integer vector containing values between `0` and `d$size`
#'   of length `n`.
#' @export
#'
random.Binomial <- function(d, n = 1L, ...) {
  rbinom(n = n, size = d$size, prob = d$p)
}

#' Evaluate the probability mass function of a Binomial distribution
#'
#' @inherit Binomial examples
#' @inheritParams random.Binomial
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Binomial <- function(d, x, ...) {
  dbinom(x = x, size = d$size, prob = d$p)
}

#' @rdname pdf.Binomial
#' @export
log_pdf.Binomial <- function(d, x, ...) {
  dbinom(x = x, size = d$size, prob = d$p, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Binomial distribution
#'
#' @inherit Binomial examples
#' @inheritParams random.Binomial
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Binomial <- function(d, x, ...) {
  pbinom(q = x, size = d$size, prob = d$p)
}

#' Determine quantiles of a Binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Binomial examples
#' @inheritParams random.Binomial
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Binomial <- function(d, p, ...) {
  qbinom(p = p, size = d$size, prob = d$p)
}

#' Fit a Binomial distribution to data
#'
#' The fit distribution will inherit the same `size` parameter as
#' the `Binomial` object passed.
#'
#' @param d A `Binomial` object.
#' @param x A vector of zeroes and ones.
#' @param ... Unused.
#'
#' @return a `Binomial` object
#' @export
fit_mle.Binomial <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Binomial(ss$trials, ss$successes / (ss$experiments * ss$trials))
}

#' Compute the sufficient statistics for the Binomial distribution from data
#'
#' @inheritParams fit_mle.Binomial
#'
#' @return A named list of the sufficient statistics of the Binomial
#'   distribution:
#'
#'   - `successes`: The total number of successful trials.
#'   - `experiments`: The number of experiments run.
#'   - `trials`: The number of trials run per experiment.
#'
#' @export
suff_stat.Binomial <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x <= d$size) & (x %% 1 == 0)
  if (any(!valid_x)) {
    stop("`x` must be an integer between zero and the size parameter of the Binomial distribution")
  }
  list(successes = sum(x), experiments = length(x), trials = d$size)
}


#' Return the support of the Binomial distribution
#'
#' @param d An `Binomial` object created by a call to [Binomial()].
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Binomial <- function(d) c(0, d$size)

