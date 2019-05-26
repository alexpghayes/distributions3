#' Create a binomial distribution
#'
#' Bernoulli distributions are used to represent situations can that can
#' thought of as `size` (often called \eqn{n} in textbooks) independent
#' coin flips, where each coin flip has probability `p` of success. The
#' [bernoulli()] distribution is a special case of binomial distribution
#' when `n = 1`.
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the binomial distribution reduces to the
#'   bernoulli distribution. Oftened called `n` in textbooks.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `binomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   The binomial distribution comes up when you are interested in the portion
#'   of people who do a thing. The binomial distribution
#'   also comes up in the sign test, sometimes called the binomial test
#'   (see [stats::binom.test()]), where you may need the binomial C.D.F. to
#'   compute p-values.
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
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
#'     P(X = k) = \binom{n}{k} p^k (1 - p)^{n-k}
#'   }{
#'     P(X = k) = choose(n, k) p^k (1 - p)^(n - k)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) = \sum_{i=0}^{\lfloor k \rfloor} \binom{n}{i} p^i (1 - p)^{n-i}
#'   }{
#'     P(X \le k) = \sum_{i=0}^k choose(n, i) p^i (1 - p)^(n-i)
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     \mathbb{E}(e^{tX}) = (1 - p + p e^t)^n
#'   }{
#'     E(e^(tX)) = (1 - p + p e^t)^n
#'   }
#'
#' @examples
#'
#' b <- binomial(10, 0.2)
#' b
#'
#' random(b, 10)
#' pdf(b, 2L)
#' log_pdf(b, 2L)
#' cdf(b, 4L)
#' quantile(b, 0.7)
#'
#' cdf(b, quantile(b, 0.7))
#' quantile(b, cdf(b, 7))
#'
binomial <- function(size, p = 0.5) {
  d <- list(size = size, p = p)
  class(d) <- "binomial"
  d
}

#' @export
print.binomial <- function(x, ...) {
  cat(glue("Binomial distribution (size = {x$size}, p = {x$p})"))
}

#' Draw a random sample from a binomial distribution
#'
#' @inherit binomial examples
#'
#' @param d A `binomial` object created by a call to [binomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return An integer vector containing values between `0` and `d$size`
#'   of length `n`.
#' @export
#'
random.binomial <- function(d, n = 1L, ...) {
  rbinom(n = n, size = d$size, prob = d$p)
}

#' Evaluate the probability mass function of a binomial distribution
#'
#' @inherit binomial examples
#' @inheritParams random.binomial
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.binomial <- function(d, x, ...) {
  dbinom(x = x, size = d$size, prob = d$p)
}

#' @rdname pdf.binomial
#' @export
#'
log_pdf.binomial <- function(d, x, ...) {
  dbinom(x = x, size = d$size, prob = d$p, log = TRUE)
}

#' Evaluate the cumulative distribution function of a binomial distribution
#'
#' @inherit binomial examples
#' @inheritParams random.binomial
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.binomial <- function(d, x, ...) {
  pbinom(q = x, size = d$size, prob = d$p)
}

#' Determine quantiles of a binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit binomial examples
#' @inheritParams random.binomial
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.binomial <- function(d, p, ...) {

  # TODO: in the documentation, more information on return and
  # how quantiles are calculated

  qbinom(p = p, size = d$size, prob = d$p)
}

#' Fit a binomial distribution to data
#'
#' The fit distribution will inherit the same `size` parameter as
#' the `binomial` object passed.
#'
#' @param d A `binomial` object.
#' @param x A vector of zeroes and ones to fit the binomial distribution to.
#'
#' @return a `binomial` object
#' @export
fit_mle.binomial <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  binomial(ss$trials, ss$successes / (ss$experiments * ss$trials))
}

#' Compute the sufficient statistics for the binomial distribution from data
#'
#' @inherit binomial
#' @export
suff_stat.binomial <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x <= d$size) & (x %% 1 == 0)
  if(any(!valid_x)) {
    stop("`x` must be an integer between zero and the size parameter of the binomial distribution")
  }
  list(successes = sum(x), experiments = length(x), trials = d$size)
}

#' Compute the likelihood of a binomial distribution given data
#'
#' @inheritParams fit_mle.binomial
#'
#' @return the likelihood
#' @export

likelihood.binomial <- function(d, x, ...) {
  prod(dbinom(x = x, size = d$size, prob = d$p))
}

#' Compute the log-likelihood of a binomial distribution given data
#'
#' @inheritParams fit_mle.binomial
#'
#' @return the log-likelihood
#' @export
log_likelihood.binomial <- function(d, x, ...) {
  sum(dbinom(x = x, size = d$size, prob = d$p, log = TRUE))
}
