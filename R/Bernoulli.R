#' Create a Bernoulli distribution
#'
#' Bernoulli distributions are used to represent events like coin flips
#' when there is single trial that is either successful or unsuccessful.
#' The Bernoulli distribution is a special case of the [Binomial()]
#' distribution with `n = 1`.
#'
#' @param p The success probability for the distribution. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `Bernoulli` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Bernoulli random variable with parameter
#'   `p` = \eqn{p}. Some textbooks also define \eqn{q = 1 - p}, or use
#'   \eqn{\pi} instead of \eqn{p}.
#'
#'   **Support**: \eqn{\{0, 1\}}{{0, 1}}
#'
#'   **Mean**: \eqn{p}
#'
#'   **Variance**: \eqn{p \cdot (1 - p) = p \cdot q}{p (1 - p)}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = x) = p^x (1 - p)^{1-x} = p^x q^{1-x}
#'   }{
#'     P(X = x) = p^x (1 - p)^(1-x)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le x) =
#'     \left \{
#'       \begin{array}{ll}
#'         0 & x < 0 \\
#'         1 - p & 0 \leq x < 1 \\
#'         1 & x \geq 1
#'       \end{array}
#'     \right.
#'   }{
#'     P(X \le x) = (1 - p) 1_{[0, 1)}(x) + 1_{1}(x)
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = (1 - p) + p e^t
#'   }{
#'     E(e^(tX)) = (1 - p) + p e^t
#'   }
#'
#' @examples
#'
#' X <- Bernoulli(0.7)
#' X
#'
#' random(X, 10)
#' pdf(X, 1)
#' log_pdf(X, 1)
#' cdf(X, 0)
#' quantile(X, 0.7)
#'
#' # TODO: make sure these are inverses of each other
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
#'
Bernoulli <- function(p = 0.5) {

  # TODO: check that 0 <= p <= 1

  d <- list(p = p)
  class(d) <- c("Bernoulli", "distribution")
  d
}

#' @export
print.Bernoulli <- function(x, ...) {
  cat(glue("Bernoulli distribution (p = {x$p})"))
}

#' Draw a random sample from a Bernoulli distribution
#'
#' @inherit Bernoulli examples
#'
#' @param d A `Bernoulli` object created by a call to [Bernoulli()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return An integer vector of zeros and ones of length `n`.
#' @export
#'
random.Bernoulli <- function(d, n = 1L, ...) {
  rbinom(n = n, size = 1, prob = d$p)
}

#' Evaluate the probability mass function of a Bernoulli distribution
#'
#' @inherit Bernoulli examples
#' @inheritParams random.Bernoulli
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Bernoulli <- function(d, x, ...) {
  dbinom(x = x, size = 1, prob = d$p)
}

#' @rdname pdf.Bernoulli
#' @export
#'
log_pdf.Bernoulli <- function(d, x, ...) {
  dbinom(x = x, size = 1, prob = d$p, log = TRUE)
}

#' Evaluate the cumulative distribution function of a Bernoulli distribution
#'
#' @inherit Bernoulli examples
#' @inheritParams random.Bernoulli
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Bernoulli <- function(d, x, ...) {
  pbinom(q = x, size = 1, prob = d$p)
}

#' Determine quantiles of a Bernoulli distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Bernoulli examples
#' @inheritParams random.Bernoulli
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Bernoulli <- function(d, p, ...) {
  qbinom(p = p, size = 1, prob = d$p)
}

#' Fit a Bernoulli distribution to data
#'
#' @param d A `Bernoulli` object.
#' @param x A vector of zeroes and ones.
#' @param ... Unused.
#'
#' @return a `Bernoulli` object
#' @export
fit_mle.Bernoulli <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Bernoulli(p = ss$successes / (ss$successes + ss$failures))
}

#' Compute the sufficient statistics for a Bernoulli distribution from data
#'
#' @inheritParams fit_mle.Bernoulli
#'
#' @return A named list of the sufficient statistics of the Bernoulli
#'   distribution:
#'
#'   - `successes`: The number of successful trials (`sum(x == 1)`)
#'   - `failures`: The number of failed trials (`sum(x == 0)`).
#'
#' @export
suff_stat.Bernoulli <- function(d, x, ...) {
  valid_x <- (x %in% c(0L, 1L))
  if(any(!valid_x)) stop("`x` contains elements other than 0 or 1")
  list(successes = sum(x == 1), failures = sum(x == 0))
}



