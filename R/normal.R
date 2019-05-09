#' Create a normal distribution
#'
#' The normal distribution is ubiquituous in statistics, partially because
#' of the central limit theorem, which states that sums of i.i.d. random
#' variables eventually become normal. Linear transformations of normal
#' random variables result in new random variables that are also normal. If
#' you are taking an intro stats course, you'll likely use the normal
#' distribution for Z-tests and in simple linear regression. The normal
#' distribution is sometimes called the gaussian distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks,
#'   which is also the mean of the distribution. Can be any real number.
#'   Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks,
#'   which is also the **standard deviation** of the distribution. Can be any
#'   positive number. Defaults to `1`. If you would like a normal
#'   distribution with **variance** \eqn{\sigma^2}, be sure to take the
#'   square root, as this is a common source of errors.
#'
#' @return A `normal` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://alexpghayes.github.io/distributions>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a normal random variable with mean
#'   `mu` = \eqn{\mu} and standard deviation `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{\mathbb{R}}{R}, the set of all real numbers
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{\sigma^2}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2}
#'   }{
#'     f(x) = 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2))
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   The cumulative distribution function has the form
#'
#'   \deqn{
#'     F(t) = \int_{-\infty}^t \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2} dx
#'   }{
#'     integral_{-\infty}^t 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2)) dx
#'   }
#'
#'   but this integral does not have a closed form solution and must be
#'   approximated numerically. The c.d.f. of a standard normal is sometimes
#'   called the "error function". The notation \eqn{\Phi(t)} also stands
#'   for the c.d.f. of a standard normal evaluated at \eqn{t}. Z-tables
#'   list the value of \eqn{\Phi(t)} for various \eqn{t}.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     \mathbb{E}(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2}
#'   }{
#'     E(e^(tX)) = e^(\mu t + \sigma^2 t^2 / 2)
#'   }
#'
#' @section Calculating p-values for Z-tests:
#'
#'   In this section, we work through an example Z-test, and point out a number
#'   of points where you might get stuck along the way.
#'
#'   Let's suppose that a student is interesting in estimating how many memes
#'   their professors know and love. So they go to class, and every time a
#'   professor uses a new meme, they write it down. After a year of classes,
#'   the student has recorded the following meme counts, where each count
#'   corresponds to a single class they took:
#'
#'   \deqn{3, 7, 11, 0, 7, 0, 4, 5, 6, 2}
#'
#'   The student talks to some other students who've done similar studies
#'   and determines that \eqn{\sigma = 2} is a reasonable value for the
#'   standard deviation of this distribution.
#'
#'   Before we can do a Z-test, we need to make check if we can reasonably
#'   treat the mean of this sample as normally distributed. This happens is
#'   the case of either of following hold:
#'
#'   1. The data comes from a normal distribution.
#'   2. We have lots of data. How much? Many textbooks use 30 data points
#'     as a rule of thumb.
#'
#'   Since we have a small sample, we let's check if the data comes from
#'   a normal distribution using a normal quantile-quantile plot.
#'
#'   ```
#'   # read in the data
#'   x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#'   ```
#'
#'   ```
#'   # make the qqplot
#'   qqnorm(x)
#'   qqline(x)
#'   ```
#'
#'   Since the data lies close the line \eqn{y = x}, and has no notable
#'   systematic deviations from line, it's safe to treat the sample as
#'   coming from a normal distribution. We can proceed with our hypothesis
#'   test.
#'
#'   Let's test the null hypothesis that, on average, professors know 3 memes.
#'   That is
#'
#'   \deqn{
#'     H_0: \mu = 3 \qquad H_A: \mu \neq 3
#'   }{
#'     H_0: \mu = 3   vs   H_A: \mu != 3
#'   }
#'
#'   First we need to calculate our Z-statistic. Let's use do this with R.
#'   Remember that the Z-statistic is defined as
#'
#'   \deqn{
#'     Z = \frac{\bar x - \mu_0}{\sigma / \sqrt{n}}
#'   }{
#'     Z = (x_bar - \mu) / (\sigma / \sqrt n)
#'   }
#'
#'   In R this looks like:
#'
#'   ```
#'   n <- length(x)
#'   ```
#'
#'   ```
#'   # calculate the z-statistic
#'   z_stat <- (mean(x) - 3) / (2 / sqrt(n))
#'   z_stat
#'   #> [1] 2.371708
#'   ```
#'
#'   To calculate a two-sided p-value, we need to find
#'
#'   \deqn{
#'     P(|Z| \ge |2.37|)
#'     = P(Z \ge 2.37) + P(Z \le -2.37)
#'     = 1 - P(Z \le 2.37) + P(Z \le -2.37)
#'     = 1 - \Phi(2.37) + \Phi(-2.37)
#'   }
#'
#'   To do this we need to c.d.f. of a standard normal
#'
#'   ```
#'   Z <- normal(0, 1)  # make a standard normal r.v.
#'   1 - cdf(Z, 2.37) + cdf(Z, -2.37)
#'   #> [1] 0.01778809
#'   ```
#'
#'   Note that we saved `z_stat` above so we could have also done
#'
#'   ```
#'   1 - cdf(Z, z_stat) + cdf(Z, -z_stat)
#'   #> [1] 0.01770607
#'   ```
#'
#'   which is slightly more accurate since there is no rounding error.
#'
#'   So our p-value is about 0.0177. You should verify this with a Z-table.
#'   Note that you should get the *same* value from `cdf(Z, 2.37)` and looking
#'   up `2.37` on a Z-table.
#'
#'   You may also have seen a different formula for the p-value of a two-sided
#'   Z-test, which makes use of the fact that the normal distribution is
#'   symmetric:
#'
#'   \deqn{
#'     P(|Z| \ge |2.37|)
#'     = 2 \cdot P(Z \le -|2.37|)
#'     = 2 \cdot \Phi(-2.37)
#'   }{
#'     P(|Z| \ge |2.37|)
#'     = 2 P(Z \le -|2.37|)
#'     = 2 \Phi(-2.37)
#'   }
#'
#'   Using this formula we get the same result:
#'
#'   ```
#'   2 * cdf(Z, -2.37)
#'   #> [1] 0.01778809
#'   ```
#'
#'   Finally, sometimes we are interest in one sided Z-tests. For the test
#'
#'   \deqn{
#'     H_0: \mu \le 3 \qquad H_A: \mu > 3
#'   }{
#'     H_0: \mu \le 3   vs   H_A: \mu > 3
#'   }
#'
#'   the p-value is given by
#'
#'   \deqn{P(Z > 2.37)}
#'
#'   which we calculate with
#'
#'   ```
#'   Z <- normal(0, 1)
#'   1 - cdf(Z, 2.37)
#'   #> 0.008894043
#'   ```
#'
#'   For the test
#'
#'   \deqn{
#'     H_0: \mu \ge 3 \qquad H_A: \mu < 3
#'   }{
#'     H_0: \mu \ge 3   vs   H_A: \mu < 3
#'   }
#'
#'   the p-value is given by
#'
#'   \deqn{P(Z < 2.37)}
#'
#'   which we calculate with
#'
#'   ```
#'   Z <- normal(0, 1)
#'   cdf(Z, 2.37)
#'   #> [1] 0.991106
#'   ```
#'
#' @section Confidence interval for a mean example:
#'
#'   The normal distribution also comes up frequently when calculating
#'   confidence intervals for sample means. Let's calculate a 88 percent
#'   confidence interval for the data used in the Z-test example:
#'
#'   \deqn{3, 7, 11, 0, 7, 0, 4, 5, 6, 2}
#'
#'   Again we will assume that \eqn{\sigma = 2}. Recall that a confidence
#'   interval for the mean based off the normal distribution is valid when:
#'
#'   1. The data comes from a normal distribution.
#'   2. We have lots of data. How much? Many textbooks use 30 data points
#'     as a rule of thumb.
#'
#'   In the Z-test example we verified that the sample seems to come from
#'   a normal distribution using a quantile-quantile plot (QQ-plot).
#'
#'   **Heads up**: It's really hard to read the following in RStudio's
#'   Help pane. Do yourself a favor and read the documentation online
#'   at <https://alexpghayes.github.io/distributions>.
#'
#'   The formula for a confidence interval with confidence coefficient
#'   \eqn{1 - \alpha} (in our case this is 0.88) is then:
#'
#'   \deqn{
#'     \left( \bar x + z_{\alpha / 2} \cdot \frac{\sigma}{\sqrt{n}},
#'     \bar x + z_{1 - \alpha / 2} \cdot \frac{\sigma}{\sqrt{n}} \right)
#'   }{
#'     (x_bar + z_{\alpha / 2} * \sigma / \sqrt n,
#'     x_bar + z_{1 - \alpha / 2} * \sigma / \sqrt n)
#'   }
#'
#'   Where \eqn{z_\alpha} stands for the alpha-ith quantile of a standard
#'   normal distribution. To get a quantile we need to take the inverse of
#'   the c.d.f. so you may also see \eqn{z_\alpha} written as
#'   \eqn{\Phi^{-1}(\alpha)}, where \eqn{\Phi} is the c.d.f. of the
#'   standard normal. Since the standard normal distribution is
#'   symmetric around zero, this is exactly equivalent to
#'
#'   \deqn{
#'     \left( \bar x - z_{1 - \alpha / 2} \cdot \frac{\sigma}{\sqrt{n}},
#'     \bar x + z_{1 - \alpha / 2} \cdot \frac{\sigma}{\sqrt{n}} \right)
#'   }{
#'     (x_bar - z_{1 - \alpha / 2} * \sigma / \sqrt n,
#'     x_bar + z_{1 - \alpha / 2} * \sigma / \sqrt n)
#'   }
#'
#'   which may look slightly more familiar. Let's go ahead and calculate
#'   this out in R. Since our confidence coefficient is 0.88 (corresponding
#'   to an 88 percent confidence interval) we have:
#'
#'   \deqn{0.88 = 1 - \alpha}
#'
#'   so that \eqn{\alpha = 0.12}. Now we can get started.
#'
#'   ```
#'   # read in the data
#'   x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#'   n <- length(x)
#'   ```
#'
#'   ```
#'   # make a standard normal r.v.
#'   Z <- normal(0, 1)
#'   ```
#'
#'   ```
#'   # first approach
#'   mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(n)
#'   > [1] 3.516675
#'   mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(n)
#'   > [1] 5.483325
#'   ```
#'
#'   So our confidence interval using the first set of equations is
#'   (3.52, 5.48). Now we use the second set of equations:
#'
#'   ```
#'   # second approach
#'   mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(n)
#'   > [1] 3.516675
#'   mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(n)
#'   > [1] 5.483325
#'   ```
#'
#'   We get the same thing! Just like we expected.
#'
#'   There's one last thing we need to address. You may not have seen either
#'   of the formulas for a Z-confidence interval that I wrote. You may have
#'   seen the formulas:
#'
#'   \deqn{
#'     \left( \bar x - z_{\alpha / 2} \cdot \frac{\sigma}{\sqrt{n}},
#'     \bar x + z_{\alpha / 2} \cdot \frac{\sigma}{\sqrt{n}} \right)
#'   }{
#'     (x_bar - z_{\alpha / 2} * \sigma / \sqrt n,
#'     x_bar + z_{\alpha / 2} * \sigma / \sqrt n)
#'   }
#'
#'   If this is you, you have my condolences, and your instructor
#'   probably hates you. Ask that they use more standard notation.
#'
#'   This looks almost like the second approach, except using
#'   \eqn{z_{\alpha / 2}} instead of \eqn{z_{1 - \alpha / 2}}. What this
#'   comes down to is whether or not \eqn{z_{\alpha / 2}} represents
#'   a *lower quantile* or an *upper quantile*. For a lower quantile,
#'   you look at the p.d.f. and start integrating from negative infinity,
#'   stop when the integral equals \eqn{\alpha}, and that take value
#'   to be the quantile. This is the only sane way to do things, although
#'   it requires being slightly more verbose so it can be inconvenient
#'   at times. Thus the upper quantile, in which case you do the same
#'   integration but start from positive infinity. See the examples
#'   section of [quantile.normal()] for a visual demonstration of the
#'   difference.
#'
#'   If you are truly unfortunate, your instructor may use
#'   \eqn{z_{\alpha / 2}} to mean lower tail quantiles at times
#'   and upper tail quantiles at other times. If this is the case,
#'   only god can help you.
#'
#' @examples
#'
#' n <- normal(5, 2)
#' n
#'
#' random(n, 10)
#' pdf(n, 2)
#' cdf(n, 4)
#' quantile(n, 0.7)
#'
#' ### example: calculating p-values for two-sided Z-test
#'
#' # here the null hypothesis is H_0: mu = 3
#' # and we assume sigma = 2
#'
#' # exactly the same as: Z <- normal(0, 1)
#' Z <- normal()
#'
#' # data to test
#' x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' nx <- length(x)
#'
#' # calculate the z-statistic
#' z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
#' z_stat
#'
#' # calculate the two-sided p-value
#' 1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#'
#' # exactly equivalent to the above
#' 2 * cdf(Z, -abs(z_stat))
#'
#' # p-value for one-sided test
#' # H_0: mu <= 3   vs   H_A: mu > 3
#' 1 - cdf(Z, z_stat)
#'
#' # p-value for one-sided test
#' # H_0: mu >= 3   vs   H_A: mu < 3
#' cdf(Z, z_stat)
#'
#' ### example: calculating a 88 percent Z CI for a mean
#'
#' # same `x` as before, still assume `sigma = 2`
#'
#' # lower-bound
#' mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # upper-bound
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # equivalent to
#' mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # also equivalent to
#' mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ### generating random samples and plugging in ks.test()
#'
#' # generate a random sample
#' normal_sample <- random(normal(3, 7), 26)
#'
#' # test if sample is normal(3, 7)
#' ks.test(normal_sample, pnorm, mean = 3, sd = 7)
#'
#' # test if sample is gamma(8, 3) using base R pgamma()
#' ks.test(normal_sample, pgamma, shape = 8, rate = 3)
#'
#'
#' ### visualizing the cdf and quantiles
#'
#' # note that the cdf() and quantile() functions are inverses
#' cdf(n, quantile(n, 0.7))
#' quantile(n, cdf(n, 7))
#'
#' library(ggplot2)
#'
#' grid <- seq(-4, 4, length.out = 300)
#' density <- pdf(Z, grid)
#' lower_trunc <- ifelse(grid <= 1.96, density, 0)
#' upper_trunc <- ifelse(grid >= 1.96, density, 0)
#'
#' ggplot(data = NULL) +
#'   geom_area(aes(grid, lower_trunc, alpha = 0.2), fill = "steelblue") +
#'   geom_line(aes(grid, density), size = 1, color = "grey") +
#'   geom_vline(xintercept = 1.96, size = 1, color = "darkgrey") +
#'   geom_text(
#'     aes(x = 3, y = 0.3, label = "z[0.95] == 1.96"),
#'     parse = TRUE,
#'     size = 4
#'   ) +
#'   labs(
#'     title = "Lower tail quantile of a standard normal",
#'     subtitle = "Integral of shaded region is 0.95",
#'     y = "Density",
#'     x = "Support"
#'   ) +
#'   theme_minimal() +
#'   theme(legend.position = "none")
#'
#' ggplot(data = NULL) +
#'   geom_area(aes(grid, upper_trunc, alpha = 0.2), fill = "steelblue") +
#'   geom_line(aes(grid, density), size = 1, color = "grey") +
#'   geom_vline(xintercept = 1.96, size = 1, color = "darkgrey") +
#'   geom_text(
#'     aes(x = 3, y = 0.3, label = "z[0.05] == 1.96"),
#'     parse = TRUE,
#'     size = 4
#'   ) +
#'   labs(
#'     title = "Upper tail quantile of a standard normal",
#'     subtitle = "Integral of shaded region is 0.05",
#'     y = "Density",
#'     x = "Support"
#'   ) +
#'   theme_minimal() +
#'   theme(legend.position = "none")
#'
normal <- function(mu = 0, sigma = 1) {
  d <- list(mu = mu, sigma = sigma)
  class(d) <- "normal"
  d
}

#' @export
print.normal <- function(d, ...) {
  cat(glue("normal distribution (mu = {d$mu}, sigma = {d$sigma})"))
}

#' Draw a random sample from a normal distribution
#'
#' Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit normal examples
#'
#' @param d A `normal` object created by a call to [normal()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family normal distribution
#'
#' @return A numeric vector of length `n`.
#' @export
#'
#'
random.normal <- function(d, n = 1L, ...) {
  rnorm(n = n, mean = d$mu, sd = d$sigma)
}

#' Evaluate the probability mass function of a normal distribution
#'
#' Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.normal <- function(d, x, ...) {
  dnorm(x = x, mean = d$mu, sd = d$sigma)
}

#' Evaluate the cumulative distribution function of a normal distribution
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family normal distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.normal <- function(d, x, ...) {
  pnorm(q = x, mean = d$mu, sd = d$sigma)
}

#' Determine quantiles of a normal distribution
#'
#' Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [normal()] for some properties
#' of the normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit normal examples
#' @inheritParams random.normal
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
#' @family normal distribution
#'
quantile.normal <- function(d, p, ...) {
  qnorm(p = p, mean = d$mu, sd = d$sigma)
}
