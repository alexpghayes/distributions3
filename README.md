
<!-- README.md is generated from README.Rmd. Please edit that file -->

# distributions

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/alexpghayes/distributions.svg?branch=master)](https://travis-ci.org/alexpghayes/distributions)
<!-- badges: end -->

`distributions`, inspired by the [eponynmous Julia
package](https://github.com/JuliaStats/Distributions.jl), provides a
generic function interface to probability distributions.

`distributions` has two goals:

1.  Replace the `rnorm()`, `pnorm()`, etc, family of functions with S3
    methods for distribution objects

2.  Be extremely well documented and friendly toward students in intro
    stat classes.

The main generics are:

  - `random()`: Draw samples from a distribution.
  - `pdf()`: Evaluate the probability density of mass at a point `x`.
  - `cdf()`: Evaluate the cumulative probability up to a point `x`.
  - `quantile()`: Determine the quantile at `p`. Inverse of `cdf()`.

## Installation

`distributions` is not yet on CRAN. You can install the development
version with:

``` r
install.packages("devtools")
devtools::install_github("alexpghayes/distributions")
```

## Basic Usage

The basic usage of `distributions` looks like:

``` r
library(distributions)

b <- bernoulli(0.1)

random(b, 10)
#>  [1] 1 0 0 0 0 0 0 0 0 0
pdf(b, 1)
#> [1] 0.1
cdf(b, 0)
#> [1] 0.9
quantile(b, 0.5)
#> [1] 0
```

Note that `quantile()` **always** returns lower tail probabilities. If
you aren’t sure what this means, please read the last several paragraphs
of `vignette("one-sample-z-confidence-interval")` and have a gander at
the plot.

## Development roadmap

The initial goal of `distributions` is to replace the base R `d`, `p`,
`q`, and `r` functions for the following families:

  - [x] `beta` - documented
  - [x] `binom` - documented
  - [x] `cauchy` - documented
  - [x] `chisq`
  - [x] `exp`
  - [x] `f`
  - [x] `gamma`
  - [x] `norm` - documented
  - [x] `pois` - documented
  - [x] `t`
  - [ ] `multinomial`
  - [ ] `unif`
  - [ ] `geom`
  - [ ] `hyper`
  - [ ] `logis`
  - [ ] `lnorm`
  - [ ] `nbinom`
  - [ ] `weibull`
  - [ ] `tukey`

We also want to write detailed vignettes describing the most common use
cases for each of these distributions. Realistically, this means we need
to document how to perform the statistical tests most common in intro
stats classes.

This documentation should live in vignettes that explicitly write down
both math and code in a cookbook style mashup. These cookbooks should
also provide code and guidance for basic assumption checking. See
`vignette("one-sample-z-test")` for an example. The vignette wishlist is
currently as follows:

**One sample tests**

  - [ ] Update the one sample z-test vignette to include rejection
    regions and power calculations
  - [ ] Update the one sample z-test for a proportion vignette to
    include rejection regions and power calculations
      - Include note about rejection region not being inverted CI due to
        different standard errors
  - [ ] Update the one sample t-test vignette to include rejection
    regions and power calculations
  - [ ] Signs test (binomial tests) using the binomial distribution
  - [ ] One sample z confidence interval for a proportion
  - [ ] Chi squared tests for counts
      - Would be really great to get help with this one

**Two sample tests**

  - [ ] Two sample T-test (independent / equal variance)
      - Comment on \(F = T^2\) relationship
  - [ ] Two sample T-test (Welch’s)
      - Show use of `power.t.test()`

**Paired tests**

  - [ ] Paired Z-test
  - [ ] Paired T-test
  - [ ] Paired sign-test

**ANOVA and linear regression**

  - [ ] One-way ANOVA using the F distribution
      - Confirm results with both `aov()` and `anova(lm())`
      - Show power calculations with `power.anova.test()`
      - Diagnostic plots / assumption checking

**Fancier tests (asymptotic versions)**

  - [ ] Likelihood ratio test
  - [ ] Wald test
  - [ ] Rao score test

After that the goal will probably to support more generics like those in
the `Distributions.jl`, with the hope that the community will contribute
additional distributions.

Candidates for the next set of generics:

  - `log_pdf()`
  - `log_cdf()`
  - `plot_pdf()`
  - `plot_cdf()`

## Contributing

The best way to get started is to look at `R/beta.R` and
`tests/testthat/test-beta.R`, copy them, and modify them for whatever
new distribution you’d like to add.

Please note that `distributions` is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.

## Related work

For a comprehensive overview of the many packages providing various
distribution related functionality see the [CRAN Task
View](https://cran.r-project.org/web/views/Distributions.html).

  - [`distr`](http://distr.r-forge.r-project.org/) is quite similar to
    `distributions`, but uses S4 objects and is less focused on
    documentation
  - [`fitdistrplus`](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    provides extensive functionality for fitting various distributions
    but does not treat distributions themselves as objects
