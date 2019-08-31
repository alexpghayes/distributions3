
<!-- README.md is generated from README.Rmd. Please edit that file -->

# distributions

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/alexpghayes/distributions.svg?branch=master)](https://travis-ci.org/alexpghayes/distributions)
[![Codecov test
coverage](https://codecov.io/gh/alexpghayes/distributions/branch/master/graph/badge.svg)](https://codecov.io/gh/alexpghayes/distributions?branch=master)
<!-- badges: end -->

`distributions3`, inspired by the [eponynmous Julia
package](https://github.com/JuliaStats/Distributions.jl), provides a
generic function interface to probability distributions. `distributions`
has two goals:

1.  Replace the `rnorm()`, `pnorm()`, etc, family of functions with S3
    methods for distribution objects

2.  Be extremely well documented and friendly for students in intro stat
    classes.

The main generics are:

  - `random()`: Draw samples from a distribution.
  - `pdf()`: Evaluate the probability density of mass at a point.
  - `cdf()`: Evaluate the cumulative probability up to a point.
  - `quantile()`: Determine the quantile for a given probability.
    Inverse of `cdf()`.

## Installation

`distributions` is not yet on CRAN. You can install the development
version with:

``` r
install.packages("devtools")
devtools::install_github("alexpghayes/distributions")
```

## Basic Usage

The basic usage of `distributions3` looks like:

``` r
library(distributions3)

X <- Bernoulli(0.1)

random(X, 10)
#>  [1] 0 0 0 0 0 0 0 0 0 0
pdf(X, 1)
#> [1] 0.1

cdf(X, 0)
#> [1] 0.9
quantile(X, 0.5)
#> [1] 0
```

Note that `quantile()` **always** returns lower tail probabilities. If
you aren’t sure what this means, please read the last several paragraphs
of `vignette("one-sample-z-confidence-interval")` and have a gander at
the plot.

## Contributing

I am very happy to review PRs and provide advice on how to add new
functionality to the package. Documentation improvements are
particularly appreciated\!

To add a new distribution, the best way to get started is to look at
`R/Beta.R` and `tests/testthat/test-Beta.R`, copy them, and modify them
for whatever new distribution you’d like to add.

Please note that `distributions` is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.

## Related work

For a comprehensive overview of the many packages providing various
distribution related functionality see the [CRAN Task
View](https://cran.r-project.org/web/views/Distributions.html).

  - [`distr`](http://distr.r-forge.r-project.org/) is quite similar to
    `distributions`, but uses S4 objects and is less focused on
    documentation.
  - [`distr6`](https://alan-turing-institute.github.io/distr6/) builds
    on `distr`, but uses R6 objects
  - [`fitdistrplus`](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    provides extensive functionality for fitting various distributions
    but does not treat distributions themselves as objects
