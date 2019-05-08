
<!-- README.md is generated from README.Rmd. Please edit that file -->

# distributions

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/alexpghayes/distributions.svg?branch=master)](https://travis-ci.org/alexpghayes/distributions)
<!-- badges: end -->

`distributions` is inspired by the [eponynmous Julia
package](https://github.com/JuliaStats/Distributions.jl), which takes a
generic function approach to probability distributions.

The key idea of `distributions` is to replace the `rnorm()`, `pnorm()`,
etc, family of functions with S3 methods for distribution objects. The
main generics are:

  - `random()`: Draw samples from a distribution.
  - `pdf()`: Evaluate the probability density of mass at a point `x`.
  - `cdf()`: Evaluate the cumulative probability up to a point `x`.
  - `quantile()`: Determine the quantile at `p`. Inverse of `cdf()`.

Another central goal of `distributions` is to use more intuitive and
explicit names, to make it easier for students in intro stats classes to
manipulate probabilities.

## Installation

`distributions` is not yet on CRAN. You can install the development
version with:

``` r
install.packages("devtools")
devtools::install_github("alexpghayes/distributions")
```

## Example

The basic usage of `distributions` looks like:

``` r
library(distributions)

b <- bernoulli(0.1)

random(b, 10)
#>  [1] 0 0 0 0 0 0 0 1 0 0
pdf(b, 1)
#> [1] 0.1
cdf(b, 0)
#> [1] 0.9
quantile(b, 0.5)
#> [1] 0
```

## Development roadmap

The initial goal of `distributions` is to replace the base R `d`, `p`,
`q`, and `r` functions for the following families:

  - [x] `beta`
  - [x] `binom`
  - [x] `cauchy`
  - [x] `chisq`
  - [x] `exp`
  - [x] `f`
  - [x] `gamma`
  - [x] `norm`
  - [x] `pois`
  - [ ] `t`
  - [ ] `multinomial`
  - [ ] `unif`
  - [ ] `geom`
  - [ ] `hyper`
  - [ ] `logis`
  - [ ] `lnorm`
  - [ ] `nbinom`
  - [ ] `weibull`
  - [ ] `tukey`

After that the goal will probably to support more generics like those in
the `Distributions.jl`, with the hope that the community will contribute
additional distributions.

Candidates for the next set of generics:

  - `log_pdf()`
  - `log_cdf()`

## Contributing

The best way to get started is to look at `R/beta.R` and
`tests/testthat/test-beta.R`, copy them, and modify them for whatever
new distribution you’d like to add.

Please note that the ‘distributions’ project is released with a
[Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By
contributing to this project, you agree to abide by its terms.

## Related work

For a comprehensive overview of the many packages providing various
distribution related functionality see the [CRAN Task
View](https://cran.r-project.org/web/views/Distributions.html).

  - [`distr`](http://distr.r-forge.r-project.org/) is quite similar to
    `distributions`, but uses S4
    objects
  - [`fitdistrplus`](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    provides extensive functionality for fitting various distributions
    but does not treat distributions themselves as objects

## Conventions

This section of the README is mostly a note to myself.

  - `quantile()` always uses lower tail probabilities.
