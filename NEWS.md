# distributions3 0.2.1

- New generics `is_discrete()` and `is_continous()` with methods for all distribution objects
  in the package. The `is_discrete()` methods return `TRUE` for every distribution that is discrete
  on the entire support and `FALSE` otherwise. Analogously, `is_continuous()` returns `TRUE` for
  every distribution that is continuous on the entire support and `FALSE` otherwise. Thus, for
  mixed discrete-continuous distributions both methods should yield `FALSE` (#90).
- New logical argument `elementwise = NULL` in `apply_dpqr()` and hence inherited in
  `cdf()`, `pdf()`, `log_pdf()`, and `quantile()`. It provides type-safety when
  applying one of the functions to a vector of distributions `d` to a numeric
  argument `x` where both `d` and `x` are of length n > 1. By setting `elementwise = TRUE`
  the function is applied element-by-element, also yielding a vector of length n.
  By setting `elementwise = FALSE` the function is applied for all combinations
  yielding an n-by-n matrix. The default `elementwise = NULL` corresponds to `FALSE`
  if `d` and `x` are of different lengths and `TRUE` if the are of the same length
  n > 1 (#87).
- Extended support for various count data distributions, now enompassing both the Poisson
  and negative binomial distributions along with various adjustments for zero counts
  (hurdle, inflation, and truncation, respectively). More details are provided in the
  following items (#86).
- New `d`/`p`/`q`/`r` functions for `hnbinom`, `zinbinom`, `ztnbinom`, and `ztpois` similar
  to the corresponding `nbinom` and `pois` functions from base R.
- New `HurdleNegativeBinomial()`, `ZINegativeBinomial()`, `ZTNegativeBinomial()`, and
  `ZTPoisson()` distribution constructors along with the corresponding S3 methods for the
  "usual" generics (except `skewness()` and `kurtosis()`).
- New `prodist()` methods for extracting the fitted/predicted probability distributions from
  models estimated by `hurdle()`, `zeroinfl()`, and `zerotrunc()` objects from either the
  `pscl` package or the `countreg` package.
- Added argument `prodist(..., sigma = "ML")` to the `lm` method for extracting the
  fitted/predicted probability distribution from a linear regression model. In the previous
  version the `prodist()` method always used the least-squares estimate of the error variance
  (= residual sum of squares divided by the residual degrees of freedom, n - k), as also
  reported by the `summary()` method. Now the default is to use the maximum-likelihood estimate
  instead (divided by the number of observations, n) which is consistent with the `logLik()`
  method. The previous behavior can be obtained by specifying `sigma = "OLS"` (#91).
- Similarly to the `lm` method the `glm` method `prodist(..., dispersion = NULL)` now, by
  default, uses the `dispersion` estimate that matches the `logLik()` output. This is based
  on the deviance divided by the number of observations, n. Alternatively,
  `dispersion = "Chisquared"` uses the estimate employed in the `summary()` method,
  based on the Chi-squared statistic divided by the residual degrees of freedom, n - k.
- Small improvements in methods for various distribution objects: Added `support()` method
  for GEV-based distributions (`GEV()`, `GP()`, `Gumbel()`, `Frechet()`). Added a
  `random()` method for the `Tukey()` distribution (using the inversion method).


# distributions3 0.2.0

- Vectorized univariate distribution objects by Moritz Lang and Achim Zeileis (#71 and #82).
  This allows representation of fitted probability distributions from regression models.
  New helper functions are provided to help setting up such distribution objects in
  a unified way. In particular, `apply_dpqr()` helps to apply the standard `d`/`p`/`q`/`r` functions
  available in base R and many packages. The accompanying manual page provides some
  worked examples and further guidance.
- New vignette (by Achim Zeileis) on using `distributions3` to go from basic probability
  theory to probabilistic regression models. Illustrated with Poisson GLMs for the
  number of goals per team in the 2018 FIFA World Cup explained by the teams' ability
  differences. (#74)
- New generic function `prodist()` to extract fitted (in-sample) or predicted (out-of-sample)
  probability distributions from model objects like `lm`, `glm`, or `arima`. (#83)
- Extended support for count data distributions (by Achim Zeileis): Alternative
  parameterization for negative binomial distribution (commonly used in regression models),
  zero-inflated Poisson, and zero-hurdle Poisson. (#80 and #81)


# distributions3 0.1.2

- Added a plotting generic for univariate distributions (@paulnorthrop, PR #56)
- Added support for the Generalised Extreme Value (GEV), Frechet, Gumbel, reversed Weibull and Generalised Pareto (GP) distributions (@paulnorthrop, PR #52)
- Added support for the Erlang distribution (@ellessenne, PR #54)
- Various minor bug fixes


# distributions3 0.1.1

- Rename to `distributions3` for CRAN


# distributions 0.1.0

- Added a `NEWS.md` file to track changes to the package.
- Initial release
