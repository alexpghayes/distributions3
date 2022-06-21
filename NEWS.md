# distributions3 0.2.0

- Vectorized univariate distribution objects by Moritz Lang and Achim Zeileis (#71 and #82).
  This allows representation of fitted probability distributions from regression models.
  New helper functions are provided to help setting up such distribution objects in
  a unified way. In particular, `apply_dpqr()` helps to apply the standard d/p/q/r functions
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
