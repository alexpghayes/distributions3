random <- function(d, ...) {
  UseMethod("random")
}

# radon-nikodym density
pdf <- function(d, ...) {
  UseMethod("pdf")
}

cdf <- function(d, ...) {
  UseMethod("cdf")
}

quantile <- function(x, ...) {
  UseMethod("quantile")
}
