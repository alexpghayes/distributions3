random <- function(x, ...) {
  UseMethod("random")
}

# radon-nikodym density
pdf <- function(x, ...) {
  UseMethod("pdf")
}

cdf <- function(x, ...) {
  UseMethod("cdf")
}

quantile <- function(x, ...) {
  UseMethod("quantile")
}
