% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Binomial.R
\name{random.Binomial}
\alias{random.Binomial}
\title{Draw a random sample from a Binomial distribution}
\usage{
\method{random}{Binomial}(x, n = 1L, drop = TRUE, ...)
}
\arguments{
\item{x}{A \code{Binomial} object created by a call to \code{\link[=Binomial]{Binomial()}}.}

\item{n}{The number of samples to draw. Defaults to \code{1L}.}

\item{drop}{logical. Should the result be simplified to a vector if possible?}

\item{...}{Unused. Unevaluated arguments will generate a warning to
catch mispellings or other possible errors.}
}
\value{
Integers containing values between \code{0} and \code{x$size}.
In case of a single distribution object or \code{n = 1}, either a numeric
vector of length \code{n} (if \code{drop = TRUE}, default) or a \code{matrix} with \code{n} columns
(if \code{drop = FALSE}).
}
\description{
Draw a random sample from a Binomial distribution
}
\examples{

set.seed(27)

X <- Binomial(10, 0.2)
X

mean(X)
variance(X)
skewness(X)
kurtosis(X)

random(X, 10)

pdf(X, 2L)
log_pdf(X, 2L)

cdf(X, 4L)
quantile(X, 0.7)

cdf(X, quantile(X, 0.7))
quantile(X, cdf(X, 7))
}
