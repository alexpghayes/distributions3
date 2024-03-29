% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gumbel.R
\name{random.Gumbel}
\alias{random.Gumbel}
\title{Draw a random sample from a Gumbel distribution}
\usage{
\method{random}{Gumbel}(x, n = 1L, drop = TRUE, ...)
}
\arguments{
\item{x}{A \code{Gumbel} object created by a call to \code{\link[=Gumbel]{Gumbel()}}.}

\item{n}{The number of samples to draw. Defaults to \code{1L}.}

\item{drop}{logical. Should the result be simplified to a vector if possible?}

\item{...}{Unused. Unevaluated arguments will generate a warning to
catch mispellings or other possible errors.}
}
\value{
In case of a single distribution object or \code{n = 1}, either a numeric
vector of length \code{n} (if \code{drop = TRUE}, default) or a \code{matrix} with \code{n} columns
(if \code{drop = FALSE}).
}
\description{
Draw a random sample from a Gumbel distribution
}
\examples{

set.seed(27)

X <- Gumbel(1, 2)
X

random(X, 10)

pdf(X, 0.7)
log_pdf(X, 0.7)

cdf(X, 0.7)
quantile(X, 0.7)

cdf(X, quantile(X, 0.7))
quantile(X, cdf(X, 0.7))
}
