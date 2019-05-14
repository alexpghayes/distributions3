% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{random}
\alias{random}
\title{Draw a random sample from a probability distribution}
\usage{
random(d, n = 1L, ...)
}
\arguments{
\item{d}{A probability distribution object such as those created by
a call to \code{\link[=bernoulli]{bernoulli()}}, \code{\link[=beta]{beta()}}, or \code{\link[=binomial]{binomial()}}.}

\item{n}{The number of samples to draw. Should be a positive
integer. Defaults to \code{1L}.}

\item{...}{Unused. Unevaluated arguments will generate a warning to
catch mispellings or other possible errors.}
}
\description{
Draw a random sample from a probability distribution
}