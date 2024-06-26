% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HyperGeometric.R
\name{HyperGeometric}
\alias{HyperGeometric}
\title{Create a HyperGeometric distribution}
\usage{
HyperGeometric(m, n, k)
}
\arguments{
\item{m}{The number of type I elements available.}

\item{n}{The number of type II elements available.}

\item{k}{The size of the sample taken.}
}
\value{
A \code{HyperGeometric} object.
}
\description{
To understand the HyperGeometric distribution, consider a set of
\eqn{r} objects, of which \eqn{m} are of the type I and
\eqn{n} are of the type II. A sample with size \eqn{k} (\eqn{k<r})
with no replacement is randomly chosen. The number of observed
type I elements observed in this sample is set to be our random
variable \eqn{X}. For example, consider that in a set of 20
car parts, there are 4 that are defective (type I).
If we take a sample of size 5 from those car parts, the
probability of finding 2 that are defective will be given by
the HyperGeometric distribution (needs double checking).
}
\details{
We recommend reading this documentation on
\url{https://alexpghayes.github.io/distributions3/}, where the math
will render with additional detail and much greater clarity.

In the following, let \eqn{X} be a HyperGeometric random variable with
success probability \code{p} = \eqn{p = m/(m+n)}.

\strong{Support}: \eqn{x \in { \{\max{(0, k-n)}, \dots, \min{(k,m)}}\}}

\strong{Mean}: \eqn{\frac{km}{n+m} = kp}

\strong{Variance}: \eqn{\frac{km(n)(n+m-k)}{(n+m)^2 (n+m-1)} =
  kp(1-p)(1 - \frac{k-1}{m+n-1})}

\strong{Probability mass function (p.m.f)}:

\deqn{
    P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
  }{
    P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
  }

\strong{Cumulative distribution function (c.d.f)}:

\deqn{
    P(X \le k) \approx \Phi\Big(\frac{x - kp}{\sqrt{kp(1-p)}}\Big)
 }
\strong{Moment generating function (m.g.f)}:

Not useful.
}
\examples{

set.seed(27)

X <- HyperGeometric(4, 5, 8)
X

random(X, 10)

pdf(X, 2)
log_pdf(X, 2)

cdf(X, 4)
quantile(X, 0.7)
}
\seealso{
Other discrete distributions: 
\code{\link{Bernoulli}()},
\code{\link{Binomial}()},
\code{\link{Categorical}()},
\code{\link{Geometric}()},
\code{\link{HurdleNegativeBinomial}()},
\code{\link{HurdlePoisson}()},
\code{\link{Multinomial}()},
\code{\link{NegativeBinomial}()},
\code{\link{Poisson}()},
\code{\link{PoissonBinomial}()},
\code{\link{ZINegativeBinomial}()},
\code{\link{ZIPoisson}()},
\code{\link{ZTNegativeBinomial}()},
\code{\link{ZTPoisson}()}
}
\concept{discrete distributions}
