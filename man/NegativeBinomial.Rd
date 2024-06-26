% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NegativeBinomial.R
\name{NegativeBinomial}
\alias{NegativeBinomial}
\title{Create a negative binomial distribution}
\usage{
NegativeBinomial(size, p = 0.5, mu = size)
}
\arguments{
\item{size}{The target number of successes (greater than \eqn{0})
until the experiment is stopped. Denoted \eqn{r} below.}

\item{p}{The success probability for a given trial. \code{p} can be any
value in \verb{[0, 1]}, and defaults to \code{0.5}.}

\item{mu}{Alternative parameterization via the non-negative mean
of the distribution (instead of the probability \code{p}), defaults to \code{size}.}
}
\value{
A \code{NegativeBinomial} object.
}
\description{
A generalization of the geometric distribution. It is the number
of failures in a sequence of i.i.d. Bernoulli trials before
a specified target number (\eqn{r}) of successes occurs.
}
\details{
We recommend reading this documentation on
\url{https://alexpghayes.github.io/distributions3/}, where the math
will render with additional detail and much greater clarity.

In the following, let \eqn{X} be a negative binomial random variable with
success probability \code{p} = \eqn{p}.

\strong{Support}: \eqn{\{0, 1, 2, 3, ...\}}

\strong{Mean}: \eqn{\frac{(1 - p) r}{p} = \mu}

\strong{Variance}: \eqn{\frac{(1 - p) r}{p^2}}

\strong{Probability mass function (p.m.f.)}:

\deqn{
     f(k) = {k + r - 1 \choose k} \cdot p^r (1-p)^k
  }{
     f(k) = (k+r-1)!/(k!(r-1)!) p^r (1-p)^k
  }

\strong{Cumulative distribution function (c.d.f.)}:

Omitted for now.

\strong{Moment generating function (m.g.f.)}:

\deqn{
     \left(\frac{p}{1 - (1 -p) e^t}\right)^r, t < -\log (1-p)
  }{
     \frac{p^r}{(1 - (1-p) e^t)^r}, t < -\log (1-p)
  }

\strong{Alternative parameterization}: Sometimes, especially when used in
regression models, the negative binomial distribution is parameterized
by its mean \eqn{\mu} (as listed above) plus the size parameter \eqn{r}.
This implies a success probability of \eqn{p = r/(r + \mu)}. This can
also be seen as a generalization of the Poisson distribution where the
assumption of equidispersion (i.e., variance equal to mean) is relaxed.
The negative binomial distribution is overdispersed (i.e., variance greater than mean)
and its variance can also be written as \eqn{\mu + 1/r \mu^2}. The Poisson
distribution is then obtained as \eqn{r} goes to infinity. Note that in this
view it is natural to also allow for non-integer \eqn{r} parameters.
The factorials in the equations above are then expressed in terms of the
gamma function.
}
\examples{

set.seed(27)

X <- NegativeBinomial(size = 5, p = 0.1)
X

random(X, 10)

pdf(X, 50)
log_pdf(X, 50)

cdf(X, 50)
quantile(X, 0.7)

## alternative parameterization of X
Y <- NegativeBinomial(mu = 45, size = 5)
Y
cdf(Y, 50)
quantile(Y, 0.7)
}
\seealso{
Other discrete distributions: 
\code{\link{Bernoulli}()},
\code{\link{Binomial}()},
\code{\link{Categorical}()},
\code{\link{Geometric}()},
\code{\link{HurdleNegativeBinomial}()},
\code{\link{HurdlePoisson}()},
\code{\link{HyperGeometric}()},
\code{\link{Multinomial}()},
\code{\link{Poisson}()},
\code{\link{PoissonBinomial}()},
\code{\link{ZINegativeBinomial}()},
\code{\link{ZIPoisson}()},
\code{\link{ZTNegativeBinomial}()},
\code{\link{ZTPoisson}()}
}
\concept{discrete distributions}
