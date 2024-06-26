% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.default.R
\name{simulate.default}
\alias{simulate.default}
\title{Simulate responses from fitted model objects}
\usage{
\method{simulate}{default}(object, nsim = 1, seed = NULL, ...)
}
\arguments{
\item{object}{An object for which a \code{\link[distributions3]{prodist}}
method is available.}

\item{nsim}{The number of response vectors to simulate. Should be a positive
integer. Defaults to 1.}

\item{seed}{An optional random seed that is to be set using \code{\link{set.seed}}
prior to drawing the random sample. The previous random seed from the global
environment (if any) is restored afterwards.}

\item{...}{Arguments passed to \code{\link[distributions3]{simulate.distribution}}.}
}
\value{
A data frame with an attribute \code{"seed"} containing the
\code{.Random.seed} from before the simulation.
}
\description{
Default method for simulating new responses from any model object
with a \code{\link[distributions3]{prodist}} method (for extracting a
probability distribution object).
}
\details{
This default method simply combines two building blocks provided in this
package: (1) \code{\link[distributions3]{prodist}} for extracting the probability
distribution from a fitted model object, (2) \code{\link[distributions3]{simulate.distribution}}
for simulating new observations from this distribution (internally calling
\code{\link[distributions3]{random}}).

Thus, this enables simulation from any fitted model object that provides a
\code{prodist} method. It waives the need to implement a dedicated
\code{\link[stats]{simulate}} method for this model class.
}
\examples{
## Poisson GLM for FIFA 2018 goals
data("FIFA2018", package = "distributions3")
m <- glm(goals ~ difference, data = FIFA2018, family = poisson)

## simulate new goals via glm method
set.seed(0)
g_glm <- simulate(m, n = 3)

## alternatively use the new default method
set.seed(0)
g_default <- simulate.default(m, n = 3)

## same results
all.equal(g_glm, g_default, check.attributes = FALSE)
}
