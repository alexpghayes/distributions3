% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{geom_auc}
\alias{geom_auc}
\title{Fill out area under the curve}
\usage{
geom_auc(mapping = NULL, data = NULL, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, from = -Inf,
  to = Inf, ...)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[=aes]{aes()}} or
\code{\link[=aes_]{aes_()}}. If specified and \code{inherit.aes = TRUE} (the
default), it is combined with the default mapping at the top level of the
plot. You must supply \code{mapping} if there is no plot mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{\link[=ggplot]{ggplot()}}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{\link[=fortify]{fortify()}} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data. A \code{function} can be created
from a \code{formula} (e.g. \code{~ head(.x, 10)}).}

\item{position}{Position adjustment, either as a string, or the result of
a call to a position adjustment function.}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[=borders]{borders()}}.}

\item{from}{Left end-point of interval}

\item{to}{right end-point of interval}

\item{...}{Other arguments passed on to \code{\link[=layer]{layer()}}. These are
often aesthetics, used to set an aesthetic to a fixed value, like
\code{colour = "red"} or \code{size = 3}. They may also be parameters
to the paired geom/stat.}
}
\description{
Fill out area under the curve
}
\examples{

X <- Normal()

#plot_pdf(X) + geom_auc(to = -0.645)
#plot_pdf(X) + geom_auc(from = -0.645, to = 0.1)
}