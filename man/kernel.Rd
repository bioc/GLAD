\name{kernelpen}
\encoding{latin1}
\title{Kernelpen function}
\alias{kernelpen}
\description{Kernel function used in the penalty term.}
\usage{
kernelpen(x, type="tricubic", param)

}
\arguments{
  \item{x}{Real Value.}

  \item{type}{Type of kernelpen to be used}

  \item{param}{a named vector.}
}
\details{The only kernel available is the "tricubic" kernel which takes
  the values \eqn{(1 - (x/d)^3)^3}. The value of d is given by
  \code{param=c(d=6)} for example.}

\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Philippe Hupé, \email{glad@curie.fr}}

\keyword{math}
