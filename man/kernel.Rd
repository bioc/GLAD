\name{kernelpen}
\title{Kernelpen function}
\alias{kernelpen}
\description{Krenel function used in the penalty term.}
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

\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}}

\keyword{}
