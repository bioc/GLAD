\name{cluster}
\alias{cluster}
\alias{cluster.default}
\title{Optimal number of cluster}
\description{
The function \code{cluster} in an internal function which allows to find
the optimal number of cluster in the data.
}
\usage{
cluster.default(CGH, Cluster, region, labels, lambda, nmin, nmax, ...)
}

\arguments{
  \item{CGH}{CGH}
  \item{Cluster}{Cluster}
  \item{labels}{labels}
  \item{lambda}{Penalty \eqn{\lambda*}.}
  \item{nmin}{Minimal number of cluster to find.}
  \item{nmax}{Maximal number of cluster to find.}
  \item{region}{Label which indicated where outliers are detected.}
  \item{...}{Parameters to be passed to the \code{\link{kernelpen}} function.}
}


\details{This function is called within \code{\link{findCluster.profileChr}}.}


\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}

\seealso{\code{\link{findCluster.profileChr}}}

\keyword{}


