\name{findCluster}
\alias{findCluster}
\alias{findCluster.profileChr}
\title{Clustering of genomic regions}
\description{
The function \code{findCluster} allows the clustering of regions.
}
\usage{
findCluster.profileChr(profileChr, region="Region", genome=TRUE,
                       lambda=10, nmin=1, nmax=10, type="tricubic", param=c(d=6), ...)


}

\arguments{
  \item{profileChr}{Object of class \code{profileChr}.}
  \item{region}{Label which indicated where outliers are detected.}
  \item{genome}{If \code{TRUE}, the "MSHR
    clustering by chromosome" clustering step is run otherwise it is the "HCSR
    clustering throughout the genome" clustering step.}
    \item{lambda}{Penalty \eqn{\lambda*}.}
      \item{nmin}{Minimal number of cluster to find.}
  \item{nmax}{Maximal number of cluster to find.}
  \item{type}{Type of kernel to use in the penalty term (this arguments is passed to the
  \code{\link{kernelpen}} function).}
  \item{param}{Parameters of the kernel (this arguments is passed to the
  \code{\link{kernelpen}} function).}
  \item{...}{Risk alpha used in \code{\link{detectOutliers}}.}
}


\value{ An object of class \code{profileChr} with the following element added
  to the data.frame \code{profileValues}:

  either:
  
   \item{ZoneChr}{indicates the regions which belongs to the same HCSR
     (if \code{genome=FALSE}).}

   
   or:

  \item{ZoneGen}{indicates HCSR which belongs to the same cluster (if \code{genome=TRUE}).}

}

\seealso{\code{\link{glad}}, \code{\link{detectOutliers}}}

\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}

\keyword{}


