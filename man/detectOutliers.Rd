\name{detectOutliers}
\alias{detectOutliers}
\alias{detectOutliers.profileChr}
\title{Outliers detection}
\description{
This function detects MAD outliers.
}
\usage{
detectOutliers.profileChr(profileChr, region="Region", alpha=0.001, ...)
}

\arguments{
  \item{profileChr}{Object of class \code{profileChr}.}
  \item{region}{Label which indicated where outliers are detected.}
  \item{alpha}{Risk alpha.}
  \item{...}{}
}


\details{The MAD outliers are detected within this function. The object
  of class \code{profileChr} has the same characteristics than the object
of class \code{\link{profileCGH}} but corresponds to the values of only
one chromosome.}

\value{
  An object of class \code{profileChr} with the following element added
  to the data.frame \code{profileValues}:

  

   \item{OutliersMad}{Each MAD outliers are flagged by -1 (if it is
    in the lower tail of the distribution) or 1 (if it is
    in the upper tail of the distribution)
    otherwise  it is 0.}

  \item{OutliersTot}{OutliersAws + OutliersMad.}

}

\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}

\seealso{\code{\link{glad.profileCGH}}}

\keyword{}



