\name{affectationGNL}
\alias{affectationGNL}
\alias{affectationGNL.profileCGH}
\title{Label Assignment}
\description{
  The function \code{affectation.profileCGH} gives a status (gained,
  normal or lost) to each BAC.
}  
\usage{
affectationGNL.profileCGH(profileCGH, alpha=0.001, region="Region", ...)
}

\arguments{
  \item{profileCGH}{Object of class \code{\link{profileCGH}}.}
  
  \item{alpha}{Outlier detected in lost or gained region that lie
    between the \eqn{\alpha/2} lower tail and the \eqn{\alpha/2} upper tail of the estimated distribution of the
    normal DNA population are considered as normal.}
  
  \item{region}{Specifies the name of the variable which has been used
    during the MSHR clustering by chromosome step.}
  \item{...}{ }
 
   
}



\value{
An object of class \code{\link{profileCGH}} with the following added information in
the data.frame attribute \code{profileValues} :
\item{ZoneGNL}{Gain is coded by 1, Loss  by -1 and Normal by 0.}
}

\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}


\seealso{\code{\link{glad}}.}

\keyword{}

