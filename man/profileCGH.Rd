\name{profileCGH}
\alias{profileCGH}
\alias{profileChr}


\title{Objects of Class profileCGH and profileChr}
\description{
Description of the objects \code{profileCGH} and \code{profileChr}. The
last object corresponds to data of only one chromosome.
}


\value{
  Objects \code{profileCGH} and \code{profileChr} are composed of a list with the first element
  \code{profileValues} which is a \code{data.frame} with the following columns
  names:
  \item{LogRatio}{Test over Reference log-ratio.}
  \item{PosOrder}{The rank position of each BAC on the genome.}
  \item{PosBase}{The base position of each BAC on the genome.}
  \item{Chromosome}{Chromosome name.}
  \item{BAC}{The name of the corresponding BAC.}
  \item{...}{Other elements can be added.}

  
  }

  \details{\code{LogRatio, Chromosome} and \code{PosOrder} are compulsory.}
  
\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}

\seealso{\code{\link{glad}}.}

\keyword{}

\examples{
  
data(snijders)
profileCGH <- list(profileValues=gm13330)
class(profileCGH) <- "profileCGH"

profileChr <- list(profileValues=gm13330[which(gm13330$Chromosome==1),])
class(profileChr) <- "profileChr"


}


