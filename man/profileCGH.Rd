\name{profileCGH}
\alias{profileCGH}
\alias{profileChr}

\encoding{latin1}
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
  \item{PosOrder}{The rank position of each clone on the genome.}
  \item{PosBase}{The base position of each clone on the genome.}
  \item{Chromosome}{Chromosome name.}
  \item{Clone}{The name of the corresponding clone.}
  \item{...}{Other elements can be added.}

  
  }

  \details{\code{LogRatio, Chromosome} and \code{PosOrder} are compulsory.}
  
\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Philippe Hupé, \email{glad@curie.fr}.}

\seealso{\code{\link{glad}}, \code{\link{as.profileCGH}}.}

\keyword{classes}

\examples{
  
data(snijders)
gm13330$Clone <- gm13330$BAC
profileCGH <- as.profileCGH(gm13330)
class(profileCGH) <- "profileCGH"

profileChr <- as.profileCGH(gm13330[which(gm13330$Chromosome==1),])
class(profileChr) <- "profileChr"


}


