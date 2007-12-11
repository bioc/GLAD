\name{as.profileCGH}
\encoding{latin1}
\title{Create an object of class profileCGH}
\alias{as.profileCGH}
\alias{as.profileCGH.data.frame}
\description{Create an object of class profileCGH.}
\usage{
as.profileCGH(object,...)
\method{as.profileCGH}{data.frame}(object, infaction=c("value","empty"), value=20, ...)
}
\arguments{
  \item{object}{A data.frame to be convert into profileCGH.}
  \item{infaction}{If "value" then the LogRatio with infinite values
  (-Inf, Inf) are replace by + or - \code{value} according to the sign. If "empty"
  then NAs are put instead.}
\item{value}{replace Inf by \code{value} if \code{infaction} is
  "value".}
\item{...}{}
}

\value{

  \item{ }{A list with the following attributes}
  \item{profileValues}{A data.frame}
  \item{profileValuesNA}{A data.frame}
}

\details{The data.frame to be convert must at least contain the
  following fields: LogRatio, PosOrder, and Chromosome. If the field
  Chromosome is of mode character,  it is automatically converted into a
numeric vector (see \code{\link{ChrNumeric}}); a field ChromosomeChar
contains the character labels. The data.frame to be converted into a
profileCGH objet is split into two data.frame: profileValuesNA contains
the rows for which there is at least a missing value for either
LogRatio, PosOrder or Chromosome; profileValues contains the remaining rows.}


\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}


\author{Philippe Hupé, \email{glad@curie.fr}}

\seealso{\code{\link{as.data.frame.profileCGH}}}

\keyword{}

\examples{

data(snijders)

### Creation of "profileCGH" object
profileCGH <- as.profileCGH(gm13330)

attributes(profileCGH)

}
