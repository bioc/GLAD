\name{ChrNumeric}
\encoding{latin1}
\title{Convert chromosome into numeric values}
\alias{ChrNumeric}
\description{Convert chromosome into numeric values.}
\usage{
ChrNumeric(Chromosome)
}
\arguments{
  \item{Chromosome}{A vector with chromosome labels.}
}
\details{For sexual chromosome, labels must contains X or Y which are
  coded by 23 and 24 respectively.}

\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Philippe Hupé, \email{glad@curie.fr}}

\keyword{}

\examples{

Chromosome <- c("1","X","Y","chr X", "ChrX", "chrX", "Chr   Y")
ChrNumeric(Chromosome)

}
