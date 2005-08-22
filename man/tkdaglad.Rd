\name{tkdaglad}
\alias{tkdaglad}
\alias{tkdaglad.default}
\alias{tkglad}
\alias{tkglad.default}


\title{Graphical interface for GLAD package}
\description{
A graphical interface to analyse array CGH data.
}


\arguments{
  \item{list}{A character vector with the array to be analysed}
  }


  
\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Philippe Hupé, \email{glad@curie.fr}.}

\seealso{\code{\link{glad}}, \code{\link{daglad}}, \code{\link{plotProfile}}.}

\keyword{}

\examples{
  
data(snijders)
array1 <- as.profileCGH(gm13330)
array2 <- as.profileCGH(gm04435)

## tkdaglad(c("array1","array2"))
## tkglad(c("array1","array2"))


}


