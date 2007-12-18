\name{arrayCGH}
\alias{arrayCGH}

\encoding{latin1}
\title{Object of Class arrayCGH}
\description{
Description of the object \code{arrayCGH}.
}


\value{
  The object \code{arrayCGH} is a list with at least a
  data.frame named \code{arrayValues} and a vector named
  \code{arrayDesign}. The data.frame \code{arrayValues} must contain the
  following fields :
  
  \item{Col}{Vector of columns coordinates.}
  \item{Row}{Vector of rows coordinates.}
  \item{...}{Other elements can be added.}

  The vector \code{arrayDesign} is composed of 4 values :  c(arrayCol,
  arrayRow, SpotCol, SpotRow). The array CGH is represented by arrayRow*arrayCol
  blocs and each bloc is composed of SpotRow*SpotCol spots.

  N.B. : Col takes the values in 1:arrayRow*SpotRow and Row takes the
  values in 1:arrayCol*SpotCol

  
  }

\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

  
\author{Philippe Hupé, \email{glad@curie.fr}.}


\seealso{\code{\link{glad}}.}

\keyword{classes}

\examples{
data(arrayCGH)

# object of class arrayCGH

array <- list(arrayValues=array2, arrayDesign=c(4,4,21,22))
class(array) <- "arrayCGH"

}


