\name{array}
\docType{data}
\alias{array1}
\alias{array2}
\alias{array3}
\title{Bladder cancer CGH data}
\description{
Bladder cancer data from 3 arrays CGH (Comparative Genomic
Hybridyzation). Arrays dimension are 4 blocs per column, 4 blocs per
row, 21 columns per bloc and 22 rows by blocs.
}
\usage{data(arrayCGH)}
\format{
A data frame composed of the following elements :
  \describe{
    \item{Log2Rat}{Log 2 ratio.}
    \item{Position}{BAC position on the genome.}
    \item{CHROMOSOME}{Chromosome.}
    \item{Col}{Column location on the array.}
    \item{Row}{Row location on the array.}
  }
}
\source{Institut Curie, \email{glad@curie.fr}.}
\examples{
data(arrayCGH)
data <- array1 #array1 to array3
}
\keyword{datasets}
