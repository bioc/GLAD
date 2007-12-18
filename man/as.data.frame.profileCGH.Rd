\name{as.data.frame.profileCGH}
\encoding{latin1}
\title{profileCGH consercion}
\alias{as.data.frame.profileCGH}
\description{Convert a profileCGH object into a data.frame.}
\usage{\method{as.data.frame}{profileCGH}(x, row.names = NULL, optional = FALSE, ...)}
\arguments{
  \item{x}{The object to converted into data.frame.}
  \item{row.names}{NULL or a character vector giving the row names for the
    data frame.  Missing values are not allowed.}

  \item{optional}{logical. If 'TRUE', setting row names and converting column
    names (to syntactic names) is optional.}

  \item{...}{ }


}

\value{A data.frame object}

\details{The attributes \code{profileValues} and \code{profileValuesNA}
  are binded into a data.frame.}


\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}


\author{Philippe Hupé, \email{glad@curie.fr}}

\seealso{\code{\link{as.profileCGH}}}

\keyword{manip}

\examples{

data(snijders)

### Creation of "profileCGH" object
profileCGH <- as.profileCGH(gm13330)



###########################################################
###
###  glad function as described in Hupé et al. (2004)
###
###########################################################


res <- glad(profileCGH, mediancenter=FALSE,
                smoothfunc="lawsglad", bandwidth=10, round=2,
                model="Gaussian", lkern="Exponential", qlambda=0.999,
                base=FALSE,
                lambdabreak=8, lambdacluster=8, lambdaclusterGen=40,
                type="tricubic", param=c(d=6),
                alpha=0.001, msize=5,
                method="centroid", nmax=8,
                verbose=FALSE)


res <- as.data.frame(res)

}
