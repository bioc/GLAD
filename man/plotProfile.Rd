\name{plotProfile}
\alias{plotProfile}
\alias{plotProfile.profileCGH}
\title{Plot genomic profile and cytogenetic banding}
\description{Plot genomic profile with breakpoints, outliers, smoothing
  line and cytogenetic banding.

}
\usage{

plotProfile.profileCGH(profileCGH, variable="LogRatio", Chromosome=NULL,
                       Smoothing=NULL, Bkp=FALSE,
                       labels=TRUE, plotband=TRUE, unit=0,
                       colDAGLAD=c("black","blue","red","green","yellow"),
                       colCytoBand=c("white","darkblue"),
                       colCentro="red", text=NULL, ...)




}

\arguments{
  \item{profileCGH}{Object of class \code{\link{profileCGH}}}
  \item{variable}{The variable to be plot.}
  \item{Chromosome}{A numeric vector with chromosome number to be
    plotted. Use 23 and 24 for chromosome X and Y respectively. If
    \code{NULL}, all the genome is plotted.}
  \item{Smoothing}{The variable used to plot the smoothing line. If
    \code{NULL}, nothing is plotted.}
  \item{Bkp}{If \code{TRUE}, the breakpoints are represented by a
    vertical red dashed line.}
  \item{labels}{If \code{TRUE}, the labels of the cytogenetic banding
    are written.}
  \item{plotband}{If \code{TRUE}, the cytogenetic banding are plotted.}
  \item{unit}{Give the unit of the PosBase. For example if \code{unit=3},
    PosBase are in Kb, if \code{unit=6}, PosBase are in Mb, ...}
  \item{colDAGLAD}{Color code to plot Deletion, Amplification, Gain,
    Lost and Normal status.}
  \item{colCytoBand}{Color code for cytogenetic banding.}
  \item{colCentro}{Color code for centromere.}
  \item{text}{A list with the parameters to be passed to the function \code{\link[graphics]{text}}.}
  \item{...}{}



  }


\details{}

\value{A plot}

\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Philippe Hupé, \email{glad@curie.fr}.}

\seealso{}

\keyword{}

\examples{

### Cytogenetic banding information
data(cytoband)

###
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


### Genomic profile on the whole genome
plotProfile(res, unit=3, Bkp=TRUE, labels=FALSE, Smoothing="Smoothing", plotband=FALSE)


### Genomic profile on the whole genome and cytogenetic banding
plotProfile(res, unit=3, Bkp=TRUE, labels=FALSE, Smoothing="Smoothing")


### Genomic profile for chromosome 1
text <- list(x=c(90000,200000),y=c(0.15,0.3),labels=c("NORMAL","GAIN"), cex=2)
plotProfile(res, unit=3, Bkp=TRUE, labels=TRUE, Chromosome=1,
Smoothing="Smoothing", plotband=FALSE, text=text)

### Genomic profile for chromosome 1 and cytogenetic banding with labels
text <- list(x=c(90000,200000),y=c(0.15,0.3),labels=c("NORMAL","GAIN"), cex=2)
plotProfile(res, unit=3, Bkp=TRUE, labels=TRUE, Chromosome=1,
Smoothing="Smoothing", text=text, main="Chromosome 1")



}

