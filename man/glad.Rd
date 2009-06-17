\name{glad}
\encoding{latin1}
\alias{glad}
\alias{glad.profileCGH}
\title{Analysis of array CGH data}
\description{
  This function allows the detection of breakpoints in genomic profiles
  obtained by array CGH technology and affects a status (gain, normal
  or lost) to each clone.
}
\usage{
\method{glad}{profileCGH}(profileCGH, mediancenter=FALSE,
                smoothfunc="lawsglad", bandwidth=10, round=1.5,
                model="Gaussian", lkern="Exponential", qlambda=0.999,
                base=FALSE, sigma,
                lambdabreak=8, lambdacluster=8, lambdaclusterGen=40,
                type="tricubic", param=c(d=6),
                alpha=0.001, msize=5,
                method="centroid", nmax=8, assignGNLOut=TRUE,
                breaksFdrQ = 0.0001, haarStartLevel = 1, haarEndLevel = 5,
                verbose=FALSE, ...)

}

\arguments{
  \item{profileCGH}{Object of class \code{\link{profileCGH}}}

  \item{mediancenter}{If \code{TRUE}, LogRatio are centered on their median.}
  
  \item{smoothfunc}{Type of algorithm used to smooth \code{LogRatio} by a
    piecewise constant function. Choose either \code{lawsglad},
    \code{haarseg}, \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}}.}

  
  \item{bandwidth}{Set the maximal bandwidth \code{hmax} in the
    \code{\link[aws]{aws}} or  \code{\link[aws]{laws}} function. For
    example, if \code{bandwidth=10} then the \code{hmax} value is set
    to 10*\eqn{X_N} where \eqn{X_N} is the position of the last clone.}

  \item{round}{The smoothing results are rounded or not depending on
    the \code{round} argument. The \code{round} value is passed to the
    argument \code{digits} of the \code{\link[base]{round}} function.}

  \item{model}{Determines the distribution type of the LogRatio. Keep
    always the model as "Gaussian" (see \code{\link[aws]{laws}}).}

  \item{lkern}{Determines the location kernel to be used (see \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}}).}

  \item{qlambda}{Determines the scale parameter for the
    stochastic penalty (see \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}})}
  
  \item{base}{If \code{TRUE}, the position of clone is the physical position on
    the chromosome, otherwise the rank position is used.}
  
  \item{sigma}{Value to be passed to either argument \code{sigma2}   
    of\code{ \link[aws]{aws}} function or \code{shape} of
    \code{\link[aws]{laws}}. If \code{NULL}, sigma is calculated from
    the data.}
  
  
  \item{lambdabreak}{Penalty term (\eqn{\lambda'}) used during the 
    \bold{Optimization of the number of breakpoints} step.}
  
  \item{lambdacluster}{Penalty term (\eqn{\lambda*}) used during the \bold{MSHR
      clustering by chromosome} step.}
  
  \item{lambdaclusterGen}{Penalty term (\eqn{\lambda*}) used during the \bold{HCSR
      clustering throughout the genome} step.}

  \item{type}{Type of kernel function used in the penalty term during the \bold{Optimization of the
      number of breakpoints} step, the \bold{MSHR
      clustering by chromosome} step and the \bold{HCSR
      clustering throughout the genome} step.}
  
  \item{param}{Parameter of kernel used in the penalty term.}
  
  \item{alpha}{Risk alpha used for the \bold{Outlier detection} step.}

  \item{msize}{The outliers MAD are calculated on regions with a
    cardinality greater or equal to msize.}
  
  \item{method}{The agglomeration method to be used during the \bold{MSHR
      clustering by chromosome} and the \bold{HCSR
      clustering throughout the genome} clustering steps.}
  
  \item{nmax}{Maximum number of clusters (N*max) allowed during
    the the \bold{MSHR
      clustering by chromosome} and the \bold{HCSR
      clustering throughout the genome} clustering steps.}


  \item{assignGNLOut}{If \code{FALSE} the status (gain/normal/loss) is
  not assigned for outliers.}  

\item{breaksFdrQ}{breaksFdrQ for HaarSeg algorithm.}
\item{haarStartLevel}{haarStartLevel for
  HaarSeg algorithm.}
\item{haarEndLevel}{ for HaarSeg algorithm.}

  \item{verbose}{If \code{TRUE} some information are printed}  
  
  \item{...}{}
}


\details{

  The function \code{glad} implements the methodology which
  is described in the article: Analysis of array CGH data: from signal
  ratio to gain and loss of DNA regions (Hupé et al., Bioinformatics, 2004).

  The principles of the GLAD algorithm:
  First, the detection of breakpoints is based on the estimation of a
  piecewise constant function with the Adaptive Weights Smoothing (AWS)
  procedure (Polzehl and Spokoiny, 2002). Alternatively, it is possible
  to use the HaarSeg algorithm (Ben-Yaacov and Eldar, Bioinformatics, 2008).
  Then, a procedure based on penalyzed maximum likelihood optimizes the number of
  breakpoints and removes the undesirable breakpoints.
  Finally, based on the regions previously identified, a two-step
  unsupervised classification (\bold{MSHR 
    clustering by chromosome} and the \bold{HCSR
    clustering throughout the genome}) with model selection criteria
  allows a status to be assigned for each region (gain, loss or
  normal).

  Main parameters to be tuned:
  \tabular{ll}{
    \code{qlambda} \tab if you want the smoothing to fit some very local effect, choose a smaller \code{qlambda}. \cr
    \code{bandwidth} \tab choose a bandwidth not to small otherwise you
    will have a lot of little discontinuities. \cr
    \code{lambdabreak} \tab The higher the parameter is, the higher the number
    of undesirable breakpoints is.\cr
    \code{lambdacluster} \tab The higher the parameter is, the higher is
  the
  number of the regions
    within a chromosome which belong to the same cluster. \cr
    \code{lambdaclusterGen} \tab More the parameter is high more the
    regions over the whole genome are supposed to belong to the same cluster.
    
    }

}

\value{



  
  \item{ }{An object of class "profileCGH" with the following attributes:}

  \item{profileValues: }{a data.frame with the following added information:

  \itemize{
    
    
    \item{\bold{Smoothing}}{The smoothing values correspond to the median
      of each \bold{MSHR (i.e. \code{Region}).}}
    
    \item{\bold{Breakpoints}}{The last position of a region with identical amount
      of DNA is flagged by 1 otherwise it is 0. Note that during the
      "Optimization of the number of breakpoints" step, removed
      breakpoints are flagged by -1.}
    

    \item{\bold{Region}}{Each position between two breakpoints are labelled the
      same way with an integer value starting from one. The label is
      incremented by one when a new breakpoint is found or when moving to
      the next chromosome. The variable \code{region} is what we call MSHR.}
    
    \item{\bold{Level}}{Each position with equal smoothing value is labelled the
      same way with an integer value starting from one. The label is
      incremented by one when a new level is found or when moving to the
      next chromosome.}
    
    \item{\bold{OutliersAws}}{Each AWS outliers are flagged -1 or 1 
      otherwise  it is 0.}
    

    \item{\bold{OutliersMad}}{Each MAD outliers are flagged -1 (if it is
      in the \eqn{\alpha/2} lower tail of the distribution) or 1 (if it is
      in the \eqn{\alpha/2} upper tail of the distribution)
      otherwise  it is 0.}

    \item{\bold{OutliersTot}}{OutliersAws + OutliersMad.}

    \item{\bold{ZoneChr}}{Clusters identified after \bold{MSHR (i.e. \code{Region})
	clustering by chromosome}.}

    \item{\bold{ZoneGen}}{Clusters identified after \bold{HCSR clustering throughout the
	genome}.}

    \item{\bold{ZoneGNL}}{Status of each clone : Gain is coded by 1, Loss by -1 and
      Normal by 0.}
    
  }
}

\item{BkpInfo: }{the data.frame attribute \code{BkpInfo} which gives
  the list of breakpoints:
\itemize{
  \item{\bold{PosOrder}}{The rank position of each clone on the genome.}
  \item{\bold{PosBase}}{The base position of each clone on the genome.}
  \item{\bold{Chromosome}}{Chromosome name.}
}
}


\item{SigmaC: }{the data.frame attribute \code{SigmaC} gives the estimation of the LogRatio standard-deviation for each chromosome:
\itemize{
  \item{\bold{Chromosome}}{Chromosome name.}
  \item{\bold{Value}}{The estimation is based on the Inter Quartile Range.}
}
}

  
}


\references{
\itemize{
  \item{Hupé et al. (Bioinformatics, 2004)}{ Analysis of array CGH data: from signal ratio to gain and loss of DNA  regions.}
  \item{Polzehl and Spokoiny (WIAS-Preprint 787, 2002)}{Local likelihood modelling by adaptive weights smoothing.}
  \item{Ben-Yaacov and Eldar (Bioinformatics, 2008)}{A fast and flexible method for the segmentation of aCGH data.}
}
}

\note{People interested in tools dealing with array CGH analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Philippe Hupé, \email{glad@curie.fr}.}

\seealso{\code{\link{profileCGH}}, \code{\link{as.profileCGH}}, \code{\link{plotProfile}}.}


\keyword{models}

\examples{

data(snijders)

### Creation of "profileCGH" object
gm13330$Clone <- gm13330$BAC
profileCGH <- as.profileCGH(gm13330)



###########################################################
###
###  glad function as described in Hupé et al. (2004)
###
###########################################################


res <- glad(profileCGH, mediancenter=FALSE,
                smoothfunc="lawsglad", bandwidth=10, round=1.5,
                model="Gaussian", lkern="Exponential", qlambda=0.999,
                base=FALSE,
                lambdabreak=8, lambdacluster=8, lambdaclusterGen=40,
                type="tricubic", param=c(d=6),
                alpha=0.001, msize=5,
                method="centroid", nmax=8,
                verbose=FALSE)

### cytoband data to plot chromosomes
data(cytoband)

### Genomic profile on the whole genome
plotProfile(res, unit=3, Bkp=TRUE, labels=FALSE, Smoothing="Smoothing",
main="Breakpoints detection: GLAD analysis", cytoband = cytoband)

###Genomic profile for chromosome 1
plotProfile(res, unit=3, Bkp=TRUE, labels=TRUE, Chromosome=1,
Smoothing="Smoothing", main="Chromosome 1: GLAD analysis", cytoband = cytoband)

### The standard-deviation of LogRatio are:
res$SigmaC

### The list of breakpoints is:
res$BkpInfo

}

