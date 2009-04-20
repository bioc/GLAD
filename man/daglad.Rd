\name{daglad}
\alias{daglad}
\alias{daglad.profileCGH}
\encoding{latin1}
\title{Analysis of array CGH data}
\description{
  This function allows the detection of breakpoints in genomic profiles
  obtained by array CGH technology and affects a status (gain, normal
  or lost) to each clone.
}
\usage{

\method{daglad}{profileCGH}(profileCGH, mediancenter=FALSE, normalrefcenter=FALSE, genomestep=FALSE,
                  smoothfunc="lawsglad", lkern="Exponential", model="Gaussian",
                  qlambda=0.999,  bandwidth=10, sigma=NULL, base=FALSE, round=1.5,
                  lambdabreak=8, lambdaclusterGen=40, param=c(d=6), alpha=0.001, msize=5,
                  method="centroid", nmin=1, nmax=8,
                  amplicon=1, deletion=-5, deltaN=0.10,  forceGL=c(-0.15,0.15), nbsigma=3,
                  MinBkpWeight=0.35, CheckBkpPos=TRUE, assignGNLOut=TRUE,
                  breaksFdrQ = 0.0001, haarStartLevel = 1, haarEndLevel = 5,
                  verbose=FALSE, ...)

}


\arguments{
  \item{profileCGH}{Object of class \code{\link{profileCGH}}}

  \item{mediancenter}{If \code{TRUE}, LogRatio are center on their median.}

  \item{genomestep}{If \code{TRUE}, a smoothing step over the whole
    genome is performed and a "clustering throughout the genome" allows
    to identify a cluster corresponding to the Normal DNA level. The threshold used in the \code{daglad}
  function (\code{deltaN, forceGL, amplicon, deletion}) and then
  compared to the median of this cluster.}
  
  \item{normalrefcenter}{If \code{TRUE}, the LogRatio are centered
    through the median of the cluster identified during the \code{genomestep}.}

  \item{smoothfunc}{Type of algorithm used to smooth \code{LogRatio} by a
    piecewise constant function. Choose either \code{lawsglad},
    \code{haarseg}, \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}}.}

  \item{lkern}{lkern determines the location kernel to be used
    (see \code{\link[aws]{laws}} for details).}

  \item{model}{model determines the distribution type of LogRatio
    (see \code{\link[aws]{laws}} for details).}

  \item{qlambda}{qlambda determines the scale parameter qlambda for the
    stochastic penalty (see \code{\link[aws]{laws}} for details).}
  
  \item{base}{If TRUE, the position of clone is the physical position onto
    the chromosome, otherwise the rank position is used.}
  
  \item{sigma}{Value to be passed to either argument \code{sigma2}   
    of \code{\link[aws]{aws}} function or \code{shape} of
    \code{\link[aws]{laws}}. If \code{NULL}, sigma is calculated from
    the data.}

  
  \item{bandwidth}{Set the maximal bandwidth \code{hmax} in the
    \code{\link[aws]{aws}} or  \code{\link[aws]{laws}} function. For
    example, if \code{bandwidth=10} then the \code{hmax} value is set
    to 10*\eqn{X_N} where \eqn{X_N} is the position of the last clone.}

  \item{round}{The smoothing results of either \code{\link[aws]{aws}}
    or \code{\link[aws]{laws}} function are rounded or not depending on
    the \code{round} argument. The \code{round} value is passed to the
    argument \code{digits} of the \code{\link[base]{round}} function.}

  
  \item{lambdabreak}{Penalty term (\eqn{\lambda'}) used during the 
    "Optimization of the number of breakpoints" step.}
  
  
  \item{lambdaclusterGen}{Penalty term (\eqn{\lambda*}) used during the "clustering throughout the genome" step.}

    
    \item{param}{Parameter of kernel used in the penalty term.}
    
    \item{alpha}{Risk alpha used for the "Outlier detection" step.}

    \item{msize}{The outliers MAD are calculated on regions with a
      cardinality greater or equal to msize.}
    
    \item{method}{The agglomeration method to be used during the "clustering throughout the genome" steps.}

    \item{nmin}{Minimum number of clusters (N*max) allowed during
      the "clustering throughout the genome" clustering step.}
    
    \item{nmax}{Maximum number of clusters (N*max) allowed during
      the "clustering throughout the genome" clustering step.}

    \item{amplicon}{Level (and outliers) with a smoothing value (log-ratio
      value) greater than this threshold are consider as amplicon. Note
      that first, the data are centered on the normal reference value
      computed during the "clustering throughout the genome" step.}

    \item{deletion}{Level (and outliers) with a smoothing value (log-ratio
      value) lower than this threshold are consider as deletion. Note
      that first, the data are centered on the normal reference value
      computed during the "clustering throughout the genome" step.}


    \item{deltaN}{Region with smoothing values in between the interval
      [-deltaN,+deltaN] are supposed to be normal.}
    
    \item{forceGL}{Level with smoothing value greater (lower) than
      \code{rangeGL[1]} (\code{rangeGL[2]}) are considered as gain
      (lost). Note that first, the data are centered on the normal reference value
      computed during the "clustering throughout the genome" step.}

    \item{nbsigma}{For each breakpoints, a weight is calculated which is
    a function of absolute value of the Gap between the smoothing values
  of the two consecutive regions. Weight = 1-
  kernelpen(abs(Gap),param=c(d=nbsigma*Sigma)) where Sigma is the
  standard deviation of the LogRatio.}

    \item{MinBkpWeight}{Breakpoints which \code{GNLchange}==0 and
      \code{Weight} less than \code{MinBkpWeight} are discarded.}

    \item{CheckBkpPos}{If \code{TRUE}, the accuracy position of each
      breakpoints is checked.}
    
  \item{assignGNLOut}{If \code{FALSE} the status (gain/normal/loss) is
  not assigned for outliers.}  


\item{breaksFdrQ}{breaksFdrQ for HaarSeg algorithm.}
\item{haarStartLevel}{haarStartLevel for
  HaarSeg algorithm.}
\item{haarEndLevel}{ for HaarSeg algorithm.}

    \item{verbose}{If \code{TRUE} some information are printed.}  
    
    \item{...}{}
}


\details{The function \code{daglad} implements a slightly modified
  version of the methodology described in the article : Analysis of array CGH data: from signal
  ratio to gain and loss of DNA regions (Hupé et al., Bioinformatics,
  2004). For smoothing, it is possible
  to use either the AWS algorithm (Polzehl and Spokoiny, 2002) or the HaarSeg algorithm (Ben-Yaacov and Eldar, Bioinformatics, 2008).
  The \code{daglad} function allows to choose some threshold to help the algorithm to identify the status of the genomic regions. The threshodls are given in the following parameters:
  \itemize{
    \item deltaN
    \item forceGL
    \item deletion
    \item amplicon
    }
}

\value{
  \item{ }{An object of class "profileCGH" with the following
    attributes:}

  \item{\bold{profileValues}}{a data.frame with the following added information:

    \itemize{
      
      \item{\bold{Smoothing}}{The smoothing values correspond to the median
	of each Level}

      
      \item{\bold{Breakpoints}}{The last position of a region with identical amount
	of DNA is flagged by 1 otherwise it is 0. Note that during the
	"Optimization of the number of breakpoints" step, removed
	breakpoints are flagged by -1.}


      \item{\bold{Level}}{Each position with equal smoothing value are labelled the
	same way with an integer value starting from one. The label is
	incremented by one when a new level occurs or when moving to the
	next chromosome.}
      
      \item{\bold{OutliersAws}}{Each AWS outliers are flagged by -1 (if it is
	in the \eqn{\alpha/2} lower tail of the distribution) or 1 (if it is
	in the \eqn{\alpha/2} upper tail of the distribution)
	otherwise  it is 0.}
      

      \item{\bold{OutliersMad}}{Each MAD outliers are flagged by -1 (if it is
	in the \eqn{\alpha/2} lower tail of the distribution) or 1 (if it is
	in the \eqn{\alpha/2} upper tail of the distribution)
	otherwise  it is 0.}

      \item{\bold{OutliersTot}}{OutliersAws + OutliersMad.}


      \item{\bold{NormalRef}}{Clusters which have been used to set the normal
	reference during the "clustering throughout the genome" step are
	code by 0. Note that if \code{genomestep=FALSE}, all the value
	are set to 0.}

      \item{\bold{ZoneGNL}}{Status of each clone: Gain is coded by 1, Loss by -1,
	Amplicon by 2, deletion by -10  and
	Normal by 0.}
      
    }

  }


  \item{\bold{BkpInfo}}{a data.frame sum up the information for each
    breakpoint:

    \itemize{
      \item{\bold{Chromosome}}{Chromosome name.}
      \item{\bold{Smoothing}}{Smoothing value for the breakpoint.}
      \item{\bold{Gap}}{absolute value of the gap between the smoothing values
  of the two consecutive regions.}
      \item{\bold{Sigma}}{The estimation of the standard-deviation of the  chromosome.}
      \item{\bold{Weight}}{1 - \code{kernelpen}(Gap, type, param=c(d=nbsigma*Sigma))}
      \item{\bold{ZoneGNL}}{Status of the level where is the breakpoint.}
      \item{\bold{GNLchange}}{Takes the value 1 if the ZoneGNL of the two
	consecutive regions are different.}
      \item{\bold{LogRatio}}{Test over Reference log-ratio.}   
    }
    

  }

  \item{\bold{NormalRef}}{If \code{genomestep=TRUE} and
    \code{normalrefcenter=FALSE}, then NormalRef is the median of the cluster which has been used to set the normal
	reference during the "clustering throughout the genome"
	step. Otherwise NormalRef is 0.}
  
  
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

\seealso{\code{\link{glad}}.}

\keyword{models}

\examples{

data(snijders)
gm13330$Clone <- gm13330$BAC
profileCGH <- as.profileCGH(gm13330)


###########################################################
###
###  daglad function
###
###########################################################


res <- daglad(profileCGH, mediancenter=FALSE, normalrefcenter=FALSE, genomestep=FALSE,
              smoothfunc="lawsglad", lkern="Exponential", model="Gaussian",
              qlambda=0.999,  bandwidth=10, base=FALSE, round=1.5,
              lambdabreak=8, lambdaclusterGen=40, param=c(d=6), alpha=0.001, msize=5,
              method="centroid", nmin=1, nmax=8,
              amplicon=1, deletion=-5, deltaN=0.10,  forceGL=c(-0.15,0.15), nbsigma=3,
              MinBkpWeight=0.35, CheckBkpPos=TRUE)


### Genomic profile on the whole genome
plotProfile(res, unit=3, Bkp=TRUE, labels=FALSE, Smoothing="Smoothing",
main="Breakpoints detection: DAGLAD analysis")



###Genomic profile for chromosome 1
plotProfile(res, unit=3, Bkp=TRUE, labels=TRUE, Chromosome=1,
Smoothing="Smoothing", main="Chromosome 1: DAGLAD analysis")


### The standard-deviation of LogRatio are:
res$SigmaC


### The list of breakpoints is:
res$BkpInfo


}

