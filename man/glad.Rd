\name{glad}
\alias{glad}
\alias{glad.profileCGH}
\title{Analysis of array CGH data}
\description{
  This function allows the detection of breakpoints in genomic profiles
  obtained by array CGH technology and affects a status (gain, normal
  or lost) to each BAC.
}
\usage{
glad.profileCGH(profileCGH, smoothfunc="aws", base=FALSE, sigma,
                   bandwidth=10, round=2, lambdabreak=8, lambdacluster=8,
                   lambdaclusterGen=40, type="tricubic",
                   param=c(d=6),alpha=0.001, method="centroid", nmax=8, verbose=FALSE, ...)
}

\arguments{
  \item{profileCGH}{Object of class \code{\link{profileCGH}}}
  \item{smoothfunc}{Type of algorithm used to smooth \code{LogRatio} by a
    piecewise constant function. Choose either \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}}.}
  
  \item{base}{If TRUE, the position of BAC is the physical position onto
    the chromosome, otherwise the rank position is used.}
    
  \item{sigma}{Value to be passed to either argument \code{sigma2}   
      of\code{ \link[aws]{aws}} function or \code{shape} of
      \code{\link[aws]{laws}}. If \code{NULL}, sigma is calculated from
      the data.}

    
   \item{bandwidth}{Set the maximal bandwidth \code{hmax} in the
     \code{\link[aws]{aws}} or  \code{\link[aws]{laws}} function. For
     example, if \code{bandwidth=10} then the \code{hmax} value is set
     to 10*\eqn{X_N} where \eqn{X_N} is the position of the last BAC.}

   \item{round}{The smoothing results of either \code{\link[aws]{aws}}
     or \code{\link[aws]{laws}} function are rounded or not depending on
     the \code{round} argument. The \code{round} value is passed to the
     argument \code{digits} of the \code{\link[base]{round}} function.}

   
  \item{lambdabreak}{Penalty term (\eqn{\lambda'}) used during the 
    "Optimization of the number of breakpoints" step.}
  
  \item{lambdacluster}{Penalty term (\eqn{\lambda*}) used during the "MSHR
    clustering by chromosome" step.}
  
  \item{lambdaclusterGen}{Penalty term (\eqn{\lambda*}) used during the "HCSR
    clustering throughout the genome" step.}

    \item{type}{Type of kernel function used in the penalty term during the "Optimization of the
      number of breakpoints" step, the "MSHR
    clustering by chromosome" step and the "HCSR
    clustering throughout the genome" step.}
    
    \item{param}{Parameter of kernel used in the penalty term.}
  
  \item{alpha}{Risk alpha used for the "Outlier detection" step.}
  
  \item{method}{The agglomeration method to be used during the "MSHR
    clustering by chromosome" and the "HCSR
    clustering throughout the genome" clustering steps.}
  
  \item{nmax}{Maximum number of clusters (N*max) allowed during
    the the "MSHR
    clustering by chromosome" and the "HCSR
    clustering throughout the genome" clustering steps.}

  \item{verbose}{If \code{TRUE} some information are printed}  
    
  \item{...}{parameters to be passed to \code{\link{chrBreakpoints}}
    function. Typically, you will have to specify the following
    arguments : \code{lkern="exponential", model="Gaussian", qlambda=0.999}. }
}


\details{The function \code{glad} implements the methodology which
  is described in the article : Analysis of array CGH data: from signal
  ratio to gain and loss of DNA regions (Hupé et al., 2004 submitted).

  First, \code{\link{chrBreakpoints}} detects breakpoints and
  \code{\link{detectOutliers}} allows the detection of MAD
  outliers. Then, the number of breakpoints is optimized with
  \code{\link{removeBreakpoints}}.  The two-step clustering ("MSHR
    clustering by chromosome" and the "HCSR
    clustering throughout the genome") is performed with
    \code{\link{findCluster}}. The function
    \code{\link{affectationGNL}} give a status to each BAC.
  
}

\value{
  
  \item{Smoothing}{Smoothing results of either \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}} function after being rounded or not
    depending on the \code{round} argument.}
  
  \item{Breakpoints}{The last position of a region with identical amount
    of DNA is flagged by 1 otherwise it is 0. Note that during the
    "Optimization of the number of breakpoints" step, removed
    breakpoints are flagged by -1.}
  

  \item{Region}{Each position between two breakpoints are labelled the
    same way with an integer value starting from one. The label is
    incremented by one when a new breakpoints occurs or when moving to
    the next chromosome. The variable \code{region} is what we call MSHR.}
  
  \item{Level}{Each position with equal smoothing value are labelled the
    same way with an integer value starting from one. The label is
    incremented by one when a new level occurs or when moving to the
    next chromosome.}
  
  \item{OutliersAws}{Each AWS outliers are flagged by -1 (if it is
    in the \eqn{\alpha/2} lower tail of the distribution) or 1 (if it is
    in the \eqn{\alpha/2} upper tail of the distribution)
    otherwise  it is 0.}
  

   \item{OutliersMad}{Each MAD outliers are flagged by -1 (if it is
    in the \eqn{\alpha/2} lower tail of the distribution) or 1 (if it is
    in the \eqn{\alpha/2} upper tail of the distribution)
    otherwise  it is 0.}

  \item{OutliersTot}{OutliersAws + OutliersMad.}

  \item{ZoneChr}{Clusters identified after MSHR (i.e. \code{Region})
    clustering by chromosome.}

  \item{ZoneGen}{Clusters identified after HCSR clustering throughout the
  genome.}

  \item{ZoneGNL}{Status of each BAC : Gain is coded by 1, Loss by -1 and Normal by 0.}

  
  
}

\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}

\seealso{\code{\link{chrBreakpoints}},
  \code{\link{removeBreakpoints}},\code{\link{detectOutliers}},
  \code{\link{findCluster}},  \code{\link{affectationGNL}}.}

\keyword{}

\examples{

data(snijders)
profileCGH <- list(profileValues=gm13330)
class(profileCGH) <- "profileCGH"

res <- glad(profileCGH, smoothfunc="laws", base=FALSE,
               bandwidth=10, round=2, lambdabreak=8, lambdacluster=8,
               lambdaclusterGen=40, alpha=0.001, method="centroid",
               nmax=8, lkern="exponential", model="Gaussian",
               qlambda=0.999)

# color code for region status

col <- rep("yellow",length(res$profileValues$PosOrder))
col[which(res$profileValues$ZoneGNL==-1)] <- "green"
col[which(res$profileValues$ZoneGNL==1)] <- "red"

# outliers

outliers <- rep(20,length(res$profileValues$PosOrder))
outliers[which(res$profileValues$OutliersTot!=0)] <- 13

plot(LogRatio ~ PosOrder, data=res$profileValues, col=col, pch=outliers)

# Limit between chromosomes

LimitChr <- unique(res$profileValues$LimitChr)+0.5
abline(v=LimitChr, col="grey", lty=2)

lines(res$profileValues$Smoothing ~ res$profileValues$PosOrder, col="black")

# Breakpoints identified

indexBP <- which(res$profileValues$Breakpoints==1)
BP <- res$profileValues$PosOrder[indexBP]+0.5
abline(v=BP, col="red", lty=2)


}

