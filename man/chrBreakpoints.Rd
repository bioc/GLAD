\name{chrBreakpoints}
\alias{chrBreakpoints}
\alias{chrBreakpoints.profileCGH}
\title{Chromosomal breakpoints detection}
\description{
The function \code{chrBreakpoints.profileCGH} find breakpoints corresponding to an abrupt change of the DNA amount along the chromosome.}
\usage{
chrBreakpoints.profileCGH (profileCGH, smoothfunc="aws", base=TRUE, sigma, bandwidth=10, round=2, ...)
}

\arguments{
  \item{profileCGH}{Object of class "profileCGH".}
  
  \item{smoothfunc}{Type of function used to smooth \code{LogRatio} by a piecewise constant function. Choose either \code{\link[aws]{aws}} or  \code{\link[aws]{laws}}.}
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
   
   \item{...}{Parameters to be passed to \code{\link[aws]{aws}} or
     \code{\link[aws]{laws}} function. Typical parameters are
     \code{qlambda}, \code{model}, \code{lkern}.}
   
}


\details{The Adaptive Weights Smoothing procedure as described by Polzehl & Spokoiny is used to fit a function piecewise constant.
Then, based on the smoothing results of either \code{\link[aws]{aws}} or \code{\link[aws]{laws}} function, a breakpoint is added when two contiguous smoothing values are different: breakpoints are flagged by 1 in the \code{Breakpoints} vector and the flag corresponds to the last position of identical amount of DNA. A specific process is implemented for outliers detection: an outlier is a position such that the smoothing value on its right and the smoothing value on its left are equal but the smoothing value in this position is different from the right and left values. If the first position (respectively the last position) is different of the second one (respectively the next to last one) then the position is considered as an outlier.
Each outliers are flagged by 1 in the \code{OutliersAws} vector.
}

\value{
An object of class "profileCGH" with the following added information in
the data.frame attribute \code{profileValues}:
\item{Smoothing}{Smoothing results of either \code{\link[aws]{aws}} or \code{\link[aws]{laws}} function after being rounded or not depending on the \code{round} argument.}
\item{Region}{Each position between two breakpoints are labelled the same way with an integer value starting from one. The label is incremented by one when a new breakpoints occurs or when moving to the next chromosome.}
\item{Level}{Each position with equal smoothing value are labelled the same way with an integer value starting from one. The label is incremented by one when a new level occurs or when moving to the next chromosome.}
\item{OutliersAws}{Each outliers detected are flagged by one otherwise it is 0.}
\item{Breakpoints}{The last position of a region with identical amount of DNA is flagged by 1 otherwise it is 0.}
}

\author{Philippe Hupé, \email{Philippe.Hupe@curie.fr}.}

\references{ 
\item{ }{Polzehl, J. and Spokoiny, S. (2002). \emph{Varying coefficient 
regression modeling 
     by adaptive weights smoothing} Manuscript}
\item{ }{Polzehl, J. and Spokoiny, S. (2002b). \emph{Local likelihood 
modelling by adaptive weights smoothing}, WIAS-Preprint 787}
\item{ }{ Polzehl, J. and Spokoiny, S. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp. 335-354}
 }

\seealso{\code{\link[aws]{aws}},\code{ \link[aws]{laws}}.}

\keyword{}

\examples{

data(snijders)
profileCGH <- list(profileValues=gm13330)
class(profileCGH) <- "profileCGH"



# Estimation of the piecewise constant function

res <- chrBreakpoints(profileCGH, smoothfunc="laws",
                      lkern="exponential", model="Gaussian",
                      qlambda=0.999, base=FALSE,bandwidth=10)

plot(LogRatio ~ PosOrder, data=res$profileValues, pch=20)



# Limit between chromosomes

LimitChr <- unique(res$profileValues$LimitChr)+0.5
abline(v=LimitChr, col="grey", lty=2)

lines(res$profileValues$Smoothing ~ res$profileValues$PosOrder, col="green")



# Breakpoints identified

indexBP <- which(res$profileValues$Breakpoints==1)
BP <- res$profileValues$PosOrder[indexBP]+0.5
abline(v=BP, col="red", lty=2)



}

