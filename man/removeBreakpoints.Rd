\name{removeBreakpoints}
\alias{removeBreakpoints}
\alias{removeBreakpoints.profileChr}
\alias{removeLevel}
\alias{removeLevel.profileChr}
\title{Remove breakpoints}
\description{
This function optimizes the number of breakpoints.
}
\usage{
removeBreakpoints.profileChr(profileChr, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, ...)
removeLevel.profileChr <- function(profileChr, lambda=10, ...)

}

\arguments{
  \item{profileChr}{Object of class \code{profileChr}.}

  \item{lambda}{Penalty term (\eqn{\lambda'}) used during the 
    "Optimization of the number of breakpoints" step.}
  
    \item{type}{Type of kernel used in the penalty term.}
    \item{param}{Parameter for the kernel used in the penalty term.}

    \item{verbose}{If \code{TRUE} some information are printed}
  \item{...}{Argument to be passed to \code{\link{detectOutliers}}, typically \code{alpha}.}
}


\details{The function \code{removeBreakpoints} uses the results of
  \code{\link{chrBreakpoints}} and optimizes the number of breakpoints
  previously detected: when a breakpoint is removed, then it is flagged
  out by -1 in the vector \code{Breakpoints} of the data.frame
  attribute of the object of class \code{\link{profileChr}}. \code{removeBreakpoints} is
  used within the global function \code{\link{glad}}. }

\value{Data frame.}

\author{Philippe Hupé, \url{Philippe.Hupe@curie.fr}.}

\seealso{\code{\link{chrBreakpoints}},\code{\link{glad}}, \code{\link{kernelpen}}.}

\keyword{}



