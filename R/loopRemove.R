### This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr


loopRemove <- function(...)
  {
    UseMethod("loopRemove")
  }



loopRemove.profileChr <- function(profileChr, sigma, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, msize=5, alpha=0.001, ...)
  {


    if (alpha>1 || alpha <0)stop("alpha must be setted between 0 and 1")
    

    alpha <- qnorm(1-alpha/2)
    
    
    l <- length(profileChr$profileValues[,1])
    res <- .C("loopRemove",
              as.double(profileChr$profileValues$LogRatio),
              Region=as.integer(profileChr$profileValues$Region),
              OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
              OutliersMad=integer(l),
              OutliersTot=integer(l),
              Breakpoints=as.integer(profileChr$profileValues$Breakpoints),
              as.integer(msize),
              as.double(alpha),
              as.double(lambda),
              as.double(param["d"]),
              as.double(sigma),
              as.integer(l),
              PACKAGE="GLAD")


    profileChr$profileValues$Region <- res$Region
    profileChr$profileValues$OutliersAws <- res$OutliersAws
    profileChr$profileValues$OutliersMad <- res$OutliersMad
    profileChr$profileValues$OutliersTot <- res$OutliersTot
    profileChr$profileValues$Breakpoints <- res$Breakpoints

    return(profileChr)
    
  }
