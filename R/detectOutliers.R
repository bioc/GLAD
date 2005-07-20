# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr



detectOutliers <- function(...)
  {
    UseMethod("detectOutliers")
  }



detectOutliers.profileChr <- function(profileChr, region="Region", msize=5, alpha=0.001, verbose=FALSE, ...)
  {


    if (verbose) print("detectOutliers: starting function")
    if (msize<1) stop("msize must be greater or equal to 1")

    if (alpha>1 || alpha <0)stop("alpha must be setted between 0 and 1")
    

    alpha <- qnorm(1-alpha/2)
    
    l <- length(profileChr$profileValues$LogRatio)
    
    res <- .C("detectOutliers",
              as.double(profileChr$profileValues$LogRatio),
              as.integer(profileChr$profileValues[,region]),
              OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
              OutliersMad=integer(l),
              OutliersTot=integer(l),
              as.integer(msize),
              as.double(alpha),
              as.integer(l),
              PACKAGE="GLAD")

    profileChr$profileValues$OutliersMad <- res$OutliersMad
    profileChr$profileValues$OutliersAws <- res$OutliersAws
    profileChr$profileValues$OutliersTot <- res$OutliersTot

    if(verbose) print("detectOutliers: ending function")

    return(profileChr)
    
    
  }
