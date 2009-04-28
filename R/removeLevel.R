### This function detects chromosomal breakpoints along genome

### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr





###############################################
#############################################
### ajout 16122004 : test sur les levels


removeLevel <- function(...)
  {
    UseMethod("removeLevel")
  }



removeLevel.profileChr <- function(profileChr, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, msize=5, alpha=0.001, BkpDetected=TRUE, ...)
  {

    if (alpha>1 || alpha <0)stop("alpha must be setted between 0 and 1")    

    if (verbose)
      {
        print("removeLevel: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }

    if(!BkpDetected)
      {
        profileChr <- detectOutliers(profileChr, region="Level", verbose=verbose, msize=msize, alpha=alpha)
        if (verbose) print("removeLevel: ending function")
        return(profileChr)
      }

    profileChr$profileValues$Region <- profileChr$profileValues$Level
    sigma <- profileChr$findClusterSigma

    
    
    l <- length(profileChr$profileValues[,1])
    res <- .C("loopRemove",
              as.double(profileChr$profileValues$LogRatio),
              Region=as.integer(profileChr$profileValues$Region),
              OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
              OutliersMad=integer(l),
              OutliersTot=integer(l),
              Breakpoints=as.integer(profileChr$profileValues$Breakpoints),
              as.integer(msize),
              as.double(qnorm(1-alpha/2)),
              as.double(lambda),
              as.double(param["d"]),
              as.double(sigma),
              as.integer(l),
              PACKAGE="GLAD")


    profileChr$profileValues[,c("Region","OutliersAws","OutliersMad", "OutliersTot", "Breakpoints")] <- res[c("Region","OutliersAws","OutliersMad", "OutliersTot", "Breakpoints")]
    


    profileChr$profileValues$Breakpoints <- 0
    profileChr$profileValues$OutliersAws <- 0
    profileChr$profileValues$NextLogRatio <- profileChr$profileValues$OutliersAws



    updateBkpRL <- .C("updateBkpRL",
                      Region=as.integer(profileChr$profileValues$Region),
                      OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
                      Breakpoints=as.integer(profileChr$profileValues$Breakpoints),
                      as.integer(profileChr$profileValues$Chromosome),
                      as.integer(profileChr$profileValues$PosOrder),
                      NextLogRatio=as.double(profileChr$profileValues$NextLogRatio),
                      as.double(profileChr$profileValues$LogRatio),
                      as.integer(l),
                      PACKAGE="GLAD")



    profileChr$profileValues[,c("Region","Breakpoints","NextLogRatio","OutliersAws")] <- updateBkpRL[c("Region","Breakpoints","NextLogRatio","OutliersAws")]

            
    profileChr$profileValues$Level <- profileChr$profileValues$Region
    profileChr <- detectOutliers(profileChr, region="Level", verbose=verbose, msize=msize, alpha=alpha)
    if (verbose) print("removeLevel: ending function")

    return(profileChr)


  }



