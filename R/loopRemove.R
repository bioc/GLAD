### This function detects chromosomal breakpoints along genome

### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hup� (Institut Curie) 2003
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


    profileChr$profileValues[,c("Region","OutliersAws","OutliersMad", "OutliersTot", "Breakpoints")] <- res[c("Region","OutliersAws","OutliersMad", "OutliersTot", "Breakpoints")]

    return(profileChr)
    
  }



removeBreakpoints <- function(...)
  {
    UseMethod("removeBreakpoints")
  }



removeBreakpoints.profileChr <- function(profileChr, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, msize=5, alpha=0.001,...)
  {


    if (verbose)
      {
        print("removeBreakpoints: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }
    
    

    sigma <- profileChr$findClusterSigma


### Appel de la fonction loopRemove
        profileChr <- loopRemove(profileChr, sigma, lambda=lambda,
                                  type=type, param=param, verbose=verbose, msize=msize, alpha=alpha)


        return(profileChr)


    if (verbose) print("removeBreakpoints: ending function")

  }



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


    ## Breakpoints et OutliersAws sont pass�s ici pour assurer
    ## la compatibili� avec la fonction R removeBreakpoint
    ## qui fait aussi appel � loopRemove mais pas � updateBkpRL
    l <- length(profileChr$profileValues[,1])
    res <- .C("loopRemove",
              as.double(profileChr$profileValues$LogRatio),
              Level=as.integer(profileChr$profileValues$Level),
              OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
              OutliersMad=integer(l),
              OutliersTot=integer(l),
              Breakpoints=as.integer(profileChr$profileValues$Breakpoints),
              as.integer(msize),
              as.double(qnorm(1-alpha/2)),
              as.double(lambda),
              as.double(param["d"]),
              as.double(profileChr$findClusterSigma),
              as.integer(l),
              PACKAGE="GLAD")


    profileChr$profileValues[,c("Level","OutliersAws","OutliersMad", "OutliersTot", "Breakpoints")] <- res[c("Level","OutliersAws","OutliersMad", "OutliersTot", "Breakpoints")]
    

    updateBkpRL <- .C("updateBkpRL",
                      Level=as.integer(profileChr$profileValues$Level),
                      OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
                      Breakpoints=as.integer(profileChr$profileValues$Breakpoints),
                      ## as.integer(profileChr$profileValues$Chromosome), on est forc�ment sur le m�me chromosome
                      as.integer(profileChr$profileValues$PosOrder),
                      NextLogRatio=as.double(profileChr$profileValues$NextLogRatio),
                      as.double(profileChr$profileValues$LogRatio),
                      as.integer(l),
                      PACKAGE="GLAD")



    profileChr$profileValues[,c("Level","Breakpoints","NextLogRatio","OutliersAws")] <- updateBkpRL[c("Level","Breakpoints","NextLogRatio","OutliersAws")]

            
    profileChr <- detectOutliers(profileChr, region="Level", verbose=verbose, msize=msize, alpha=alpha)


    if (verbose) print("removeLevel: ending function")

    return(profileChr)


  }


