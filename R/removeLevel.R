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


### Appel de la fonction loopRemove

    
### 03062005: à quoi cela sert de passer profileChr$profileValues à
###      la fonction loopRemove
    profileChr <- loopRemove(profileChr, sigma, lambda=lambda,
                             type=type, param=param, verbose=verbose, msize=msize, alpha=alpha)


    
    profileChr$profileValues <- profileChr$profileValues[order(profileChr$profileValues$PosOrder),]
    profileChr$profileValues$Breakpoints <- 0
    profileChr$profileValues$OutliersAws <- 0
    profileChr$profileValues$NextLogRatio <- profileChr$profileValues$OutliersAws



    nb <- length(profileChr$profileValues[,1])
    updateBkpRL <- .C("updateBkpRL",
                      Region=as.integer(profileChr$profileValues$Region),
                      OutliersAws=as.integer(profileChr$profileValues$OutliersAws),
                      Breakpoints=as.integer(profileChr$profileValues$Breakpoints),
                      as.integer(profileChr$profileValues$Chromosome),
                      as.integer(profileChr$profileValues$PosOrder),
                      NextLogRatio=as.double(profileChr$profileValues$NextLogRatio),
                      as.double(profileChr$profileValues$LogRatio),
                      as.integer(nb),
                      PACKAGE="GLAD")

    ##     profileChr$profileValues$Region <- updateBkpRL$Region
    ##     profileChr$profileValues$Breakpoints <- updateBkpRL$Breakpoints
    ##     profileChr$profileValues$NextLogRatio <- updateBkpRL$NextLogRatio
    ##     profileChr$profileValues$OutliersAws <- updateBkpRL$OutliersAws

    profileChr$profileValues[,c("Region","Breakpoints","NextLogRatio","OutliersAws")] <- updateBkpRL[c("Region","Breakpoints","NextLogRatio","OutliersAws")]


            
    profileChr$profileValues$Level <- profileChr$profileValues$Region
    profileChr <- detectOutliers(profileChr, region="Level", verbose=verbose, msize=msize, alpha=alpha)
    if (verbose) print("removeLevel: ending function")

    return(profileChr)


  }



