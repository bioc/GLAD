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
    

    CGH <- profileChr$profileValues
	
        Median <- CGH$LogRatio
        MAD <- rep(0, length(CGH$LogRatio))
        CGH$OutliersMad <- MAD	
        labelRegion <- unique(CGH[,region])	
        
        alpha <- qnorm(1-alpha/2)
        

        
        for (i in 1: length(labelRegion))
          {                
            indexRegion <- which(CGH[,region]==labelRegion[i]&CGH$OutliersAws==0)
            
            if(length(indexRegion)>=msize)
              {
                MAD[indexRegion] <- mad(CGH$LogRatio[indexRegion])
                Median[indexRegion] <- median(CGH$LogRatio[indexRegion])	
              }              
          }


        OutPlus <- CGH$LogRatio >(Median + alpha*MAD)
        OutMoins <- CGH$LogRatio <(Median - alpha*MAD)
        CGH$OutliersMad[which(OutPlus==TRUE)] <- 1
        CGH$OutliersMad[which(OutMoins==TRUE)] <- -1
        ### vérification des ouliers AWS: dans le nouveau glad
        ### certains ouliers AWS doivent être supprimés
        OutMad <- OutPlus + OutMoins
        OutAwsToDel <- which(OutMad==0 & CGH$OutliersAws!=0)
        CGH$OutliersAws[OutAwsToDel] <- 0
        
        CGH$OutliersTot <- CGH$OutliersAws + CGH$OutliersMad

        profileChr$profileValues <- CGH
        if(verbose) print("detectOutliers: ending function")
        return(profileChr)
        

    

  }
