# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: bioinfo-staff@curie.fr
# It is strictly forbidden to transfer, use or re-use this code 
# or part of it without explicit written authorization from Institut Curie.

detectOutliers <- function(...)
{
	UseMethod("detectOutliers")
}

detectOutliers.profileChr <- function(profileChr, region="Region", alpha=0.001, ...)
{


  if (alpha>1 || alpha <0)stop("alpha must be setted between 0 and 1")
	

  CGH <- profileChr$profileValues
  indexNA <- attr(na.omit(CGH[,c("Chromosome","LogRatio","PosOrder")]),"na.action")
  if (length(indexNA)<length(CGH[,1]))
    {
      CGHna <- CGH[indexNA,]

      if (!is.null(indexNA))CGH <- CGH[-indexNA,]
	
      Median <- CGH$LogRatio
      MAD <- rep(0, length(CGH$LogRatio))
      CGH$OutliersMad <- MAD	
      #regionAux <- MAD
      labelRegion <- unique(CGH[,region])	
      #MAD[which(CGH$OutliersAws!=0)] <- 0
			
      alpha <- qnorm(1-alpha/2)
	      

      if (region=="Level")
        {          
          # attention ici le terme region désigne le Level

          
          for (i in 1: length(labelRegion))
            {
              # recherche des OutliersAws qui ne sont pas seuls dans leur Level
              indexL <- which(CGH[,region]==labelRegion[i])
              indexLOut <- which(CGH[,region]==labelRegion[i]&CGH$OutliersAws!=0)

              if (length(indexL) > length(indexLOut))
                {
                  indexRegion <- indexL
                 
                  if(length(indexRegion)>1)
                    {
                      MAD[indexRegion] <- mad(CGH$LogRatio[indexRegion])
                      Median[indexRegion] <- median(CGH$LogRatio[indexRegion])	
                    } 
                }                                 
            }


          OutPlus <- CGH$LogRatio >(Median + alpha*MAD)
          OutMoins <- CGH$LogRatio <(Median - alpha*MAD)
          CGH$OutliersMad[which(OutPlus==TRUE)] <- 1
          CGH$OutliersMad[which(OutMoins==TRUE)] <- -1

          # il faut regarder si certains OutliersAws passent
          # OutliersMad : ces outliers là deviennent OutliersMad

          indexOutAwsMad <- which(CGH$OutliersMad!=0&CGH$OutliersAws!=0)
          if (!is.null(indexOutAwsMad))
            {
              CGH$OutliersAws[indexOutAwsMad] <- 0
            }
          
          
          CGH$OutliersTot <- CGH$OutliersAws + CGH$OutliersMad
          CGHna$OutliersMad <- rep(NA, length(CGHna$PosOrder))
          CGHna$OutliersTot <- CGHna$OutliersMad
          CGH <- rbind(CGH, CGHna)

          profileChr$profileValues <- CGH

          return(profileChr)
        }

      else
        {
          for (i in 1: length(labelRegion))
            {                
              indexRegion <- which(CGH[,region]==labelRegion[i]&CGH$OutliersAws==0)
              
              if(length(indexRegion)>1)
                {
                  MAD[indexRegion] <- mad(CGH$LogRatio[indexRegion])
                  Median[indexRegion] <- median(CGH$LogRatio[indexRegion])	
                }              
            }


          OutPlus <- CGH$LogRatio >(Median + alpha*MAD)
          OutMoins <- CGH$LogRatio <(Median - alpha*MAD)
          CGH$OutliersMad[which(OutPlus==TRUE)] <- 1
          CGH$OutliersMad[which(OutMoins==TRUE)] <- -1
          
          CGH$OutliersTot <- CGH$OutliersAws + CGH$OutliersMad
          CGHna$OutliersMad <- rep(NA, length(CGHna$PosOrder))
          CGHna$OutliersTot <- CGHna$OutliersMad
          CGH <- rbind(CGH, CGHna)

          profileChr$profileValues <- CGH

          return(profileChr)
        }

    }


  else
    {
      CGH$OutliersMad <- CGH$OutliersTot <- rep(NA,length(CGH[,1]))
      profileChr$profileValues <- CGH
      return(profileChr)
    }    


        
	 
}
