# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

affectationGNL <- function(...)
{
  UseMethod("affectationGNL")
}


affectationGNL.profileCGH <- function(profileCGH, alpha=0.001, verbose=FALSE, assignGNLOut=TRUE, ...)
{
  if (verbose) print("affectationGNL: starting function")
  CGH <- profileCGH$profileValues


  indexout <- which(CGH$OutliersTot==0)
  CGHaux <- CGH[indexout,]
  aggZone <- aggregate(CGHaux$LogRatio, list(ZoneGen=CGHaux$ZoneGen), median)
  names(aggZone) <- c("ZoneGen","Median")
  aggZone$ZoneGen <- as.numeric(as.character(aggZone$ZoneGen))

  aggZone$MedianSquare <- aggZone$Median*aggZone$Median

  MinMedian <- aggZone$Median[which(min(aggZone$MedianSquare)==aggZone$MedianSquare)][1]
  aggZone$MinMedian <- aggZone$Median
  aggZone$MinMedian <- MinMedian


  aggZone$ZoneGNL <- aggZone$ZoneGen

  aggZone$ZoneGNL <- 0

  aggZone$ZoneGNL[aggZone$Median>aggZone$MinMedian] <- 1

  aggZone$ZoneGNL[aggZone$Median<aggZone$MinMedian] <- -1


  lengthDest <- length(CGH[,1])
  lengthSrc <- length(aggZone$ZoneGNL)
  myzoneGNL <- .C("my_merge_int",
                  as.integer(CGH$ZoneGen),
                  ZoneGNL=integer(lengthDest),
                  as.integer(aggZone$ZoneGen),
                  as.integer(aggZone$ZoneGNL),
                  as.integer(lengthDest),
                  as.integer(lengthSrc),
                  PACKAGE="GLAD")

  CGH$ZoneGNL <- myzoneGNL$ZoneGNL

##  CGH <- merge(CGH,aggZone[,c("ZoneGen","ZoneGNL")], by="ZoneGen")



###################################################################
###
###      traitement particulier pour les outliers
###
###################################################################

  if(assignGNLOut==FALSE) print("GNL will no be assigned for outliers")
  if(assignGNLOut)
    {
      indexout <- which(CGH$OutliersTot!=0)



      if (length(indexout)>0)
        {


          
          CGHnotout <- CGH[-indexout,]

          sagg <- split(CGHnotout$LogRatio, CGHnotout$ZoneGNL)

          Mean <- sapply(sagg,mean)
          std <- sapply(sagg,var)
          std <- std^0.5
          ZoneGNL <- as.numeric(as.character(names(Mean)))
          
          statnormal <- data.frame(ZoneGNL,Mean,std)

          
          
          statnormal <- statnormal[which(statnormal$ZoneGNL==0),]
          

          
          CGHout <- CGH[indexout,]
          ZoneGNLaux <- CGHout$ZoneGNL      
          CGHout$ZoneGNL <- CGHout$OutliersTot

          


### intervalle de confiance pour les régions normales
          seuilinf <- statnormal$Mean - statnormal$std*qnorm(1-alpha/2)
          seuilsup <- statnormal$Mean + statnormal$std*qnorm(1-alpha/2)



###################################################################################
###
###      cas des outliers situés dans les régions gagnées et outlier -1
###
###################################################################################



### ceux qui sont perdus
          indexgainmoins <-  which(ZoneGNLaux==1&CGHout$OutliersTot==-1&CGHout$LogRatio<seuilinf)

          if (length(indexgainmoins)>0)
            {
              CGHout[indexgainmoins,"ZoneGNL"] <- -1

            }
          

          

### ceux qui sont normaux
          indexgainmoins <-  which(ZoneGNLaux==1&CGHout$OutliersTot==-1&CGHout$LogRatio>=seuilinf&CGHout$LogRatio<=seuilsup)

          if (length(indexgainmoins)>0)
            {
              CGHout[indexgainmoins,"ZoneGNL"] <- 0

            }


          
### ceux qui sont outliers mais quand meme gagnés
          indexgainmoins <-  which(ZoneGNLaux==1&CGHout$OutliersTot==-1&CGHout$LogRatio>seuilsup)

          if (length(indexgainmoins)>0)
            {
              CGHout[indexgainmoins,"ZoneGNL"] <- 1

            }
          

###################################################################################
###
###     cas des outliers situés dans les régions perdues et outlier +1
###
###################################################################################

### ceux qui sont gagnés
          indexlostplus <-  which(ZoneGNLaux==-1&CGHout$OutliersTot==1&CGHout$LogRatio>seuilsup)

          if (length(indexlostplus)>0)
            {
              CGHout[indexlostplus,"ZoneGNL"] <- 1

            }

### ceux qui sont normaux
          indexlostplus <-  which(ZoneGNLaux==-1&CGHout$OutliersTot==1&CGHout$LogRatio<=seuilsup&CGHout$LogRatio>=seuilinf)

          if (length(indexlostplus)>0)
            {
              CGHout[indexlostplus,"ZoneGNL"] <- 0

            }
          

### ceux qui sont des outliers mais perdus quand meme
          indexlostplus <-  which(ZoneGNLaux==-1&CGHout$OutliersTot==1&CGHout$LogRatio<seuilinf)

          if (length(indexlostplus)>0)
            {
              CGHout[indexlostplus,"ZoneGNL"] <- -1

            }
          
          CGH[indexout,] <- CGHout
          
          
        }

    }


  profileCGH$profileValues <- CGH
  if (verbose) print("affectationGNL: ending function")
  return(profileCGH)
  
}
