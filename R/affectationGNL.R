### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr

affectationGNL <- function(...)
{
  UseMethod("affectationGNL")
}


affectationGNL.profileCGH <- function(profileCGH, alpha=0.001, verbose=FALSE, assignGNLOut=TRUE, ...)
{
  if (verbose) print("affectationGNL: starting function")
  CGH <- profileCGH$profileValues

  indexout <- which(CGH[["OutliersTot"]] == 0)
                                        #  CGHaux <- CGH[indexout,]
  aggZone <- aggregate(CGH[["LogRatio"]][indexout], list(ZoneGen=CGH[["ZoneGen"]][indexout]), median)
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


  lengthDest <- length(CGH[[1]])
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

  if(assignGNLOut == FALSE) print("GNL will no be assigned for outliers")
  if(assignGNLOut)
    {

      indexout <- which(CGH[["OutliersTot"]] != 0)


      if (length(indexout) > 0)
        {
          
                                        #          CGHnotout <- CGH[-indexout,]

          sagg <- split(CGH[["LogRatio"]][-indexout], CGH[["ZoneGNL"]][-indexout])

          Mean <- sapply(sagg,mean)
          std <- sapply(sagg,var)
          std <- std^0.5
          ZoneGNL <- as.numeric(as.character(names(Mean)))
          
          statnormal <- data.frame(ZoneGNL,Mean,std)

          
          
          statnormal <- statnormal[which(statnormal$ZoneGNL==0),]
          

          
                                        #          CGHout <- CGH[indexout,]
          ZoneGNLaux <- CGH[["ZoneGNL"]][indexout]      
          CGH[["ZoneGNL"]][indexout] <- CGH[["OutliersTot"]][indexout]

          


### intervalle de confiance pour les régions normales
          seuilinf <- statnormal$Mean - statnormal$std*qnorm(1-alpha/2)
          seuilsup <- statnormal$Mean + statnormal$std*qnorm(1-alpha/2)


###################################################################################
###
###      cas des outliers situés dans les régions gagnées et outlier -1
###
###################################################################################



### ceux qui sont perdus
          indexgainmoins <-  which(ZoneGNLaux == 1 &
                                   CGH[["OutliersTot"]][indexout] == -1 &
                                   CGH[["LogRatio"]][indexout] < seuilinf)

          if (length(indexgainmoins) > 0)
            {
              CGH[["ZoneGNL"]][indexout][indexgainmoins] <- -1

            }
          

          

### ceux qui sont normaux
          indexgainmoins <-  which(ZoneGNLaux == 1 &
                                   CGH[["OutliersTot"]][indexout] == -1 &
                                   CGH[["LogRatio"]][indexout] >= seuilinf &
                                   CGH[["LogRatio"]][indexout] <= seuilsup)

          if (length(indexgainmoins) > 0)
            {
              CGH[["ZoneGNL"]][indexout][indexgainmoins] <- 0

            }


          
### ceux qui sont outliers mais quand meme gagnés
          indexgainmoins <-  which(ZoneGNLaux == 1 &
                                   CGH[["OutliersTot"]][indexout] == -1 &
                                   CGH[["LogRatio"]][indexout] > seuilsup)

          if (length(indexgainmoins) > 0)
            {
              CGH[["ZoneGNL"]][indexout][indexgainmoins] <- 1

            }



          
###################################################################################
###
###     cas des outliers situés dans les régions perdues et outlier +1
###
###################################################################################

### ceux qui sont gagnés
          indexlostplus <-  which(ZoneGNLaux == -1 &
                                  CGH[["OutliersTot"]][indexout] == 1 &
                                  CGH[["LogRatio"]][indexout] > seuilsup)

          if (length(indexlostplus) > 0)
            {
              CGH[["ZoneGNL"]][indexout][indexlostplus] <- 1

            }

          
### ceux qui sont normaux
          indexlostplus <-  which(ZoneGNLaux== -1 &
                                  CGH[["OutliersTot"]][indexout] == 1 &
                                  CGH[["LogRatio"]][indexout] <= seuilsup &
                                  CGH[["LogRatio"]][indexout] >= seuilinf)

          if (length(indexlostplus)>0)
            {
              CGH[["ZoneGNL"]][indexout][indexlostplus] <- 0

            }
          

          
### ceux qui sont des outliers mais perdus quand meme
          indexlostplus <-  which(ZoneGNLaux == -1 &
                                  CGH[["OutliersTot"]][indexout] == 1 &
                                  CGH[["LogRatio"]][indexout] < seuilinf)

          if (length(indexlostplus) > 0)
            {
              CGH[["ZoneGNL"]][indexout][indexlostplus] <- -1

            }
          
        }

    }


  profileCGH$profileValues <- CGH
  if (verbose) print("affectationGNL: ending function")
  return(profileCGH)
  
}
