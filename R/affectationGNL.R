# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

affectationGNL <- function(...)
{
  UseMethod("affectationGNL")
}


affectationGNL.profileCGH <- function(profileCGH, alpha=0.001, region="Region", verbose=FALSE, ...)
{
  if (verbose) print("affectationGNL: starting function")
  CGH <- profileCGH$profileValues
  indexna <- attr(na.omit(CGH[,c("Chromosome","LogRatio","PosOrder","ZoneGen")]),"na.action")
  if (!is.null(indexna))
    {
      CGHna <- CGH[indexna,]
      CGH <- CGH[-indexna,]
    }

  
  labelZone <- unique(CGH$ZoneGen)
  medianZone <- rep(NULL, length(labelZone))
  GNL <- medianZone


#i <- 1	
  for (i in 1:length(labelZone))	
    {	
      indexZone <- which(CGH$ZoneGen==labelZone[i]&CGH$OutliersTot==0)	
      medianZone[i] <- median(CGH$LogRatio[indexZone])	
    }	

  medianZone2 <- medianZone*medianZone
  normal <- medianZone[which(medianZone2==min(medianZone2))]


#i <- 1
  for (i in 1:length(labelZone))	
    {
      if (medianZone[i]==normal)
        {
          GNL[i] <- 0
        }	
      
      else
        {
          if (medianZone[i]>normal)
            {
              GNL[i] <- 1
            }
          
          else
            {
              GNL[i] <- -1
            }
        }
    }

  zoneGNL <- rep(NULL, length(CGH$ZoneGen))
#i <- 1
  for (i in 1:length(CGH$ZoneGen))
    {
      indexZone <- which(CGH$ZoneGen==labelZone[i])
      zoneGNL[indexZone] <- GNL[i]
    }
  
  CGH$ZoneGNL <- zoneGNL

###################################################################
#
#      traitement particulier pour les outliers
#
###################################################################


  indexout <- which(CGH$OutliersTot!=0)

# calcul de la moyenne des ratios pour chacun des statuts (G,N,L)
#CGH[-indexout,"LogRatio"]
#CGH[-indexout,"ZoneGNL"]

  
  if (length(indexout)>0)
    {
      meanzone <- aggregate(CGH[-indexout,"LogRatio"],list(ZoneGNL=CGH[-indexout,"ZoneGNL"]),mean)
      colnames(meanzone) <- c("ZoneGNL","mean")


      stdzone <- aggregate(CGH[-indexout,"LogRatio"],list(ZoneGNL=CGH[-indexout,"ZoneGNL"]),var)
      stdzone$x <- stdzone$x^.5
      colnames(stdzone) <- c("ZoneGNL","std")

      statnormal <- merge(meanzone,stdzone)
      statnormal$ZoneGNL <- as.numeric(as.character(statnormal$ZoneGNL))
      statnormal <- statnormal[which(statnormal$ZoneGNL==0),]
      


      

      CGHout <- CGH[indexout,]
      ZoneGNLaux <- CGHout$ZoneGNL


      
      CGHout$ZoneGNL <- CGHout$OutliersTot

# intervalle de confiance pour les régions normales
      seuilinf <- statnormal$mean - statnormal$std*qnorm(1-alpha/2)
      seuilsup <- statnormal$mean + statnormal$std*qnorm(1-alpha/2)



###################################################################################
#
#      cas des outliers situés dans les régions gagnées et outlier -1
#
###################################################################################



# ceux qui sont perdus
      indexgainmoins <-  which(ZoneGNLaux==1&CGHout$OutliersTot==-1&CGHout$LogRatio<seuilinf)

      if (length(indexgainmoins)>0)
        {
          CGHout[indexgainmoins,"ZoneGNL"] <- -1

        }
      

      

# ceux qui sont normaux
      indexgainmoins <-  which(ZoneGNLaux==1&CGHout$OutliersTot==-1&CGHout$LogRatio>=seuilinf&CGHout$LogRatio<=seuilsup)

      if (length(indexgainmoins)>0)
        {
          CGHout[indexgainmoins,"ZoneGNL"] <- 0

        }


      
# ceux qui sont outliers mais quand meme gagnés
      indexgainmoins <-  which(ZoneGNLaux==1&CGHout$OutliersTot==-1&CGHout$LogRatio>seuilsup)

      if (length(indexgainmoins)>0)
        {
          CGHout[indexgainmoins,"ZoneGNL"] <- 1

        }
      

###################################################################################
#
#     cas des outliers situés dans les régions perdues et outlier +1
#
###################################################################################

# ceux qui sont gagnés
      indexlostplus <-  which(ZoneGNLaux==-1&CGHout$OutliersTot==1&CGHout$LogRatio>seuilsup)

      if (length(indexlostplus)>0)
        {
          CGHout[indexlostplus,"ZoneGNL"] <- 1

        }

# ceux qui sont normaux
      indexlostplus <-  which(ZoneGNLaux==-1&CGHout$OutliersTot==1&CGHout$LogRatio<=seuilsup&CGHout$LogRatio>=seuilinf)

      if (length(indexlostplus)>0)
        {
          CGHout[indexlostplus,"ZoneGNL"] <- 0

        }
      

# ceux qui sont des outliers mais perdus quand meme
      indexlostplus <-  which(ZoneGNLaux==-1&CGHout$OutliersTot==1&CGHout$LogRatio<seuilinf)

      if (length(indexlostplus)>0)
        {
          CGHout[indexlostplus,"ZoneGNL"] <- -1

        }
      
      CGH[indexout,] <- CGHout
      
      
    }





#######################################################################
#
# il faut maintenant vérifier le statuts des OutliersAws car il peut
# y avoir des incohérences : en effet, des OutliersAws situés dans une
# région de GAIN mais dont le Level correspond à celui d'une région
# NORMALE seront considérés comme PERDUS
# 
#######################################################################

  if (region=="Level")
    {
      

      labelLevel <- sort(unique(CGH$Level))
      for (i in 1:length(labelLevel))
        {
          indexLevel <- which(CGH$Level==labelLevel[i])
          subset <- CGH[indexLevel,]
          indexoutAws <- which(subset$OutliersAws!=0)
          indexoutTot <- which(subset$OutliersTot!=0)
          if (length(indexoutAws)>0)
            {
              statutGNL <- unique(subset[-indexoutTot,"ZoneGNL"])
              subset[indexoutAws,"ZoneGNL"] <- statutGNL
              CGH[indexLevel,] <- subset
            }
        }
    }
  if (!is.null(indexna))
    {
      CGHna <- data.frame(CGHna,ZoneGNL=rep(NA,length(CGHna$LogRatio)))
      CGH <- rbind(CGH, CGHna)
    }
  profileCGH$profileValues <- CGH
  if (verbose) print("affectationGNL: ending function")
  return(profileCGH)
  
}
