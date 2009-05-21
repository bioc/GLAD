# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

## OutliersGNL <- function(...)
## {
##   UseMethod("OutliersGNL")
## }


## OutliersGNL.profileCGH <- function(profileCGH, alpha=0.001, sigma, NormalRef, amplicon, deletion, verbose=FALSE, assignGNLOut=TRUE, ...)
## {
##   if (verbose) print("OutliersGNL: starting function")
##   if(!assignGNLOut)
##     {
##       print("GNL will not be assigned for outliers")
##       return(profileCGH)
##     }

##   ### seuils de détection pour les outliers
##   seuilsup <- NormalRef + sigma*qnorm(1-alpha/2)
##   seuilinf <- NormalRef - sigma*qnorm(1-alpha/2)
  

##   ### Tous les outliers
##   indexout <- which(profileCGH$profileValues$OutliersTot!=0)
##   profileCGH$profileValues$ZoneGNL[indexout] <- 0
  
##   ### Outliers Gain
##   indexoutGain <- which(profileCGH$profileValues$OutliersTot!=0 & (profileCGH$profileValues$LogRatio - NormalRef)>seuilsup)
##   profileCGH$profileValues$ZoneGNL[indexoutGain] <- 1

##   ### Outliers Perte
##   indexoutLost <- which(profileCGH$profileValues$OutliersTot!=0 & (profileCGH$profileValues$LogRatio - NormalRef)<seuilinf)
##   profileCGH$profileValues$ZoneGNL[indexoutLost] <- -1

##   ### Outliers Amplification
##   indexOutAmp <- which(profileCGH$profileValues$OutliersTot!=0 & (profileCGH$profileValues$LogRatio - NormalRef) > amplicon)
##   profileCGH$profileValues$ZoneGNL[indexOutAmp] <- 2

##   ### Outliers Deletion
##   indexOutDel <- which(profileCGH$profileValues$OutliersTot!=0 & (profileCGH$profileValues$LogRatio - NormalRef) < deletion)
##   profileCGH$profileValues$ZoneGNL[indexOutDel] <- -10    

##   checkGain <- FALSE
##   checkLost <- FALSE
##   checkNormal <- FALSE
##   indexNormal <- which(profileCGH$profileValues$ZoneGNL==0 & profileCGH$profileValues$OutliersTot==0)
##   if (length(indexNormal)>0)
##     {
##       maxNormal <- max(profileCGH$profileValues$Smoothing[indexNormal])
##       minNormal <- min(profileCGH$profileValues$Smoothing[indexNormal])
## ###      indMinN <- which(profileCGH$profileValues$Smoothing[indexNormal]==minNormal)
##       checkNormal <- TRUE
##     }


##   ### Il faut vérifier que des Outliers normaux n'aient pas un logratio supérieur au minimum du smoothing
##   ### des régions gagnés (même principe pour les pertes)
##   ###GNL <- unique(profileCGH$profileValues$ZoneGNL)
##   ### A-t-on détecté des gains
##   indexGain <- which(profileCGH$profileValues$ZoneGNL==1 & profileCGH$profileValues$OutliersTot==0)
##   if (length(indexGain)>0)
##     {
##       checkGain <- TRUE
##       minGain <- min(profileCGH$profileValues$Smoothing[indexGain])
##       indexOutGain <- which(profileCGH$profileValues$OutliersTot!=0 & profileCGH$profileValues$LogRatio>minGain & profileCGH$profileValues$ZoneGNL==0)
##       profileCGH$profileValues$ZoneGNL[indexOutGain] <- 1
##     }

  
##   ### A-t-on détecté des pertes
##   indexLost <- which(profileCGH$profileValues$ZoneGNL==-1 & profileCGH$profileValues$OutliersTot==0)
##   if (length(indexLost)>0)
##     {
##       checkLost <- TRUE
##       maxLost <- max(profileCGH$profileValues$Smoothing[indexLost])
##       indexOutLost <- which(profileCGH$profileValues$OutliersTot!=0 & profileCGH$profileValues$LogRatio<maxLost & profileCGH$profileValues$ZoneGNL==0)
##       profileCGH$profileValues$ZoneGNL[indexOutLost] <- -1
##     }


##   print(paste("minNormal=%f\n",minNormal));
##   print(paste("maxNormal=%f\n",maxNormal));
##   print(paste("minGain=%f\n",minGain));
##   print(paste("maxLost=%f\n",maxLost));

  
##   if (checkGain & checkLost)
##     {
##       checkAlert <- FALSE
##       if (maxLost>minGain)
##         {
##           checkAlert <- TRUE
##         }
##       if (checkNormal)
##         {
##           if (maxLost>minNormal)
##             {
##               ind <- which(profileCGH$profileValues$ZoneGNL==0 & profileCGH$profileValues$OutliersTot==0 & profileCGH$profileValues$Smoothing<=maxLost)
##               profileCGH$profileValues$ZoneGNL[ind] <- -1
##               checkAlert <- TRUE
##             }

##           if (minGain<maxNormal)
##             {
##               ind <- which(profileCGH$profileValues$ZoneGNL==0 & profileCGH$profileValues$OutliersTot==0 & profileCGH$profileValues$Smoothing>=minGain)
##               profileCGH$profileValues$ZoneGNL[ind] <- 1
##               checkAlert <- TRUE
##             }
##         }

##       if (checkAlert)
##         {
##           print("In function OutliersGNL: Inconsistency for smoothing values vs. GNL status has been corrected)")
##         }
      
##     }
  
  
##   if (verbose) print("OutliersGNL: ending function")
##   return(profileCGH)
  
## }


OutliersGNL <- function(...)
{
  UseMethod("OutliersGNL")
}

OutliersGNL.profileCGH <- function(profileCGH, alpha=0.001, sigma, NormalRef, amplicon, deletion, verbose=FALSE, assignGNLOut=TRUE, ...)
{

  if (verbose) print("OutliersGNL: starting function")
  if(!assignGNLOut)
    {
      print("GNL will not be assigned for outliers")
      return(profileCGH)
    }

  ### seuils de détection pour les outliers
  seuilsup <- NormalRef + sigma*qnorm(1-alpha/2)
  seuilinf <- NormalRef - sigma*qnorm(1-alpha/2)

  myOutliersGNL <- .C("OutliersGNL",
                      OutliersTot = as.integer(profileCGH$profileValues[,"OutliersTot"]),
                      ZoneGNL = as.integer(profileCGH$profileValues[,"ZoneGNL"]),
                      as.double(profileCGH$profileValues[,"LogRatio"]),
                      as.double(profileCGH$profileValues[,"Smoothing"]),
                      as.double(seuilsup),
                      as.double(seuilinf),
                      as.double(amplicon),
                      as.double(deletion),
                      as.double(NormalRef),
                      as.integer(length(profileCGH$profileValues[,"Smoothing"])),
                      PACKAGE = "GLAD")

  profileCGH$profileValues[,c("OutliersTot", "ZoneGNL")] <- unlist(myOutliersGNL[c("OutliersTot", "ZoneGNL")])
  
  
  if (verbose) print("OutliersGNL: ending function")
  return(profileCGH)
  
}
