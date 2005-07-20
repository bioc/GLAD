# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

OutliersGNL <- function(...)
{
  UseMethod("OutliersGNL")
}


OutliersGNL.profileCGH <- function(profileCGH, alpha=0.001, sigma, NormalRef, amplicon, deletion, verbose=FALSE, ...)
{
  if (verbose) print("OutliersGNL: starting function")

  CGH <- profileCGH$profileValues

  
  ### dire qu'il faut la BaseLine
  indexout <- which(CGH$OutliersTot!=0)
  CGH$ZoneGNL[indexout] <- 0
  seuilsup <- NormalRef + sigma*qnorm(1-alpha/2)
  seuilinf <- NormalRef - sigma*qnorm(1-alpha/2)
  indexoutGain <- which(CGH$OutliersTot!=0 & CGH$LogRatio>seuilsup)
  CGH$ZoneGNL[indexoutGain] <- 1
  indexoutLost <- which(CGH$OutliersTot!=0 & CGH$LogRatio<seuilinf)
  CGH$ZoneGNL[indexoutLost] <- -1

  indexOutAmp <- which(CGH$OutliersTot!=0 & (CGH$LogRatio - NormalRef) > amplicon)
  CGH$ZoneGNL[indexOutAmp] <- 2

  ###print("dans OutliersGNL")
  ###print(CGH[indexOutAmp,])

  indexOutDel <- which(CGH$OutliersTot!=0 & (CGH$LogRatio - NormalRef) < deletion)
  CGH$ZoneGNL[indexOutDel] <- -10    

  checkGain <- FALSE
  checkLost <- FALSE


  ### Il faut vérifier que des Outliers normaux n'aient pas un logratio supérieur au minimum du smoothing
  ### des régions gagnés (même principe pour les pertes)
  GNL <- unique(CGH$ZoneGNL)
  ### A-t-on détecté des gains
  indexGain <- which(CGH$ZoneGNL==1 & CGH$OutliersTot==0)
  if (length(indexGain)>0)
    {
      checkGain <- TRUE
      minGain <- min(CGH$Smoothing[indexGain])
      #print("G")
      #print(minGain)
      #ind1 <- which(CGH$Smoothing==min(CGH$Smoothing[indexGain]) & CGH$ZoneGNL==1)
      #print(CGH[ind1,c("PosOrder","Smoothing","Level","Breakpoints","ZoneGNL","OutliersTot")])
      indexOutGain <- which(CGH$OutliersTot!=0 & CGH$LogRatio>minGain & CGH$ZoneGNL==0)
      CGH$ZoneGNL[indexOutGain] <- 1
      ###print("pour les gains")
      ###print(CGH[indexOutGain,])
    }

  
  ### A-t-on détecté des pertes
  indexLost <- which(CGH$ZoneGNL==-1 & CGH$OutliersTot==0)
  if (length(indexLost)>0)
    {
      checkLost <- TRUE
      maxLost <- max(CGH$Smoothing[indexLost])
      #print("L")
      #print(maxLost)
      #ind2 <- which(CGH$Smoothing==max(CGH$Smoothing[indexLost]) & CGH$ZoneGNL==-1)
      #print(CGH[ind2,c("PosOrder","Smoothing","Level","Breakpoints","ZoneGNL","OutliersTot")])
      indexOutLost <- which(CGH$OutliersTot!=0 & CGH$LogRatio<maxLost & CGH$ZoneGNL==0)
      CGH$ZoneGNL[indexOutLost] <- -1
      ###print("pour les pertes")
      ###print(CGH[indexOutLost,])
    }

  if (checkGain & checkLost)
    {
      checkAlert <- FALSE
      if (maxLost>minGain)
        {
          checkAlert <- TRUE
        }
      if (maxLost>NormalRef)
        {
          checkAlert <- TRUE
        }

      if (minGain<NormalRef)
        {
          checkAlert <- TRUE
        }

      if (checkAlert)
        {
          print("In function OutliersGNL: Alert (Inconsistency for smoothing values vs. Gain and Lost Status)")
        }
      
    }
  
  
  profileCGH$profileValues <- CGH
  if (verbose) print("OutliersGNL: ending function")
  return(profileCGH)
  
}
