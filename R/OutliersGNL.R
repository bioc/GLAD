# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

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

###  CGH <- profileCGH$profileValues

  
  ### dire qu'il faut la BaseLine
  indexout <- which(profileCGH$profileValues$OutliersTot!=0)
  profileCGH$profileValues$ZoneGNL[indexout] <- 0
  seuilsup <- NormalRef + sigma*qnorm(1-alpha/2)
  seuilinf <- NormalRef - sigma*qnorm(1-alpha/2)
  indexoutGain <- which(profileCGH$profileValues$OutliersTot!=0 & profileCGH$profileValues$LogRatio>seuilsup)
  profileCGH$profileValues$ZoneGNL[indexoutGain] <- 1
  indexoutLost <- which(profileCGH$profileValues$OutliersTot!=0 & profileCGH$profileValues$LogRatio<seuilinf)
  profileCGH$profileValues$ZoneGNL[indexoutLost] <- -1

  indexOutAmp <- which(profileCGH$profileValues$OutliersTot!=0 & (profileCGH$profileValues$LogRatio - NormalRef) > amplicon)
  profileCGH$profileValues$ZoneGNL[indexOutAmp] <- 2

  ###print("dans OutliersGNL")
  ###print(profileCGH$profileValues[indexOutAmp,])

  indexOutDel <- which(profileCGH$profileValues$OutliersTot!=0 & (profileCGH$profileValues$LogRatio - NormalRef) < deletion)
  profileCGH$profileValues$ZoneGNL[indexOutDel] <- -10    

  checkGain <- FALSE
  checkLost <- FALSE
  checkNormal <- FALSE
  indexNormal <- which(profileCGH$profileValues$ZoneGNL==0 & profileCGH$profileValues$OutliersTot==0)
  if (length(indexNormal)>0)
    {

      maxNormal <- max(profileCGH$profileValues$Smoothing[indexNormal])
      minNormal <- min(profileCGH$profileValues$Smoothing[indexNormal])
      indMinN <- which(profileCGH$profileValues$Smoothing[indexNormal]==minNormal)
      #print(indMinN)
      ###print(profileCGH$profileValues[indexNormal,][indMinN,])
      checkNormal <- TRUE

      ###print("maxNormal")
      ###print(maxNormal)
      ###print("minNormal")
      ###print(minNormal)
    }


  ### Il faut vérifier que des Outliers normaux n'aient pas un logratio supérieur au minimum du smoothing
  ### des régions gagnés (même principe pour les pertes)
  ###GNL <- unique(profileCGH$profileValues$ZoneGNL)
  ### A-t-on détecté des gains
  indexGain <- which(profileCGH$profileValues$ZoneGNL==1 & profileCGH$profileValues$OutliersTot==0)
  if (length(indexGain)>0)
    {
      checkGain <- TRUE
      minGain <- min(profileCGH$profileValues$Smoothing[indexGain])
      ###print("G")
      ###print(minGain)
      #ind1 <- which(profileCGH$profileValues$Smoothing==min(profileCGH$profileValues$Smoothing[indexGain]) & profileCGH$profileValues$ZoneGNL==1)
      #print(profileCGH$profileValues[ind1,c("PosOrder","Smoothing","Level","Breakpoints","ZoneGNL","OutliersTot")])
      indexOutGain <- which(profileCGH$profileValues$OutliersTot!=0 & profileCGH$profileValues$LogRatio>minGain & profileCGH$profileValues$ZoneGNL==0)
      profileCGH$profileValues$ZoneGNL[indexOutGain] <- 1
      ###print("pour les gains")
      ###print(profileCGH$profileValues[indexOutGain,])
    }

  
  ### A-t-on détecté des pertes
  indexLost <- which(profileCGH$profileValues$ZoneGNL==-1 & profileCGH$profileValues$OutliersTot==0)
  if (length(indexLost)>0)
    {
      checkLost <- TRUE
      maxLost <- max(profileCGH$profileValues$Smoothing[indexLost])
      ###print("L")
      ###print(maxLost)
      ###ind2 <- which(profileCGH$profileValues$Smoothing==max(profileCGH$profileValues$Smoothing[indexLost]) & profileCGH$profileValues$ZoneGNL==-1)
      ###print(profileCGH$profileValues[ind2,c("PosOrder","Smoothing","Level","Breakpoints","ZoneGNL","OutliersTot", "Chromosome")])
      indexOutLost <- which(profileCGH$profileValues$OutliersTot!=0 & profileCGH$profileValues$LogRatio<maxLost & profileCGH$profileValues$ZoneGNL==0)
      profileCGH$profileValues$ZoneGNL[indexOutLost] <- -1
      ###print("pour les pertes")
      ###print(profileCGH$profileValues[indexOutLost,])
    }

  if (checkGain & checkLost)
    {
      checkAlert <- FALSE
      if (maxLost>minGain)
        {
          checkAlert <- TRUE
###          print("A1")
        }
      if (checkNormal)
        {
          if (maxLost>minNormal)
            {
              ind <- which(profileCGH$profileValues$ZoneGNL==0 & profileCGH$profileValues$OutliersTot==0 & profileCGH$profileValues$Smoothing<=maxLost)
              profileCGH$profileValues$ZoneGNL[ind] <- -1
              checkAlert <- TRUE
###              print("A2")
            }

          if (minGain<maxNormal)
            {

              ind <- which(profileCGH$profileValues$ZoneGNL==0 & profileCGH$profileValues$OutliersTot==0 & profileCGH$profileValues$Smoothing>=minGain)
              profileCGH$profileValues$ZoneGNL[ind] <- 1
              checkAlert <- TRUE
###              print("A3")
            }
        }

      if (checkAlert)
        {
          print("In function OutliersGNL: Inconsistency for smoothing values vs. GNL status has been corrected)")
        }
      
    }
  
  
###  profileCGH$profileValues <- CGH
  if (verbose) print("OutliersGNL: ending function")
  return(profileCGH)
  
}
