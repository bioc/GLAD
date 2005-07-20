### This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr



testBkpToMove <- function(...)
  {
    UseMethod("testBkpToMove")
  }


testBkpToMove.profileCGH <- function(profileCGH, ...)
  {
    BkpInfo <- profileCGH$BkpInfo
    testSingle <- function(LogRatio, NextLogRatio, Smoothing, SmoothingNext)
      {
        moveBkp <- 0
### le créneau est plus bas à droite qu'à gauche
        if (Smoothing > SmoothingNext)
          {
            if (SmoothingNext <= LogRatio & LogRatio <= Smoothing)
              {
                if((LogRatio-SmoothingNext) < (Smoothing-LogRatio))
                  {
### il faut déplacer le Bkp vers la gauche
                    moveBkp <- -1
                  }                        
              }
            
            if (SmoothingNext <= NextLogRatio & NextLogRatio <= Smoothing)
              {
                if ( (NextLogRatio-SmoothingNext)>(Smoothing-NextLogRatio))
                  {
### il faut déplacer le Bkp vers la droite
                    moveBkp <- 1
                  }
                
              }

            if (LogRatio <= SmoothingNext)
              {
                moveBkp <- -1
              }

            if (NextLogRatio>=Smoothing)
              {
                moveBkp <- 1
              }
          }
### le créneau est plus bas à gauche qu'à droite
        else
          {
            if (SmoothingNext >= LogRatio & LogRatio >= Smoothing)
              {
                if ((SmoothingNext-LogRatio) < (LogRatio - Smoothing))
                  {
### il faut déplacer le Bkp vers la gauche
                    moveBkp <- -1
                  }
              }

            if (SmoothingNext >= NextLogRatio & NextLogRatio >= Smoothing)
              {
                if ((SmoothingNext-NextLogRatio) > (NextLogRatio-Smoothing))
                  {
### il faut déplacer le Bkp vers la droite
                    moveBkp <- 1
                  }
              }

            if (LogRatio>=SmoothingNext)
              {
                moveBkp <- -1
              }

            if (NextLogRatio<=Smoothing)
              {
                moveBkp <- 1
              }
          }
        return(moveBkp)

      }

    for (NbBkp in 1:length(profileCGH$BkpInfo[,1]))
      {
        BkpInfo$MoveBkp[NbBkp] <- testSingle(BkpInfo$LogRatio[NbBkp], BkpInfo$NextLogRatio[NbBkp],BkpInfo$Smoothing[NbBkp], BkpInfo$SmoothingNext[NbBkp])
      }
    profileCGH$BkpInfo <- BkpInfo


    profileCGH$BkpInfo$NextPosOrder <- profileCGH$BkpInfo$PosOrder+1
    profileCGH$BkpInfo$BeforePosOrder <- profileCGH$BkpInfo$PosOrder-1
### on vérifie qu'on ne déplace pas des Bkp au niveau des extrémités
    indexRight <- which(profileCGH$BkpInfo$MoveBkp==1 & profileCGH$BkpInfo$NextPosOrder==profileCGH$BkpInfo$MaxPosOrder)
    profileCGH$BkpInfo$MoveBkp[indexRight] <- 0
    indexLeft <- which(profileCGH$BkpInfo$MoveBkp==-1 & profileCGH$BkpInfo$BeforePosOrder==profileCGH$BkpInfo$MinPosOrder)
    profileCGH$BkpInfo$MoveBkp[indexLeft] <- 0

    profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-NextPosOrder)
    profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-BeforePosOrder)
    return(profileCGH)
    
  }


