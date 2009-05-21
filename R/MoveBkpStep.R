### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr



MoveBkpStep <- function(...)
  {
    UseMethod("MoveBkpStep")
  }

MoveBkpStep.profileCGH <- function(profileCGH, assignGNLOut=TRUE,...)
  {
    
    fin <- 0
    maxiter <- 2
    nbiter <- 0

    print("Check Breakpoints Position")
    if (is.data.frame(profileCGH$BkpInfo))
      {

        FieldsBkp <- names(profileCGH$BkpInfo)

        profileCGH$BkpInfo$MoveBkp <- 0


        profileCGH <- testBkpToMove(profileCGH)
        indexMoveBkp <- which(profileCGH$BkpInfo$MoveBkp != 0)
       

        if (length(indexMoveBkp)>0)
          {
            fin <- 0
          }
        else
          {
            fin <- 1
          }

        while (fin!=1)
          {
            nbiter <- nbiter + 1
            profileCGH <- MoveBkp(profileCGH, assignGNLOut = assignGNLOut)

            if (is.data.frame(profileCGH$BkpInfo))
              {
                profileCGH$BkpInfo$MoveBkp <- 0
                profileCGH <- testBkpToMove(profileCGH)
                indexMoveBkp <- which(profileCGH$BkpInfo$MoveBkp != 0)
               

                if (length(indexMoveBkp)>0)
                  {
                    fin <- 0
                  }
                else
                  {
                    fin <- 1
                  }
                if(nbiter>maxiter)
                  {
                    fin <- 1
                  }
              }
            else
              {
                fin <- 1
              }
          }
        
        profileCGH$BkpInfo <- profileCGH$BkpInfo[,FieldsBkp]
        
      }

    return(profileCGH)
    
  }



testBkpToMove <- function(...)
  {
    UseMethod("testBkpToMove")
  }


testBkpToMove.profileCGH <- function(profileCGH, ...)
  {
##     BkpInfo <- profileCGH$BkpInfo
##     testSingle <- function(LogRatio, NextLogRatio, Smoothing, SmoothingNext)
##       {
##         moveBkp <- 0
## ### le créneau est plus bas à droite qu'à gauche
##         if (Smoothing > SmoothingNext)
##           {
##             if (SmoothingNext <= LogRatio & LogRatio <= Smoothing)
##               {
##                 if((LogRatio-SmoothingNext) < (Smoothing-LogRatio))
##                   {
## ### il faut déplacer le Bkp vers la gauche
##                     moveBkp <- -1
##                   }                        
##               }
            
##             if (SmoothingNext <= NextLogRatio & NextLogRatio <= Smoothing)
##               {
##                 if ( (NextLogRatio-SmoothingNext)>(Smoothing-NextLogRatio))
##                   {
## ### il faut déplacer le Bkp vers la droite
##                     moveBkp <- 1
##                   }
                
##               }

##             if (LogRatio <= SmoothingNext)
##               {
##                 moveBkp <- -1
##               }

##             if (NextLogRatio>=Smoothing)
##               {
##                 moveBkp <- 1
##               }
##           }
## ### le créneau est plus bas à gauche qu'à droite
##         else
##           {
##             if (SmoothingNext >= LogRatio & LogRatio >= Smoothing)
##               {
##                 if ((SmoothingNext-LogRatio) < (LogRatio - Smoothing))
##                   {
## ### il faut déplacer le Bkp vers la gauche
##                     moveBkp <- -1
##                   }
##               }

##             if (SmoothingNext >= NextLogRatio & NextLogRatio >= Smoothing)
##               {
##                 if ((SmoothingNext-NextLogRatio) > (NextLogRatio-Smoothing))
##                   {
## ### il faut déplacer le Bkp vers la droite
##                     moveBkp <- 1
##                   }
##               }

##             if (LogRatio>=SmoothingNext)
##               {
##                 moveBkp <- -1
##               }

##             if (NextLogRatio<=Smoothing)
##               {
##                 moveBkp <- 1
##               }
##           }
##         return(moveBkp)

##       }

##     for (NbBkp in 1:length(profileCGH$BkpInfo[,1]))
##       {
##         BkpInfo$MoveBkp[NbBkp] <- testSingle(BkpInfo$LogRatio[NbBkp], BkpInfo$NextLogRatio[NbBkp],BkpInfo$Smoothing[NbBkp], BkpInfo$SmoothingNext[NbBkp])
##       }
##     profileCGH$BkpInfo <- BkpInfo


    NbBkp <- length(profileCGH$BkpInfo[,1])
    myMoveBkp <- .C("loopTestBkpToMove",
                    as.double(profileCGH$BkpInfo$LogRatio),
                    as.double(profileCGH$BkpInfo$NextLogRatio),
                    as.double(profileCGH$BkpInfo$Smoothing),
                    as.double(profileCGH$BkpInfo$SmoothingNext),
                    as.integer(profileCGH$BkpInfo$PosOrder),
                    as.integer(profileCGH$BkpInfo$MaxPosOrder),
                    as.integer(profileCGH$BkpInfo$MinPosOrder),                                                            
                    MoveBkp=integer(NbBkp),
                    as.integer(NbBkp),
                    PACKAGE="GLAD")
    
    profileCGH$BkpInfo$MoveBkp <- myMoveBkp$MoveBkp    

    return(profileCGH)
    
  }


