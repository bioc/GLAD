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


