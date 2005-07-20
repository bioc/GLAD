### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr



MoveBkpStep <- function(...)
  {
    UseMethod("MoveBkpStep")
  }

MoveBkpStep.profileCGH <- function(profileCGH, ...)
  {
    
    fin <- 0
    maxiter <- 3
    nbiter <- 0

    print("deb MoveBkp")
    if (is.data.frame(profileCGH$BkpInfo))
      {

        FieldsBkp <- names(profileCGH$BkpInfo)

        profileCGH$BkpInfo$MoveBkp <- rep(0,length(profileCGH$BkpInfo[,1]))


        profileCGH <- testBkpToMove(profileCGH)
        indexMoveBkp <- which(profileCGH$BkpInfo$MoveBkp!=0)
       

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
            profileCGH <- MoveBkp(profileCGH)

            if (is.data.frame(profileCGH$BkpInfo))
              {
                profileCGH$BkpInfo$MoveBkp <- rep(0,length(profileCGH$BkpInfo[,1]))
                profileCGH <- testBkpToMove(profileCGH)
                BkpInfo <- profileCGH$BkpInfo
                indexMoveBkp <- which(BkpInfo$MoveBkp!=0)
               

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
