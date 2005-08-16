### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr



filterBkpStep <- function(...)
  {
    UseMethod("filterBkpStep")
  }



filterBkpStep.profileCGH <- function(profileCGH, MinBkpWeight=0.35, ...)
  {


    profileCGH <- filterBkp(profileCGH, MinBkpWeight=MinBkpWeight)
     
    fin <- 0
    maxiter <- 100
    nbiter <- 0
    


    if (is.data.frame(profileCGH$BkpInfo))
      {
        BkpInfo <- profileCGH$BkpInfo
        indexWeightToSmall <- which(BkpInfo$Weight<MinBkpWeight & BkpInfo$GNLchange==0 & BkpInfo$ZoneGNL!=2)
        indexWeightToSmall <- c(indexWeightToSmall, which(BkpInfo$Weight==0))
        indexWeightZero <- which(BkpInfo$Weight==0 & BkpInfo$GNLchange==1)        
        if (length(indexWeightToSmall)>0 || length(indexWeightZero)>0)
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
            profileCGH <- filterBkp(profileCGH, MinBkpWeight=MinBkpWeight)
            if (is.data.frame(profileCGH$BkpInfo))
              {
                BkpInfo <- profileCGH$BkpInfo
                indexWeightToSmall <- which(BkpInfo$Weight<MinBkpWeight & BkpInfo$GNLchange==0 & BkpInfo$ZoneGNL!=2)
                indexWeightZero <- which(BkpInfo$Weight==0 & BkpInfo$GNLchange==1)
                if (length(indexWeightToSmall)>0 || length(indexWeightZero)>0)
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
                    print("There is something wrong in the loop calling filterBkp")
                  }
              }
            else
              {
                fin <- 1
              }
          }

      }

    return(profileCGH)
    
  }
