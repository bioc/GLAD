BkpInfo <- function(...)
  {
    UseMethod("BkpInfo")
  }

BkpInfo.profileCGH <- function(profileCGH, order=TRUE, ...)
  {

    if (order)
      {
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
      }
    
    indexBP <- which(profileCGH$profileValues$Breakpoints==1)
    if (length(indexBP)>0)
      {
        
        indexBPplusun <- indexBP + 1

        nomchamp <- c("PosOrder","PosBase",
        "Smoothing","Chromosome","ZoneGNL", "LogRatio",
        "NextLogRatio", "MinPosOrder", "MaxPosOrder", "Clone")

        nomchamptot <- names(profileCGH$profileValues)

        champinter <- intersect(nomchamp,nomchamptot)
        

        BP <- profileCGH$profileValues[indexBP,champinter]
        
        BPlag <- profileCGH$profileValues[indexBPplusun,c("Smoothing","ZoneGNL")]
        names(BPlag) <- c("Next","ZoneGNLnext")
        BP <- data.frame(BP,BPlag)
        ### cette jointure ne prend pas trop de temps (0.02s)
        BP <- merge(BP, profileCGH$SigmaC)
        BP$Gap <- abs(BP$Smoothing-BP$Next)
        BP$Weight <- BP$Gap
        BP$GNLchange <- 0


        makeBkpInfo <- .C("make_BkpInfo",
                          as.double(BP$Gap),
                          GNLchange=as.integer(BP$GNLchange),
                          as.double(BP$Value),
                          Weight=as.double (BP$Weight),
                          as.integer(BP$ZoneGNL),
                          as.integer(BP$ZoneGNLnext),
                          as.integer(length(BP[,1])),
                          as.double(profileCGH$nbsigma),
                          PACKAGE="GLAD")

        BP$Weight <- makeBkpInfo$Weight
        BP$GNLchange <- makeBkpInfo$GNLchange
                
        
##         for (i in 1:length(BP[,1]))
##           {
##             BP$Weight[i] <- 1 - kernelpen(BP$Gap[i], param=c(d=profileCGH$nbsigma*BP$Value[i]))
##             if (BP$ZoneGNL[i]==BP$ZoneGNLnext[i])
##               {
##                 BP$GNLchange[i] <- 0
##               }
##             else
##               {
##                 BP$GNLchange[i] <- 1
##               }
##           }


        BP$SmoothingNext <- BP$Next
        BP$Sigma <- BP$Value

        
        nomchamp <- c("Clone","PosOrder","PosBase","Chromosome",
        "Weight","GNLchange",  "Smoothing", "SmoothingNext", "Gap",
        "Sigma", "LogRatio", "NextLogRatio", "ZoneGNL", "ZoneGNLnext",
        "MinPosOrder", "MaxPosOrder")
        

        nomchamptot <- names(BP)

        champinter <- intersect(nomchamp,nomchamptot)

##         if (length(intersect(names(profileCGH$profileValues),"PosBase"))>=1)
##           {
##             profileCGH$BkpInfo <- BP[,c("PosOrder","PosBase", "Chromosome", "Weight","GNLchange", "Smoothing", "Next", "Gap", "Value", "LogRatio", "NextLogRatio", "ZoneGNL", "ZoneGNLnext", "MinPosOrder", "MaxPosOrder")]
##             names(profileCGH$BkpInfo) <- c("PosOrder","PosBase","Chromosome", "Weight","GNLchange",  "Smoothing", "SmoothingNext", "Gap", "Sigma", "LogRatio", "NextLogRatio", "ZoneGNL", "ZoneGNLnext", "MinPosOrder", "MaxPosOrder")
##           }
##         else
##           {
##             profileCGH$BkpInfo <- BP[,c("PosOrder","Chromosome", "Weight","GNLchange", "Smoothing", "Next", "Gap", "Value", "LogRatio", "NextLogRatio", "ZoneGNL", "ZoneGNLnext", "MinPosOrder", "MaxPosOrder")]
##             names(profileCGH$BkpInfo) <- c("PosOrder","Chromosome", "Weight","GNLchange",  "Smoothing", "SmoothingNext", "Gap", "Sigma", "LogRatio", "NextLogRatio", "ZoneGNL", "ZoneGNLnext", "MinPosOrder", "MaxPosOrder")  
##           }

         profileCGH$BkpInfo <- BP[,champinter]
        
      }
    else
      {
        profileCGH$BkpInfo <- NA
      }

    return(profileCGH$BkpInfo)

  }
