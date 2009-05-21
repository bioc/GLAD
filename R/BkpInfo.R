BkpInfo <- function(...)
  {
    UseMethod("BkpInfo")
  }

BkpInfo.profileCGH <- function(profileCGH, order=TRUE, ...)
  {

    if (order)
      {
##        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
      }

    indexBP <- which(profileCGH$profileValues[,"Breakpoints"] == 1)
    if (length(indexBP)>0)
      {
        
        indexBPplusun <- indexBP + 1

        nomchamp <- c("PosOrder","PosBase",
                      "Smoothing","Chromosome","ZoneGNL", "LogRatio",
                      "NextLogRatio", "MinPosOrder", "MaxPosOrder", "Clone")

        nomchamptot <- colnames(profileCGH$profileValues)

        champinter <- intersect(nomchamp,nomchamptot)
        
        BP <- profileCGH$profileValues[indexBP,champinter]

        
        BPlag <- profileCGH$profileValues[indexBPplusun,c("Smoothing","ZoneGNL")]
        
        colnames(BPlag) <- c("Next","ZoneGNLnext")
        BP <- data.frame(BP,BPlag)
### cette jointure ne prend pas trop de temps (0.02s)
        t1 <- system.time(BP <- merge(BP, profileCGH$SigmaC))
##         print("jointure BkpInfo")
##         print(t1)
        BP$Gap <- abs(BP$Smoothing - BP$Next)
        BP$Weight <- BP$Gap
        BP$GNLchange <- 0


        makeBkpInfo <- .C("make_BkpInfo",
                          as.double(BP$Gap),
                          GNLchange = as.integer(BP$GNLchange),
                          as.double(BP$Value),
                          Weight = as.double (BP$Weight),
                          as.integer(BP$ZoneGNL),
                          as.integer(BP$ZoneGNLnext),
                          as.integer(length(BP[,1])),
                          as.double(profileCGH$nbsigma),
                          PACKAGE = "GLAD")


        
        BP$Weight <- makeBkpInfo$Weight
        BP$GNLchange <- makeBkpInfo$GNLchange


        BP$SmoothingNext <- BP$Next
        BP$Sigma <- BP$Value

        
        nomchamp <- c("Clone","PosOrder","PosBase","Chromosome",
                      "Weight","GNLchange",  "Smoothing", "SmoothingNext", "Gap",
                      "Sigma", "LogRatio", "NextLogRatio", "ZoneGNL", "ZoneGNLnext",
                      "MinPosOrder", "MaxPosOrder")
        

        nomchamptot <- names(BP)

        champinter <- intersect(nomchamp,nomchamptot)


        profileCGH$BkpInfo <- BP[,champinter]
        
      }
    else
      {
        profileCGH$BkpInfo <- NA
      }

    return(profileCGH$BkpInfo)

  }
