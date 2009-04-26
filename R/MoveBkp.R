### This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr

MoveBkp <- function(profileCGH, ...)
  {
    UseMethod("MoveBkp")
  }

MoveBkp.profileCGH <- function(profileCGH, region="Level", assignGNLOut=TRUE,...)
  {
    

    if (is.data.frame(profileCGH$BkpInfo))
      {

        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
        RecomputeBkpInfo <- FALSE

        ZoneGNLAux <- profileCGH$profileValues$ZoneGNL

        l <- length(profileCGH$profileValues[,1])

        updateGNL <- .C("updateGNL",
                        ZoneGNL=as.integer(profileCGH$profileValues$ZoneGNL),
                        as.double(profileCGH$profileValues$Smoothing),
                        as.integer(profileCGH$profileValues$OutliersTot),
                        as.integer(l),
                        PACKAGE="GLAD")
        

        profileCGH$profileValues$ZoneGNL <- updateGNL$ZoneGNL




################################################################################
###
### Déplacement des Bkp avec MoveBkp!=0
###
################################################################################        


        indexMoveBkp <- which(profileCGH$BkpInfo$MoveBkp!=0)

        if (length(indexMoveBkp)>0)
          {

            RecomputeBkpInfo <- TRUE
            subBkpInfo <- profileCGH$BkpInfo[indexMoveBkp,]

            l <- length(subBkpInfo[,1])


            
            res <- .C("loopMoveBkp",
                      as.integer(subBkpInfo$MoveBkp),
                      as.integer(subBkpInfo$PosOrder),
                      Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                      OutliersTot=as.integer(profileCGH$profileValues$OutliersTot),
                      OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                      OutliersMad=as.integer(profileCGH$profileValues$OutliersMad),
                      Level=as.integer(profileCGH$profileValues[,region]),
                      Smoothing=as.double(profileCGH$profileValues$Smoothing),
                      ZoneGNL=as.integer(profileCGH$profileValues$ZoneGNL),
                      as.integer(l),
                      PACKAGE="GLAD")

##             profileCGH$profileValues$Breakpoints <- res$Breakpoints
##             profileCGH$profileValues$OutliersTot <- res$OutliersTot
##             profileCGH$profileValues$OutliersAws <- res$OutliersAws
##             profileCGH$profileValues$OutliersMad <- res$OutliersMad
##             profileCGH$profileValues[,region] <- res$Level
##             profileCGH$profileValues$Smoothing <- res$Smoothing
##             profileCGH$profileValues$ZoneGNL <- res$ZoneGNL

            
            profileCGH$profileValues[,c("Breakpoints",
                                        "OutliersTot",
                                        "OutliersAws",
                                        "OutliersMad",
                                        region,
                                        "Smoothing",
                                        "ZoneGNL")] <- res[c("Breakpoints",
                                                             "OutliersTot",
                                                             "OutliersAws",
                                                             "OutliersMad",
                                                             "Level",
                                                             "Smoothing",
                                                             "ZoneGNL")]

          }



        if (RecomputeBkpInfo)
          {
        

            l <- length(profileCGH$profileValues[,1])

            updateLevel <- .C("updateLevel",
                              as.integer(profileCGH$profileValues$Chromosome),
                              Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                              Level=as.integer(profileCGH$profileValues[,region]),
                              as.integer(profileCGH$profileValues$PosOrder),
                              NextLogRatio=as.double(rep(0,l)),
                              as.double(profileCGH$profileValues$LogRatio),
                              as.integer(max(profileCGH$profileValues$Level)),
                              as.integer(l),
                              PACKAGE="GLAD")


            profileCGH$profileValues[,region] <- updateLevel$Level
            profileCGH$profileValues$Breakpoints <- updateLevel$Breakpoints
            profileCGH$profileValues$NextLogRatio <- updateLevel$NextLogRatio

            

            
            updateOutliers <- .C("updateOutliersMoveBkp",
                                 OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                                 as.integer(profileCGH$profileValues$OutliersTot),
                                 Level=as.integer(profileCGH$profileValues[,region]),
                                 Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                                 Smoothing=as.double(profileCGH$profileValues$Smoothing),
                                 ZoneGNL=as.integer(profileCGH$profileValues$ZoneGNL),
                                 as.integer(l),
                                 PACKAGE="GLAD")
            
            profileCGH$profileValues[,region] <- updateOutliers$Level
            profileCGH$profileValues$Breakpoints <- updateOutliers$Breakpoints
            profileCGH$profileValues$OutliersAws <- updateOutliers$OutliersAws
            profileCGH$profileValues$Smoothing <- updateOutliers$Smoothing
            profileCGH$profileValues$ZoneGNL <- updateOutliers$ZoneGNL

            
            
            profileCGH$BkpInfo <- BkpInfo(profileCGH, order=FALSE)

            class(profileCGH) <- "profileChr"
            profileCGH <- detectOutliers(profileCGH, region=region, alpha=profileCGH$alpha, msize=profileCGH$msize)

### Mise à jour du GNL des Outliers
            class(profileCGH) <- "profileCGH"
            if(assignGNLOut)
              {
                profileCGH <- OutliersGNL(profileCGH, alpha=profileCGH$alpha, sigma=profileCGH$SigmaG$Value, NormalRef=profileCGH$NormalRef, amplicon=profileCGH$amplicon, deletion=profileCGH$deletion)

              }
          }

        else
          {
            profileCGH$profileValues$ZoneGNL <- ZoneGNLAux
          }
        

      }




    return(profileCGH)

    
  }


