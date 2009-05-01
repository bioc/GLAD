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
        updateGNL <- .C("MoveBkp_updateGNL",
                        ZoneGNL=as.integer(profileCGH$profileValues$ZoneGNL),
                        as.double(profileCGH$profileValues$Smoothing),
                        as.integer(profileCGH$profileValues$OutliersTot),
                        as.integer(l),
                        PACKAGE="GLAD")
        

        profileCGH$profileValues$ZoneGNL <- updateGNL$ZoneGNL


        ## ###################################
        ##
        ## Déplacement des Bkp avec MoveBkp!=0
        ##
        ## ###################################

        indexMoveBkp <- which(profileCGH$BkpInfo$MoveBkp!=0)

        if (length(indexMoveBkp)>0)
          {

            
            RecomputeBkpInfo <- TRUE
            subBkpInfo <- profileCGH$BkpInfo[indexMoveBkp,]



            if(region == "Region")
              warning("Level and Region are the same vector: the code may not works")
            
            l <- length(subBkpInfo[,1])            
            res <- .C("MoveBkp_Delete_Bkp",
                      as.integer(subBkpInfo$MoveBkp),
                      as.integer(subBkpInfo$PosOrder),
                      Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                      OutliersTot=as.integer(profileCGH$profileValues$OutliersTot),
                      OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                      OutliersMad=as.integer(profileCGH$profileValues$OutliersMad),
                      Level=as.integer(profileCGH$profileValues[,region]),
                      Region=as.integer(profileCGH$profileValues$Region),
                      Smoothing=as.double(profileCGH$profileValues$Smoothing),
                      ZoneGNL=as.integer(profileCGH$profileValues$ZoneGNL),
                      as.integer(l),
                      PACKAGE="GLAD")

            
            profileCGH$profileValues[,c("Breakpoints",
                                        "OutliersTot",
                                        "OutliersAws",
                                        "OutliersMad",
                                        region,
                                        "Region",
                                        "Smoothing",
                                        "ZoneGNL")] <- res[c("Breakpoints",
                                                             "OutliersTot",
                                                             "OutliersAws",
                                                             "OutliersMad",
                                                             "Level",
                                                             "Region",
                                                             "Smoothing",
                                                             "ZoneGNL")]

            

          }



        if (RecomputeBkpInfo)
          {

            lengthDest <- length(profileCGH$profileValues$Level)
            l <- lengthDest             
            ## recalcul de la smoothing line
            ##             agg <- aggregate(profileCGH$profileValues$LogRatio, list(Level=profileCGH$profileValues$Level), median)
            ##             agg$Level <- as.numeric(as.character(agg$Level))
            ##             names(agg) <- c("Level","Smoothing")


            ##             lengthDest <- length(profileCGH$profileValues$Level)
            ##             lengthSrc <- length(agg$Level)
            ##             mySmoothing <- .C("my_merge",
            ##                               as.integer(profileCGH$profileValues$Level),
            ##                               Smoothing=double(lengthDest),
            ##                               as.integer(agg$Level),
            ##                               as.double(agg$Smoothing),
            ##                               as.integer(lengthDest),
            ##                               as.integer(lengthSrc),
            ##                               PACKAGE="GLAD")

            mySmoothing <- .C("compute_median_smoothing",
                              as.double(profileCGH$profileValues$LogRatio),
                              as.integer(profileCGH$profileValues$Level),
                              Smoothing=double(l),
                              as.integer(l),
                              PACKAGE="GLAD")

            
            profileCGH$profileValues$Smoothing <- mySmoothing$Smoothing

            

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

            ## ICI on peut optimiser pour récupérer les valeurs
            profileCGH$profileValues[,region] <- updateLevel$Level
            profileCGH$profileValues$Breakpoints <- updateLevel$Breakpoints
            profileCGH$profileValues$NextLogRatio <- updateLevel$NextLogRatio

            
            ## ##################
            ## Mise à jour du GNL
            ## ##################

            profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGNL")]
            
            ## on prend comme référence ceux qui sont compris entre certaines valeurs            
            NormalRange <- .C("compute_NormalRange",
                              as.double(profileCGH$profileValues$Smoothing),
                              as.double(profileCGH$NormalRef),
                              as.integer(profileCGH$profileValues$Level),
                              NormalRange=integer(l),
                              as.double(profileCGH$deltaN),
                              as.integer(l),
                              PACKAGE="GLAD")

            profileCGH$profileValues$NormalRange <- NormalRange$NormalRange
            
            
            ## le clustering est fait sur les niveaux NormalRange
            class(profileCGH) <- "profileChr"
            profileCGH <- findCluster(profileCGH, region="NormalRange", method=profileCGH$method, genome=TRUE,
                                      lambda=profileCGH$lambdaclusterGen,
                                      nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)

            ## le cluster correspondant au normal est celui qui comprend
            ## le NormalRange 0
            ##             indexNormalRange <- which(profileCGH$profileValues$NormalRange==0)
            ##             NormalCluster <- unique(profileCGH$profileValues$ZoneGen[indexNormalRange])
            ##             MedianCluster <- aggregate(profileCGH$profileValues$LogRatio, list(ZoneGen=profileCGH$profileValues$ZoneGen),median,na.rm=TRUE)
            ##             MedianCluster$ZoneGen <- as.numeric(as.character(MedianCluster$ZoneGen))
            ##             names(MedianCluster) <- c("ZoneGen","Median")
            ##             RefNorm <- MedianCluster$Median[which(MedianCluster$ZoneGen==NormalCluster)]
            ##             MedianCluster$ZoneGNL <- 0
            ##             indexClusterGain <- which(MedianCluster$Median>RefNorm)
            ##             MedianCluster$ZoneGNL[indexClusterGain] <- 1
            ##             indexClusterLost <- which(MedianCluster$Median<RefNorm)
            ##             MedianCluster$ZoneGNL[indexClusterLost] <- -1

            ##             lengthSrc <- length(MedianCluster$ZoneGen)
            ##             myZoneGNL <- .C("my_merge_int_forceGL",
            ##                             as.integer(profileCGH$profileValues$ZoneGen),
            ##                             ZoneGNL=integer(lengthDest),
            ##                             as.integer(MedianCluster$ZoneGen),
            ##                             as.integer(MedianCluster$ZoneGNL),
            ##                             as.integer(lengthDest),
            ##                             as.integer(lengthSrc),
            ##                             as.double(profileCGH$profileValues$Smoothing),
            ##                             as.double(profileCGH$forceGL[1]),
            ##                             as.double(profileCGH$forceGL[2]),
            ##                             as.double(profileCGH$NormalRef),
            ##                             as.double(profileCGH$amplicon),
            ##                             as.double(profileCGH$deletion),                                                                                    
            ##                             PACKAGE="GLAD")

            lengthDest <- length(profileCGH$profileValues$ZoneGen)
            myZoneGNL <- .C("compute_cluster_LossNormalGain",
                            ## variables pour la jointure
                            as.integer(profileCGH$profileValues$ZoneGen),
                            ZoneGNL=integer(lengthDest),
                            as.integer(lengthDest),
                            as.double(profileCGH$profileValues$Smoothing),
                            as.double(profileCGH$forceGL[1]),
                            as.double(profileCGH$forceGL[2]),
                            as.double(profileCGH$NormalRef),
                            as.double(profileCGH$amplicon),
                            as.double(profileCGH$deletion),                                                                                    
                            ## variables pour le calcul de la médiane par cluster
                            as.double(profileCGH$profileValues$LogRatio),
                            as.integer(profileCGH$profileValues$NormalRange),
                            PACKAGE="GLAD")

            

            profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL

            
            ## ########################
            ## Mise à jour des outliers
            ## ########################            
            updateOutliers <- .C("MoveBkp_updateOutliers",
                                 OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                                 OutliersTot=as.integer(profileCGH$profileValues$OutliersTot),
                                 Level=as.integer(profileCGH$profileValues[,region]),
                                 Region=as.integer(profileCGH$profileValues$Region),
                                 Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                                 Smoothing=as.double(profileCGH$profileValues$Smoothing),
                                 ZoneGNL=as.integer(profileCGH$profileValues$ZoneGNL),
                                 as.integer(l),
                                 PACKAGE="GLAD")


            ## ici on peut optimiser
            profileCGH$profileValues[,region] <- updateOutliers$Level
            profileCGH$profileValues$Region <- updateOutliers$Region
            profileCGH$profileValues$Breakpoints <- updateOutliers$Breakpoints
            profileCGH$profileValues$OutliersAws <- updateOutliers$OutliersAws
            profileCGH$profileValues$OutliersTot <- updateOutliers$OutliersTot            
            profileCGH$profileValues$Smoothing <- updateOutliers$Smoothing
            profileCGH$profileValues$ZoneGNL <- updateOutliers$ZoneGNL



            ## recalcul de la smoothing line
            ##             agg <- aggregate(profileCGH$profileValues$LogRatio, list(Level=profileCGH$profileValues$Level), median)
            ##             agg$Level <- as.numeric(as.character(agg$Level))
            ##             names(agg) <- c("Level","Smoothing")


            ##             lengthDest <- length(profileCGH$profileValues$Level)
            ##             lengthSrc <- length(agg$Level)
            ##             mySmoothing <- .C("my_merge",
            ##                              as.integer(profileCGH$profileValues$Level),
            ##                              Smoothing=double(lengthDest),
            ##                              as.integer(agg$Level),
            ##                              as.double(agg$Smoothing),
            ##                              as.integer(lengthDest),
            ##                              as.integer(lengthSrc),
            ##                              PACKAGE="GLAD")

            mySmoothing <- .C("compute_median_smoothing",
                              as.double(profileCGH$profileValues$LogRatio),
                              as.integer(profileCGH$profileValues$Level),
                              Smoothing=double(l),
                              as.integer(l),
                              PACKAGE="GLAD")


            profileCGH$profileValues$Smoothing <- mySmoothing$Smoothing
            
            
            ## ##################
            ## Mise à jour du GNL
            ## ##################
            profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGNL")]
            
            ## on prend comme référence ceux qui sont compris entre certaines valeurs            
            NormalRange <- .C("compute_NormalRange",
                              as.double(profileCGH$profileValues$Smoothing),
                              as.double(profileCGH$NormalRef),
                              as.integer(profileCGH$profileValues$Level),
                              NormalRange=integer(l),
                              as.double(profileCGH$deltaN),
                              as.integer(l),
                              PACKAGE="GLAD")

            profileCGH$profileValues$NormalRange <- NormalRange$NormalRange
            
            
            
            
            ## le clustering est fait sur les niveaux NormalRange
            class(profileCGH) <- "profileChr"
            profileCGH <- findCluster(profileCGH, region="NormalRange", method=profileCGH$method, genome=TRUE,
                                      lambda=profileCGH$lambdaclusterGen,
                                      nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)

            ## le cluster correspondant au normal est celui qui comprend
            ## le NormalRange 0
            ##             indexNormalRange <- which(profileCGH$profileValues$NormalRange==0)
            ##             NormalCluster <- unique(profileCGH$profileValues$ZoneGen[indexNormalRange])
            ##             MedianCluster <- aggregate(profileCGH$profileValues$LogRatio, list(ZoneGen=profileCGH$profileValues$ZoneGen),median,na.rm=TRUE)
            ##             MedianCluster$ZoneGen <- as.numeric(as.character(MedianCluster$ZoneGen))
            ##             names(MedianCluster) <- c("ZoneGen","Median")
            ##             RefNorm <- MedianCluster$Median[which(MedianCluster$ZoneGen==NormalCluster)]
            ##             MedianCluster$ZoneGNL <- 0
            ##             indexClusterGain <- which(MedianCluster$Median>RefNorm)
            ##             MedianCluster$ZoneGNL[indexClusterGain] <- 1
            ##             indexClusterLost <- which(MedianCluster$Median<RefNorm)
            ##             MedianCluster$ZoneGNL[indexClusterLost] <- -1

            ##             lengthSrc <- length(MedianCluster$ZoneGen)
            ##             myZoneGNL <- .C("my_merge_int_forceGL",
            ##                             as.integer(profileCGH$profileValues$ZoneGen),
            ##                             ZoneGNL=integer(lengthDest),
            ##                             as.integer(MedianCluster$ZoneGen),
            ##                             as.integer(MedianCluster$ZoneGNL),
            ##                             as.integer(lengthDest),
            ##                             as.integer(lengthSrc),
            ##                             as.double(profileCGH$profileValues$Smoothing),
            ##                             as.double(profileCGH$forceGL[1]),
            ##                             as.double(profileCGH$forceGL[2]),
            ##                             as.double(profileCGH$NormalRef),
            ##                             as.double(profileCGH$amplicon),
            ##                             as.double(profileCGH$deletion),                                                                                    
            ##                             PACKAGE="GLAD")

            lengthDest <- length(profileCGH$profileValues$ZoneGen)
            myZoneGNL <- .C("compute_cluster_LossNormalGain",
                            ## variables pour la jointure
                            as.integer(profileCGH$profileValues$ZoneGen),
                            ZoneGNL=integer(lengthDest),
                            as.integer(lengthDest),
                            as.double(profileCGH$profileValues$Smoothing),
                            as.double(profileCGH$forceGL[1]),
                            as.double(profileCGH$forceGL[2]),
                            as.double(profileCGH$NormalRef),
                            as.double(profileCGH$amplicon),
                            as.double(profileCGH$deletion),                                                                                    
                            ## variables pour le calcul de la médiane par cluster
                            as.double(profileCGH$profileValues$LogRatio),
                            as.integer(profileCGH$profileValues$NormalRange),
                            PACKAGE="GLAD")
            
            profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL
            
            class(profileCGH) <- "profileCGH"
            profileCGH$BkpInfo <- BkpInfo(profileCGH, order=FALSE)

            
            class(profileCGH) <- "profileChr"
            profileCGH <- detectOutliers(profileCGH, region=region, alpha=profileCGH$alpha, msize=profileCGH$msize)

            ## Mise à jour du GNL des Outliers
            class(profileCGH) <- "profileCGH"
            if(assignGNLOut)
              {
                profileCGH <- OutliersGNL(profileCGH, alpha=profileCGH$alpha, sigma=profileCGH$SigmaG$Value, NormalRef=profileCGH$NormalRef, amplicon=profileCGH$amplicon, deletion=profileCGH$deletion)

              }

            profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("ZoneGen","NormalRange"))]            
          }

        else
          {
            profileCGH$profileValues$ZoneGNL <- ZoneGNLAux
          }
        

      }

    return(profileCGH)

    
  }


