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

##        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
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


        ## ######################################
        ##
        ## Déplacement des Bkp avec MoveBkp != 0
        ##
        ## ######################################

        indexMoveBkp <- which(profileCGH$BkpInfo$MoveBkp != 0)

        if (length(indexMoveBkp) > 0)
          {
            
            subBkpInfo <- profileCGH$BkpInfo[indexMoveBkp,]


            if(region == "Region")
              {
                warning("Level and Region are the same vector: the code may not works")
              }
            
            lensubBkp <- length(subBkpInfo[,1])
            lengthDest <- length(profileCGH$profileValues$Level)
            l <- lengthDest             
            
                       
##             res <- .C("MoveBkp_Delete_Bkp",
##                       as.integer(subBkpInfo$MoveBkp),
##                       as.integer(subBkpInfo$PosOrder),
##                       Breakpoints = as.integer(profileCGH$profileValues$Breakpoints),
##                       OutliersTot = as.integer(profileCGH$profileValues$OutliersTot),
##                       OutliersAws = as.integer(profileCGH$profileValues$OutliersAws),
##                       OutliersMad = as.integer(profileCGH$profileValues$OutliersMad),
##                       Level = as.integer(profileCGH$profileValues[,region]),
##                       Region = as.integer(profileCGH$profileValues$Region),
##                       Smoothing = as.double(profileCGH$profileValues$Smoothing),
##                       ZoneGNL = as.integer(profileCGH$profileValues$ZoneGNL),
##                       as.integer(lensubBkp),
##                       PACKAGE="GLAD")

            
            

##             ## recalcul de la smoothing line
##             mySmoothing <- .C("compute_median_smoothing",
##                               as.double(profileCGH$profileValues$LogRatio),
##                               as.integer(profileCGH$profileValues$Level),
##                               Smoothing = double(l),
##                               as.integer(l),
##                               PACKAGE="GLAD")

            
##             profileCGH$profileValues$Smoothing <- mySmoothing$Smoothing

            
##             updateLevel <- .C("updateLevel",
##                               as.integer(profileCGH$profileValues$Chromosome),
##                               Breakpoints = as.integer(profileCGH$profileValues$Breakpoints),
##                               Level = as.integer(profileCGH$profileValues[,region]),
##                               as.integer(profileCGH$profileValues$PosOrder),
##                               NextLogRatio = double(l),
##                               as.double(profileCGH$profileValues$LogRatio),
##                               as.integer(max(profileCGH$profileValues$Level)),
##                               as.integer(l),
##                               PACKAGE="GLAD")

##             profileCGH$profileValues[,c(region, "Breakpoints", "NextLogRatio")] <- updateLevel[c("Level", "Breakpoints", "NextLogRatio")]
            
            ## ##################
            ## Mise à jour du GNL
            ## ##################

##            profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGNL")]
            
##             ## on prend comme référence ceux qui sont compris entre certaines valeurs            
##             NormalRange <- .C("compute_NormalRange",
##                               as.double(profileCGH$profileValues$Smoothing),
##                               as.double(profileCGH$NormalRef),
##                               as.integer(profileCGH$profileValues$Level),
##                               NormalRange = integer(l),
##                               as.double(profileCGH$deltaN),
##                               as.integer(l),
##                               PACKAGE="GLAD")

##             profileCGH$profileValues$NormalRange <- NormalRange$NormalRange


            #####################################
            ## assemblage des différentes étapes
            res <- .C("MoveBkp_Step1",
                      as.integer(subBkpInfo$MoveBkp),
                      as.integer(subBkpInfo$PosOrder),
                      as.double(profileCGH$profileValues$LogRatio),
                      NextLogRatio = double(l),
                      as.integer(profileCGH$profileValues$Chromosome),                      
                      as.integer(profileCGH$profileValues$PosOrder),
                      Breakpoints = as.integer(profileCGH$profileValues$Breakpoints),
                      OutliersTot = as.integer(profileCGH$profileValues$OutliersTot),
                      OutliersAws = as.integer(profileCGH$profileValues$OutliersAws),
                      OutliersMad = as.integer(profileCGH$profileValues$OutliersMad),
                      Level = as.integer(profileCGH$profileValues[,region]),
                      Region = as.integer(profileCGH$profileValues$Region),
                      Smoothing = as.double(profileCGH$profileValues$Smoothing),
                      ZoneGNL = as.integer(profileCGH$profileValues$ZoneGNL),
                      NormalRange = integer(l),
                      ## seuils
                      as.double(profileCGH$NormalRef),
                      as.double(profileCGH$deltaN),
                      as.integer(lensubBkp),
                      as.integer(l),
                      PACKAGE="GLAD")
            ## fin de l'assemblage
            ######################################


            ## #########################################            
            ## Récuparation des données après assemblage
            profileCGH$profileValues[,c("Breakpoints",
                                        "OutliersTot",
                                        "OutliersAws",
                                        "OutliersMad",
                                        region,
                                        "Region",
                                        "Smoothing",
                                        "ZoneGNL",
                                        "NextLogRatio",
                                        "NormalRange")] <- res[c("Breakpoints",
                                                                 "OutliersTot",
                                                                 "OutliersAws",
                                                                 "OutliersMad",
                                                                 "Level",
                                                                 "Region",
                                                                 "Smoothing",
                                                                 "ZoneGNL",
                                                                 "NextLogRatio",
                                                                 "NormalRange")]
            
            ## fin de la récupération des données
            ## #########################################                        

            
            ## le clustering est fait sur les niveaux NormalRange
            class(profileCGH) <- "profileChr"
            profileCGH <- findCluster(profileCGH, region = "NormalRange", method = profileCGH$method, genome = TRUE,
                                      lambda = profileCGH$lambdaclusterGen,
                                      nmin = profileCGH$NbClusterOpt, nmax = profileCGH$NbClusterOpt)


            ## ###################################
            ## assemblage des différentes étapes

            res <- .C("MoveBkp_Step2",
                      OutliersAws = as.integer(profileCGH$profileValues$OutliersAws),
                      OutliersTot = as.integer(profileCGH$profileValues$OutliersTot),
                      Level = as.integer(profileCGH$profileValues[,region]),
                      Region = as.integer(profileCGH$profileValues$Region),
                      Breakpoints = as.integer(profileCGH$profileValues$Breakpoints),                      
                      ## variables pour la jointure
                      as.integer(profileCGH$profileValues$ZoneGen),
                      ZoneGNL = integer(lengthDest),
                      as.integer(lengthDest),
                      Smoothing = as.double(profileCGH$profileValues$Smoothing),
                      as.double(profileCGH$forceGL[1]),
                      as.double(profileCGH$forceGL[2]),
                      as.double(profileCGH$NormalRef),
                      as.double(profileCGH$amplicon),
                      as.double(profileCGH$deletion),
                      as.double(profileCGH$deltaN),
                      ## variables pour le calcul de la médiane par cluster
                      as.double(profileCGH$profileValues$LogRatio),
                      NormalRange = as.integer(profileCGH$profileValues$NormalRange),
                      PACKAGE="GLAD")

            ## fin de l'assemblage
            ## ###################################


            ## ##################################
            ## récupération des données
            profileCGH$profileValues[,c(region,"Region","Breakpoints","OutliersAws","OutliersTot", "Smoothing", "ZoneGNL", "NormalRange")] <-  res[c("Level","Region","Breakpoints","OutliersAws","OutliersTot", "Smoothing", "ZoneGNL", "NormalRange")]

## ##            lengthDest <- length(profileCGH$profileValues$ZoneGen)
##             myZoneGNL <- .C("compute_cluster_LossNormalGain",
##                             ## variables pour la jointure
##                             as.integer(profileCGH$profileValues$ZoneGen),
##                             ZoneGNL = integer(lengthDest),
##                             as.integer(lengthDest),
##                             as.double(profileCGH$profileValues$Smoothing),
##                             as.double(profileCGH$forceGL[1]),
##                             as.double(profileCGH$forceGL[2]),
##                             as.double(profileCGH$NormalRef),
##                             as.double(profileCGH$amplicon),
##                             as.double(profileCGH$deletion),
##                             ## variables pour le calcul de la médiane par cluster
##                             as.double(profileCGH$profileValues$LogRatio),
##                             as.integer(profileCGH$profileValues$NormalRange),
##                             PACKAGE="GLAD")
            

##             profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL

            
##             ## ##########################
##             ## Mise à jour des outliers
##             ## ##########################            
##             updateOutliers <- .C("MoveBkp_updateOutliers",
##                                  OutliersAws = as.integer(profileCGH$profileValues$OutliersAws),
##                                  OutliersTot = as.integer(profileCGH$profileValues$OutliersTot),
##                                  Level = as.integer(profileCGH$profileValues[,region]),
##                                  Region = as.integer(profileCGH$profileValues$Region),
##                                  Breakpoints = as.integer(profileCGH$profileValues$Breakpoints),
##                                  Smoothing = as.double(profileCGH$profileValues$Smoothing),
##                                  ZoneGNL = as.integer(profileCGH$profileValues$ZoneGNL),
##                                  as.integer(l),
##                                  PACKAGE="GLAD")




##             ## recalcul de la smoothing line
##             mySmoothing <- .C("compute_median_smoothing",
##                               as.double(profileCGH$profileValues$LogRatio),
##                               as.integer(profileCGH$profileValues$Level),
##                               Smoothing = double(l),
##                               as.integer(l),
##                               PACKAGE="GLAD")


##             profileCGH$profileValues$Smoothing <- mySmoothing$Smoothing
            
            
            ## ##################
            ## Mise à jour du GNL
            ## ##################
  ##          profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGNL")]
            
##             ## on prend comme référence ceux qui sont compris entre certaines valeurs            
##             NormalRange <- .C("compute_NormalRange",
##                               as.double(profileCGH$profileValues$Smoothing),
##                               as.double(profileCGH$NormalRef),
##                               as.integer(profileCGH$profileValues$Level),
##                               NormalRange = integer(l),
##                               as.double(profileCGH$deltaN),
##                               as.integer(l),
##                               PACKAGE="GLAD")

##             profileCGH$profileValues$NormalRange <- NormalRange$NormalRange                        
            
            
            ## le clustering est fait sur les niveaux NormalRange
            class(profileCGH) <- "profileChr"
            profileCGH <- findCluster(profileCGH, region="NormalRange", method=profileCGH$method, genome=TRUE,
                                      lambda=profileCGH$lambdaclusterGen,
                                      nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)


##            lengthDest <- length(profileCGH$profileValues$ZoneGen)
            myZoneGNL <- .C("compute_cluster_LossNormalGain",
                            ## variables pour la jointure
                            as.integer(profileCGH$profileValues$ZoneGen),
                            ZoneGNL = integer(lengthDest),
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


