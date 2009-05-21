## Copyright (C) 2005 Institut Curie
## Author(s): Philippe Hupé (Institut Curie) 2005
## Contact: glad@curie.fr

filterBkp <- function(...)
  {
    UseMethod("filterBkp")
  }

filterBkp.profileCGH <- function(profileCGH, MinBkpWeight=0.25, assignGNLOut=TRUE, verbose=FALSE, ...)
  {
    
    if (verbose) print("filterBkp: starting function")
    
    if (is.data.frame(profileCGH$BkpInfo))
      {
        
        
        RecomputeGNL <- FALSE


        ## ################################################################################
        ## On supprime les Breakpoints qui sont situés au sein des régions amplifiées
        ## ################################################################################
        


        if (verbose) print("filterBkp: Breakpoints in amplified regions are removed")
        
        indexBkpToDel <- which(profileCGH$BkpInfo["GNLchange"] == 0 & profileCGH$BkpInfo["ZoneGNL"] == 2)
        if (length(indexBkpToDel) > 0)
          {
            RecomputeGNL <- TRUE
            profileCGH$profileValues$Breakpoints[profileCGH$BkpInfo$PosOrder[indexBkpToDel]] <- -1
            profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexBkpToDel,]
          }

        

        ## ################################################################################
        ## On Déplace les Bkp qui sont aussi Outliers et dont
        ## le GNL correspond à celui du BAC d'après
        ## On déplace également les Bkp après lequel il y a un outlier
        ## correspondant au statut du Bkp        
        ## ################################################################################
        
        if (verbose) print("filterBkp: move breakpoints which are outliers")        

        nb <- length(profileCGH$profileValues[,1]) - 1        
        moveBkp <- .C("filterBkp_moveBkp_Outliers",
                      as.integer(profileCGH$profileValues[,"ZoneGNL"]),
                      Level = as.integer(profileCGH$profileValues[,"Level"]),
                      Breakpoints = as.integer(profileCGH$profileValues[,"Breakpoints"]),
                      OutliersTot = as.integer(profileCGH$profileValues[,"OutliersTot"]),
                      OutliersAws = as.integer(profileCGH$profileValues[,"OutliersAws"]),
                      as.integer(profileCGH$profileValues[,"Chromosome"]),
                      RecomputeSmt = as.integer(0),
                      as.integer(nb),
                      PACKAGE = "GLAD")

        if (moveBkp$RecomputeSmt == 1)
          {            
            rownames(profileCGH$profileValues) <- 0:(length(profileCGH$profileValues[,1])-1)
            RecomputeGNL <- TRUE

            profileCGH$profileValues[,c("Level", "Breakpoints", "OutliersTot", "OutliersAws")] <- moveBkp[c("Level", "Breakpoints", "OutliersTot", "OutliersAws")]
          }
        


        ## ################################################################################
        ## Suppression des Bkp dont est poids est inférieur à un seuil
        ## et qui ne correspondent pas à un changement de GNL
        ## ################################################################################        


        if (verbose) print("filterBkp: delete breakpoint with small weight")
        indexWeightToSmall <- which(profileCGH$BkpInfo["Weight"] < MinBkpWeight &
                                    profileCGH$BkpInfo["GNLchange"] == 0 &
                                    profileCGH$BkpInfo["ZoneGNL"] != 2)
        if (length(indexWeightToSmall) > 0)
          {
            RecomputeGNL <- TRUE
            indexPos <- profileCGH$BkpInfo[,"PosOrder"][indexWeightToSmall]
            profileCGH$profileValues[,"Breakpoints"][indexPos] <- -1            
            if (length(indexWeightToSmall) == length(profileCGH$BkpInfo[,1]))
              {
                profileCGH$BkpInfo <- NA
              }
            else
              {
                profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexWeightToSmall,]

              }
          }



        ## ################################################################################
        ## Suppression des Bkp dont le poids vaut 0
        ## et qui correspondent à un changement de GNL
        ## cette situation peut arriver après élimination des Bkp
        ## dont le poids est inférieur à un seuil
        ## ################################################################################

        if (verbose) print("filterBkp: delete breakpoint with null weight")
        indexWeightZero <- which(profileCGH$BkpInfo["Weight"] == 0 &
                                 profileCGH$BkpInfo["GNLchange"] == 1)
        if (length(indexWeightZero) > 0)
          {
            RecomputeGNL <- TRUE
            indexPos <- profileCGH$BkpInfo$PosOrder[indexWeightZero]
            profileCGH$profileValues[,"Breakpoints"][indexPos] <- -1            
            
            if (length(indexWeightZero) == length(profileCGH$BkpInfo[,1]))
              {
                profileCGH$BkpInfo <- NA
              }
            else
              {
                profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexWeightZero,]
              }
          }


        ## Quand je vais recalculer les Outliers, il faut le NormalRef
        ## Attention à ce que celui-ci soit bien transmis
        ## Normalement NormalRef vaut 0 puisqu'en sortie de gladLA
        ## les log-ratios sont centrés sur NormalRef
        

        if (RecomputeGNL)
          {
            if (verbose) print("filterBkp: recomputeGNL")        

            ## La détection des Outliers va être faite directement dans la fonction C
            alpha <- profileCGH$alpha
            msize <- profileCGH$msize
            if (msize<1) stop("msize must be greater or equal to 1")
            if (alpha>1 || alpha <0)stop("alpha must be setted between 0 and 1")
            alpha <- qnorm(1-alpha/2)            
            
            l <- length(profileCGH$profileValues[,1])
            updateFilterBkp <- .C("updateFilterBkp",
                                  as.integer(profileCGH$profileValues[,"Chromosome"]),
                                  Breakpoints = as.integer(profileCGH$profileValues[,"Breakpoints"]),  ## valeur de sortie
                                  Level = as.integer(profileCGH$profileValues[,"Level"]),              ## valeur de sortie
                                  as.integer(profileCGH$profileValues[,"PosOrder"]),
                                  NextLogRatio  = as.double(profileCGH$profileValues[,"NextLogRatio"]), ## valeur de sortie
                                  as.double(profileCGH$profileValues[,"LogRatio"]),
                                  as.integer(max(profileCGH$profileValues[,"Level"])),
                                  ## ajout des variables pour updateOutliers
                                  OutliersAws = as.integer(profileCGH$profileValues[,"OutliersAws"]),  ## valeur de sortie
                                  Smoothing = as.double(profileCGH$profileValues[,"Smoothing"]),       ## valeur de sortie
                                  ## ajout des variables pour detectOutliers
                                  OutliersMad = integer(l),                                        ## valeur de sortie
                                  OutliersTot = integer(l),                                        ## valeur de sortie
                                  as.integer(msize),
                                  as.double(alpha),
                                  as.integer(l),
                                  as.double(profileCGH$NormalRef),
                                  as.double(profileCGH$deltaN),
                                  NormalRange = integer(l),
                                  PACKAGE = "GLAD")


            
            profileCGH$profileValues[,c("Level", "NextLogRatio", "Breakpoints", "OutliersAws", "Smoothing", "OutliersTot", "OutliersMad", "NormalRange")] <- updateFilterBkp[c("Level", "NextLogRatio", "Breakpoints", "OutliersAws", "Smoothing", "OutliersTot", "OutliersMad", "NormalRange")]

            
            ## le clustering est fait sur les niveaux NormalRange
            ## il faut prendre nmin comme le min(nmin,profileCGH$NbClusterOpt)
            class(profileCGH) <- "profileChr"
            profileCGH <- findCluster(profileCGH, region = "NormalRange", method = profileCGH$method, genome = TRUE,
                                      lambda = profileCGH$lambdaclusterGen,
                                      nmin = min(profileCGH$nmin, profileCGH$NbClusterOpt), nmax = profileCGH$NbClusterOpt, param = profileCGH$param)
            class(profileCGH) <- "profileCGH"


            lengthDest <- length(profileCGH$profileValues[,"ZoneGen"])
            myZoneGNL <- .C("compute_cluster_LossNormalGain",
                            ## variables pour la jointure
                            as.integer(profileCGH$profileValues[,"ZoneGen"]),
                            ZoneGNL = integer(lengthDest),
                            as.integer(lengthDest),
                            as.double(profileCGH$profileValues[,"Smoothing"]),
                            as.double(profileCGH$forceGL[1]),
                            as.double(profileCGH$forceGL[2]),
                            as.double(profileCGH$NormalRef),
                            as.double(profileCGH$amplicon),
                            as.double(profileCGH$deletion),                                                                                    
                            ## variables pour le calcul de la médiane par cluster
                            as.double(profileCGH$profileValues[,"LogRatio"]),
                            as.integer(profileCGH$profileValues[,"NormalRange"]),
                            PACKAGE = "GLAD")



            profileCGH$profileValues[,"ZoneGNL"] <- myZoneGNL$ZoneGNL
            
            ## ################################
            ## Mise des infos sur les Bkp
            ## ################################            


            profileCGH$BkpInfo <- BkpInfo(profileCGH)

            ## ################################
            ## Mise à jour du GNL des Outliers
            ## ################################            

            if(assignGNLOut)
              {
                profileCGH <- OutliersGNL(profileCGH, alpha = profileCGH$alpha,
                                          sigma = profileCGH$SigmaG$Value, NormalRef = profileCGH$NormalRef,
                                          amplicon = profileCGH$amplicon, deletion = profileCGH$deletion)
              }

          }
        
      }


    if (verbose) print("filterBkp: ending function")                

    return(profileCGH)

    
  }

