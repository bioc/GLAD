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

        
        ## choix de la méthode de clustering
        METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                     "median", "centroid")
        method <- pmatch(profileCGH$method, METHODS)


        if (is.na(method)) 
          stop("invalid clustering method")
        if (method == -1) 
          stop("ambiguous clustering method")
        

        RecomputeBkpInfo <- FALSE

        ZoneGNLAux <- profileCGH$profileValues[,"ZoneGNL"]

        l <- length(profileCGH$profileValues[,1])
        updateGNL <- .C("MoveBkp_updateGNL",
                        ZoneGNL = as.integer(profileCGH$profileValues[,"ZoneGNL"]),
                        as.double(profileCGH$profileValues[,"Smoothing"]),
                        as.integer(profileCGH$profileValues[,"OutliersTot"]),
                        as.integer(l),
                        PACKAGE = "GLAD")
        

        profileCGH$profileValues[,"ZoneGNL"] <- updateGNL$ZoneGNL


        ## ######################################
        ## Déplacement des Bkp avec MoveBkp != 0
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
            lengthDest <- length(profileCGH$profileValues[,"Level"])
            l <- lengthDest             
            
            res <- .C("MoveBkp_StepAll",
                      as.integer(subBkpInfo$MoveBkp),
                      as.integer(subBkpInfo$PosOrder),
                      as.double(profileCGH$profileValues[,"LogRatio"]),
                      NextLogRatio = double(l),
                      as.integer(profileCGH$profileValues[,"Chromosome"]),                      
                      as.integer(profileCGH$profileValues[,"PosOrder"]),
                      Breakpoints = as.integer(profileCGH$profileValues[,"Breakpoints"]),
                      OutliersTot = as.integer(profileCGH$profileValues[,"OutliersTot"]),
                      OutliersAws = as.integer(profileCGH$profileValues[,"OutliersAws"]),
                      OutliersMad = as.integer(profileCGH$profileValues[,"OutliersMad"]),
                      Level = as.integer(profileCGH$profileValues[,region]),
                      Region = as.integer(profileCGH$profileValues[,"Region"]),
                      Smoothing = as.double(profileCGH$profileValues[,"Smoothing"]),
                      ZoneGNL = as.integer(profileCGH$profileValues[,"ZoneGNL"]),
                      NormalRange = integer(l),
                      ## seuils
                      as.double(profileCGH$NormalRef),
                      as.double(profileCGH$deltaN),
                      as.double(profileCGH$forceGL[1]),
                      as.double(profileCGH$forceGL[2]),
                      as.double(profileCGH$amplicon),
                      as.double(profileCGH$deletion),            
                      ## paramètres pour findCluster
                      as.integer(method),
                      as.double(profileCGH$findClusterSigma),
                      as.double(profileCGH$param["d"]),
                      as.double(profileCGH$lambdaclusterGen),
                      as.integer(profileCGH$NbClusterOpt),
                      as.integer(profileCGH$NbClusterOpt),
                      nbclasses = integer(1), ## valeur à récupérer
                      as.integer(lensubBkp),
                      as.integer(l),
                      PACKAGE="GLAD")


            ## ########################
            ## Récuparation des données 
            ## ########################
            profileCGH$profileValues[,c("Breakpoints",
                                        "OutliersTot",
                                        "OutliersAws",
                                        "OutliersMad",
                                        region,
                                        "Region",
                                        "Smoothing",
                                        "ZoneGNL",
                                        "NextLogRatio")] <- res[c("Breakpoints",
                                                                 "OutliersTot",
                                                                 "OutliersAws",
                                                                 "OutliersMad",
                                                                 "Level",
                                                                 "Region",
                                                                 "Smoothing",
                                                                 "ZoneGNL",
                                                                 "NextLogRatio")]

            profileCGH$NbClusterOpt <-  res$nbclasses


            ## Mise à jour du BkpInfo
            class(profileCGH) <- "profileCGH"
            profileCGH$BkpInfo <- BkpInfo(profileCGH, order=FALSE)

            
            class(profileCGH) <- "profileChr"
            profileCGH <- detectOutliers(profileCGH, region = region, alpha = profileCGH$alpha, msize = profileCGH$msize)

            ## Mise à jour du GNL des Outliers
            class(profileCGH) <- "profileCGH"
            if(assignGNLOut)
              {
                profileCGH <- OutliersGNL(profileCGH, alpha = profileCGH$alpha, sigma = profileCGH$SigmaG$Value,
                                          NormalRef = profileCGH$NormalRef, amplicon = profileCGH$amplicon,
                                          deletion = profileCGH$deletion)
              }

          }

        else
          {
            profileCGH$profileValues[,"ZoneGNL"] <- ZoneGNLAux
          }        

      }

    return(profileCGH)

    
  }


