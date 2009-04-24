filterBkp <- function(...)
  {
    UseMethod("filterBkp")
  }

filterBkp.profileCGH <- function(profileCGH, MinBkpWeight=0.25, assignGNLOut=TRUE, verbose=FALSE, ...)
  {


    profilage <- FALSE
    if (verbose) print("filterBkp: starting function")            
    if (is.data.frame(profileCGH$BkpInfo))
      {
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
        RecomputeGNL <- FALSE


##################################################################################
###
### On supprime les Breakpoints qui sont situés au sein des régions amplifiées
###
##################################################################################
        


        if (verbose) print("filterBkp: Breakpoints in amplified regions are removed")
        indexBkpToDel <- which(profileCGH$BkpInfo["GNLchange"]==0 & profileCGH$BkpInfo["ZoneGNL"]==2)
        if (length(indexBkpToDel)>0)
          {
            RecomputeGNL <- TRUE
            profileCGH$profileValues$Breakpoints[profileCGH$BkpInfo$PosOrder[indexBkpToDel]] <- -1
            profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexBkpToDel,]
          }


##################################################################################
###
### On Déplace les Bkp qui sont aussi Outliers et dont
### le GNL correspond à celui du BAC d'après
### On déplace également les Bkp après lequel il y a un outlier
### correspondant au statut du Bkp        
###
##################################################################################
        
        if (verbose) print("filterBkp: move breakpoints which are outliers")        
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]

        nb <- length(profileCGH$profileValues[,1])-1
        

        moveBkp <- .C("moveBkp",
                      as.integer(profileCGH$profileValues$ZoneGNL),
                      Level=as.integer(profileCGH$profileValues$Level),
                      Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                      OutliersTot=as.integer(profileCGH$profileValues$OutliersTot),
                      OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                      as.integer(profileCGH$profileValues$Chromosome),
                      RecomputeSmt=as.integer(0),
                      as.integer(nb),
                      PACKAGE="GLAD")

        if (moveBkp$RecomputeSmt==1)
          {            
            rownames(profileCGH$profileValues) <- 0:(length(profileCGH$profileValues[,1])-1)
            RecomputeGNL <- TRUE
            profileCGH$profileValues$Level <- moveBkp$Level
            profileCGH$profileValues$Breakpoints <- moveBkp$Breakpoints
            profileCGH$profileValues$OutliersTot <- moveBkp$OutliersTot
            profileCGH$profileValues$OutliersAws <- moveBkp$OutliersAws
          }
        

##################################################################################
###        
### Suppression des Bkp dont est poids est inférieur à un seuil
### et qui ne correspondent pas à un changement de GNL
###
##################################################################################        


        if (verbose) print("filterBkp: delete breakpoint with small weight")
        indexWeightToSmall <- which(profileCGH$BkpInfo["Weight"]<MinBkpWeight & profileCGH$BkpInfo["GNLchange"]==0 & profileCGH$BkpInfo["ZoneGNL"]!=2)
        if (length(indexWeightToSmall)>0)
          {
            RecomputeGNL <- TRUE
            indexPos <- profileCGH$BkpInfo$PosOrder[indexWeightToSmall]
            profileCGH$profileValues$Breakpoints[indexPos] <- -1            
            ##            print(profileCGH$BkpInfo$PosOrder[indexWeightToSmall])
            ##             for (PosBkp in profileCGH$BkpInfo$PosOrder[indexWeightToSmall])
            ##               {
            ##                                         #indexPos <- which(profileCGH$profileValues$PosOrder==PosBkp)
            ##                                         #profileCGH$profileValues$Breakpoints[indexPos] <- -1
            ##                 profileCGH$profileValues$Breakpoints[PosBkp] <- -1
            ##               }
            if (length(indexWeightToSmall)==length(profileCGH$BkpInfo[,1]))
              {
                profileCGH$BkpInfo <- NA
              }
            else
              {
                profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexWeightToSmall,]
              }

          }


##################################################################################
###        
### Suppression des Bkp dont le poids vaut 0
### et qui correspondent à un changement de GNL
### cette situation peut arriver après élimination des Bkp
### dont le poids est inférieur à un seuil
###
##################################################################################

        if (verbose) print("filterBkp: delete breakpoint with null weight")
        indexWeightZero <- which(profileCGH$BkpInfo["Weight"]==0 & profileCGH$BkpInfo["GNLchange"]==1)
        if (length(indexWeightZero)>0)
          {

            RecomputeGNL <- TRUE
            indexPos <- profileCGH$BkpInfo$PosOrder[indexWeightZero]
            profileCGH$profileValues$Breakpoints[indexPos] <- -1            
            ##             for (PosBkp in profileCGH$BkpInfo$PosOrder[indexWeightZero])
            ##               {
            ##                                         #indexPos <- which(profileCGH$profileValues$PosOrder==PosBkp)
            ##                                         #profileCGH$profileValues$Breakpoints[indexPos] <- -1
            ##                 profileCGH$profileValues$Breakpoints[PosBkp] <- -1
            ##               }
            
            if (length(indexWeightZero)==length(profileCGH$BkpInfo[,1]))
              {
                profileCGH$BkpInfo <- NA
              }
            else
              {
                profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexWeightZero,]
              }

          }


### Quand je vais recalculer les Outliers, il faut le NormalRef
### Attention à ce que celui-ci soit bien transmis
### Normalement NormalRef vaut 0 puisqu'en sortie de gladLA
### les log-ratios sont centrés sur NormalRef
        

        if (RecomputeGNL)
          {
            if(profilage) Rprof("/tmp/Step00FilterBkp.dat")
            if (verbose) print("filterBkp: recomputeGNL")        
            profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
            l <- length(profileCGH$profileValues[,1])
            updateLevel <- .C("updateLevel",
                              as.integer(profileCGH$profileValues$Chromosome),
                              Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                              Level=as.integer(profileCGH$profileValues$Level),
                              as.integer(profileCGH$profileValues$PosOrder),
                              NextLogRatio=as.double(profileCGH$profileValues$NextLogRatio),
                              as.double(profileCGH$profileValues$LogRatio),
                              as.integer(max(profileCGH$profileValues$Level)),
                              as.integer(l),
                              PACKAGE="GLAD")

            profileCGH$profileValues$Level <- updateLevel$Level
            profileCGH$profileValues$Breakpoints <- updateLevel$Breakpoints
            profileCGH$profileValues$NextLogRatio <- updateLevel$NextLogRatio

            updateOutliers <- .C("updateOutliers",
                                 OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                                 Level=as.integer(profileCGH$profileValues$Level),
                                 Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                                 Smoothing=as.double(profileCGH$profileValues$Smoothing),
                                 as.integer(l),
                                 PACKAGE="GLAD")
            profileCGH$profileValues$Level <- updateOutliers$Level
            profileCGH$profileValues$Breakpoints <- updateOutliers$Breakpoints
            profileCGH$profileValues$OutliersAws <- updateOutliers$OutliersAws
            profileCGH$profileValues$Smoothing <- updateOutliers$Smoothing

            if(profilage) Rprof(NULL)
            
            if(profilage) Rprof("/tmp/Step01FilterBkpOutliers.dat")            
### Recalcul des Outliers
            class(profileCGH) <- "profileChr"
            profileCGH <- detectOutliers(profileCGH, region="Level", alpha=profileCGH$alpha, msize=profileCGH$msize)
            
### recalcul de la smoothing line
            agg <- aggregate(profileCGH$profileValues$LogRatio, list(Level=profileCGH$profileValues$Level), median)
            agg$Level <- as.numeric(as.character(agg$Level))
            names(agg) <- c("Level","Smoothing")

### ne sert plus à rien car on utilise la fonction my_merge à la place de la fonction merge de R
###            profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Smoothing"))

            lengthDest <- length(profileCGH$profileValues$Level)
            lengthSrc <- length(agg$Level)
            mySmooting <- .C("my_merge",
                             as.integer(profileCGH$profileValues$Level),
                             Smoothing=double(lengthDest),
                             as.integer(agg$Level),
                             as.double(agg$Smoothing),
                             as.integer(lengthDest),
                             as.integer(lengthSrc),
                             PACKAGE="GLAD")

            profileCGH$profileValues$Smoothing <- mySmooting$Smoothing
            
           
###             t1 <- system.time(profileCGH$profileValues <- merge(profileCGH$profileValues, agg, by="Level", all=TRUE))


####################
            profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGNL")]
            indexNormalLevel <- which(abs(profileCGH$profileValues$Smoothing-profileCGH$NormalRef)<=profileCGH$deltaN)
            profileCGH$profileValues$NormalRange <- profileCGH$profileValues$Level
            profileCGH$profileValues$NormalRange[indexNormalLevel] <- 0

### le clustering est fait sur les niveaux NormalRange
            profileCGH <- findCluster(profileCGH, region="NormalRange", method=profileCGH$method, genome=TRUE,
                                      lambda=profileCGH$lambdaclusterGen,
                                      nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)

### le cluster correspondant au normal est celui qui comprend
### le NormalRange 0
            indexNormalRange <- which(profileCGH$profileValues$NormalRange==0)
            NormalCluster <- unique(profileCGH$profileValues$ZoneGen[indexNormalRange])
            MedianCluster <- aggregate(profileCGH$profileValues$LogRatio, list(ZoneGen=profileCGH$profileValues$ZoneGen),median,na.rm=TRUE)
            MedianCluster$ZoneGen <- as.numeric(as.character(MedianCluster$ZoneGen))
            names(MedianCluster) <- c("ZoneGen","Median")
            RefNorm <- MedianCluster$Median[which(MedianCluster$ZoneGen==NormalCluster)]
            MedianCluster$ZoneGNL <- 0
            indexClusterGain <- which(MedianCluster$Median>RefNorm)
            MedianCluster$ZoneGNL[indexClusterGain] <- 1
            indexClusterLost <- which(MedianCluster$Median<RefNorm)
            MedianCluster$ZoneGNL[indexClusterLost] <- -1

            lengthDest <- length(profileCGH$profileValues$ZoneGen)
            lengthSrc <- length(MedianCluster$ZoneGen)
            myZoneGNL <- .C("my_merge_int_forceGL",
                            as.integer(profileCGH$profileValues$ZoneGen),
                            ZoneGNL=integer(lengthDest),
                            as.integer(MedianCluster$ZoneGen),
                            as.integer(MedianCluster$ZoneGNL),
                            as.integer(lengthDest),
                            as.integer(lengthSrc),
                            as.double(profileCGH$profileValues$Smoothing),
                            as.double(profileCGH$forceGL[1]),
                            as.double(profileCGH$forceGL[2]),
                            as.double(profileCGH$NormalRef),
                            as.double(profileCGH$amplicon),
                            as.double(profileCGH$deletion),                                                                                    
                            PACKAGE="GLAD")

            profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL


            
###             t2 <- system.time(profileCGH$profileValues <- merge(profileCGH$profileValues, MedianCluster[,c("ZoneGen","ZoneGNL")], all=TRUE, by="ZoneGen"))


## ### on force les gains et les pertes pour certaines valeur de smoothing
##             indexForceGain <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) >= profileCGH$forceGL[2])
##             profileCGH$profileValues$ZoneGNL[indexForceGain] <- 1
##             indexForceLost <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) <= profileCGH$forceGL[1])
##             profileCGH$profileValues$ZoneGNL[indexForceLost] <- -1

## ### Amplicon et deletion
##             indexAmp <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) >= profileCGH$amplicon)
##             profileCGH$profileValues$ZoneGNL[indexAmp] <- 2
##             indexDel <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) <= profileCGH$deletion)
##             profileCGH$profileValues$ZoneGNL[indexDel] <- -10
            


            profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("ZoneGen","NormalRange"))]
##             profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGen"))
##             profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"NormalRange"))

            

            class(profileCGH) <- "profileCGH"

            profileCGH$BkpInfo <- BkpInfo(profileCGH)
            
### Mise à jour du GNL des Outliers
            class(profileCGH) <- "profileCGH"

            if(assignGNLOut)
              {
                profileCGH <- OutliersGNL(profileCGH, alpha=profileCGH$alpha,
                                          sigma=profileCGH$SigmaG$Value, NormalRef=profileCGH$NormalRef,
                                          amplicon=profileCGH$amplicon, deletion=profileCGH$deletion)
              }

            if(profilage) Rprof(NULL)
          }
        
      }


    if (verbose) print("filterBkp: ending function")                

    return(profileCGH)

    
  }


