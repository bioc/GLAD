# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

findCluster <- function(...)
  {

    UseMethod("findCluster")
  }


findCluster.profileChr <- function(profileChr, region="Region", genome=TRUE, lambda=10, nmin=1, nmax=10, type="tricubic", param=c(d=6), verbose=FALSE, ...)
  {

    if (verbose) print("findCluster: starting function")


    subsetdata <- profileChr$profileValues[which(profileChr$profileValues$OutliersTot==0),]

    
### vérifier le comportement pour les clusters de cardinalité 1

    sagg <- split(subsetdata$LogRatio,subsetdata[,region])
    Mean <- sapply(sagg,mean)
    Card <- sapply(sagg,NROW)
    Var <- sapply(sagg,var)
    Region <- as.numeric(as.character(names(Card)))
    clusterRegion <- data.frame(Region, Card, Var, Mean)
    clusterRegion$Var <- clusterRegion$Var*((clusterRegion$Card-1)/clusterRegion$Card)
    clusterRegion$VarLike <- clusterRegion$Var
    indexSingle <- which(clusterRegion$Card==1)
    clusterRegion$Var[indexSingle] <- 0
    clusterRegion$VarLike[indexSingle] <- 1
    
    
###    zone <- rep(1, length(profileChr$profileValues$PosOrder))


    if (length(clusterRegion[,1])==1) 	
      {	
        nbclasses <- 1
        profileChr$NbClusterOpt <- nbclasses
###        profileChr$profileValues$zone <- zone
        if (genome==FALSE)
          {
            profileChr$profileValues$ZoneChr <- 1
          }
        else
          {
            profileChr$profileValues$ZoneGen <- 1
          }
        
        profileChr$profileValues$zone <- 1        
      }	
    
    else	
      {

        sigma <- profileChr$findClusterSigma
        dist <- dist(clusterRegion$Mean)
        cluster.res <- hclustglad(dist, members=clusterRegion$Card, ...)
        nbclasses <- clusterglad(cluster.res, region, clusterRegion, lambda, nmin, nmax, sigma, type, param)
        classes <- cutree(cluster.res, k=nbclasses)
        profileChr$NbClusterOpt <- nbclasses
        clusterRegion <- data.frame(clusterRegion, zone=classes)

        lengthDest <- length(profileChr$profileValues[,region])
        lengthSrc <- length(clusterRegion$Region)
        myzone <- .C("my_merge_int",
                     as.integer(profileChr$profileValues[,region]),
                     zone=integer(lengthDest),
                     as.integer(clusterRegion$Region),
                     as.integer(clusterRegion$zone),
                     as.integer(lengthDest),
                     as.integer(lengthSrc),
                     PACKAGE="GLAD")

###        profileChr$profileValues$zone <- myzone$zone


        if (genome==FALSE)
          {
            profileChr$profileValues$ZoneChr <- myzone$zone
          }
        else
          {
            profileChr$profileValues$ZoneGen <- myzone$zone
          }

      }	

    
    
    if (verbose) print("findCluster: ending function")
    return(profileChr)

  }

