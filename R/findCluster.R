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
    
    
    zone <- rep(1, length(profileChr$profileValues$PosOrder))


    if (length(clusterRegion[,1])==1) 	
      {	
        nbclasses <- 1
        profileChr$NbClusterOpt <- nbclasses
        profileChr$profileValues$zone <- zone
      }	
    
    else	
      {

        sigma <- profileChr$findClusterSigma
        dist <- dist(clusterRegion$Mean)
        cluster <- hclustglad(dist, members=clusterRegion$Card, ...)
        nbclasses <- cluster(cluster, region, clusterRegion, lambda, nmin, nmax, sigma, type, param)
        classes <- cutree(cluster, k=nbclasses)
        profileChr$NbClusterOpt <- nbclasses
        clusterRegion <- data.frame(clusterRegion, zone=classes)
### jointure à optimiser
        profileChr$profileValues <- merge(profileChr$profileValues,
                                          clusterRegion[,c("Region","zone")],
                                          by.x=region, by.y="Region")

      }	

    
    if (genome==FALSE)
      {
        profileChr$profileValues$ZoneChr <- profileChr$profileValues$zone
        profileChr$profileValues <- subset(profileChr$profileValues, select=-zone)
      }
    else
      {
        profileChr$profileValues$ZoneGen <- profileChr$profileValues$zone
        profileChr$profileValues <- subset(profileChr$profileValues, select=-zone)

      }

    
    if (verbose) print("findCluster: ending function")
    return(profileChr)

  }

