## Copyright (C) 2003 Institut Curie
## Author(s): Philippe Hupé (Institut Curie) 2003
## Contact: glad@curie.fr

findCluster <- function(...)
  {

    UseMethod("findCluster")
  }


findCluster.profileChr <- function(profileChr, region="Region", genome=TRUE, lambda=10, nmin=1, nmax=10, type="tricubic", param=c(d=6), verbose=FALSE, ...)
  {

    if (verbose) print("findCluster: starting function")


    t.start <- Sys.time()
    
    t1 <- system.time(subsetdata <- profileChr$profileValues[which(profileChr$profileValues$OutliersTot == 0),c(region,"LogRatio")])
    print("subsetdata")
    print(t1)

    
### vérifier le comportement pour les clusters de cardinalité 1

    t1 <- Sys.time()

    sagg <- split(subsetdata$LogRatio,subsetdata[,region])
    Mean <- sapply(sagg,mean)
    Card <- sapply(sagg,NROW)
    Var <- sapply(sagg,var)
    Region <- as.numeric(as.character(names(Card)))
    clusterRegion <- data.frame(Region, Card, Var, Mean)
    clusterRegion$Var <- clusterRegion$Var*((clusterRegion$Card-1)/clusterRegion$Card)
    clusterRegion$VarLike <- clusterRegion$Var
    indexSingle <- which(clusterRegion$Card == 1)
    clusterRegion$Var[indexSingle] <- 0
    clusterRegion$VarLike[indexSingle] <- 1


    t2 <- Sys.time()
    print("diff")
    print(t2 - t1)
    
    
    if (length(clusterRegion[,1]) == 1) 	
      {	
        nbclasses <- 1
        profileChr$NbClusterOpt <- nbclasses

        if (genome == FALSE)
          {
            profileChr$profileValues$ZoneChr <- 1
          }
        else
          {
            profileChr$profileValues$ZoneGen <- 1
          }
      }	
    
    else	
      {

        sigma <- profileChr$findClusterSigma
        dist <- dist(clusterRegion$Mean)
        t3 <- system.time(cluster.res <- hclustglad(dist, members=clusterRegion$Card, ...))
        print("hclustglad")
        print(t3)
        t3 <- system.time(nbclasses <- clusterglad(cluster.res, region, clusterRegion, lambda, nmin, nmax, sigma, type, param))
        print("clusterglad")
        print(t3)
        classes <- cutree(cluster.res, k=nbclasses)
        profileChr$NbClusterOpt <- nbclasses
        clusterRegion <- data.frame(clusterRegion, zone=classes)


        lengthDest <- length(profileChr$profileValues[,region])
        lengthSrc <- length(clusterRegion$Region)
        t1 <- system.time(myzone <- .C("my_merge_int",
                     as.integer(profileChr$profileValues[,region]),
                     zone=integer(lengthDest),
                     as.integer(clusterRegion$Region),
                     as.integer(clusterRegion$zone),
                     as.integer(lengthDest),
                     as.integer(lengthSrc),
                     PACKAGE="GLAD"))
        print("merge")
        print(t1)


        if (genome == FALSE)
          {
            profileChr$profileValues[,"ZoneChr"] <- myzone$zone
          }
        else
          {
            t1 <- system.time(profileChr$profileValues[,"ZoneGen"] <- myzone$zone)
            print("zone")
            print(t1)
          }

      }	

    t.end <- Sys.time()

    print(paste("Temps findCluster:", t.end - t.start))
    
    if (verbose) print("findCluster: ending function")
    return(profileChr)

  }

