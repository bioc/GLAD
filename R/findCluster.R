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

    ## choix de la méthode de clustering
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method <- pmatch(profileChr$method, METHODS)

#    print(profileChr$method)
    if (is.na(method)) 
      stop("invalid clustering method")
    if (method == -1) 
      stop("ambiguous clustering method")
    


#    print(match.call())

    nbregion <- length(unique(profileChr$profileValues[,region]))
    
    t0.start <- Sys.time()
    


    if (nbregion == 1) 	
      {	
        profileChr$NbClusterOpt <- 1

        if (genome == FALSE)
          {
            profileChr$profileValues$ZoneChr <- 1
          }
        else
          {
            profileChr$profileValues$ZoneGen <- 1
          }

        t5.merge <- Sys.time()
      }	
    
    else	
      {

        subsetdata <- profileChr$profileValues[which(profileChr$profileValues$OutliersTot == 0),c(region, "LogRatio")]

        
### vérifier le comportement pour les clusters de cardinalité 1

        t1.subset <- Sys.time()

#        print("subset")
#        print(t1.subset - t0.start)

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




##         ## test  aggregation ###
## #        print(clusterRegion)
##         l <- length(profileChr$profileValues$LogRatio)
##         res <- .C("findCluster",
##                   as.double(profileChr$profileValues$LogRatio),
##                   as.integer(profileChr$profileValues[,region]),              
##                   as.integer(profileChr$profileValues$OutliersTot),
##                   ## paramètres pour clusterglad
##                   as.double(profileChr$findClusterSigma),
##                   as.integer(nbregion),
##                   as.integer(l),
##                   PACKAGE = "GLAD")

##         ## fin test ###

        
        t2.agg <- Sys.time()
        
##         print("aggregation")
##         print(t2.agg - t1.subset)

        
        sigma <- profileChr$findClusterSigma
        dist <- dist(clusterRegion$Mean)
        cluster.res <- hclustglad(dist, members = clusterRegion$Card, ...)

        t3.hclust <- Sys.time()
##         print("hclust")
##         print(t3.hclust - t2.agg)


        nbclasses <- clusterglad(Cluster = cluster.res, clusterRegion = clusterRegion, lambda = lambda, nmin = nmin, nmax = nmax, sigma = sigma, type = type, param = param)

        t4.cluster <- Sys.time()
##         print("cluster")
##         print(t4.cluster - t3.hclust)
        
        classes <- cutree(cluster.res, k=nbclasses)

##         ## ###########
##         ## test cutree
##         nbelt <- length(clusterRegion$Mean)
##         res.test <- .C("R_cutree",
##                        as.integer(cluster.res$merge),
##                        as.integer(nbclasses),
##                        classes = integer(nbelt), ## valeur de sortie
##                        as.integer(nbelt),
##                        PACKAGE = "GLAD")
##         print("VERIF")
##         print(which(res.test[["classes"]] != classes))
##         print("END VERIF")
##         ## fin test
##         ## ###########
        
        profileChr$NbClusterOpt <- nbclasses
        clusterRegion <- data.frame(clusterRegion, zone = classes)


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


        if (genome == FALSE)
          {
            profileChr$profileValues[,"ZoneChr"] <- myzone$zone
          }
        else
          {
            profileChr$profileValues[,"ZoneGen"] <- myzone$zone
          }

        t5.merge <- Sys.time()
##         print("merge")
##         print(t5.merge - t4.cluster)

      }	

    t.end <- Sys.time()
##     print("clustering")
##     print(t.end - t5.merge)

##     print(paste("Temps findCluster:", t.end - t0.start))
    
    if (verbose) print("findCluster: ending function")
    return(profileChr)

  }

