## Copyright (C) 2003 Institut Curie
## Author(s): Philippe Hupé (Institut Curie) 2003
## Contact: glad@curie.fr

findCluster <- function(...)
  {

    UseMethod("findCluster")
  }


findCluster.profileChr <- function(profileChr, region="Region", genome = TRUE,
                                   lambda = 10, nmin = 1, nmax = 10,
                                   type = "tricubic", param = c(d = 6), verbose = FALSE, ...)
  {

    print("param")
    print(param)
    
    doinR <- 0
    
    if (verbose) print("findCluster: starting function")

    ## choix de la méthode de clustering
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method <- pmatch(profileChr$method, METHODS)


    if (is.na(method)) 
      stop("invalid clustering method")
    if (method == -1) 
      stop("ambiguous clustering method")
    
    

    if(doinR)
      {

        nbregion <- length(unique(profileChr$profileValues[,region]))
        
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
          }
        else
          {
            

            subsetdata <- profileChr$profileValues[which(profileChr$profileValues$OutliersTot == 0),c(region, "LogRatio")]
            
            ## vérifier le comportement pour les clusters de cardinalité 1

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

            
            sigma <- profileChr$findClusterSigma
            dist <- dist(clusterRegion$Mean)
            cluster.res <- hclustglad(dist, members = clusterRegion$Card, ...)


            nbclasses <- clusterglad(Cluster = cluster.res, clusterRegion = clusterRegion, lambda = lambda, nmin = nmin, nmax = nmax, sigma = sigma, type = type, param = param)

            
            classes <- cutree(cluster.res, k=nbclasses)

            
            profileChr$NbClusterOpt <- nbclasses
            clusterRegion <- data.frame(clusterRegion, zone = classes)


            lengthDest <- length(profileChr$profileValues[,region])
            lengthSrc <- length(clusterRegion$Region)
            myzone <- .C("my_merge_int",
                         as.integer(profileChr$profileValues[,region]),
                         zone = integer(lengthDest),
                         as.integer(clusterRegion$Region),
                         as.integer(clusterRegion$zone),
                         as.integer(lengthDest),
                         as.integer(lengthSrc),
                         PACKAGE="GLAD")
          }
      }
    else
      {

        l <- length(profileChr$profileValues$LogRatio)
        myzone <- .C("findCluster",
                     as.double(profileChr$profileValues$LogRatio),
                     as.integer(profileChr$profileValues[,region]),              
                     as.integer(profileChr$profileValues$OutliersTot),
                     zone = integer(l),
                     as.integer(method),
                     ## paramètres pour clusterglad
                     as.double(profileChr$findClusterSigma),
                     as.double(param["d"]),
                     as.double(lambda),
                     as.integer(nmin),
                     as.integer(nmax),
                     nbclasses = integer(1),
                     as.integer(l),
                     PACKAGE = "GLAD")


        profileChr$NbClusterOpt <- myzone$nbclasses

      }
    
    if (genome == FALSE)
      {
        profileChr$profileValues[,"ZoneChr"] <- myzone$zone
      }
    else
      {
        profileChr$profileValues[,"ZoneGen"] <- myzone$zone
      }



    
    if (verbose) print("findCluster: ending function")
    return(profileChr)

  }

