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

    if (verbose) print("findCluster: starting function")

    ## choix de la méthode de clustering
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method <- pmatch(profileChr$method, METHODS)


    if (is.na(method)) 
      stop("invalid clustering method")
    if (method == -1) 
      stop("ambiguous clustering method")
    
    

    l <- length(profileChr$profileValues[["LogRatio"]])
    myzone <- .C("findCluster",
                 as.double(profileChr$profileValues[["LogRatio"]]),
                 as.integer(profileChr$profileValues[[region]]),              
                 as.integer(profileChr$profileValues[["OutliersTot"]]),
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


    
    if (genome == FALSE)
      {
        profileChr$profileValues[["ZoneChr"]] <- myzone$zone
      }
    else
      {
        profileChr$profileValues[["ZoneGen"]] <- myzone$zone
      }

    
    if (verbose) print("findCluster: ending function")
    return(profileChr)

  }



clusterglad <- function(...)
{
  UseMethod("clusterglad")
}




clusterglad.hclust <- function(Cluster = NULL, clusterRegion = NULL, lambda = NULL, nmin = NULL, nmax = NULL, sigma = NULL, ...)
{

  if(class(Cluster) != "hclust")stop("wrong type for Cluster")
  if(nmin > nmax)stop("nmin greater than nmax")

  if (nmin == nmax)
    {
      nmin <- min(length(Cluster$order), nmin)
      return(nmin)
    }

  
  
  mergeLike <- function(data)
    {
      nbobs <- sum(data$Card)
      barycentre <- sum(data$Mean * data$Card)
      barycentre <- barycentre / nbobs
      within <- sum(data$Card * data$Var)
      within <- within / nbobs
      between <- sum(data$Card * (data$Mean - barycentre)^2)
      between <- between / nbobs
      variance <- within + between


      if (nbobs == 1)
        {
          res <- data.frame(logVar = 0,Mean = barycentre)
        }
      else
        {
          logVar <- nbobs * (log(variance) + (1+log(2 * pi)))
          res <- data.frame(logVar, Mean = barycentre)
        }
      
      return(res)
      
    }

  
  NbTotObs <- sum(clusterRegion$Card)

  

  if (nmax>length(clusterRegion[,1]))
    {
      nmax <- length(clusterRegion[,1])
    }

  
  logLike <- rep(0, nmax - nmin + 1)

  for (i in nmin:nmax)
    {
      
      Classe <- data.frame(Classe = cutree(Cluster, k = i), Region = clusterRegion$Region)
      newtab <- merge(x = clusterRegion, y = Classe, by = "Region")
      newtab <- by(newtab, newtab$Classe, mergeLike)
      Aux <- rep(0,attr(newtab, "dim"))
      newtabAux <- data.frame(logVar = Aux, Mean = Aux)

      for (l in 1:i)
        {
          newtabAux[l,] <- newtab[[l]]
        }

      newtabAux <- newtabAux[order(newtabAux$Mean),]

      deltaoversigma <- abs(diff(newtabAux$Mean) / sigma)

      if (i == 1)
        {
          sumkernelpen <- 1
        }
      else
        {
          sumkernelpen <- sum(kernelpen(deltaoversigma, ...))
        }

      logLike[i - nmin + 1] <- sum(newtabAux$logVar) + lambda * sumkernelpen * log(NbTotObs)
      
    }


  return(nmin + which(logLike == min(logLike))[1] - 1)

}
