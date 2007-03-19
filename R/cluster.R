# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr


clusterglad <- function(...)
{
  UseMethod("clusterglad")
}




clusterglad <- function(Cluster, region, clusterRegion, lambda, nmin, nmax, sigma, ...)
{

  if(class(Cluster)!="hclust")stop("wrong type for Cluster")
  if(nmin>nmax)stop("nmin greater than nmax")

  if (nmin==nmax)
    {
      nmin <- min(length(Cluster$order),nmin)
      return(nmin)
    }

  
  
  mergeLike <- function(data)
    {
      nbobs <- sum(data$Card)
      barycentre <- sum(data$Mean*data$Card)
      barycentre <- barycentre/nbobs
      within <- sum(data$Card*data$Var)
      within <- within/nbobs
      between <- sum(data$Card*(data$Mean-barycentre)^2)
      between <- between/nbobs
      variance <- within + between

      if (nbobs==1)
        {
          res <- data.frame(logVar=0,Mean=barycentre)
        }
      else
        {
          logVar <- nbobs*(log(variance) + (1+log(2*pi)))
          res <- data.frame(logVar,Mean=barycentre)
        }
      
      return(res)
      
    }

  
  NbTotObs <- sum(clusterRegion$Card)

  

  if (nmax>length(clusterRegion[,1]))
    {
      nmax <- length(clusterRegion[,1])
    }

  
  logLike <- rep(0,nmax-nmin+1)

  for (i in nmin:nmax)
    {
      
      Classe <- data.frame(Classe=cutree(Cluster, k=i), Region=clusterRegion$Region)
      newtab <- merge(x=clusterRegion, y=Classe, by="Region")
      newtab <- by(newtab,newtab$Classe,mergeLike)
      Aux <- rep(0,attr(newtab,"dim"))
      newtabAux <- data.frame(logVar=Aux,Mean=Aux)

      for (l in 1:i)
        {
          newtabAux[l,] <- newtab[[l]]
        }
      
      newtabAux <- newtabAux[order(newtabAux$Mean),]

      deltaoversigma <- abs(diff(newtabAux$Mean)/sigma)

      logLike[i-nmin+1] <- sum(newtabAux$logVar) + lambda*sum(kernelpen(deltaoversigma, ...))*log(NbTotObs)
      
    }


  return(nmin+which(logLike==min(logLike))[1]-1)

}
