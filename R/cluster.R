# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: bioinfo-staff@curie.fr
# It is strictly forbidden to transfer, use or re-use this code 
# or part of it without explicit written authorization from Institut Curie.


cluster <- function(...)
{
	UseMethod("cluster")
}




cluster.default <- function(CGH, Cluster, region, labels, lambda, nmin, nmax, ...)
{

	#if (class(CGH)!="CGH")stop("wrong type for CGH")
	if(class(Cluster)!="hclust")stop("wrong type for Cluster")
	if(nmin>nmax)stop("nmin greater than nmax")
	#data <- CGH.data.frame(CGH)
	data <- CGH
	data <- data[which(data$OutliersTot==0),]

	#param <- c(d=6)
		
	IQRdiff <- function(y) IQR(diff(y))/1.908
	
	sigma <- IQRdiff(data$LogRatio)
	

	if (nmax>length(labels))
	{
		nmax <- length(labels)
	}

	
	logLike <- rep(NULL,nmax-nmin+1)

	i <- nmin
	for (i in nmin:nmax)
	{

		Classe <- data.frame(Classe=cutree(Cluster, k=i), Region=labels)
		tab <- merge(x=data,y=Classe, by.x=region, by.y="Region")
		aggregation <- aggregate(tab["LogRatio"],list(Classe=tab$Classe),mean)
		aggregation <- aggregation[order(aggregation$LogRatio),]
		deltaoversigma <- abs(diff(aggregation$LogRatio)/sigma)

	
		j <- 1
		minus2logL <- 0
		for (j in 1:i)	
		{
			indexClasse <- which(tab$Classe==j)
			if (length(indexClasse) > 1)
			{
				minus2logL <- minus2logL + length(indexClasse)*log(var(tab$LogRatio[indexClasse])*(length(indexClasse)-1)/length(indexClasse),exp(1))
				minus2logL <- minus2logL + length(indexClasse)*(1+log(2*pi,exp(1)))
			}
	
			#logLike[i] <- minus2logL + i*2*log(length(tab$LogRatio), exp(1))
			logLike[i-nmin+1] <- minus2logL + lambda*sum(kernelpen(deltaoversigma, ...))*log(length(data$LogRatio), exp(1))
		}
	
	}

	return(nmin+which(logLike==min(logLike))[1]-1)

}
