# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: bioinfo-staff@curie.fr
# It is strictly forbidden to transfer, use or re-use this code 
# or part of it without explicit written authorization from Institut Curie.

findCluster <- function(...)
{
	UseMethod("findCluster")
}

findCluster.profileChr <- function(profileChr, region="Region", genome=TRUE, lambda=10, nmin=1, nmax=10, type="tricubic", param=c(d=6), ...)
{

		
	CGH <- profileChr$profileValues
	#CGHna <- CGH[attr(na.omit(CGH),"na.action"),]


        indexNA <- attr(na.omit(CGH),"na.action")
        if (length(indexNA)<length(CGH[,1]))
          {
	CGHna <- CGH[indexNA,]
	if (!is.null(indexNA))CGH <- CGH[-indexNA,]

        

	subset <- CGH[which(CGH$OutliersTot==0),]
	
	j <- 1	
	labelRegion <- unique(subset[,region])	
	freq <- rep(0, length(labelRegion))	
	std <- freq	
	reg <- mean <- freq	
	
	for (j in 1:length(labelRegion))	
	{	
		indexRegion <- which(subset[,region]==labelRegion[j])	
		freq[j] <- length(indexRegion)	
		if (length(indexRegion)>1)	
		{	
			std[j] <- var(subset$LogRatio[indexRegion])^0.5	
		}	
		reg[j] <- labelRegion[j]	
		mean[j] <- mean(subset$LogRatio[indexRegion])	
	}	
	
	clusterRegion <- data.frame(reg, freq, std, mean)
	class(clusterRegion) <- "Cluster"
	zone <- rep(1, length(CGH$PosOrder))
				
	if (length(labelRegion)==1) 	
	{	
		nbclasses <- 1	
	}	
	
	else	
	{	
		dist <- dist(clusterRegion$mean)	
		cluster <- hclustglad(dist, members=clusterRegion$freq, ...)
		#if(genome==TRUE)
		#	plclust(cluster,  labels=clusterRegion$reg, hang=-1, main=paste("Genome clustering on",region))	
		
		#subset <- data.frame.CGH(subset)
		nbclasses <- cluster(subset, cluster, region, clusterRegion$reg, lambda, nmin, nmax, type, param)	
	
		classes <- cutree(cluster, k=nbclasses)

		j <- 1	
		for (j in 1:length(classes))	
		{	
			zone[which(CGH[,region]==clusterRegion$reg[j])] <- classes[j]	
		}	
	}	

	if (genome==FALSE)
	{
		CGH$ZoneChr <- zone
		CGHna$ZoneChr <- rep(NA,length(CGHna$PosOrder))
	}
	else
	{
		CGH$ZoneGen <- zone
		CGHna$ZoneGen <- rep(NA,length(CGHna$PosOrder))

	}
	
		
	CGH <- rbind(CGH, CGHna)
	profileChr$profileValues <- CGH
	return(profileChr)
      }

        else
          {

 	if (genome==FALSE)
	{
		CGH$ZoneChr <- rep(NA,length(CGH$PosOrder))
	}
	else
	{
		CGH$ZoneGen <- rep(NA,length(CGH$PosOrder))

	}
        
	profileChr$profileValues <- CGH
        return(profileChr)
	           
          }
	 
}
