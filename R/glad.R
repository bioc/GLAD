# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: bioinfo-staff@curie.fr
# It is strictly forbidden to transfer, use or re-use this code 
# or part of it without explicit written authorization from Institut Curie.



glad <- function(...)
{
  
  UseMethod("glad")
  
}


glad.profileCGH <- function(profileCGH, smoothfunc="aws", base=FALSE, sigma, bandwidth=10, round=2, lambdabreak=8, lambdacluster=8, lambdaclusterGen=40, type="tricubic", param=c(d=6), alpha=0.001, method="centroid", nmax=8, ...)
  {
 
	

    #profileCGH <- list(profileValues=data)	
    #class(profileCGH) <- "profileCGH"	
	
    # Breakpoints detection
    #profileCGH <- chrBreakpoints(profileCGH, smoothfunc="laws", lkern="Triangle", model="Gaussian", qlambda=0.999, base=FALSE, bandwidth=1, ...)		
    
    profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma, bandwidth=bandwidth, round=round, ...)
    
    # LogRatio are median-centered	
    median <- median(na.omit(profileCGH$profileValues$LogRatio)) 
    profileCGH$profileValues$LogRatio <- profileCGH$profileValues$LogRatio - median	
    profileCGH$profileValues$Smoothing <- profileCGH$profileValues$Smoothing - median	
	
	
    profileAux <- NULL	
    # profile by chromosome	
    nbzonetot <- 0 #total number of zones that have been previously identify	

    labelChr <- sort(unique(profileCGH$profileValues$Chromosome))	
    for (i in 1:length(labelChr))
      {
        indexChr <- which(profileCGH$profileValues$Chromosome==labelChr[i])	
        subset <- profileCGH$profileValues[indexChr,]	
	
        profileChr <- list(profileValues=subset)	
        class(profileChr) <- "profileChr"	
	
        profileChr <- removeBreakpoints(profileChr, lambda=lambdabreak, alpha=alpha, type=type, param=param)	
        profileChr <- detectOutliers(profileChr, region="Region", alpha=alpha)	
	
	# ça ne doit pas servir : à vérifier
        nmin <- 1	
        if (length(which(profileChr$profileValues$Breakpoints==1))>=1)
          {
            nmin <- 2
          }

        profileChr <- findCluster(profileChr, method=method, genome=FALSE, lambda=lambdacluster, nmin=1, nmax=nmax,type=type, param=param)
        profileChr <- detectOutliers(profileChr, region="ZoneChr", alpha=alpha)	
 	
        nbzone <- length(unique((profileChr$profileValues$ZoneChr[which(profileChr$profileValues$ZoneChr!=0)])))	
        profileChr$profileValues$ZoneChr <- profileChr$profileValues$ZoneChr + nbzonetot	
        nbzonetot <- nbzonetot + nbzone	
	
        profileAux <- rbind(profileAux, profileChr$profileValues)
      }
	

	
    profileCGH$profileValues <- profileAux	
	
    class(profileCGH) <- "profileChr"	

    profileCGH <- findCluster(profileCGH, region="ZoneChr", method=method, genome=TRUE, lambda=lambdaclusterGen, nmin=1, nmax=nmax, type=type, param=param)	
    class(profileCGH) <- "profileCGH"
    profileCGH <- affectationGNL(profileCGH)


    return(profileCGH)





}
