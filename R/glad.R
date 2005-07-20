### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr
### http://bioinfo.curie.fr



glad <- function(...)
  {   
    UseMethod("glad")    
  }


glad.profileCGH <- function(profileCGH, mediancenter=FALSE,
                            smoothfunc="lawsglad", bandwidth=10, round=1.5,
                            model="Gaussian", lkern="Exponential", qlambda=0.999,
                            base=FALSE, sigma=NULL,
                            lambdabreak=8, lambdacluster=8, lambdaclusterGen=40,
                            type="tricubic", param=c(d=6),
                            alpha=0.001, msize=5,
                            method="centroid", nmax=8,
                            verbose=FALSE, ...)
  {
    

### champs tels qu'ils sont en entrée
    fieldinput <- names(profileCGH$profileValues)
    
    if (verbose)
      {
        print("GLAD: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }


    if(smoothfunc=="lawsglad" & model!="Gaussian")
      stop("Error in glad: for lawsglad, only Gaussian model is available")

    if(smoothfunc=="lawsglad" & base==TRUE)
      stop("Error in glad: for lawsglad it is not possible to use the option base=TRUE. Choose laws smoothfunc instead.")

    
    if (smoothfunc!="lawsglad")
      {
        print(paste("You have chosen smoothfunc=", smoothfunc))
        print(paste("Choose smoothfunc=lawsglad if you want the process runs faster"))
      }

    ### Méthode d'estimation du sigma

    profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma, bandwidth=bandwidth, round=round, verbose=verbose, model=model, lkern=lkern, qlambda=qlambda)
### LogRatio are median-centered
    if (mediancenter)
      {
        med <- median(profileCGH$profileValues$LogRatio) 
        profileCGH$profileValues$LogRatio <- profileCGH$profileValues$LogRatio - med
        profileCGH$profileValues$Smoothing <- profileCGH$profileValues$Smoothing - med
      }
    
    
    
### profile by chromosome	
    nbzonetot <- 0 #total number of zones that have been previously identify	



    indice <- 1:length(profileCGH$profileValues[,1])
    ChrIndice <- split(indice,profileCGH$profileValues$Chromosome)
    NbChr <- length(names(ChrIndice))

    Init <- indice
    Init <- 0

    profileCGH$profileValues <- data.frame(profileCGH$profileValues, NextLogRatio=Init, OutliersMad=Init, OutliersTot=Init, ZoneChr=Init)
    FieldOrder <- names(profileCGH$profileValues)
    
    for (i in 1:NbChr)
      {
        indexChr <- ChrIndice[[i]]
        subset <- profileCGH$profileValues[indexChr,]	
	
        profileChr <- list(profileValues=subset)	
        class(profileChr) <- "profileChr"	

        profileChr$findClusterSigma <- profileCGH$Sigma$Value[i]

        profileChr <- removeBreakpoints(profileChr, lambda=lambdabreak, alpha=alpha, msize=msize,
                                        type=type, param=param, verbose=verbose)
        
        #profileChr <- detectOutliers(profileChr, region="Region", alpha=alpha, msize=msize, verbose=verbose)	
	
### ça ne doit pas servir : à vérifier
        nmin <- 1	
        if (length(which(profileChr$profileValues$Breakpoints==1))>=1)
          {
            nmin <- 2
          }

        
        profileChr <- findCluster(profileChr, method=method, genome=FALSE,
                                  lambda=lambdacluster, nmin=1, nmax=nmax,type=type,
                                  param=param, verbose=verbose)
        
        profileChr <- detectOutliers(profileChr, region="ZoneChr", alpha=alpha, msize=msize)	
 	
        nbzone <- length(unique((profileChr$profileValues$ZoneChr[which(profileChr$profileValues$ZoneChr!=0)])))	
        profileChr$profileValues$ZoneChr <- profileChr$profileValues$ZoneChr + nbzonetot	
        nbzonetot <- nbzonetot + nbzone

        profileCGH$profileValues[indexChr,] <- profileChr$profileValues[,FieldOrder]

      }
    

    
    class(profileCGH) <- "profileChr"

    if (is.null(sigma))
      {

        IQRdiff <- function(y) IQR(diff(y))/1.908
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$Chromosome, profileCGH$profileValues$PosOrder),]
        profileCGH$findClusterSigma <- IQRdiff(profileCGH$profileValues$LogRatio)
                
        
      }
    
    else
      {
        profileCGH$findClusterSigma <- sigma^0.5
      }

    if (verbose) print("GLAD: starting clustering for whole genome")
    profileCGH <- findCluster(profileCGH, region="ZoneChr", method=method, genome=TRUE,
                              lambda=lambdaclusterGen, nmin=1, nmax=nmax, type=type,
                              param=param, verbose=verbose)
    if (verbose)
      {
        print("GLAD: ending clustering for whole genome")
        print("GLAD: starting affectationGNL")      
      }
    
    class(profileCGH) <- "profileCGH"
    profileCGH <- affectationGNL(profileCGH, verbose=verbose)

    if (verbose)
      {
        print("GLAD: ending affectationGNL")
        print("GLAD: ending function")
      }

    ### recalcul du Smoothing:
    ### on prend la médiane par région

    agg <- aggregate(profileCGH$profileValues$LogRatio, list(Region=profileCGH$profileValues$Region), median)
    names(agg) <- c("Region","Smoothing")
    agg$Region <- as.numeric(as.character(agg$Region))
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Smoothing)

    profileCGH$profileValues <- merge(profileCGH$profileValues, agg)

    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]

### champs ajoutés par glad
    fieldglad <- c("Smoothing","Region","Level","OutliersAws","Breakpoints","OutliersMad","OutliersTot","ZoneChr","ZoneGen","ZoneGNL")

### ordre des champs
    fieldorder <- c(fieldinput, fieldglad)
    profileCGH$profileValues <- profileCGH$profileValues[, fieldorder]


    indexBkp <- which(profileCGH$profileValues$Breakpoints==1)
    if (length(indexBkp)>0)
      {
        nomchamp <- c("PosOrder","PosBase","Clone","Chromosome")
        BkpInfo <- profileCGH$profileValues[indexBkp,intersect(nomchamp,names(profileCGH$profileValues))]
        profileCGH$BkpInfo <- BkpInfo
      }

    nomliste <- names(profileCGH)
    nomliste[which(nomliste=="Sigma")] <- "SigmaC"
    names(profileCGH) <- nomliste
    

    return(profileCGH)

  }






