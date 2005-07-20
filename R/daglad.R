### This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hup� (Institut Curie) 2005
### Contact: glad@curie.fr




daglad <- function(...)
  {  
    UseMethod("daglad")
  }



daglad.profileCGH <- function(profileCGH, mediancenter=FALSE, normalrefcenter=FALSE, genomestep=FALSE,
                              smoothfunc="lawsglad", lkern="Exponential", model="Gaussian",
                              qlambda=0.999,  bandwidth=10, sigma=NULL, base=FALSE, round=1.5,
                              lambdabreak=8, lambdaclusterGen=40, param=c(d=6), alpha=0.001, msize=5,
                              method="centroid", nmin=1, nmax=8,
                              amplicon=1, deletion=-5, deltaN=0.10,  forceGL=c(-0.15,0.15), nbsigma=3,
                              MinBkpWeight=0.35, CheckBkpPos=TRUE,
                              verbose=FALSE, ...)
  {


    
    if(smoothfunc=="lawsglad" & model!="Gaussian")
      stop("Error in daglad: for lawsglad, only Gaussian model is available")

    if(smoothfunc=="lawsglad" & base==TRUE)
      stop("Error in daglad: for lawsglad it is not possible to use the option base=TRUE. Choose laws smoothfunc instead.")

    
    if (smoothfunc!="lawsglad")
      {
        print(paste("You have chosen smoothfunc=", smoothfunc))
        print(paste("Choose smoothfunc=lawsglad if you want the process runs faster"))
      }


    ### Il faut des PosOrder uniques
    profileCGH$profileValues$NewPosOrder <- profileCGH$profileValues$PosOrder
    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$Chromosome,profileCGH$profileValues$PosOrder),]
    profileCGH$profileValues$PosOrder <- 1:length(profileCGH$profileValues[,1])
    


##############################################################################
###
### D�but de l'analyse daglad
###
##############################################################################

#######################################################
### AWS smoothing chromosome par chromosome
#######################################################    

    print("Smoothing for each Chromosome")
    profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma,
                                 bandwidth=bandwidth, lkern=lkern, model=model, qlambda=qlambda,
                                 round=round, verbose=verbose)


    
    PosOrderRange <- profileCGH$PosOrderRange


### Estimation des �cart-type pour chaque chromosome    
    IQRinfoC <- profileCGH$Sigma


### LogRatio are median-centered
### attention: ne pas appeler la variable median sinon
### l'utilisation de aggregate(..., ..., median)
### renvoie une erreur
    
    if (mediancenter)
      {
        med <- median(profileCGH$profileValues$LogRatio)
        profileCGH$profileValues$LogRatio <- profileCGH$profileValues$LogRatio - med
        profileCGH$profileValues$Smoothing <- profileCGH$profileValues$Smoothing - med
      }



#######################################################
### AWS smoothing sur l'ensemble du g�nome
#######################################################   

### il y a d�j� les noms Smoothing Level Region Breakpoints OutliersAws
### ces champs sont renomm�s pour qu'ils soient sp�cifiques
### de l'�tape chromosome par chromosome:
### ils portent � la fin la lettre C

    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Smoothing)

    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile=="Level")] <- "LevelC"
    namesprofile[which(namesprofile=="Region")] <- "RegionC"
    namesprofile[which(namesprofile=="Breakpoints")] <- "BreakpointsC"
    namesprofile[which(namesprofile=="OutliersAws")] <- "OutliersAwsC"
    namesprofile[which(namesprofile=="MinPosOrder")] <- "MinPosOrderC"
    namesprofile[which(namesprofile=="MaxPosOrder")] <- "MaxPosOrderC"
    names(profileCGH$profileValues) <- namesprofile


    if (genomestep)
      {

### pour le smoothing sur l'ensemble du g�nome, il ne faut qu'un seul chromosome
        profileCGH$profileValues$ChromosomeTrue <- profileCGH$profileValues$Chromosome
        profileCGH$profileValues$Chromosome <- 0

### AWS smoothing sur l'ensemble du g�nome
### ici, on impose un qlamba �lev� (l'utilisateur n'a pas le choix du param�tre)
### pour d�gager la tendance des donn�es

        print("Smoothing over the genome")
        profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=FALSE, sigma=sigma,
                                     bandwidth=bandwidth, round=round, verbose=verbose,
                                     lkern=lkern, model=model, qlambda=0.9999999)


        profileCGH$PosOrderRange <- PosOrderRange
        

### estimation de l'�cart-type sur l'ensemble du g�nome    
        IQRinfoG <- profileCGH$Sigma
        profileCGH$SigmaC <- IQRinfoC
        profileCGH$SigmaG <- IQRinfoG

### maintenant, le champ ChromosomeTrue ne sert plus � rien
        profileCGH$profileValues$Chromosome <- profileCGH$profileValues$ChromosomeTrue
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-ChromosomeTrue)
        
### les champs Smoothing Level Region Breakpoints OutliersAws
### sont renomm�s pour qu'ils soient sp�cifiques
### de l'�tape sur l'ensemble du g�nome:
### ils portent � la fin la lettre G


        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Smoothing)
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Region)
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Breakpoints)
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-MinPosOrder)
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-MaxPosOrder)

        namesprofile <- names(profileCGH$profileValues)
        namesprofile[which(namesprofile=="Level")] <- "LevelG"
        namesprofile[which(namesprofile=="MinPosOrderC")] <- "MinPosOrder"
        namesprofile[which(namesprofile=="MaxPosOrderC")] <- "MaxPosOrder"
        names(profileCGH$profileValues) <- namesprofile

      
        
########################################################
###
### on cherche la ligne de base � partir du
### smoothing sur l'ensemble du g�nome
###
##########################################################

### l'object doit etre de type profileChr
### � l'avenir, il faudra supprimer cette classe
### (attention � faire les modifs �galement dans la fonction glad)
        class(profileCGH) <- "profileChr"

### d�tection des outliers
        profileCGH <- detectOutliers(profileCGH, region="LevelG", alpha=alpha, msize=msize)

        
### le clustering est fait sur les niveaux LevelG

        print("Find cluster over the Genome")
        profileCGH$findClusterSigma <- profileCGH$SigmaG$Value[1]
        profileCGH <- findCluster(profileCGH, region="LevelG", method=method, genome=TRUE,
                                  lambda=lambdaclusterGen, nmin=nmin, nmax=nmax)
        class(profileCGH) <- "profileCGH"


### On a besoin d'un champ OutliersAws:
### on ne consid�re comme Outliers AWS que ceux d�tect�s        
### � l'�tape chromosome
        
        profileCGH$profileValues$OutliersAws <- profileCGH$profileValues$OutliersAwsC
                

        print("Find the Normal Baseline")
### on trouve le normal (ligne de base) � partir de la fonction affectationGNL
        profileCGH <- affectationGNL(profileCGH)

       
        
### On renomme le champ ZoneGNL en ZoneGNLGen
        namesprofile <- names(profileCGH$profileValues)
        namesprofile[which(namesprofile=="ZoneGNL")] <- "ZoneGNLGen"   
        names(profileCGH$profileValues) <- namesprofile
        
### notre niveau de r�f�rence pour le normal correspond
### � la m�diane des log-ratios correspondant � une ZoneGNLGen de 0
        NormalRef <- median(profileCGH$profileValues[which(profileCGH$profileValues$ZoneGNLGen==0),"LogRatio"], na.rm=TRUE)

        profileCGH$NormalRef <- NormalRef

        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-LevelG)
              

      }

    else
      {

        IQRdiff <- function(y) IQR(diff(y))/1.908
        profileCGH$profileValues$OutliersAws <- profileCGH$profileValues$OutliersAwsC
        profileCGH$profileValues$OutliersMad <- profileCGH$profileValues$OutliersAws
        profileCGH$profileValues$OutliersTot <- profileCGH$profileValues$OutliersAws
        profileCGH$profileValues$MinPosOrder <- profileCGH$profileValues$MinPosOrderC
        profileCGH$profileValues$MaxPosOrder <- profileCGH$profileValues$MaxPosOrderC
        
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-MinPosOrderC)
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-MaxPosOrderC)
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-OutliersAwsC)

        NormalRef <- 0
        profileCGH$NormalRef <- NormalRef
        profileCGH$SigmaC <- IQRinfoC
        profileCGH$profileValues$ZoneGen <- profileCGH$profileValues$PosOrder


        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]


        if (is.null(sigma))
          {
            IQRinfoG <- IQRdiff(profileCGH$profileValues$LogRatio)
          }
        else
          {
            IQRinfoG <- sigma
          }
        profileCGH$SigmaG <- data.frame(Chromosome=0,Value=IQRinfoG)
        profileCGH$findClusterSigma <- profileCGH$SigmaG$Value[1]
      }

    
##################################################
#### fin de l'�tape sur le g�nome
##############################################

       
    
### On va appliquer la fonction removeLevel
### l'objet doit donc �tre de type profileChr
    class(profileCGH) <- "profileChr"


    profileCGH$profileValues$Level <- profileCGH$profileValues$LevelC    
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-LevelC)

### il faut le champ Breakpoints et Region
### mais on peut supprimer son utilisation dans removeLevel
    profileCGH$profileValues$Breakpoints <- profileCGH$profileValues$BreakpointsC
    profileCGH$profileValues$Region <- profileCGH$profileValues$RegionC
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-RegionC)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-BreakpointsC)
      

    indice <- 1:length(profileCGH$profileValues[,1])
    ChrIndice <- split(indice,profileCGH$profileValues$Chromosome)
    NbChr <- length(names(ChrIndice))

    Init <- 0       

    profileCGH$profileValues <- data.frame(profileCGH$profileValues, NextLogRatio=Init)
    FieldOrder <- names(profileCGH$profileValues)
    
    print("Optimization of the Breakpoints")
    for (i in 1:NbChr)
      {
        indexChr <- ChrIndice[[i]]
        subset <- profileCGH$profileValues[indexChr,]	
        profileChr <- list(profileValues=subset)	
        class(profileChr) <- "profileChr"
	profileChr$findClusterSigma <- profileCGH$SigmaC$Value[i]
        #profileChr <- detectOutliers(profileChr, region="Level", alpha=alpha, msize=msize)       
        profileChr <- removeLevel(profileChr, lambda=lambdabreak, param=param, alpha=alpha, msize=msize, verbose=verbose)
        #profileChr <- detectOutliers(profileChr, region="Level", alpha=alpha, msize=msize)       
        profileCGH$profileValues[indexChr,] <- profileChr$profileValues[,FieldOrder]
      }


 
    agg <- aggregate(profileCGH$profileValues$LogRatio, list(Level=profileCGH$profileValues$Level), median)
    agg$Level <- as.numeric(as.character(agg$Level))
    names(agg) <- c("Level","Smoothing")
    profileCGH$profileValues <- merge(profileCGH$profileValues, agg, by="Level", all=TRUE)    
    class(profileCGH) <- "profileCGH"
    
    

### pour le calcul des poids, l'�cart-type utilis�
### est celui calcul� sur l'ensemble de g�nome    
    Sigma <- profileCGH$SigmaG$Value[1]    
    profileCGH$profileValues <- data.frame(profileCGH$profileValues, NormalRef=NormalRef, Sigma=Sigma)
    profileCGH$profileValues$DiffBase <- profileCGH$profileValues$Smoothing - profileCGH$profileValues$NormalRef

    
### on prend comme r�f�rence ceux qui sont compris entre certaines valeurs
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-ZoneGen)
    indexNormalLevel <- which(abs(profileCGH$profileValues$DiffBase)<=deltaN)
    profileCGH$profileValues$NormalRange <- profileCGH$profileValues$Level
    profileCGH$profileValues$NormalRange[indexNormalLevel] <- 0


    
### le clustering est fait sur les niveaux NormalRange
    class(profileCGH) <- "profileChr"
    profileCGH <- findCluster(profileCGH, region="NormalRange", method=method, genome=TRUE,
                              lambda=lambdaclusterGen, nmin=nmin, nmax=nmax, verbose=verbose)
    class(profileCGH) <- "profileCGH"


### le cluster correspondant au normal est celui qui comprend
### le NormalRange 0
    indexNormalRange <- which(profileCGH$profileValues$NormalRange==0)
    NormalCluster <- unique(profileCGH$profileValues$ZoneGen[indexNormalRange])
    MedianCluster <- aggregate(profileCGH$profileValues$LogRatio, list(ZoneGen=profileCGH$profileValues$ZoneGen),median,na.rm=TRUE)
    MedianCluster$ZoneGen <- as.numeric(as.character(MedianCluster$ZoneGen))
    names(MedianCluster) <- c("ZoneGen","Median")
    RefNorm <- MedianCluster$Median[which(MedianCluster$ZoneGen==NormalCluster)]
    MedianCluster$ZoneGNL <- rep(0,length(MedianCluster[,1]))
    indexClusterGain <- which(MedianCluster$Median>RefNorm)
    MedianCluster$ZoneGNL[indexClusterGain] <- 1
    indexClusterLost <- which(MedianCluster$Median<RefNorm)
    MedianCluster$ZoneGNL[indexClusterLost] <- -1
    profileCGH$profileValues <- merge(profileCGH$profileValues, MedianCluster[,c("ZoneGen","ZoneGNL")], all=TRUE, by="ZoneGen")


   
    
### on force les gains et les pertes pour certaines valeur de smoothing
    indexForceGain <- which(profileCGH$profileValues$Smoothing - NormalRef>=forceGL[2])
    profileCGH$profileValues$ZoneGNL[indexForceGain] <- 1
    indexForceLost <- which(profileCGH$profileValues$Smoothing - NormalRef<=forceGL[1])
    profileCGH$profileValues$ZoneGNL[indexForceLost] <- -1

### Amplicon et deletion
    indexAmp <- which(profileCGH$profileValues$Smoothing - NormalRef>= amplicon)
    profileCGH$profileValues$ZoneGNL[indexAmp] <- 2
    indexDel <- which(profileCGH$profileValues$Smoothing - NormalRef<= deletion)
    profileCGH$profileValues$ZoneGNL[indexDel] <- -10
    



### Calcul d'un poids pour les Breakpoints
### Attention: comme on calcul une variable GNLchange
### il ne faut pas le statut des outliers pour d�finir le changement de statut
### car un outlier peut �tre aussi un breakpoints


    profileCGH$nbsigma <- nbsigma
    profileCGH$BkpInfo <- BkpInfo(profileCGH)
    

    
    RecomputeGNL <- FALSE
### Epuration des Breakpoints cons�cutifs
    if (is.data.frame(profileCGH$BkpInfo))
      {
        profileCGH$BkpInfo <- profileCGH$BkpInfo[order(profileCGH$BkpInfo$PosOrder),]
        profileCGH$BkpInfo$BkpToDel <- rep(0, length(profileCGH$BkpInfo[,1]))
        profileCGH$BkpInfo$Side <- profileCGH$BkpInfo$BkpToDel
        profileCGH$BkpInfo$NextPosOrder <- profileCGH$BkpInfo$PosOrder + 1

### on regarde que si il y a au moins 2 Bkp
        if (length(profileCGH$BkpInfo[,1])>1)
          {
            for (i in 2:length(profileCGH$BkpInfo[,1]))
              {
                if (profileCGH$BkpInfo$PosOrder[i]==profileCGH$BkpInfo$NextPosOrder[i-1] & profileCGH$BkpInfo$BkpToDel[i-1]==0)
                  {
                    SigChr <- profileCGH$BkpInfo$Sigma[i]
### on regarde d'abord � gauche
                    UpLeft <- profileCGH$BkpInfo$Smoothing[i-1] + 3*SigChr
                    LowLeft <- profileCGH$BkpInfo$Smoothing[i-1] - 3*SigChr

### On regarde ce qui se passe � droite
                    UpRight <- profileCGH$BkpInfo$SmoothingNext[i] + 3*SigChr
                    LowRight <- profileCGH$BkpInfo$SmoothingNext[i] - 3*SigChr

                    LRV <- profileCGH$BkpInfo$LogRatio[i]
                    
                    if ((LRV > LowLeft & LRV < UpLeft) || (LRV > LowRight & LRV < UpRight))
                      {
                                        #print("des Bkp vont �tre supprim�s")
                        RecomputeGNL <- TRUE
                        if ((LRV > LowLeft & LRV < UpLeft) & (LRV > LowRight & LRV < UpRight))
                          {
### attention, lors de la suppression d'un Bkp, le gap n'est plus bon
### d'o� recalcul du weight
### je ne vois pas o� il est recalcul�!!!!
                            DiffLeft <- abs(LRV - profileCGH$BkpInfo$Smoothing[i-1])
                            DiffRight <- abs(LRV - profileCGH$BkpInfo$SmoothingNext[i])
                            if (DiffRight<DiffLeft)
                              {
### On fusionne le Bkp avec la r�gion � droite
                                profileCGH$BkpInfo$BkpToDel[i] <- 1
                                profileCGH$BkpInfo$Side[i] <- 1
                                profileCGH$BkpInfo$Gap[i-1] <- abs(profileCGH$BkpInfo$Smoothing[i-1]-profileCGH$BkpInfo$SmoothingNext[i])
                                profileCGH$BkpInfo$Weight[i-1] <- 1 - kernelpen(profileCGH$BkpInfo$Gap[i-1], param=c(d=nbsigma*profileCGH$BkpInfo$Sigma[i-1]))
                                
                              }
                            else
                              {
### On fusionne le Bkp avec la r�gion � gauche
                                profileCGH$BkpInfo$BkpToDel[i-1] <- 1
                                profileCGH$BkpInfo$Side[i-1] <- 0
                                profileCGH$BkpInfo$Gap[i] <- abs(profileCGH$BkpInfo$Smoothing[i-1]-profileCGH$BkpInfo$SmoothingNext[i])
                                profileCGH$BkpInfo$Weight[i] <- 1 - kernelpen(profileCGH$BkpInfo$Gap[i-1], param=c(d=nbsigma*profileCGH$BkpInfo$Sigma[i-1]))
                                

                              }
                          }
                        else
                          {
                            if ((LRV > LowLeft & LRV < UpLeft))
                              {
### On fusionne le Bkp avec la r�gion � gauche
                                profileCGH$BkpInfo$BkpToDel[i-1] <- 1
                                profileCGH$BkpInfo$Side[i-1] <- 0
                                profileCGH$BkpInfo$Gap[i] <- abs(profileCGH$BkpInfo$Smoothing[i-1]-profileCGH$BkpInfo$SmoothingNext[i])
                                profileCGH$BkpInfo$Weight[i] <- 1 - kernelpen(profileCGH$BkpInfo$Gap[i-1], param=c(d=nbsigma*profileCGH$BkpInfo$Sigma[i-1]))
                                
                              }
                            else
                              {
### On fusionne le Bkp avec la r�gion � droite
                                profileCGH$BkpInfo$BkpToDel[i] <- 1
                                profileCGH$BkpInfo$Side[i] <- 1
                                profileCGH$BkpInfo$Gap[i-1] <- abs(profileCGH$BkpInfo$Smoothing[i-1]-profileCGH$BkpInfo$SmoothingNext[i])
                                profileCGH$BkpInfo$Weight[i-1] <- 1 - kernelpen(profileCGH$BkpInfo$Gap[i-1], param=c(d=nbsigma*profileCGH$BkpInfo$Sigma[i-1]))                            
                              }
                          }                                        
                      }
                  }
              }
          }


        indexBPtoDel <- which(profileCGH$BkpInfo$BkpToDel==1)
        if (length(indexBPtoDel)>0)
          {
                                       
            profile <- profileCGH$profileValues
            profile <- profile[order(profile$PosOrder),]

            for (ind in indexBPtoDel)
              {
                
                profile <- profile[order(profile$PosOrder),]
                #indextochange <- which(profile$PosOrder==profileCGH$BkpInfo$PosOrder[ind])
                indextochange <- profileCGH$BkpInfo$PosOrder[ind]
                profile$Breakpoints[indextochange] <- 0
                
                if (profileCGH$BkpInfo$Side[ind]==0)
                  {
                    profile$Level[indextochange+1] <- profile$Level[indextochange]
                    profile$Smoothing[indextochange+1] <- profile$Smoothing[indextochange]
                  }
                else
                  {
                    profile$Level[indextochange] <- profile$Level[indextochange+1]
                    profile$Smoothing[indextochange] <- profile$Smoothing[indextochange+1] 
                  }
                
              }


            profileCGH$profileValues <- profile

            profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexBPtoDel,]
          }

        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-Side)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-BkpToDel)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-NextPosOrder)

      }

    

    if (RecomputeGNL)
      {
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-ZoneGNL)

        profileCGH$profileValues$DiffBase <- profileCGH$profileValues$Smoothing - profileCGH$profileValues$NormalRef
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-ZoneGen)
        
        indexNormalLevel <- which(abs(profileCGH$profileValues$DiffBase)<=deltaN)
        profileCGH$profileValues$NormalRange <- profileCGH$profileValues$Level
        profileCGH$profileValues$NormalRange[indexNormalLevel] <- 0

        
### le clustering est fait sur les niveaux NormalRange
        class(profileCGH) <- "profileChr"
        profileCGH <- findCluster(profileCGH, region="NormalRange", method=method, genome=TRUE,
                                  lambda=lambdaclusterGen, nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)
        class(profileCGH) <- "profileCGH"

### le cluster correspondant au normal est celui qui comprend
### le NormalRange 0
        indexNormalRange <- which(profileCGH$profileValues$NormalRange==0)
        NormalCluster <- unique(profileCGH$profileValues$ZoneGen[indexNormalRange])
        MedianCluster <- aggregate(profileCGH$profileValues$LogRatio, list(ZoneGen=profileCGH$profileValues$ZoneGen),median,na.rm=TRUE)
        MedianCluster$ZoneGen <- as.numeric(as.character(MedianCluster$ZoneGen))
        names(MedianCluster) <- c("ZoneGen","Median")
        RefNorm <- MedianCluster$Median[which(MedianCluster$ZoneGen==NormalCluster)]
        MedianCluster$ZoneGNL <- rep(0,length(MedianCluster[,1]))
        indexClusterGain <- which(MedianCluster$Median>RefNorm)
        MedianCluster$ZoneGNL[indexClusterGain] <- 1
        indexClusterLost <- which(MedianCluster$Median<RefNorm)
        MedianCluster$ZoneGNL[indexClusterLost] <- -1
        profileCGH$profileValues <- merge(profileCGH$profileValues, MedianCluster[,c("ZoneGen","ZoneGNL")], all=TRUE, by="ZoneGen")

### on force les gains et les pertes pour certaines valeur de smoothing
        indexForceGain <- which(profileCGH$profileValues$Smoothing - NormalRef>=forceGL[2])
        profileCGH$profileValues$ZoneGNL[indexForceGain] <- 1
        indexForceLost <- which(profileCGH$profileValues$Smoothing - NormalRef<=forceGL[1])
        profileCGH$profileValues$ZoneGNL[indexForceLost] <- -1

### Amplicon et deletion
        indexAmp <- which(profileCGH$profileValues$Smoothing - NormalRef>= amplicon)
        profileCGH$profileValues$ZoneGNL[indexAmp] <- 2
        indexDel <- which(profileCGH$profileValues$Smoothing - NormalRef<= deletion)
        profileCGH$profileValues$ZoneGNL[indexDel] <- -10

      }


    

### Statut des Outliers
    profileCGH <- OutliersGNL(profileCGH, alpha=alpha, sigma=Sigma, NormalRef=NormalRef, amplicon=amplicon, deletion=deletion)    
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-NormalRange)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-ZoneGen)
    

### R�cup�ration des informations sur l'analyse
    profileCGH <- profileCGH[-which(names(profileCGH)=="Sigma")]
### Les donn�es sont centr�es sur NormalRef

    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
    rownames(profileCGH$profileValues) <- 1:length(profileCGH$profileValues$LogRatio)

    class(profileCGH) <- "profileCGH"

    profileCGH$alpha <- alpha
    profileCGH$msize <- msize
    profileCGH$amplicon <- amplicon
    profileCGH$deletion <- deletion
    profileCGH$deltaN <- deltaN
    profileCGH$method <- method
    profileCGH$lambdaclusterGen <- lambdaclusterGen 
    profileCGH$nmax <- nmax
    profileCGH$nmin <- nmin
    profileCGH$forceGL <- forceGL


    
#################################################################################
###
###  Bkp filter
###
#################################################################################        

    profileCGH <- filterBkpStep(profileCGH, MinBkpWeight=MinBkpWeight)

#################################################################################
###
###  D�placement des Bkp
###
#################################################################################    

    if (CheckBkpPos)
      {
        profileCGH <- MoveBkpStep(profileCGH)        
      }


### fin comment
### suppression des champs inutiles

    if (genomestep)
      {
        profileCGH$profileValues$NormalRef <- profileCGH$profileValues$ZoneGNLGen
        profileCGH$profileValues <- subset(profileCGH$profileValues, select=-ZoneGNLGen)
      }

    
    if (normalrefcenter)
      {
        profileCGH$profileValues$Smoothing <- profileCGH$profileValues$Smoothing - NormalRef
        profileCGH$profileValues$LogRatio <- profileCGH$profileValues$LogRatio - NormalRef
        
        
        if (is.data.frame(profileCGH$BkpInfo))
          {
            profileCGH$BkpInfo$LogRatio <- profileCGH$BkpInfo$LogRatio - NormalRef
            profileCGH$BkpInfo$Smoothing <- profileCGH$BkpInfo$Smoothing - NormalRef
            profileCGH$BkpInfo$SmoothingNext <- profileCGH$BkpInfo$SmoothingNext - NormalRef
          }

        profileCGH$NormalRef <- 0
      }


    if (is.data.frame(profileCGH$BkpInfo))
      {
        profileCGH$BkpInfo <- merge(profileCGH$BkpInfo, profileCGH$profileValues[,c("PosOrder","NewPosOrder")],
                                    by="PosOrder",all.x=TRUE)
        PosAux <- profileCGH$BkpInfo$PosOrder
        profileCGH$BkpInfo$PosOrder <- profileCGH$BkpInfo$NewPosOrder
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-NewPosOrder)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-SmoothingNext)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-NextLogRatio)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-MaxPosOrder)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-MinPosOrder)
        profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=-ZoneGNLnext)
      }

### il faut supprimer le champ r�gion que ne sert plus � rien
### ou bien s'assurer qu'il est � jour

## il faut le champ ZoneGen ou un �quivalent
    ### remettre les champs dans l'ordre
    
    PosAux <- profileCGH$profileValues$PosOrder
    profileCGH$profileValues$PosOrder <- profileCGH$profileValues$NewPosOrder
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-NewPosOrder)
#    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-NormalRef)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Sigma)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-DiffBase)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-Region)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-NextLogRatio)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-MaxPosOrder)
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=-MinPosOrder)






    at <- setdiff(attributes(profileCGH)$names,c("PosOrderRange","findClusterSigma","NbClusterOpt"))
    profileCGH <- profileCGH[at]
    class(profileCGH) <- "profileCGH"
    
        
    return (profileCGH)
    
    
  }

