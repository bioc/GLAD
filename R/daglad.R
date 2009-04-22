## This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
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
                              MinBkpWeight=0.35, CheckBkpPos=TRUE, assignGNLOut=TRUE,
                              breaksFdrQ = 0.0001, haarStartLevel = 1, haarEndLevel = 5,
                              verbose=FALSE, ...)
  {


    profilage <- FALSE
##############################################################################
###
### Préparation des données
###
##############################################################################

    if(profilage) Rprof("/tmp/Step00DataPreparation.dat")
### vérification de l'object donné en entrée

    if (verbose) print("daglad - step CheckData")
    CheckData(profileCGH, bandwidth=bandwidth, smoothfunc=smoothfunc)

    if(base==TRUE)
      {
        if(!require(aws))
          {
            stop("Error in daglad: the aws package is required to use these function. The aws package can be installed from http://www.r-project.org")
          }
      }

    
    if(smoothfunc=="lawsglad" & model!="Gaussian")
      stop("Error in daglad: for lawsglad, only Gaussian model is available")

    if(smoothfunc=="lawsglad" & base==TRUE)
      stop("Error in daglad: for lawsglad it is not possible to use the option base=TRUE. Choose laws smoothfunc instead.")

    
    if (smoothfunc!="lawsglad" && smoothfunc!="haarseg")      
      {
        print(paste("You have chosen smoothfunc=", smoothfunc))
        print(paste("Choose smoothfunc=lawsglad or smoothfunc=haarseg if you want the process runs faster"))
      }

    inputfields <- names(profileCGH$profileValues)
    excdudefields <- c("Level", "OutliersAws", "OutliersMad",
                       "OutliersTot", "Breakpoints", "Smoothing",
                       "NormalRef", "ZoneGNL")
    fieldstodel <- intersect(inputfields, excdudefields)
    inputfields <- setdiff(inputfields,fieldstodel)
    profileCGH$profileValues <- profileCGH$profileValues[,inputfields]

### Il faut des PosOrder uniques
    profileCGH$profileValues$NewPosOrder <- profileCGH$profileValues$PosOrder
    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$Chromosome,profileCGH$profileValues$PosOrder),]
    profileCGH$profileValues$PosOrder <- 1:length(profileCGH$profileValues[,1])

    if(profilage)Rprof(NULL)

##############################################################################
###
### Début de l'analyse daglad
###
##############################################################################

#######################################################
### AWS smoothing chromosome par chromosome
#######################################################    

    print("Smoothing for each Chromosome")
    if(profilage) Rprof("/tmp/Step01DataSegmentation.dat")    
    profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma,
                                 bandwidth=bandwidth, lkern=lkern, model=model, qlambda=qlambda,
                                 round=round, verbose=verbose,
                                 breaksFdrQ=breaksFdrQ, haarStartLevel=haarStartLevel , haarEndLevel=haarEndLevel)


    if(profilage)Rprof(NULL)
    
    PosOrderRange <- profileCGH$PosOrderRange


### Estimation des écart-type pour chaque chromosome    
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
### AWS smoothing sur l'ensemble du génome
#######################################################   

### il y a déjà les noms Smoothing Level Region Breakpoints OutliersAws
### ces champs sont renommés pour qu'ils soient spécifiques
### de l'étape chromosome par chromosome:
### ils portent à la fin la lettre C

    if(profilage) Rprof("/tmp/Step02DataPreparation.dat")
    profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"Smoothing")]
###    profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Smoothing"))

    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile=="Level")] <- "LevelC"
    namesprofile[which(namesprofile=="Region")] <- "RegionC"
    namesprofile[which(namesprofile=="Breakpoints")] <- "BreakpointsC"
    namesprofile[which(namesprofile=="OutliersAws")] <- "OutliersAwsC"
    namesprofile[which(namesprofile=="MinPosOrder")] <- "MinPosOrderC"
    namesprofile[which(namesprofile=="MaxPosOrder")] <- "MaxPosOrderC"
    names(profileCGH$profileValues) <- namesprofile

    if(profilage) Rprof(NULL)    


    if (genomestep)
      {

        if (verbose) print("daglad - step genomestep")        
### pour le smoothing sur l'ensemble du génome, il ne faut qu'un seul chromosome
        profileCGH$profileValues$ChromosomeTrue <- profileCGH$profileValues$Chromosome
        profileCGH$profileValues$Chromosome <- 0

### AWS smoothing sur l'ensemble du génome
### ici, on impose un qlamba élevé (l'utilisateur n'a pas le choix du paramètre)
### pour dégager la tendance des données

        print("Smoothing over the genome")
        profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=FALSE, sigma=sigma,
                                     bandwidth=bandwidth, round=round, verbose=verbose,
                                     lkern=lkern, model=model, qlambda=0.9999999,
                                     breaksFdrQ=breaksFdrQ, haarStartLevel=haarStartLevel , haarEndLevel=haarEndLevel)


        profileCGH$PosOrderRange <- PosOrderRange
        

### estimation de l'écart-type sur l'ensemble du génome    
        IQRinfoG <- profileCGH$Sigma
        profileCGH$SigmaC <- IQRinfoC
        profileCGH$SigmaG <- IQRinfoG

### maintenant, le champ ChromosomeTrue ne sert plus à rien
###        profileCGH$profileValues$Chromosome <- profileCGH$profileValues$ChromosomeTrue
        namesprofile <- names(profileCGH$profileValues)
        namesprofile[which(namesprofile=="ChromosomeTrue")] <- "Chromosome"
        names(profileCGH$profileValues) <- namesprofile        
###        profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ChromosomeTrue"))
        
### les champs Smoothing Level Region Breakpoints OutliersAws
### sont renommés pour qu'ils soient spécifiques
### de l'étape sur l'ensemble du génome:
### ils portent à la fin la lettre G

        profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("Smoothing","Region","Breakpoints","MinPosOrder","MaxPosOrder"))]
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Smoothing"))
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Region"))
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Breakpoints"))
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MinPosOrder"))
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MaxPosOrder"))

        namesprofile <- names(profileCGH$profileValues)
        namesprofile[which(namesprofile=="Level")] <- "LevelG"
        namesprofile[which(namesprofile=="MinPosOrderC")] <- "MinPosOrder"
        namesprofile[which(namesprofile=="MaxPosOrderC")] <- "MaxPosOrder"
        names(profileCGH$profileValues) <- namesprofile

      
        
########################################################
###
### on cherche la ligne de base à partir du
### smoothing sur l'ensemble du génome
###
##########################################################

### l'object doit etre de type profileChr
### à l'avenir, il faudra supprimer cette classe
### (attention à faire les modifs également dans la fonction glad)
        class(profileCGH) <- "profileChr"

### détection des outliers
        profileCGH <- detectOutliers(profileCGH, region="LevelG", alpha=alpha, msize=msize)

        
### le clustering est fait sur les niveaux LevelG

        print("Find cluster over the Genome")
        profileCGH$findClusterSigma <- profileCGH$SigmaG$Value[1]
        profileCGH <- findCluster(profileCGH, region="LevelG", method=method, genome=TRUE,
                                  lambda=lambdaclusterGen, nmin=nmin, nmax=nmax)
        class(profileCGH) <- "profileCGH"


### On a besoin d'un champ OutliersAws:
### on ne considère comme Outliers AWS que ceux détectés        
### à l'étape chromosome
        
        profileCGH$profileValues$OutliersAws <- profileCGH$profileValues$OutliersAwsC
                

        print("Find the Normal Baseline")
### on trouve le normal (ligne de base) à partir de la fonction affectationGNL
        profileCGH <- affectationGNL(profileCGH)


        
### On renomme le champ ZoneGNL en ZoneGNLGen
        namesprofile <- names(profileCGH$profileValues)
        namesprofile[which(namesprofile=="ZoneGNL")] <- "ZoneGNLGen"   
        names(profileCGH$profileValues) <- namesprofile
        
### notre niveau de référence pour le normal correspond
### à la médiane des log-ratios correspondant à une ZoneGNLGen de 0
        NormalRef <- median(profileCGH$profileValues[which(profileCGH$profileValues$ZoneGNLGen==0),"LogRatio"], na.rm=TRUE)

        profileCGH$NormalRef <- NormalRef

        profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"LevelG")]
###        profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"LevelG"))


      }

    else
      {


        IQRdiff <- function(y) IQR(diff(y))/1.908


        if(profilage) Rprof("/tmp/Step03DataPreparation.dat")
        
##         profileCGH$profileValues$OutliersAws <- profileCGH$profileValues$OutliersAwsC
##         profileCGH$profileValues$MinPosOrder <- profileCGH$profileValues$MinPosOrderC
##         profileCGH$profileValues$MaxPosOrder <- profileCGH$profileValues$MaxPosOrderC

        namesprofile <- names(profileCGH$profileValues)
        namesprofile[which(namesprofile=="OutliersAwsC")] <- "OutliersAws"
        namesprofile[which(namesprofile=="MinPosOrderC")] <- "MinPosOrder"
        namesprofile[which(namesprofile=="MaxPosOrderC")] <- "MaxPosOrder"
        names(profileCGH$profileValues) <- namesprofile
        
        profileCGH$profileValues$OutliersMad <- profileCGH$profileValues$OutliersAws
        profileCGH$profileValues$OutliersTot <- profileCGH$profileValues$OutliersAws        
        
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MinPosOrderC"))
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MaxPosOrderC"))
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"OutliersAwsC"))

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

        if(profilage)Rprof(NULL)
      }

    
##################################################
#### fin de l'étape sur le génome
##############################################

       

    if(profilage) Rprof("/tmp/Step04DataPreparation.dat")    
### On va appliquer la fonction removeLevel
### l'objet doit donc être de type profileChr
    class(profileCGH) <- "profileChr"


    

### il faut le champ Breakpoints et Region
### mais on peut supprimer son utilisation dans removeLevel

##     profileCGH$profileValues$Level <- profileCGH$profileValues$LevelC    
##     profileCGH$profileValues$Breakpoints <- profileCGH$profileValues$BreakpointsC
##     profileCGH$profileValues$Region <- profileCGH$profileValues$RegionC

    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile=="LevelC")] <- "Level"
    namesprofile[which(namesprofile=="BreakpointsC")] <- "Breakpoints"
    namesprofile[which(namesprofile=="RegionC")] <- "Region"
    names(profileCGH$profileValues) <- namesprofile    
    
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"RegionC"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"BreakpointsC"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"LevelC"))    
      

    indice <- 1:length(profileCGH$profileValues[,1])
    ChrIndice <- split(indice,profileCGH$profileValues$Chromosome)
    NbChr <- length(names(ChrIndice))

    Init <- 0       

    profileCGH$profileValues <- data.frame(profileCGH$profileValues, NextLogRatio=Init)
    FieldOrder <- names(profileCGH$profileValues)

    if(profilage) Rprof(NULL)
    print("Optimization of the Breakpoints")

    if(profilage) Rprof("/tmp/Step05BkpOptimization.dat")        
    for (i in 1:NbChr)
      {
        indexChr <- ChrIndice[[i]]
        subset <- profileCGH$profileValues[indexChr,]	
        profileChr <- list(profileValues=subset)	
        class(profileChr) <- "profileChr"
	profileChr$findClusterSigma <- profileCGH$SigmaC$Value[i]
        profileChr <- removeLevel(profileChr, lambda=lambdabreak, param=param, alpha=alpha, msize=msize, verbose=verbose)
        profileCGH$profileValues[indexChr,] <- profileChr$profileValues[,FieldOrder]
      }

    if(profilage) Rprof(NULL)



    if(profilage) Rprof("/tmp/Step06DataPreparation.dat")        
### Ajout le 23062006
    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
##     profileCGH$profileValues$Level <- .C("makeRegion",
##                                          as.integer(profileCGH$profileValues$Level),
##                                          as.integer(profileCGH$profileValues$Chromosome),
##                                          ResLevel = as.integer(profileCGH$profileValues$Level),
##                                          as.integer(length(profileCGH$profileValues[,1])),
##                                          PACKAGE="GLAD")$ResLevel
    
 
    agg <- aggregate(profileCGH$profileValues$LogRatio, list(Level=profileCGH$profileValues$Level), median)
    agg$Level <- as.numeric(as.character(agg$Level))
    names(agg) <- c("Level","Smoothing")

    lengthDest <- length(profileCGH$profileValues$Level)
    lengthSrc <- length(agg$Level)
    mySmooting <- .C("my_merge",
                     as.integer(profileCGH$profileValues$Level),
                     Smoothing=double(lengthDest),
                     as.integer(agg$Level),
                     as.double(agg$Smoothing),
                     as.integer(lengthDest),
                     as.integer(lengthSrc),
                     PACKAGE="GLAD")

    profileCGH$profileValues$Smoothing <- mySmooting$Smoothing
    
    ##     t1 <- system.time(profileCGH$profileValues <- merge(profileCGH$profileValues, agg, by="Level", all=TRUE))
    ##     print(t1)
    class(profileCGH) <- "profileCGH"

            


    
### pour le calcul des poids, l'écart-type utilisé
### est celui calculé sur l'ensemble de génome    
    Sigma <- profileCGH$SigmaG$Value[1]    
    profileCGH$profileValues <- data.frame(profileCGH$profileValues, NormalRef=NormalRef, Sigma=Sigma)
    profileCGH$profileValues$DiffBase <- profileCGH$profileValues$Smoothing - profileCGH$profileValues$NormalRef

    
### on prend comme référence ceux qui sont compris entre certaines valeurs
    profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGen")]
#    profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGen"))
    indexNormalLevel <- which(abs(profileCGH$profileValues$DiffBase)<=deltaN)
    profileCGH$profileValues$NormalRange <- profileCGH$profileValues$Level
    profileCGH$profileValues$NormalRange[indexNormalLevel] <- 0


    if(profilage)Rprof(NULL)
    
### le clustering est fait sur les niveaux NormalRange
    print("DNA copy number calling")
    if(profilage) Rprof("/tmp/Step07DataCalling.dat")            
    class(profileCGH) <- "profileChr"
    profileCGH <- findCluster(profileCGH, region="NormalRange", method=method, genome=TRUE,
                              lambda=lambdaclusterGen, nmin=nmin, nmax=nmax, verbose=verbose)
    class(profileCGH) <- "profileCGH"
    if(profilage)Rprof(NULL)

### le cluster correspondant au normal est celui qui comprend
### le NormalRange 0
    if(profilage) Rprof("/tmp/Step08DataPreparation.dat")                
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


    lengthDest <- length(profileCGH$profileValues$ZoneGen)
    lengthSrc <- length(MedianCluster$ZoneGen)
    myZoneGNL <- .C("my_merge_int",
                    as.integer(profileCGH$profileValues$ZoneGen),
                    ZoneGNL=integer(lengthDest),
                    as.integer(MedianCluster$ZoneGen),
                    as.integer(MedianCluster$ZoneGNL),
                    as.integer(lengthDest),
                    as.integer(lengthSrc),
                    PACKAGE="GLAD")

    profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL
    
##     t1 <- system.time(profileCGH$profileValues <- merge(profileCGH$profileValues, MedianCluster[,c("ZoneGen","ZoneGNL")], all=TRUE, by="ZoneGen"))
##     print(t1)

    
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
    
    if(profilage)Rprof(NULL)


### Calcul d'un poids pour les Breakpoints
### Attention: comme on calcul une variable GNLchange
### il ne faut pas le statut des outliers pour définir le changement de statut
### car un outlier peut être aussi un breakpoints

    if (verbose) print("daglad - step BkpInfo")
    if(profilage) Rprof("/tmp/Step09BkpInfo.dat")                    
    profileCGH$nbsigma <- nbsigma
    profileCGH$BkpInfo <- BkpInfo(profileCGH)
    if(profilage)Rprof(NULL)

    
    RecomputeGNL <- FALSE
### Epuration des Breakpoints consécutifs
    if (is.data.frame(profileCGH$BkpInfo))
      {

        if(profilage) Rprof("/tmp/Step10BkpInfo.dat")                            
        profileCGH$BkpInfo <- profileCGH$BkpInfo[order(profileCGH$BkpInfo$PosOrder),]
        profileCGH$BkpInfo$BkpToDel <- rep(0, length(profileCGH$BkpInfo[,1]))
        profileCGH$BkpInfo$Side <- profileCGH$BkpInfo$BkpToDel
        profileCGH$BkpInfo$NextPosOrder <- profileCGH$BkpInfo$PosOrder + 1

### on regarde que si il y a au moins 2 Bkp
        if (length(profileCGH$BkpInfo[,1])>1)
          {

            deleteContiguousBkp <- .C("delete_contiguous_bkp",
                                      BkpToDel=as.integer(profileCGH$BkpInfo$BkpToDel),
                                      Gap=as.double(profileCGH$BkpInfo$Gap),
                                      LogRatio=as.double(profileCGH$BkpInfo$LogRatio),
                                      NextPosOrder=as.integer(profileCGH$BkpInfo$NextPosOrder),
                                      PosOrder=as.integer(profileCGH$BkpInfo$PosOrder),
                                      Side=as.integer(profileCGH$BkpInfo$Side),
                                      Sigma=as.double(profileCGH$BkpInfo$Sigma),
                                      Smoothing=as.double(profileCGH$BkpInfo$Smoothing),
                                      SmoothingNext=as.double(profileCGH$BkpInfo$SmoothingNext),
                                      Weight=as.double(profileCGH$BkpInfo$Weight),
                                      as.integer(length(profileCGH$BkpInfo[,1])),
                                      RecomputeGNL=as.integer(0),
                                      as.integer(nbsigma),                                      
                                      PACKAGE="GLAD")


            profileCGH$BkpInfo$BkpToDel <- deleteContiguousBkp$BkpToDel
            profileCGH$BkpInfo$Gap <- deleteContiguousBkp$Gap
            profileCGH$BkpInfo$LogRatio <- deleteContiguousBkp$LogRatio
            profileCGH$BkpInfo$NextPosOrder <- deleteContiguousBkp$NextPosOrder
            profileCGH$BkpInfo$PosOrder <- deleteContiguousBkp$PosOrder
            profileCGH$BkpInfo$Side <- deleteContiguousBkp$Side
            profileCGH$BkpInfo$Sigma <- deleteContiguousBkp$Sigma
            profileCGH$BkpInfo$Smoothing <- deleteContiguousBkp$Smoothing
            profileCGH$BkpInfo$SmoothingNext <- deleteContiguousBkp$SmoothingNext
            profileCGH$BkpInfo$Weight <- deleteContiguousBkp$Weight
            
            RecomputeGNL <- deleteContiguousBkp$RecomputeGNL

            if(profilage)Rprof(NULL)
          }


### ICI on peut encore optimiser en supprimant la boucle for
        indexBPtoDel <- which(profileCGH$BkpInfo$BkpToDel==1)
        if (length(indexBPtoDel)>0)
          {

            if(profilage) Rprof("/tmp/Step11BkpInfo.dat")                                                                   
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

            if(profilage)Rprof(NULL)            
          }


        profileCGH$BkpInfo <- profileCGH$BkpInfo[,setdiff(names(profileCGH$BkpInfo),c("Side","BkpToDel","NextPosOrder"))]
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"Side"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"BkpToDel"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"NextPosOrder"))

      }



    if (verbose) print("daglad - step recomputeGNL")
    if (RecomputeGNL)
      {
        if(profilage) Rprof("/tmp/Step12RecomputeGNL.dat")                                                                           
        profileCGH$profileValues$DiffBase <- profileCGH$profileValues$Smoothing - profileCGH$profileValues$NormalRef


        profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("ZoneGNL","ZoneGen"))]
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGNL"))        
##         profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGen"))
        
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

        lengthDest <- length(profileCGH$profileValues$ZoneGen)
        lengthSrc <- length(MedianCluster$ZoneGen)
        myZoneGNL <- .C("my_merge_int",
                        as.integer(profileCGH$profileValues$ZoneGen),
                        ZoneGNL=integer(lengthDest),
                        as.integer(MedianCluster$ZoneGen),
                        as.integer(MedianCluster$ZoneGNL),
                        as.integer(lengthDest),
                        as.integer(lengthSrc),
                        PACKAGE="GLAD")

        profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL
        
##         t1 <- system.time(profileCGH$profileValues <- merge(profileCGH$profileValues, MedianCluster[,c("ZoneGen","ZoneGNL")], all=TRUE, by="ZoneGen"))
##         print(t1)

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

        if(profilage)Rprof(NULL)                    
      }



### Statut des Outliers

    if (verbose) print("daglad - step OutliersGNL")

    if(profilage) Rprof("/tmp/Step13OutliersGNL.dat")                                                                               
    profileCGH <- OutliersGNL(profileCGH, alpha=alpha, sigma=Sigma, NormalRef=NormalRef, amplicon=amplicon, deletion=deletion, assignGNLOut=assignGNLOut)


    profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("NormalRange","ZoneGen"))]
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"NormalRange"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGen"))


### Récupération des informations sur l'analyse
    profileCGH <- profileCGH[-which(names(profileCGH)=="Sigma")]
### Les données sont centrées sur NormalRef

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

    if(profilage)Rprof(NULL)                    

    
#################################################################################
###
###  Bkp filter
###
#################################################################################        

    if (verbose) print("daglad - step filterBkpStep")
    if(profilage) Rprof("/tmp/Step14FilterBkp.dat")                                                                                   
    profileCGH <- filterBkpStep(profileCGH, MinBkpWeight=MinBkpWeight, assignGNLOut=assignGNLOut, verbose=verbose)
    if(profilage)Rprof(NULL)                        

#################################################################################
###
###  Déplacement des Bkp
###
#################################################################################    

    if (CheckBkpPos)
      {
        if (verbose) print("daglad - step MoveBkpStep")
        if(profilage) Rprof("/tmp/Step15MoveBkp.dat")                                                                                           
        profileCGH <- MoveBkpStep(profileCGH, assignGNLOut=assignGNLOut)
        if(profilage)Rprof(NULL)                        
      }


    print("Results Preparation")
    
### fin comment
### suppression des champs inutiles


    if(profilage) Rprof("/tmp/Step15Results.dat")                                                                                               
    if (genomestep)
      {
        profileCGH$profileValues$NormalRef <- profileCGH$profileValues$ZoneGNLGen
        profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),"ZoneGNLGen")]
##        profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGNLGen"))
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
###        PosAux <- profileCGH$BkpInfo$PosOrder
        namesprofile <- names(profileCGH$BkpInfo)
        namesprofile[which(namesprofile=="NewPosOrder")] <- "PosOrder"
        names(profileCGH$BkpInfo) <- namesprofile        
        
###        profileCGH$BkpInfo$PosOrder <- profileCGH$BkpInfo$NewPosOrder
        profileCGH$BkpInfo <- profileCGH$BkpInfo[,setdiff(names(profileCGH$BkpInfo),c("SmoothingNext","NextLogRatio","MaxPosOrder","MinPosOrder","ZoneGNLnext"))]
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"NewPosOrder"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"SmoothingNext"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"NextLogRatio"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"MaxPosOrder"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"MinPosOrder"))
##         profileCGH$BkpInfo <- subset(profileCGH$BkpInfo, select=setdiff(names(profileCGH$BkpInfo),"ZoneGNLnext"))
      }

### il faut supprimer le champ région que ne sert plus à rien
### ou bien s'assurer qu'il est à jour

## il faut le champ ZoneGen ou un équivalent
    ### remettre les champs dans l'ordre
    
###    PosAux <- profileCGH$profileValues$PosOrder
    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile=="NewPosOrder")] <- "PosOrder"
    names(profileCGH$profileValues) <- namesprofile        
    
###    profileCGH$profileValues$PosOrder <- profileCGH$profileValues$NewPosOrder

    profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("Sigma","DiffBase","Region","NextLogRatio","MaxPosOrder","MinPosOrder"))]
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"-NewPosOrder"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Sigma"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"DiffBase"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Region"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"NextLogRatio"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MaxPosOrder"))
##     profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MinPosOrder"))



    outputfields <- setdiff(names(profileCGH$profileValues),inputfields)
    at <- setdiff(attributes(profileCGH)$names,c("PosOrderRange","findClusterSigma","NbClusterOpt"))
    profileCGH <- profileCGH[at]
    profileCGH$profileValues <- profileCGH$profileValues[,c(inputfields,outputfields)]
    class(profileCGH) <- "profileCGH"
    

    if(profilage)Rprof(NULL)
    
    return (profileCGH)
    
    
  }


