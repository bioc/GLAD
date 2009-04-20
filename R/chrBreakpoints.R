### This function detects chromosomal breakpoints along genome

### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr




chrBreakpoints <- function(...)
  {    
    UseMethod("chrBreakpoints")
  }



chrBreakpoints.profileCGH <- function(profileCGH, smoothfunc="lawsglad", base=FALSE, sigma=NULL,
                                      model="Gaussian", bandwidth=10, round=1.5, verbose=FALSE,
                                      breaksFdrQ = 0.0001, haarStartLevel = 1, haarEndLevel = 5, ...)
  {
 

    if (verbose)
      {
        print("chrBreakpoints: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }

    if (smoothfunc!="laws" && smoothfunc!="aws" && smoothfunc!="lawsglad" && smoothfunc!="haarseg")stop("Choose either aws, laws or haarseg for smoothfunc")

    if (base==TRUE)
      {
        if(smoothfunc!="lawsglad" && smoothfunc!="haarseg")stop("Choose either aws, or laws when base=TRUE")
      }
    if (is.null(sigma))
      resetsigma <- TRUE

    else
      resetsigma <- FALSE


### le round est fait a de nombreux endroits différents: cela peut
### surement être changé



#### local functions
    IQRdiff <- function(y) IQR(diff(y))/1.908

    roundglad <- function(x, digits=2)
      {
        dec <- (digits-trunc(digits))
        if (dec==0)
          dec <- 1
        
        r <- 10^(-trunc(digits))*dec
        x <- r*round(x/r)
        return(x)
      }


    ### on ordonne par position et chromosome
    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$Chromosome,profileCGH$profileValues$PosOrder),]
    indice <- 1:length(profileCGH$profileValues[,1])

    ChrIndice <- split(indice,profileCGH$profileValues$Chromosome)
    ChrName <- names(ChrIndice)
    if (is.numeric(profileCGH$profileValues$Chromosome[1]))
      {
        labelChr <- as.numeric(ChrName)
      }
    else
      {
        labelChr <- as.character(ChrName)
      }


    NbChr <- length(labelChr)
    lg <- rep(0,NbChr)
    PosOrderRange <- data.frame(Chromosome=lg, MinPosOrder=lg, MaxPosOrder=lg)
  

    
### ajout des informations
    NbRow <- length(profileCGH$profileValues[,1])
    Smoothing <- rep(NA, NbRow)
    OutliersAws <- rep(0, NbRow)
    Region <- Level <- Breakpoints <- MinPosOrder <- MaxPosOrder <- Smoothing
    profileCGH$profileValues <- data.frame(profileCGH$profileValues, Smoothing,
                                           OutliersAws, Region, Level, Breakpoints,
                                           MinPosOrder, MaxPosOrder)
     


### initialization of region number to 0
    nbregion <- 0

### initialization of level number to 0
    nblevel <- 0


    ### Il ne faut pas les champs MedianLevel et LevelNewOrder
    ### pour la jointure entre subsetdata et MedianLevel
    ### par contre on en a besoin au momment du remplissage des données
    FieldInit <- names(profileCGH$profileValues)
    profileCGH$profileValues <- data.frame(profileCGH$profileValues,
                                           MedianLevel=profileCGH$profileValues$Region,
                                           LevelNewOrder=profileCGH$profileValues$Region)

    FieldOrder <- names(profileCGH$profileValues)

    IQRvalue <- IQRChr <- NULL
    for (i in 1:NbChr)
      {
        if (verbose) print(paste("chrBreakpoints: starting chromosome", labelChr[i]))


### location  of data related to each chromosome

        indexChr <- ChrIndice[[i]]
        subsetdata <- profileCGH$profileValues[indexChr,FieldInit]

### information sur les bornes d'un chromosome
        PosOrderRange$Chromosome[i] <- labelChr[i]
        rangePos <- range(subsetdata$PosOrder)
        subsetdata$MinPosOrder <- PosOrderRange$MinPosOrder[i] <- rangePos[1]
        subsetdata$MaxPosOrder <- PosOrderRange$MaxPosOrder[i] <- rangePos[2]

### les données doivent être ordonnées par position

        if (length(indexChr)>1)
          {
            if (resetsigma)
              {

                sigma <- IQRdiff(subsetdata$LogRatio)^2
                IQRvalue[i] <- sigma^(0.5)
                IQRChr[i] <- labelChr[i]                

              }
            else
              {
                IQRvalue[i] <- sigma^0.5
                IQRChr[i] <- labelChr[i]
              }
            


            if (base==TRUE)
              {
                x <- subsetdata$PosBase
                datarange <- range(x)
                hmax <- diff(datarange)*bandwidth
                hinit <- median(diff(x))
                
                
### smoothing
                if (smoothfunc=="laws")
                  {
                    dim(x) <- c(1,length(x)) #à supprimer dans la nouvelle version du package AWS
                    awsres <- laws(y=subsetdata$LogRatio, x=x,
                                   hinit=hinit, hmax=hmax, shape=sigma, NN=FALSE,
                                   symmetric=TRUE, model=model, ...)$theta

                    if (is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }		
                  }

                if (smoothfunc=="aws")
                  {
                    awsres <- aws(y=subsetdata$LogRatio, x=x, hinit=hinit, hmax=hmax, sigma2=sigma, NN=FALSE, symmetric=TRUE, ...)$theta

                    if (is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }
                  }
                
              }
            
            else          
              {
                hinit <- 1
                hmax <- length(subsetdata$PosOrder)*bandwidth #si on laisse hmax à la valeur de length(data$PosOrder[indexChr]) ce n'est pas assez et les créneaux ne sont pas bien fittés

                
                if (smoothfunc=="lawsglad")
                  {
                    awsres <- lawsglad(y=subsetdata$LogRatio,
                                       hinit=hinit, hmax=hmax, shape=sigma, model=model, ...)

                    if(is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }                        
                  }                
                
                if (smoothfunc=="laws")
                  {
                    awsres <- laws(y=subsetdata$LogRatio, hinit=hinit,
                                   hmax=hmax, shape=sigma, symmetric=TRUE, model=model,
                                   ...)$theta

                    if(is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }                        
                  }

                if (smoothfunc=="aws")
                  {
                    awsres <- aws(y=subsetdata$LogRatio, hinit=hinit, hmax=hmax, sigma2=sigma, symmetric=TRUE, ...)$theta
                    
                    if(is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }                    
                    
                 }

                if (smoothfunc=="haarseg")
                  {
                    awsres <- HaarSeg(subsetdata$LogRatio, breaksFdrQ=breaksFdrQ, haarStartLevel=haarStartLevel , haarEndLevel=haarEndLevel)$Segmented

                    if (is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }
                  }
                
              }


            nbregion <- nbregion + 1
            
                                        #vector with label zone (initialized for each chromosome)
            regionChr <- NULL
            rupture <- NULL


                                        #set label zone for each position along chromosome
            regionChr <- c(regionChr, nbregion)  #first BAC corresponds to nbregion
            
            
                                        #level assignation

            #j <- 1
            

            ### 29062005: cette partie peut-être optimisée:

            
            #labelLevel <- unique(roundglad(subsetdata$Smoothing,round))


            #########################################
            ########################################
            ########################################

            
##             labelLevel <- unique(subsetdata$Smoothing)
  
##             for (j in 1:length(labelLevel))
##               {
##                 nblevel <- nblevel + 1
##                 ###subsetdata$Level[which(roundglad(subsetdata$Smoothing, round)==labelLevel[j])] <- nblevel	
##                 subsetdata$Level[which(subsetdata$Smoothing==labelLevel[j])] <- nblevel	
##               }

            l <- length(subsetdata$Smoothing)
            putLevel <- .C("putLevel",
                           as.double(subsetdata$Smoothing),
                           Level=integer(l),
                           nblevel=as.integer(nblevel),
                           as.integer(l),
                           PACKAGE="GLAD")

            subsetdata$Level <- putLevel$Level
            nblevel <- putLevel$nblevel
            
            ###########################################"
            ###########################################
            ############################################
            
            ### 02062005: les levels sont ordonnées par ordre croissant
            ### ceci est important dans le cas de l'utilisation de la fonction removeLevel
            MedianLevel <- aggregate(subsetdata$LogRatio, list(Level=subsetdata$Level),median)
            names(MedianLevel) <- c("Level","MedianLevel")
            MedianLevel <- MedianLevel[order(MedianLevel$MedianLevel),]
            MedianLevel$Level <- as.numeric(as.character(MedianLevel$Level))
            MedianLevel$LevelNewOrder <- min(MedianLevel$Level):max(MedianLevel$Level)
            subsetdata <- merge(subsetdata, MedianLevel, by="Level",all=TRUE)




###breakpoints detection
            rupture <- 0


### il faut remettre les données dans l'ordre par position
            subsetdata <- subsetdata[order(subsetdata$PosOrder),]


            intl <- length(subsetdata$LogRatio)
            awsBkp <- .C("awsBkp",
                         #as.double(roundglad(subsetdata$Smoothing,round)),
                         as.double(subsetdata$Smoothing,round),
                         OutliersAws=as.integer(subsetdata$OutliersAws),
                         Level=as.integer(subsetdata$LevelNewOrder),
                         nbregion=as.integer(nbregion),
                         regionChr=as.integer(c(nbregion,rep(0,intl-1))),
                         rupture=as.integer(rep(0,intl)),
                         as.integer(intl),
                         PACKAGE="GLAD")

            
            subsetdata$Breakpoints <- c(awsBkp$rupture[2:intl],0)
            subsetdata$Region <- awsBkp$regionChr
            subsetdata$Level <- awsBkp$Level
            subsetdata$OutliersAws <- awsBkp$OutliersAws
            nbregion <- awsBkp$nbregion                        

          }

        else
          {
            subsetdata$Region <- -1
            subsetdata$Level <- -1
            subsetdata$Breakpoints <- 0
            subsetdata$LevelNewOrder <- subsetdata$Breakpoints
            subsetdata$MedianLevel <- subsetdata$Breakpoints
            IQRvalue[i] <- 0
            IQRChr[i] <- labelChr[i]    
          }


        profileCGH$profileValues[indexChr,] <- subsetdata[,FieldOrder]
        
        
        if (verbose) print(paste("chrBreakpoints: ending chromosome", labelChr[i]))

      }


    profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"MedianLevel"))
    profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"LevelNewOrder"))
    nomdata <- names(profileCGH$profileValues)
    nomdata <- nomdata[order(nomdata)]
    profileCGH$profileValues <- profileCGH$profileValues[,nomdata]


    
### permutation des Champs Level et LevelTrue pour faciliter l'utilisation dans les test sur removeLevel


    profileCGH$Sigma <- data.frame(Chromosome=IQRChr,Value=IQRvalue)

    profileCGH$PosOrderRange <- PosOrderRange
    
    if (verbose) print("chrBreakpoints: ending function")

    return(profileCGH)
  }



