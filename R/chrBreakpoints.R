### This function detects chromosomal breakpoints along genome

### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr




chrBreakpoints <- function(...)
  {    
    UseMethod("chrBreakpoints")
  }



chrBreakpoints.profileCGH <- function(profileCGH, smoothfunc="laws", base=TRUE, sigma=NULL, bandwidth=10, round=2, verbose=FALSE, ...)
  {
 

    if (verbose)
      {
        print("chrBreakpoints: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }

    if (smoothfunc!="laws" && smoothfunc!="aws" && smoothfunc!="lawsglad")stop("Choose either aws or laws for smoothfunc")

    if (is.null(sigma))
      resetsigma <- TRUE

    else
      resetsigma <- FALSE




#### local function
    IQRdiff <- function(y) IQR(diff(y))/1.908


####data <- data.frame(LogRatio=LogRatio, Position=Position, Chromosome=Chromosome)

    data <- profileCGH$profileValues

    ### on ordonne par position et chromosome
    data <- data[order(data$Chromosome,data$PosOrder),]
    indice <- 1:length(data[,1])

    ChrIndice <- split(indice,data$Chromosome)
    ChrName <- names(ChrIndice)
    if (is.numeric(data$Chromosome[1]))
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
    Smoothing <- rep(NA, length(data[,1]))
    OutliersAws <- rep(0, length(data[,1]))
    Region <- Level <- Breakpoints <- MinPosOrder <- MaxPosOrder <- Smoothing
    data <- data.frame(data, Smoothing, OutliersAws, Region, Level, Breakpoints, MinPosOrder, MaxPosOrder)
    


    #labelChr <- unique(data$Chromosome)
 


### initialization of region number to 0
    nbregion <- 0

### initialization of level number to 0
    nblevel <- 0


    dataaux <- data.frame(data, MeanLevel=data$Region, LevelNewOrder=data$Region)
    #print(class(dataaux$BAC))
    FieldOrder <- names(dataaux)
    #dataaux <- NULL
    #print(names(dataaux))
    IQRvalue <- IQRChr <- NULL
    for (i in 1:NbChr)
      {
        if (verbose) print(paste("chrBreakpoints: starting chromosome", labelChr[i]))


### location  of data related to each chromosome

        indexChr <- ChrIndice[[i]]
        #indexChr <- which(data$Chromosome==labelChr[i])
        subsetdata <- data[indexChr,]

### information sur les bornes d'un chromosome
        PosOrderRange$Chromosome[i] <- labelChr[i]
        rangePos <- range(subsetdata$PosOrder)
        subsetdata$MinPosOrder <- PosOrderRange$MinPosOrder[i] <- rangePos[1]
        subsetdata$MaxPosOrder <- PosOrderRange$MaxPosOrder[i] <- rangePos[2]

### les données doivent être ordonnées par position
        #subsetdata <- subsetdata[order(subsetdata$PosOrder),]

        if (length(indexChr)>1)
          {
            if (resetsigma)
              {
                #print(length(subsetdata$LogRatio))
                #print(paste("in chrBkp:", IQRdiff(subsetdata$LogRatio)))
                sigma <- IQRdiff(subsetdata$LogRatio)^2
                IQRvalue[i] <- sigma^(0.5)
                IQRChr[i] <- labelChr[i]
                #write.table(subsetdata[,c("LogRatio","PosOrder")],file="~/tmp/Bkp.txt", row.names=FALSE, sep="\t")

              }
            else
              {
                IQRvalue[i] <- sigma
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
                    awsres <- laws(y=subsetdata$LogRatio, x=x, hinit=hinit, hmax=hmax, shape=sigma, NN=FALSE, symmetric=TRUE, ...)$theta

                    if (is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- round(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }		
                  }

                if (smoothfunc=="aws")
                  {
                    awsres <- aws(y=subsetdata$LogRatio, x=x, hinit=hinit, hmax=hmax, sigma2=sigma, NN=FALSE, symmetric=TRUE, ...)$theta
                    print(timeBkp)

                    if (is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- round(awsres, round)
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
                    awsres <- lawsglad(y=subsetdata$LogRatio, hinit=hinit, hmax=hmax, shape=sigma, ...)

                    if(is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- round(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }                        
                  }                
                
                if (smoothfunc=="laws")
                  {
                    awsres <- laws(y=subsetdata$LogRatio, hinit=hinit, hmax=hmax, shape=sigma, symmetric=TRUE, ...)$theta

                    if(is.null(awsres)==FALSE)
                      {
                        subsetdata$Smoothing <- round(awsres, round)
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
                        subsetdata$Smoothing <- round(awsres, round)
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

            j <- 1
            

            labelLevel <- unique(round(subsetdata$Smoothing,round))
            
            for (j in 1:length(labelLevel))
              {
                nblevel <- nblevel + 1
                subsetdata$Level[which(round(subsetdata$Smoothing, round)==labelLevel[j])] <- nblevel	
              }
            
            MeanLevel <- aggregate(subsetdata$LogRatio, list(Level=subsetdata$Level),mean)
            names(MeanLevel) <- c("Level","MeanLevel")
            MeanLevel <- MeanLevel[order(MeanLevel$MeanLevel),]
            MeanLevel$Level <- as.numeric(as.character(MeanLevel$Level))
            MeanLevel$LevelNewOrder <- min(MeanLevel$Level):max(MeanLevel$Level)
            subsetdata <- merge(subsetdata, MeanLevel, by="Level",all=TRUE)




###breakpoints detection
            rupture <- 0


### il faut remettre les données dans l'ordre par position
            subsetdata <- subsetdata[order(subsetdata$PosOrder),]

                                        #           ### debut
                                        #           for (j in 2:length(subsetdata$LogRatio))
                                        #             {
            
                                        # ####Outliers detection
            
                                        # ### not Outliers
                                        # ### Pour les outliers, leur Level correspond à celui de la région
                                        # ### de laquelle ils sont issus. Cela permet de faciliter l'affectation
                                        # ### de leur statut dans la fonction affectationGNL              


            intl <- length(subsetdata$LogRatio)
            awsBkp <- .C("awsBkp",
                         as.double(round(subsetdata$Smoothing,round)),
                         OutliersAws=as.integer(subsetdata$OutliersAws),
                         Level=as.integer(subsetdata$Level),
                         nbregion=as.integer(nbregion),
                         regionChr=as.integer(c(nbregion,rep(0,intl-1))),
                         rupture=as.integer(rep(0,intl)),
                         as.integer(intl),
                         PACKAGE="GLAD")

            
            subsetdata$Breakpoints <- c(awsBkp$rupture[2:intl],0)
            subsetdata$Region <- awsBkp$regionChr
            subsetdata$Level <- awsBkp$regionChr
            subsetdata$OutliersAws <- awsBkp$OutliersAws
            nbregion <- awsBkp$nbregion
            



            

          }

        else
          {
            subsetdata$Region <- -1
            subsetdata$Level <- -1
            subsetdata$Breakpoints <- 0
            subsetdata$LevelNewOrder <- subsetdata$Breakpoints
            subsetdata$MeanLevel <- subsetdata$Breakpoints
          }

#        dataaux <- rbind(dataaux, subsetdata)
#        ns <- print(names(subsetdata))
#         print(ns==nd)
#         print("ns")
#         print(ns)
#         print("nd")
#         print(nd)
#         print("sd1")
#         print(setdiff(ns,nd))
#         print("sd2")
#         print(setdiff(nd,ns))
#         print(class(subsetdata$BAC))
        dataaux[indexChr,] <- subsetdata[,FieldOrder]
        
        
        if (verbose) print(paste("chrBreakpoints: ending chromosome", labelChr[i]))

      }

    #print(dataaux[1:10,])

    dataaux$Level <- dataaux$LevelNewOrder
    data <- subset(dataaux, select=-MeanLevel)
    data <- subset(data, select=-LevelNewOrder)
    nomdata <- names(data)
    nomdata <- nomdata[order(nomdata)]
    data <- data[,nomdata]


    
### permutation des Champs Level et LevelTrue pour faciliter l'utilisation dans les test sur removeLevel

    profileCGH$profileValues <- data

    profileCGH$Sigma <- data.frame(Chromosome=IQRChr,Value=IQRvalue)

    profileCGH$PosOrderRange <- PosOrderRange
    
    if (verbose) print("chrBreakpoints: ending function")

    return(profileCGH)
  }



