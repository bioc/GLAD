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

    if (smoothfunc != "laws" &&
        smoothfunc != "aws" &&
        smoothfunc != "lawsglad" &&
        smoothfunc != "haarseg")
      {
        stop("Choose either aws, laws or haarseg for smoothfunc")
      }

    if (base == TRUE)
      {
        if(smoothfunc != "lawsglad" && smoothfunc != "haarseg")
          {
            stop("Choose either aws, or laws when base=TRUE")
          }
      }
    
    if (is.null(sigma))
      {
        resetsigma <- TRUE
      }

    else
      {
        resetsigma <- FALSE
      }


### le round est fait a de nombreux endroits différents: cela peut
### surement être changé



#### local functions
    IQRdiff <- function(y) IQR(diff(y)) / 1.908

    roundglad <- function(x, digits=2)
      {
        ## en C        return floor(n * pow(10., d) + .5) / pow(10., d);
        dec <- (digits-trunc(digits))
        if (dec==0)
          dec <- 1
        
        r <- 10^(-trunc(digits))*dec
        x <- r*round(x/r)
        return(x)
      }

    ##     roundtest <- function(x, d=2)
    ##       {
    ##         return(floor(x * 10^d + .5) / 10^d)
    ##       }


    indice <- 1:length(profileCGH$profileValues[,1])

    ChrIndice <- split(indice,profileCGH$profileValues$Chromosome)
    ChrName <- names(ChrIndice)
###    print(ChrName)
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
    ##     NbRow <- length(profileCGH$profileValues[,1])
    ##     Smoothing <- rep(NA, NbRow)
    ##     OutliersAws <- rep(0, NbRow)
    ##     Region <- Level <- Breakpoints <- MinPosOrder <- MaxPosOrder <- Smoothing

    profileCGH$profileValues <- data.frame(profileCGH$profileValues, Smoothing = NA,
                                           OutliersAws = 0, Region = 0, Level = 0, Breakpoints = 0,
                                           MinPosOrder = 0, MaxPosOrder = 0)
    


### initialization of region number to 0
    nbregion <- 0

### initialization of level number to 0
    nblevel <- 0


### Il ne faut pas les champs MedianLevel et LevelNewOrder
### pour la jointure entre subsetdata et MedianLevel
### par contre on en a besoin au momment du remplissage des données
    ##    FieldInit <- names(profileCGH$profileValues)
    ##     profileCGH$profileValues <- data.frame(profileCGH$profileValues,
    ##                                            MedianLevel=profileCGH$profileValues$Region,
    ##                                            LevelNewOrder=profileCGH$profileValues$Region)

    ##    FieldOrder <- names(profileCGH$profileValues)

    
    profileCGH$BkpDetected <- data.frame(Chromosome = as.integer(labelChr), BkpDetected = 0)

    
    IQRvalue <- IQRChr <- NULL
    for (i in 1:NbChr)
      {
        if (verbose) print(paste("chrBreakpoints: starting chromosome", labelChr[i]))


        ## location  of data related to each chromosome
        indexChr <- ChrIndice[[i]]
        subsetdata <- profileCGH$profileValues[indexChr,]

        ## information sur les bornes d'un chromosome
        PosOrderRange$Chromosome[i] <- labelChr[i]
        rangePos <- range(subsetdata$PosOrder)
        subsetdata$MinPosOrder <- PosOrderRange$MinPosOrder[i] <- rangePos[1]
        subsetdata$MaxPosOrder <- PosOrderRange$MaxPosOrder[i] <- rangePos[2]


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
            

            if (base == TRUE)
              {
                x <- subsetdata$PosBase
                datarange <- range(x)
                hmax <- diff(datarange)*bandwidth
                hinit <- median(diff(x))
                
                
                ## smoothing
                if (smoothfunc == "laws")
                  {
                    dim(x) <- c(1,length(x)) #à supprimer dans la nouvelle version du package AWS
                    awsres <- laws(y=subsetdata$LogRatio, x=x,
                                   hinit=hinit, hmax=hmax, shape=sigma, NN=FALSE,
                                   symmetric=TRUE, model=model, ...)$theta

                    if (is.null(awsres) == FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }		
                  }

                if (smoothfunc == "aws")
                  {
                    awsres <- aws(y=subsetdata$LogRatio, x=x, hinit=hinit, hmax=hmax, sigma2=sigma, NN=FALSE, symmetric=TRUE, ...)$theta

                    if (is.null(awsres) == FALSE)
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

                
                if (smoothfunc == "lawsglad")
                  {
                    awsres <- lawsglad(y=subsetdata$LogRatio,
                                       hinit=hinit, hmax=hmax, shape=sigma, model=model, ...)

                    if(is.null(awsres) == FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }                        
                  }                
                
                if (smoothfunc == "laws")
                  {
                    awsres <- laws(y=subsetdata$LogRatio, hinit=hinit,
                                   hmax=hmax, shape=sigma, symmetric=TRUE, model=model,
                                   ...)$theta

                    if(is.null(awsres) == FALSE)
                      {
                        subsetdata$Smoothing <- roundglad(awsres, round)
                      }

                    else
                      {
                        subsetdata$Smoothing <- 99999
                      }                        
                  }

                if (smoothfunc == "aws")
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

                if (smoothfunc == "haarseg")
                  {

                    awsres <- HaarSegGLADCPP(subsetdata$LogRatio, breaksFdrQ=breaksFdrQ, haarStartLevel=haarStartLevel , haarEndLevel=haarEndLevel)                    

                    if (is.null(awsres) == FALSE)
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
            

            l <- length(subsetdata$Smoothing)            
            putLevel <- .C("putLevel",
                           Smoothing = as.double(subsetdata$Smoothing), ## valeur de sortie
                           as.double(subsetdata$LogRatio),                              
                           Level = integer(l),                          ## valeur de sortie
                           nblevel = as.integer(nblevel),               ## valeur de sortie
                           as.integer(l),
                           PACKAGE="GLAD")

            subsetdata[,c("Smoothing", "Level")] <- putLevel[c("Smoothing", "Level")]
            nblevel <- putLevel$nblevel

            
            awsBkp <- .C("awsBkp",
                         as.double(subsetdata$Smoothing),
                         OutliersAws = as.integer(subsetdata$OutliersAws), ## valeur de sortie
                         Level = as.integer(subsetdata$Level),             ## valeur de sortie
                         nbregion = as.integer(nbregion),                  ## valeur de sortie
                         regionChr = integer(l),                           ## valeur de sortie
                         Breakpoints = integer(l),                         ## valeur de sortie
                         BkpDetected = integer(1),                         ## valeur de sortie
                         as.integer(l),
                         PACKAGE="GLAD")


            profileCGH$BkpDetected$BkpDetected[i] <- awsBkp$BkpDetected
            
            
            subsetdata[,c("Region","Level","OutliersAws", "Breakpoints")] <- awsBkp[c("regionChr","Level","OutliersAws", "Breakpoints")]

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


        profileCGH$profileValues[indexChr,] <- subsetdata

        
        if (verbose) print(paste("chrBreakpoints: ending chromosome", labelChr[i]))

      } ## fin de la boucle par chromosome


    ##    profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("MedianLevel","LevelNewOrder"))]    
    ##     nomdata <- names(profileCGH$profileValues)
    ##     nomdata <- nomdata[order(nomdata)]
    ##     profileCGH$profileValues <- profileCGH$profileValues[,nomdata]

    

    profileCGH$Sigma <- data.frame(Chromosome=IQRChr,Value=IQRvalue)

    profileCGH$PosOrderRange <- PosOrderRange
    
    if (verbose) print("chrBreakpoints: ending function")

    return(profileCGH)
  }



