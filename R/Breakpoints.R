#This function detects chromosomal breakpoints along genome

# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr


chrBreakpoints <- function(...)
{
  
  UseMethod("chrBreakpoints")
  
}




chrBreakpoints.profileCGH <- function(profileCGH, smoothfunc="aws", base=TRUE, sigma, bandwidth=10, round=2, verbose=FALSE, ...)
{

  # ce n'est plus la peine, c'est vérifier au chargement du package
  #if (!require("aws")) stop("Could not load required package aws")
  if (verbose)
    {
      print("chrBreakpoints: starting function")
      call <- match.call()
      print(paste("Call function:", call))
    }

  if (smoothfunc!="laws" && smoothfunc!="aws")stop("Choose either aws or laws for smoothfunc")

  if (missing(sigma))
    resetsigma <- TRUE

  else
    resetsigma <- FALSE




#local function
  IQRdiff <- function(y) IQR(diff(y))/1.908


#data <- data.frame(LogRatio=LogRatio, Position=Position, Chromosome=Chromosome)

  data <- profileCGH$profileValues
  indexna <- attr(na.omit(data[,c("Chromosome","LogRatio","PosOrder")]),"na.action")
  if(!is.null(indexna))
    {    
      datana <- data[indexna,]      
      data <- data[-indexna,]
    }
  else
    {
      datana <- NULL     
    }
  
  data <- data[order(data$Chromosome,data$PosOrder),]
  Smoothing <- rep(NA, length(data$PosOrder))
  OutliersAws <- rep(0, length(data$PosOrder))
  Region <- rep(NA, length(data$Chromosome)) #label zone
  Level <- Region #label level
  Breakpoints <- Level
  labelChr <- unique(data$Chromosome)



#initialization of region number to 0
  nbregion <- 0

#initialization of level number to 0
  nblevel <- 0


  for (i in 1:length(labelChr))
    {
      if (verbose) print(paste("chrBreakpoints: starting chromosome", labelChr[i]))

#location  of data related to each chromosome
      indexChr <- which(data$Chromosome==labelChr[i])


      if (length(indexChr)>1)
        {                 
          if (resetsigma)
            sigma <- IQRdiff(data$LogRatio[indexChr])^2


          if (base==TRUE)
            {
              x <- data$PosBase[indexChr]
              datarange <- range(x)
              hmax <- diff(datarange)*bandwidth
              hinit <- mean(diff(x))
              
              
#smoothing
              if (smoothfunc=="laws")
                {
                  dim(x) <- c(1,length(x)) #à supprimer dans la nouvelle version du package AWS
                  awsres <- laws(y=data$LogRatio[indexChr], x=x, hinit=hinit, hmax=hmax, shape=sigma, NN=FALSE, ...)$theta

                  if (is.null(awsres)==FALSE)
                    {
                      Smoothing[indexChr] <- round(awsres, round)
                    }

                  else
                    {
                      Smoothing[indexChr] <- 99999
                    }		
                }

              if (smoothfunc=="aws")
                {
                  awsres <- aws(y=data$LogRatio[indexChr], x=x, hinit=hinit, hmax=hmax, sigma2=sigma, NN=FALSE, ...)$theta

                  if (is.null(awsres)==FALSE)
                    {
                      Smoothing[indexChr] <- round(awsres, round)
                    }

                  else
                    {
                      Smoothing[indexChr] <- 99999
                    }
                }
            }
          
          else          
            {
              hinit <- 1
              hmax <- length(data$PosOrder[indexChr])*bandwidth #si on laisse hmax à la valeur de length(data$PosOrder[indexChr]) ce n'est pas assez et les créneaux ne sont pas bien fittés
              
              if (smoothfunc=="laws")
                {
                  awsres <- laws(y=data$LogRatio[indexChr], hinit=hinit, hmax=hmax, shape=sigma, ...)$theta

                  if(is.null(awsres)==FALSE)
                    {
                      Smoothing[indexChr] <- round(awsres, round)
                    }

                  else
                    {
                      Smoothing[indexChr] <- 99999
                    }                        
                }

              if (smoothfunc=="aws")
                {
                  awsres <- aws(y=data$LogRatio[indexChr], hinit=hinit, hmax=hmax, sigma2=sigma, ...)$theta
                  
                  if(is.null(awsres)==FALSE)
                    {
                      Smoothing[indexChr] <- round(awsres, round)
                    }

                  else
                    {
                      Smoothing[indexChr] <- 99999
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
          

          labelLevel <- unique(round(Smoothing[indexChr],round))
          LevelChr <- rep(NULL, length(indexChr))
          
          for (j in 1:length(labelLevel))
            {
              nblevel <- nblevel + 1
              LevelChr[which(round(Smoothing[indexChr], round)==labelLevel[j])] <- nblevel	
            }
          
          Level[indexChr] <- LevelChr

#breakpoints detection
          rupture <- 0


          for (j in 2:length(data$LogRatio[indexChr]))
            {
#Outliers detection

              if (
                  (round(Smoothing[indexChr][j],round) != round(Smoothing[indexChr][j-1],round))
                  &
                  (round(Smoothing[indexChr][j+1],round) != round(Smoothing[indexChr][j],round))
                  &
                  (round(Smoothing[indexChr][j+1],round) == round(Smoothing[indexChr][j-1],round))
                  &
                  (j<length(data$LogRatio[indexChr]))
                  )
                {
                  if (OutliersAws[indexChr][j-1]==0)
                    {
                      if (Smoothing[indexChr][j]>Smoothing[indexChr][j+1])
                        OutliersAws[indexChr][j] <- 1

                      else
			OutliersAws[indexChr][j] <- -1
                    }
                  
                  regionChr <- c(regionChr, nbregion) #set label zone for position j
                  rupture <- c(rupture,0)
                }
              
#not Outliers

              else
                {
                  if (
                      (round(Smoothing[indexChr][j],round) != round(Smoothing[indexChr][j-1],round))
                      &
                      (OutliersAws[indexChr][j-1] ==0)
                      )
                    {

                      if (j==2 || j==length(data$LogRatio[indexChr]))
                        {
                          regionChr <- c(regionChr, nbregion)
                          rupture <- c(rupture,0)

                          if (j==2)
                            {
                              if (Smoothing[indexChr][1]>Smoothing[indexChr][2])
                                OutliersAws[indexChr][1] <- 1

                              else
                                OutliersAws[indexChr][1] <- -1
                            }
                          
                          else
                            {
                              if (Smoothing[indexChr][j]>Smoothing[indexChr][j-1])
                                OutliersAws[indexChr][j] <- 1

                              else
                                OutliersAws[indexChr][j] <- -1
                            }
                        }
                      
                      else
                        {
                          nbregion <- nbregion + 1
                          regionChr <- c(regionChr, nbregion) #set label zone for position j
                          rupture <- c(rupture, 1)
                        }
                    }

                  else
                    {
                      regionChr <- c(regionChr, nbregion) #set label zone for position j
                      rupture <- c(rupture, 0)
                    }
                }
            }

          Region[indexChr] <- regionChr
          rupture <- rep(0,length(indexChr)) + c(rupture[2:length(indexChr)],0)
          Breakpoints[indexChr] <- rupture

        }

      else
        {
          Region[indexChr] <- -1
          Level[indexChr] <- -1
          Breakpoints[indexChr] <- 0
        }
      
      if (verbose) print(paste("chrBreakpoints: ending chromosome", labelChr[i]))

    }


  data <- data.frame(data, Smoothing=Smoothing, Region=Region, Level=Level, OutliersAws=OutliersAws, Breakpoints=Breakpoints)

  if (!is.null(indexna))
    {
      Smoothing <- rep(NA,length(datana$LogRatio))
      Region <- Smoothing
      Level <- Smoothing
      OutliersAws <- Smoothing
      Breakpoints <- Smoothing
      datana <- data.frame(datana, Smoothing=Smoothing, Region=Region, Level=Level, OutliersAws=OutliersAws, Breakpoints=Breakpoints)
      data <- rbind(data, datana)
    }
#data$Smoothing[which(data$Smoothing==99999)] <- NA
  profileCGH$profileValues <- data
  
  if (verbose) print("chrBreakpoints: ending function")

  return(profileCGH)

}


removeBreakpoints <- function(...)
  {
    UseMethod("removeBreakpoints")
  }



removeBreakpoints.profileChr <- function(profileChr, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, ...)
{

  if (verbose)
    {
      print("removeBreakpoints: starting function")
      call <- match.call()
      print(paste("Call function:", call))
    }
  
  CGH <- profileChr$profileValues
  CGH$LikeliRatio <- rep(0,length(CGH[,1]))

###########################################################################
#
# Suppression des valeurs manquantes
#
##############################################################################

# Les autres champs peuvent être NA, cela ne pose pas de problème
  indexNA <- attr(na.omit(CGH[,c("Chromosome","LogRatio","PosOrder")]),"na.action")


  
#############################################################################
#
#  Il faut tester si, pour certains chromosomes, toutes les valeurs sont NA
#
#############################################################################

  if (length(indexNA)<length(CGH[,1]))
    {
      CGHna <- CGH[indexNA,]
      if (!is.null(indexNA))
        {
          CGH <- CGH[-indexNA,]
        }
      
      subset <- list(profileValues=CGH)
      class(subset) <- "profileChr"
      

#########################################################################################
#
#  la fonction IQR permet d'estimer l'écart-type qui sera utilisé dans la fonction kernelpen
#
##########################################################################################

      
      IQRdiff <- function(y) IQR(diff(y))/1.908

# estimation de l'écart-type :
# les données doivent etre ordonnées à cause du diff
      subset$profileValues[order(subset$profileValues$PosOrder),]
      sigma <- IQRdiff(subset$profileValues$LogRatio)
      

      
      stop <- 0
      
      while (stop!=1)
	{
	  subset <- detectOutliers(subset, region="Region", verbose=verbose, ...)
	  indexOutliersTot <- which(subset$profileValues$OutliersTot==0)
          profileValues <- subset$profileValues[indexOutliersTot,]
	  aggregation <- aggregate(profileValues["LogRatio"],list(Region=profileValues$Region),mean)
	  aggregation <- aggregation[order(aggregation$Region),]
	  deltaoversigma <- abs(diff(aggregation$LogRatio)/sigma)

          
          if (verbose)
            {
              print("")
            }

          
	  labelRegion <- sort(unique(profileValues$Region))

	  if (length(labelRegion)>1)
            {
              logsigma <- 0
              
	      for (i in 1:length(labelRegion))
                {
		  indexRegion <- which(profileValues$Region==labelRegion[i])

                  
                  
		  if (length(indexRegion)>1)
                    {
                      logsigmaregion <- length(indexRegion)*log(var(profileValues$LogRatio[indexRegion])*((length(indexRegion)-1)/length(indexRegion)), exp(1))
		      logsigma <- logsigma + logsigmaregion

                      if (verbose)
                        {
                          print(paste("removeBreakpoints: logsigma of region", labelRegion[i], "=",logsigmaregion, "(Nb Obs = ",length(indexRegion), ")"))
                        }
                    }

                  else
                    {
                      
                      if (verbose)
                        {
                          print(paste("removeBreakpoints: logsigma of region", labelRegion[i], "= 0", "( Nb Obs = ",length(indexRegion),")"))
                        }
                    }
                }

             
              
              
              
#likelihoodGLOBAL <- logsigma + lambda*(length(labelRegion)-1)*log(length(profileValues$LogRatio),exp(1))
              likelihoodGLOBAL <- logsigma + lambda*sum(kernelpen(deltaoversigma,type=type,param=param))*log(length(profileValues$LogRatio),exp(1))
              if (verbose)
                {
                  print(paste("removeBreakpoints: likelihood for the profile = ", logsigma))
                  print(paste("removeBreakpoints: penalised likelihood for the profile = ", likelihoodGLOBAL))
                }

              likelihood <- rep(NULL, length(labelRegion)-1)
              
              for (i in 1:(length(labelRegion)-1))
                {
                  if (verbose)
                    {
                      print(paste("removeBreakpoints: suppression of region", labelRegion[i+1]))
                    }
                  RegionAux <- profileValues$Region
                  indexRegion <- which(profileValues$Region==labelRegion[i+1])
                  RegionAux[indexRegion] <- labelRegion[i]
                  indexRegion <- which(RegionAux==labelRegion[i])
                  
                  aggregation <- aggregate(profileValues["LogRatio"],list(Region=RegionAux),mean)
                  aggregation <- aggregation[order(aggregation$Region),]
                  deltaoversigma <- abs(diff(aggregation$LogRatio)/sigma)
                  
                  logsigma <-0
                  labelRegionAux <- unique(RegionAux)
                   
                  
                  for (j in 1:length(labelRegionAux))
                    {
                      indexRegionAux <- which(RegionAux==labelRegionAux[j])
                      
                      if (length(indexRegionAux)>1)
                        {

                          logsigmaregion <- length(indexRegionAux)*log(var(profileValues$LogRatio[indexRegionAux])*((length(indexRegionAux)-1)/length(indexRegionAux)), exp(1))
                          logsigma <- logsigma + logsigmaregion

                          if (verbose)
                            {
                              print(paste("removeBreakpoints: logsigma of region", labelRegionAux[j], "=",logsigmaregion, "( Nb Obs =",length(indexRegionAux), ")"))
                            }
                          
                        }

                      else
                        {
                          
                          if (verbose)
                            {
                              print(paste("removeBreakpoints: logsigma of region", labelRegionAux[j], "= 0", "( Nb Obs = ",length(indexRegionAux), ")"))
                            }
                          
                        }
                    }
                  
                  likelihood[i] <- logsigma + lambda*sum(kernelpen(deltaoversigma,type=type,param=param))*log(length(profileValues$LogRatio),exp(1))

                  if (verbose)
                    {
                       print(paste("removeBreakpoints: likelihood for the profile = ", logsigma))
                  print(paste("removeBreakpoints: penalised likelihood for the profile = ", likelihood[i]))
                    }
                }
              
              
              
              if (min(likelihood)<likelihoodGLOBAL)
                {
                  regionRemoved <- which(likelihood==min(likelihood))[1]  #au cas ou il y ait plus d'une région
                  subset$profileValues$Breakpoints[which(subset$profileValues$Region==labelRegion[regionRemoved]&subset$profileValues$Breakpoints==1)] <- -1

                  if (verbose)
                    {
                      print("The following breakpoint has been removed:")
                      print(CGH[which(CGH$Region==labelRegion[regionRemoved]&CGH$Breakpoints==1),])
                      print(paste("Global -LogLibelihood:", likelihoodGLOBAL))
                      print(paste("-LogLikelihood without the breakpoint:", min(likelihood))[1])
                    }
                  CGH$Breakpoints[which(CGH$Region==labelRegion[regionRemoved]&CGH$Breakpoints==1)] <- -1
                  
                  indexRegionRemoved <- which(subset$profileValues$Region==labelRegion[regionRemoved+1])
                  subset$profileValues$Region[indexRegionRemoved] <- labelRegion[regionRemoved]
                  indexRegionRemoved <- which(CGH$Region==labelRegion[regionRemoved+1])
                  CGH$Region[indexRegionRemoved] <- labelRegion[regionRemoved]
                }

              else
                {
                  stop <- 1
                  LikeliRatio <- likelihood/likelihoodGLOBAL
                  infoBP <- data.frame(Region=labelRegion[1:(length(labelRegion)-1)], LikeliRatioAux=LikeliRatio)

                  if (length(LikeliRatio) >=1)
                    {
                      indexBP <- which(CGH$Breakpoints==1)
                      CGHBP <- CGH[indexBP,]
                      CGHWBP <- CGH[-indexBP,]
                      CGHBP <- merge(CGHBP,infoBP)
                      CGHBP$LikeliRatio <- CGHBP$LikeliRatioAux
                      indexcol <- which(names(CGHBP)=="LikeliRatioAux")
                      CGHBP <- CGHBP[,-indexcol]
                      CGH <- rbind(CGHBP,CGHWBP)
                      
                    }
                  
                  
                }
            }
          
          
          else
            {
              stop <- 1
            }
        }

      indexLikeliRatio <- which(names(CGH)=="LikeliRatio")      
      CGH <- rbind(CGH, CGHna)
      CGH <- CGH[,-indexLikeliRatio] 
      profileChr$profileValues <- CGH
      return(profileChr)
    }
  else
    {
      return(profileChr)
    }

  if (verbose) print("removeBreakpoints: ending function")

}




removeLevel <- function(...)
  {
    UseMethod("removeLevel")
  }



removeLevel.profileChr <- function(profileChr, lambda=10, type="tricubic", param=c(d=6), ...)
{
  CGH <- profileChr$profileValues
  indexNA <- attr(na.omit(CGH[,c("Chromosome","LogRatio","PosOrder")]),"na.action")

#############################################################################
#
#  Il faut tester si, pour certains chromosomes, toutes les valeurs sont NA
#
#############################################################################

  if (length(indexNA)<length(CGH[,1]))
    {
      CGHna <- CGH[indexNA,]
      if (!is.null(indexNA))
        {
          CGH <- CGH[-indexNA,]
        }
      
      subset <- list(profileValues=CGH)
      class(subset) <- "profileChr"
      

#########################################################################################
#
#  la fonction IQR permet d'estimer l'écart-type qui sera utiliser dans la fonction kernelpen
#
##########################################################################################

      
      IQRdiff <- function(y) IQR(diff(y))/1.908

# estimation de l'écart-type :
# les données doivent etre ordonnées à cause du diff
      subset$profileValues[order(subset$profileValues$PosOrder),]
      sigma <- IQRdiff(subset$profileValues$LogRatio)
      

# détection des outliers pour chacun des Levels
      subset <- detectOutliers(subset, region="Level", ...)
      indexOutliersTot <- which(subset$profileValues$OutliersTot==0)
      #profileValues <- subset$profileValues[indexOutliersTot,]
      profileValues <- subset$profileValues

# moyenne des LogRatios par Level
      aggregation <- aggregate(profileValues["LogRatio"],list(Level=profileValues$Level),mean)
      aggregation <- aggregation[order(aggregation$LogRatio),]

# les Levels sont ordonnés de manière croissante
      LevelOrder <- 1:length(unique(profileValues$Level))
      aggregation <- data.frame(aggregation, LevelOrder=LevelOrder)


# la valeur du champ région est comparable à celle du champ Level
      CGH <- merge(CGH, aggregation[,c("Level","LevelOrder")], by.x="Level", by.y="Level", all=TRUE)
      CGH <- CGH[order(CGH$LevelOrder),]

# sauvegarde des régions telles qu'elles sont en entrée
      CGH <- data.frame(CGH, RegionCGH=CGH$Region)
# on donne la valeur des Levels à la variable CGH$Region
# pour utiliser le meme script que dans la fonction removeBreakpoints
      CGH$Region <- CGH$LevelOrder
      CGH <- CGH[,-which(names(CGH)=="LevelOrder")]

      subset <- list(profileValues=CGH)
      class(subset) <- "profileChr"

      
      stop <- 0
      
      while (stop!=1)
	{
          subset <- detectOutliers(subset, region="Region", ...)
	  indexOutliersTot <- which(subset$profileValues$OutliersTot==0)          
	  profileValues <- subset$profileValues[indexOutliersTot,]         
	  aggregation <- aggregate(profileValues["LogRatio"],list(Region=profileValues$Region),mean)
	  aggregation <- aggregation[order(aggregation$Region),]
# attention : est-ce-que le diff veut tjs dire qqch comme on travaille maintenant sur les Levels?
	  deltaoversigma <- abs(diff(aggregation$LogRatio)/sigma)
	  
	  labelRegion <- sort(unique(profileValues$Region))

	  if (length(labelRegion)>1)
            {
              logsigma <- 0
              
	      for (i in 1:length(labelRegion))
                {
		  indexRegion <- which(profileValues$Region==labelRegion[i])
                  
		  if (length(indexRegion)>1)
                    {
		      logsigma <- logsigma + length(indexRegion)*log(var(profileValues$LogRatio[indexRegion])*((length(indexRegion)-1)/length(indexRegion)), exp(1))
                    }
                }
              
              
              
              
#likelihoodGLOBAL <- logsigma + lambda*(length(labelRegion)-1)*log(length(profileValues$LogRatio),exp(1))
              likelihoodGLOBAL <- logsigma + lambda*sum(kernelpen(deltaoversigma,type=type,param=param))*log(length(profileValues$LogRatio),exp(1))

              likelihood <- rep(NULL, length(labelRegion)-1)
              
              for (i in 1:(length(labelRegion)-1))
                {
                  RegionAux <- profileValues$Region
                  indexRegion <- which(profileValues$Region==labelRegion[i+1])
                  RegionAux[indexRegion] <- labelRegion[i]
                  indexRegion <- which(RegionAux==labelRegion[i])
                  
                  aggregation <- aggregate(profileValues["LogRatio"],list(Region=RegionAux),mean)
                  aggregation <- aggregation[order(aggregation$Region),]
                  deltaoversigma <- abs(diff(aggregation$LogRatio)/sigma)
                  
                  logsigma <-0
                  labelRegionAux <- unique(RegionAux)
                  
                  
                  for (j in 1:length(labelRegionAux))
                    {
                      indexRegionAux <- which(RegionAux==labelRegionAux[j])
                      
                      if (length(indexRegionAux)>1)
                        {
                          logsigma <- logsigma + length(indexRegionAux)*log(var(profileValues$LogRatio[indexRegionAux])*((length(indexRegionAux)-1)/length(indexRegionAux)), exp(1))
                        }
                    }
                  
                  likelihood[i] <- logsigma + lambda*sum(kernelpen(deltaoversigma,type=type,param=param))*log(length(profileValues$LogRatio),exp(1))
                }
              
              
              
              if (min(likelihood)<likelihoodGLOBAL)
                {
                  regionRemoved <- which(likelihood==min(likelihood))[1]  #au cas ou il y ait plus d'une région
#subset$profileValues$Breakpoints[which(subset$profileValues$Region==labelRegion[regionRemoved]&subset$profileValues$Breakpoints==1)] <- -1
#CGH$Breakpoints[which(CGH$Region==labelRegion[regionRemoved]&CGH$Breakpoints==1)] <- -1
                  indexRegionRemoved <- which(subset$profileValues$Region==labelRegion[regionRemoved+1])
                  subset$profileValues$Region[indexRegionRemoved] <- labelRegion[regionRemoved]
                  indexRegionRemoved <- which(CGH$Region==labelRegion[regionRemoved+1])
                  CGH$Region[indexRegionRemoved] <- labelRegion[regionRemoved]
                }

              else
                {
                  stop <- 1
                }
            }
          
          
          else
            {
              stop <- 1
            }
        }


# on n'a pas besoin de comparer si on est sur le meme chromosome car cette function est appliqué séparémment pour chaque chromosome
      CGH <- CGH[order(CGH$PosOrder),]

      for (i in 1:(length(CGH$PosOrder)-1))
        {
          if (CGH$Region[i]==CGH$Region[i+1])
            {
              CGH$Breakpoints[i] <- 0
            }
          else
            {
              CGH$Breakpoints[i] <- 1
            }
        }

      CGH$Level <- CGH$Region
      CGH$Region <- CGH$RegionCGH
      CGH <- CGH[,-which(names(CGH)=="RegionCGH")]


# renumérotation des Levels
      labelLevel <- unique(CGH$Level) 
      newnumber <- data.frame(Level=labelLevel, LevelAux=1:length(labelLevel))
      CGH <- merge(CGH, newnumber, by.x="Level", by.y="Level")
      CGH$Level <- CGH$LevelAux
      CGH <- CGH[,-which(names(CGH)=="LevelAux")]

# récupération des données manquantes
      CGH <- rbind(CGH, CGHna)
      profileChr$profileValues <- CGH
      return(profileChr)
    }
  else
    {
      return(profileChr)
    }
  
}
