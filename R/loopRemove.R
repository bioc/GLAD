
loopRemove <- function(...)
  {
    UseMethod("loopRemove")
  }



loopRemove.profileChr <- function(subset, CGH, sigma, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, ...)
  {
### attention de passer les paramètres pour detecter outliers
### et le paramètre verbose
### plus les paramètres du noyau
### le paramètre sigma claculé avec IQR.
### debut boucle

    stop <- 0
    
    while (stop!=1)
      {
        subset <- detectOutliers(subset, region="Region", verbose=verbose, ...)
        indexOutliersTot <- which(subset$profileValues$OutliersTot==0)
        profileValues <- subset$profileValues[indexOutliersTot,]

        ### utilisation de split et sapply
        sagg <- split(profileValues$LogRatio, profileValues$Region)
        Mean <- sapply(sagg,mean)
        #print(Mean)
        #print(names(Mean))
        Card <- sapply(sagg,NROW)
        Var <- sapply(sagg,var)
        Region <- as.numeric(as.character(names(Mean)))
        aggTot <- data.frame(Region, Card, Mean, Var)

        ##############################################################
        
#         aggregation <- aggregate(profileValues["LogRatio"],list(Region=profileValues$Region),mean)
#         names(aggregation) <- c("Region","Mean")
#         aggregation$Region <- as.numeric(as.character(aggregation$Region))
        aggTot <- aggTot[order(aggTot$Region),]
        #print("aggregation")
        #print(aggregation)
        deltaoversigma <- abs(diff(aggTot$Mean)/sigma)
        #print("deltaoversigma")
        #print(deltaoversigma)

        
        if (verbose)
          {
            print("")
          }


        
# ### Effectif
#         aggCard <- aggregate(profileValues$LogRatio, list(Region=profileValues$Region),NROW)
#         aggCard$Region <- as.numeric(as.character(aggCard$Region))
#         names(aggCard) <- c("Region","Card")
#         #print("aggCard vaut:")
#         #print(aggCard)

# ### Variance
#         aggVar <- aggregate(profileValues$LogRatio, list(Region=profileValues$Region),var)
#         aggVar$Region <- as.numeric(as.character(aggVar$Region))
#         names(aggVar) <- c("Region","Var")
#         #print("aggVar vaut:")
#         #print(aggVar)

# ### jointure entre les tables
#         aggTot <- merge(aggregation, aggCard, by="Region")
#         aggTot <- merge(aggTot, aggVar, by="Region")


        #print("aggTot vaut:")
        #print(aggTot)
        
### variance de l'estimateur du maximum de vraisemblance
        aggTot$Var <- aggTot$Var*((aggTot$Card-1)/aggTot$Card)
        aggTot$VarLike <- aggTot$Var
### les enregistrements pour lesquels le cardinal vaut 1
        indexSingle <- which(aggTot$Card==1)
        aggTot$Var[indexSingle] <- 0
        aggTot$VarLike[indexSingle] <- 1

        #print("aggTot vaut:")
        #print(aggTot)

        computeLike <- function(agg, lambda, sumkernelpen)
          {
### sumkernelpen = sum(kernelpen(deltaoversigma,type=type,param=param))
            logsigma <- agg$Card*log(agg$VarLike,exp(1))
            logsigma <- sum(logsigma)
            nbdata <- log(sum(agg$Card),exp(1))
            likelihoodGLOBAL <- logsigma + lambda*sumkernelpen*nbdata
            return(likelihoodGLOBAL)
          }

        

        nbregion <- dim(aggTot)[1]
        


        if (nbregion>1)
          {
                                        #             {
                                        #               logsigma <- aggTot$Card*log(aggTot$Var,exp(1))
                                        #               logsigma <- sum(logsigma)
                                        #               nbdata <- log(length(profileValues$LogRatio),exp(1)) 
                                        #               likelihoodGLOBAL <- logsigma + lambda*sum(kernelpen(deltaoversigma,type=type,param=param))*nbdata
            


            skpen <- sum(kernelpen(deltaoversigma,type=type,param=param))
            likelihoodGLOBAL <- computeLike(aggTot, lambda=lambda, sumkernelpen=skpen)
            #print("likelihoodGLOBAL")
            #print(likelihoodGLOBAL)

            likelihood <- rep(0, nbregion-1)

            for (i in 1:(nbregion-1))
              {
                indRegion <- i+1
                aggTotAux <- aggTot[-indRegion,]
                aggTotAux$Card[i] <- aggTot$Card[i] + aggTot$Card[indRegion]
                barycentre <- aggTot$Card[i]*aggTot$Mean[i] + aggTot$Card[indRegion]*aggTot$Mean[indRegion]
                barycentre <- barycentre/aggTotAux$Card[i]
                aggTotAux$Mean[i] <- barycentre
                within <- aggTot$Card[i]*aggTot$Var[i] + aggTot$Card[indRegion]*aggTot$Var[indRegion]
                within <- within/aggTotAux$Card[i]
                between <- aggTot$Card[i]*(aggTot$Mean[i]-barycentre)^2 + aggTot$Card[indRegion]*(aggTot$Mean[indRegion]-barycentre)^2
                between <- between/aggTotAux$Card[i]
                aggTotAux$Var[i] <- within + between
                aggTotAux$VarLike[i] <- aggTotAux$Var[i]

                #print("dans boucle for : aggTotAux")
                #print(aggTotAux)

                deltaoversigma <- abs(diff(aggTotAux$Mean)/sigma)
                skpen <- sum(kernelpen(deltaoversigma,type=type,param=param))

                likelihood[i] <- computeLike(aggTotAux, lambda=lambda, sumkernelpen=skpen)
                
              }

            #print("likelihood")
            #print(likelihood)
            
            if (min(likelihood)<likelihoodGLOBAL)
              {
                regionRemoved <- which(likelihood==min(likelihood))[1]  #au cas ou il y ait plus d'une région
                subset$profileValues$Breakpoints[which(subset$profileValues$Region==aggTot$Region[regionRemoved]&subset$profileValues$Breakpoints==1)] <- -1

                if (verbose)
                  {
                    print("The following Region has been removed:")
                    print(aggTot$Region[regionRemoved])
                    print(paste("Global -LogLibelihood:", likelihoodGLOBAL))
                    print(paste("-LogLikelihood without the breakpoint:", min(likelihood))[1])
                  }
                
                CGH$Breakpoints[which(CGH$Region==aggTot$Region[regionRemoved]&CGH$Breakpoints==1)] <- -1
                
                indexRegionRemoved <- which(subset$profileValues$Region==aggTot$Region[regionRemoved+1])
                subset$profileValues$Region[indexRegionRemoved] <- aggTot$Region[regionRemoved]
                indexRegionRemoved <- which(CGH$Region==aggTot$Region[regionRemoved+1])
                CGH$Region[indexRegionRemoved] <- aggTot$Region[regionRemoved]
                
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

    return(CGH)
    
  }
