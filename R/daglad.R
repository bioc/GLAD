### This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr


daglad <- function(...)
  {  
    UseMethod("daglad")
  }



daglad.profileCGH <- function(profileCGH, mediancenter = FALSE, normalrefcenter = FALSE, genomestep = FALSE,
                              OnlySmoothing = FALSE, OnlyOptimCall = FALSE, 
                              smoothfunc = "lawsglad", lkern = "Exponential", model = "Gaussian",
                              qlambda = 0.999,  bandwidth = 10, sigma = NULL, base = FALSE, round = 2,
                              lambdabreak = 8, lambdaclusterGen = 40, param = c(d = 6), alpha = 0.001, msize = 2,
                              method = "centroid", nmin = 1, nmax = 8, region.size = 2,
                              amplicon = 1, deletion = -5, deltaN = 0.10,  forceGL = c(-0.15,0.15), nbsigma = 3,
                              MinBkpWeight = 0.35, DelBkpInAmp=TRUE, CheckBkpPos = TRUE, assignGNLOut = TRUE,
                              breaksFdrQ = 0.0001, haarStartLevel = 1, haarEndLevel = 5, weights.name = NULL,
                              verbose = FALSE, ...)
  {

    IQRdiff <- function(y) IQR(diff(y))/1.908
    
    ## Récupération des paramètres de la fonction    
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
    profileCGH$nbsigma <- nbsigma    
    profileCGH$smoothfunc <- smoothfunc
    profileCGH$lambdabreak <- lambdabreak
    profileCGH$param <- param
    profileCGH$NbProbes <- length(profileCGH$profileValues[["PosOrder"]])
    profileCGH$TooSmall <- FALSE    

    ## ############################################################################
    ## Préparation des données
    ## ############################################################################

    ## vérification de l'object donné en entrée

    if (verbose) print("daglad - step CheckData")
    CheckData(profileCGH, bandwidth = bandwidth, smoothfunc = smoothfunc, weights.name = weights.name, OnlyOptimCall = OnlyOptimCall)

    if(base == TRUE)
      {
        if(!require(aws))
          {
            stop("Error in daglad: the aws package is required to use these function. The aws package can be installed from http://www.r-project.org")
          }
      }

    if((msize > region.size) & (region.size != 0))
      stop("Error in daglad: msize must be lower than region.size")

### cette condition a été ajouté pour les cas ou par exemple une régions de 3 clones comporte un outliers
### compte-tenu du fait que le clustering se fait sur des moyennes et non pas sur des médianes
### il peut y avoir une forte différence entre les deux valeurs ce qui peut conduire à une incohérence sur le GNL
    
    if(smoothfunc == "lawsglad" & model != "Gaussian")
      stop("Error in daglad: for lawsglad, only Gaussian model is available")

    if(smoothfunc == "lawsglad" & base == TRUE)
      stop("Error in daglad: for lawsglad it is not possible to use the option base=TRUE. Choose laws smoothfunc instead.")

    
    if (smoothfunc != "lawsglad" && smoothfunc != "haarseg")      
      {
        print(paste("You have chosen smoothfunc=", smoothfunc))
        print(paste("Choose smoothfunc=lawsglad or smoothfunc=haarseg if you want the process runs faster"))
      }

    inputfields <- names(profileCGH$profileValues)

    if (!OnlyOptimCall)
      {
        
        excdudefields <- c("Level", "OutliersAws", "OutliersMad",
                           "OutliersTot", "Breakpoints", "Smoothing",
                           "NormalRef", "ZoneGNL")
      }
    else
      {
        excdudefields <- c("Level", "OutliersAws", "OutliersMad",
                           "OutliersTot", "Breakpoints", "Region",
                           "NormalRef", "ZoneGNL")
          
      }
    
    fieldstodel <- intersect(inputfields, excdudefields)
    if (length(fieldstodel) > 0)
      {
        print("Error in daglad: the following fields must be removed from profileValues before starting the function:")
        print(fieldstodel)
        stop()
      }
    

    ## on trie les données par chromosome et position
    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues[["Chromosome"]], profileCGH$profileValues[["PosOrder"]]),]
    
    ## ajout des champs nécessaires à la procédure
    ## et transformation des données en liste

    if (OnlyOptimCall)
      {
        new.fields <- c("NewPosOrder",
                        "OutliersAws", "Region",
                        "Level", "Breakpoints",
                        "MinPosOrder", "MaxPosOrder",
                        "OutliersMad",  "OutliersTot",
                        "NextLogRatio", "NormalRange",
                        "ZoneGen", "ZoneGNL")

        ## on force à haarseg de telle sorte que daglad définisse les points de cassure
        ## et les régions
        smoothfunc <- "haarseg"
      }
    else
      {
        
        new.fields <- c("NewPosOrder", "Smoothing",
                        "OutliersAws", "Region",
                        "Level", "Breakpoints",
                        "MinPosOrder", "MaxPosOrder",
                        "OutliersMad",  "OutliersTot",
                        "NextLogRatio", "NormalRange",
                        "ZoneGen", "ZoneGNL")
      }

    nb.new.fields <- length(new.fields)

    if (smoothfunc == "haarseg")
      {
        profileCGH$profileValues <- as.list(profileCGH$profileValues)
        profileCGH$profileValues[new.fields] <- lapply(as.list(rep(profileCGH$NbProbes, nb.new.fields)), numeric)
      }
    else
      {
        profileCGH$profileValues[new.fields] <- 0
      }

    
    ## Il faut des PosOrder uniques
    profileCGH$profileValues[["NewPosOrder"]] <- profileCGH$profileValues[["PosOrder"]]
    profileCGH$profileValues[["PosOrder"]] <- 1:profileCGH$NbProbes
    


    ## LogRatio are median-centered
    if (mediancenter)
      {
        med <- median(profileCGH$profileValues[["LogRatio"]])
        profileCGH$profileValues[["LogRatio"]] <- profileCGH$profileValues[["LogRatio"]] - med
      }


    ## ############################################################################
    ##
    ## Début de l'analyse daglad
    ##
    ## ############################################################################

    ## #####################################################
    ## AWS smoothing chromosome par chromosome
    ## #####################################################    

    print("Smoothing for each Chromosome")
    profileCGH <- chrBreakpoints(profileCGH, smoothfunc = smoothfunc, OnlyOptimCall = OnlyOptimCall, base = base, sigma = sigma,
                                 bandwidth=bandwidth, lkern=lkern, model=model, qlambda=qlambda,
                                 round=round, verbose=verbose,
                                 breaksFdrQ=breaksFdrQ, haarStartLevel=haarStartLevel , haarEndLevel=haarEndLevel, weights.name = weights.name)

    ## estimation de l'écart-type chromosome par chromosome
    profileCGH$SigmaC <- profileCGH$Sigma
    profileCGH$Sigma <- NULL

    if (OnlySmoothing)
      {

        junk.fields <- c("NewPosOrder",
                        "OutliersAws",
                        "MinPosOrder", "MaxPosOrder",
                        "OutliersMad",  "OutliersTot",
                        "NextLogRatio", "NormalRange",
                        "ZoneGen", "ZoneGNL")
        

        ## conversion en data.frame
        profileCGH$profileValues <- as.data.frame(profileCGH$profileValues[setdiff(names(profileCGH$profileValues),junk.fields)])
        
        return(profileCGH)
      }
    
    ## #####################################################
    ##  AWS smoothing sur l'ensemble du génome
    ## #####################################################   

    if (genomestep)
      {

        
        profileCGH <- dogenomestep(profileCGH, nb.new.fields = nb.new.fields, new.fields = new.fields,
                                   smoothfunc = smoothfunc, lkern = lkern, model = model,
                                   qlambda = qlambda,  bandwidth = bandwidth, sigma = sigma, base = FALSE, round = round,
                                   lambdabreak = lambdabreak, lambdaclusterGen = lambdaclusterGen, param = param, alpha = alpha, msize = msize,
                                   method = method, nmin = nmin, nmax = nmax,
                                   amplicon = amplicon, deletion = deletion, deltaN = deltaN,  forceGL = forceGL, nbsigma = nbsigma,
                                   MinBkpWeight = MinBkpWeight, DelBkpInAmp=DelBkpInAmp, CheckBkpPos = CheckBkpPos, assignGNLOut = assignGNLOut,
                                   breaksFdrQ = breaksFdrQ, haarStartLevel = haarStartLevel, haarEndLevel = haarEndLevel, weights.name = weights.name,
                                   verbose = verbose)

        
      } ## fin (if) de l'étape sur le génome
    else
      {
        profileCGH$NormalRef <- 0

        if (is.null(sigma))
          {
            IQRinfoG <- IQRdiff(profileCGH$profileValues[["LogRatio"]])
          }
        else
          {
            IQRinfoG <- sigma
          }
        profileCGH$SigmaG <- data.frame(Chromosome = 0,Value = IQRinfoG)
        profileCGH$findClusterSigma <- profileCGH$SigmaG$Value[1]


      } ### fin (else) de l'étape sur le génome


    
    FieldOrder <- names(profileCGH$profileValues)

    ## position des points de cassure absolus (extrémités de chromosomes)
    profileCGH$AbsoluteBkp <- data.frame(Chromosome = c(0, profileCGH$PosOrderRange$Chromosome),
                                         PosOrder = c(0, profileCGH$PosOrderRange$MaxPosOrder),
                                         AbsoluteBkp = 1)


    print("Optimization of the Breakpoints and DNA copy number calling")

    ## Optimization of the Breakpoints and DNA copy number calling
    profileCGH <- OptimBkpFindCluster(profileCGH)
 
    
    ## suppression des points de cassure qui délimitent une région trop petite (il fau faire les deux étapes)
    profileCGH <- DelRegionTooSmall(profileCGH, region.size = region.size) ### STEP 1
    profileCGH <- DelRegionTooSmall(profileCGH, region.size = region.size) ### STEP 2
                 
 
    
    ## Calcul d'un poids pour les Breakpoints
    ## Attention: comme on calcul une variable GNLchange
    ## il ne faut pas le statut des outliers pour définir le changement de statut
    ## car un outlier peut être aussi un breakpoints

    if (verbose) print("daglad - step BkpInfo")

    profileCGH$BkpInfo <- BkpInfo(profileCGH)


    
    if (verbose) print("daglad - step OutliersGNL")


    ## ###############################################################################
    ##  GNL des Outliers
    ## ###############################################################################        
        
    
    ## on peut récupérer directement les paramètre passés à OutliersGNL depuis l'objet profileCGH
    profileCGH <- OutliersGNL(profileCGH, alpha = alpha, sigma = profileCGH$SigmaG$Value[1], NormalRef = profileCGH$NormalRef,
                              amplicon = amplicon, deletion = deletion, assignGNLOut = assignGNLOut)


    


    ## ###############################################################################
    ##  Bkp filter
    ## ###############################################################################        


    
    
    if (verbose) print("daglad - step filterBkpStep (pass 1)")
    profileCGH <- filterBkpStep(profileCGH, MinBkpWeight=MinBkpWeight, DelBkpInAmp=DelBkpInAmp, assignGNLOut=assignGNLOut, verbose=verbose)


    
    ## ###############################################################################
    ##  Déplacement des Bkp
    ## ###############################################################################    

    if (CheckBkpPos)
      {
        if (verbose) print("daglad - step MoveBkpStep")
        profileCGH <- MoveBkpStep(profileCGH, assignGNLOut=assignGNLOut)
      }



    if (verbose) print("daglad - step filterBkpStep (pass 2)")    
    profileCGH <- filterBkpStep(profileCGH, MinBkpWeight=MinBkpWeight, DelBkpInAmp=DelBkpInAmp, assignGNLOut=assignGNLOut, verbose=verbose)    
  
    
    print("Results Preparation")
    
    ###################################
    ## suppression des champs inutiles
    ###################################

    profileCGH <- prepare.output.daglad(profileCGH = profileCGH, genomestep = genomestep, normalrefcenter = normalrefcenter, inputfields = inputfields)    
    
    return (profileCGH)
    
    
  }


prepare.output.daglad <- function(profileCGH = NULL, genomestep = NULL, normalrefcenter = NULL, inputfields = NULL)
  {

    
    if (genomestep)
      {
        profileCGH$profileValues[["NormalRef"]] <- profileCGH$profileValues[["ZoneGNLGen"]]
        profileCGH$profileValues <- profileCGH$profileValues[setdiff(names(profileCGH$profileValues),"ZoneGNLGen")]
      }

    if (normalrefcenter)
      {
        profileCGH$profileValues[["Smoothing"]] <- profileCGH$profileValues[["Smoothing"]] - profileCGH$NormalRef
        profileCGH$profileValues[["LogRatio"]] <- profileCGH$profileValues[["LogRatio"]] - profileCGH$NormalRef
        
        
        if (is.data.frame(profileCGH$BkpInfo))
          {
            profileCGH$BkpInfo$LogRatio <- profileCGH$BkpInfo$LogRatio - profileCGH$NormalRef
            profileCGH$BkpInfo$Smoothing <- profileCGH$BkpInfo$Smoothing - profileCGH$NormalRef
            profileCGH$BkpInfo$SmoothingNext <- profileCGH$BkpInfo$SmoothingNext - profileCGH$NormalRef
          }

        profileCGH$NormalRef <- 0
      }

    ## conversion en data.frame
    profileCGH$profileValues <- as.data.frame(profileCGH$profileValues)
    
    
    if (is.data.frame(profileCGH$BkpInfo))
      {

        profileCGH$BkpInfo <- merge(profileCGH$BkpInfo, profileCGH$profileValues[c("PosOrder","NewPosOrder")], by="PosOrder",all.x=TRUE)


        profileCGH$BkpInfo <- profileCGH$BkpInfo[setdiff(names(profileCGH$BkpInfo),c("SmoothingNext","NextLogRatio","MaxPosOrder","MinPosOrder","ZoneGNLnext","PosOrder"))]

        namesprofile <- names(profileCGH$BkpInfo)
        namesprofile[which(namesprofile == "NewPosOrder")] <- "PosOrder"
        names(profileCGH$BkpInfo) <- namesprofile

      }


    profileCGH$profileValues <- profileCGH$profileValues[setdiff(names(profileCGH$profileValues),c("Sigma","DiffBase","Region","NextLogRatio","MaxPosOrder","MinPosOrder","PosOrder", "NormalRange", "ZoneGen"))]

    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile=="NewPosOrder")] <- "PosOrder"
    names(profileCGH$profileValues) <- namesprofile        
    

    outputfields <- setdiff(names(profileCGH$profileValues),inputfields)
    at <- setdiff(attributes(profileCGH)$names,c("PosOrderRange","findClusterSigma","NbClusterOpt"))
    profileCGH <- profileCGH[at]
    profileCGH$profileValues <- profileCGH$profileValues[c(inputfields,outputfields)]
    class(profileCGH) <- "profileCGH"
    

    return(profileCGH)
    
  }




DelRegionTooSmall <- function(profileCGH, region.size = 0, verbose = FALSE)
{


#  print("suis dans la fonction DelRegionTooSmall")
  
  if(region.size == 0)
    return(profileCGH)

  
  BkpInfoTmp <- BkpInfo(profileCGH)
  if(!is.data.frame(BkpInfoTmp))
    return(profileCGH)    

  RegionSize <- data.frame(BkpInfoTmp[c("Chromosome", "PosOrder")], AbsoluteBkp = 0)
  RegionSize <- rbind(profileCGH$AbsoluteBkp, RegionSize)
  RegionSize <- RegionSize[order(RegionSize$Chromosome, RegionSize$PosOrder),]
  RegionSize$indice <- 1:length(RegionSize[,1])
  RegionSize$Size <- c(0,diff(RegionSize$PosOrder))
  RegionSize$contig <- 0
  RegionSize$contig.first <- 0
  RegionSize$contig.last <- 0    


  ## récupération des régions trop petites
  ind.region <- which((RegionSize$Size <= region.size) & (RegionSize$Chromosome != 0))

  if(length(ind.region) > 0)
    {
      RegionWithSmallSize <- RegionSize[ind.region,]
      

      ## a-t-on des régions contigues      
#      print("regions contigues")
      RegionWithSmallSize <- RegionWithSmallSize[order(RegionWithSmallSize$Chromosome, RegionWithSmallSize$PosOrder),]     
#      print("taille de petites regions")
#      print(RegionWithSmallSize[,-c(2:3)])      
      ind.contig <- which(diff(RegionWithSmallSize$indice) == 1)
#     print(ind.contig)
      RegionWithSmallSize$contig[ind.contig] <- 1
      ind.contig.first <- which(diff(RegionWithSmallSize$contig) == 1)
      ind.contig.last <- which(diff(RegionWithSmallSize$contig) == -1)      
      RegionWithSmallSize$contig.first[1+ind.contig.first] <- 1
      RegionWithSmallSize$contig.last[ind.contig.last] <- 1            
#     print(RegionWithSmallSize[,-c(2:3)])      
      
      BkpToDel <- NULL

      ## breakpoint dans les chromosomes        
 #      ind.notabsolute <- which((RegionWithSmallSize$AbsoluteBkp == 0) & (RegionWithSmallSize$contig.first != 1) & (RegionWithSmallSize$contig.last != 1) & (RegionWithSmallSize$contig != 1))
      ind.notabsolute <- which((RegionWithSmallSize$AbsoluteBkp == 0) & (RegionWithSmallSize$contig.first != 1) & (RegionWithSmallSize$contig.last != 1) )


#      ind.notabsolute <- which(RegionWithSmallSize$AbsoluteBkp == 0)

      
      if(length(ind.notabsolute) > 0)
        {
          BkpToDel <- RegionWithSmallSize$indice[ind.notabsolute]
        }

      ## breakpoint aux extrémités des chromosomes
      ind.absolute <- which(RegionWithSmallSize$AbsoluteBkp == 1)

      
      if(length(ind.absolute) > 0)
        {
          ## BkpToDel <- c(BkpToDel, RegionWithSmallSize$indice[ind.absolute] + 1)
            BkpToDel <- c(BkpToDel, RegionWithSmallSize$indice[ind.absolute] - 1)  ## fix PG
        }

      if(length(BkpToDel) > 0)
        {
          PosOrderToDel <- RegionSize$PosOrder[BkpToDel]
        }

      
      ind.bkp <- which(profileCGH$profileValues[["PosOrder"]] %in% PosOrderToDel)


      profileCGH$profileValues$Breakpoints[ind.bkp] <- -1



      l <- length(profileCGH$profileValues[[1]])
      res <- .C("updateLevel",
                as.integer(profileCGH$profileValues[["Chromosome"]]),
                as.integer(profileCGH$profileValues[["Breakpoints"]]),
                Level = as.integer(profileCGH$profileValues[["Level"]]),
                as.integer(profileCGH$profileValues[["PosOrder"]]),
                NextLogRatio = as.double(profileCGH$profileValues[["NextLogRatio"]]),
                as.double(profileCGH$profileValues[["LogRatio"]]),                                
                as.integer(max(profileCGH$profileValues[["Level"]])),
#                Smoothing = as.double(profileCGH$profileValues[["Smoothing"]]),
                as.integer(l),
                PACKAGE = "GLAD")


#      profileCGH$profileValues[c("Level", "NextLogRatio", "Smoothing")] <- res[c("Level", "NextLogRatio", "Smoothing")]
      profileCGH$profileValues[c("Level", "NextLogRatio")] <- res[c("Level", "NextLogRatio")]      

      profileCGH$TooSmall <- TRUE      
      
    }


    ## #########################################
    ## On optimise à nouveau les points de cassure si des petites régions ont été supprimées
    ## #########################################

    
    if(profileCGH$TooSmall)
      {
        if (verbose)
          print("daglad - Small regions detected")

        ## Optimization of the Breakpoints and DNA copy number calling
        profileCGH <- OptimBkpFindCluster(profileCGH)
                
      }
    
  return(profileCGH)

}



dogenomestep <- function(profileCGH, nb.new.fields = NULL, new.fields = NULL,
                         smoothfunc = "lawsglad", lkern = "Exponential", model = "Gaussian",
                         qlambda = 0.999,  bandwidth = 10, sigma = NULL, base = FALSE, round = 2,
                         lambdabreak = 8, lambdaclusterGen = 40, param = c(d = 6), alpha = 0.001, msize = 5,
                         method = "centroid", nmin = 1, nmax = 8,
                         amplicon = 1, deletion = -5, deltaN = 0.10,  forceGL = c(-0.15,0.15), nbsigma = 3,
                         MinBkpWeight = 0.35, DelBkpInAmp=TRUE, CheckBkpPos = TRUE, assignGNLOut = TRUE,
                         breaksFdrQ = 0.0001, haarStartLevel = 1, haarEndLevel = 5, weights.name = NULL,
                         verbose = FALSE, ...)
  {
    ## il y a déjà les noms Smoothing Level Region Breakpoints OutliersAws
    ## ces champs sont renommés pour qu'ils soient spécifiques
    ## de l'étape chromosome par chromosome:
    ## ils portent à la fin la lettre C

    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile == "Smoothing")] <- "SmoothingC"    
    namesprofile[which(namesprofile == "Level")] <- "LevelC"
    namesprofile[which(namesprofile == "Region")] <- "RegionC"
    namesprofile[which(namesprofile == "Breakpoints")] <- "BreakpointsC"
    namesprofile[which(namesprofile == "OutliersAws")] <- "OutliersAwsC"
    namesprofile[which(namesprofile == "MinPosOrder")] <- "MinPosOrderC"
    namesprofile[which(namesprofile == "MaxPosOrder")] <- "MaxPosOrderC"
    names(profileCGH$profileValues) <- namesprofile
    
    if (verbose) print("daglad - step genomestep")        
    ## pour le smoothing sur l'ensemble du génome, il ne faut qu'un seul chromosome
    profileCGH$profileValues[["ChromosomeTrue"]] <- profileCGH$profileValues[["Chromosome"]]
    profileCGH$profileValues[["Chromosome"]] <- integer(profileCGH$NbProbes)

    ## Indicateur des Bkp détéctés chromosome par chromosome
    BkpDetectedAux <- profileCGH$BkpDetected

    ## Range des PosOrder chromosome par chromosome
    PosOrderRangeAux <- profileCGH$PosOrderRange        

    ## AWS smoothing sur l'ensemble du génome
    ## ici, on impose un qlamba élevé (l'utilisateur n'a pas le choix du paramètre)
    ## pour dégager la tendance des données


    ## ajout des champs nécessaires
    new.fields <- setdiff(new.fields,"NewPosOrder")
    nb.new.fields <- length(nb.new.fields)
    
    if (smoothfunc == "haarseg")
      {
        profileCGH$profileValues[new.fields] <- lapply(as.list(rep(profileCGH$NbProbes, nb.new.fields)), numeric)
      }
    else
      {
        profileCGH$profileValues[new.fields] <- 0
      }

    
    
    print("Smoothing over the genome")
    profileCGH <- chrBreakpoints(profileCGH, smoothfunc = smoothfunc, base = FALSE, sigma = sigma,
                                 bandwidth = bandwidth, round = round, verbose = verbose,
                                 lkern = lkern, model = model, qlambda = 0.9999999,
                                 breaksFdrQ = breaksFdrQ, haarStartLevel = haarStartLevel , haarEndLevel = haarEndLevel, weights.name = weights.name)

    ## estimation de l'écart-type sur l'ensemble du génome    
    profileCGH$SigmaG <- profileCGH$Sigma
    profileCGH$Sigma <- NULL
    profileCGH$AbsoluteBkp <- NULL


    ## réinitialisation de BkpDetected aux valeurs détéctées chromosome par chromosome
    profileCGH$BkpDetected <- BkpDetectedAux

    ## réinitialisation de PosOrderRange aux valeurs détéctées chromosome par chromosome        
    profileCGH$PosOrderRange <- PosOrderRangeAux
    

    ## maintenant, le champ ChromosomeTrue ne sert plus à rien
    ## on récupère le nom des chromosomes de départ
    profileCGH$profileValues <- profileCGH$profileValues[setdiff(names(profileCGH$profileValues),c("Chromosome"))]        
    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile == "ChromosomeTrue")] <- "Chromosome"
    names(profileCGH$profileValues) <- namesprofile        
    
    ## les champs Smoothing Level Region Breakpoints MinPosOrder MaxPosOrder sont supprimés
    ## on conserve le champ Level qui est renommé en LevelG
    ## on réinitialise les Min et MaxposOrder aux valeurs chromosome par chromosome
    fields.kept <- setdiff(names(profileCGH$profileValues),c("Smoothing","Region","Breakpoints","MinPosOrder","MaxPosOrder"))
    profileCGH$profileValues <- profileCGH$profileValues[fields.kept]

    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile == "Level")] <- "LevelG"
    namesprofile[which(namesprofile == "MinPosOrderC")] <- "MinPosOrder"
    namesprofile[which(namesprofile == "MaxPosOrderC")] <- "MaxPosOrder"
    names(profileCGH$profileValues) <- namesprofile

    
    ## #####################################################
    ## on cherche la ligne de base à partir du
    ## smoothing sur l'ensemble du génome
    ## #######################################################

    ## l'object doit etre de type profileChr
    ## à l'avenir, il faudra supprimer cette classe
    ## (attention à faire les modifs également dans la fonction glad)
    class(profileCGH) <- "profileChr"

    ## détection des outliers
    profileCGH <- detectOutliers(profileCGH, region="LevelG", alpha=alpha, msize=msize)

    
    ## le clustering est fait sur les niveaux LevelG
    print("Find cluster over the Genome")
    profileCGH$findClusterSigma <- profileCGH$SigmaG$Value[1]
    profileCGH <- findCluster(profileCGH, region = "LevelG", method = method, genome = TRUE,
                              lambda = lambdaclusterGen, nmin = nmin, nmax = nmax, param = profileCGH$param)
    class(profileCGH) <- "profileCGH"


    ## On a besoin d'un champ OutliersAws:
    ## on ne considère comme Outliers AWS que ceux détectés        
    ## à l'étape chromosome
    profileCGH$profileValues[["OutliersAws"]] <- profileCGH$profileValues[["OutliersAwsC"]]

    
    print("Find the Normal Baseline")
    ## on trouve le normal (ligne de base) à partir de la fonction affectationGNL
    profileCGH <- affectationGNL(profileCGH, verbose=verbose)

    
    
    ## notre niveau de référence pour le normal correspond
    ## à la médiane des log-ratios correspondant à une ZoneGNLGen de 0
    profileCGH$NormalRef <- median(profileCGH$profileValues[["LogRatio"]][which(profileCGH$profileValues[["ZoneGNL"]] == 0)], na.rm=TRUE)

    ## on supprime les champs LevelG et ZoneGNL
    fields.kept <- setdiff(names(profileCGH$profileValues),c("LevelG", "OutliersAws", "ZoneGNL"))
    profileCGH$profileValues <- profileCGH$profileValues[fields.kept]

    ## on récupère les champs calculés chromosome par chromosome
    namesprofile <- names(profileCGH$profileValues)
    namesprofile[which(namesprofile == "SmoothingC")] <- "Smoothing"    
    namesprofile[which(namesprofile == "LevelC")] <- "Level"
    namesprofile[which(namesprofile == "RegionC")] <- "Region"
    namesprofile[which(namesprofile == "BreakpointsC")] <- "Breakpoints"
    namesprofile[which(namesprofile == "OutliersAwsC")] <- "OutliersAws"
    names(profileCGH$profileValues) <- namesprofile


    ## ajout des champs nécessaires
    profileCGH$profileValues[setdiff(new.fields, names(profileCGH$profileValues))] <- 0

    return(profileCGH)
    
  }

### 01052009 Suppression de cette étape
### Les Breakpoints consécutifs ont été supprimés
### dans les fonctions précédentes

#####################        
### DEBUT SUPPRESSION
#####################            

##     RecomputeGNL <- FALSE
##     ## Epuration des Breakpoints consécutifs
##     if (is.data.frame(profileCGH$BkpInfo))
##       {
##         profileCGH$BkpInfo <- profileCGH$BkpInfo[order(profileCGH$BkpInfo$PosOrder),]
##         profileCGH$BkpInfo$BkpToDel <- rep(0, length(profileCGH$BkpInfo[,1]))
##         profileCGH$BkpInfo$Side <- profileCGH$BkpInfo$BkpToDel
##         profileCGH$BkpInfo$NextPosOrder <- profileCGH$BkpInfo$PosOrder + 1

##         ## on regarde que si il y a au moins 2 Bkp
##         if (length(profileCGH$BkpInfo[,1])>1)
##           {

##             deleteContiguousBkp <- .C("delete_contiguous_bkp",
##                                       BkpToDel=as.integer(profileCGH$BkpInfo$BkpToDel),
##                                       Gap=as.double(profileCGH$BkpInfo$Gap),
##                                       LogRatio=as.double(profileCGH$BkpInfo$LogRatio),
##                                       NextPosOrder=as.integer(profileCGH$BkpInfo$NextPosOrder),
##                                       PosOrder=as.integer(profileCGH$BkpInfo$PosOrder),
##                                       Side=as.integer(profileCGH$BkpInfo$Side),
##                                       Sigma=as.double(profileCGH$BkpInfo$Sigma),
##                                       Smoothing=as.double(profileCGH$BkpInfo$Smoothing),
##                                       SmoothingNext=as.double(profileCGH$BkpInfo$SmoothingNext),
##                                       Weight=as.double(profileCGH$BkpInfo$Weight),
##                                       as.integer(length(profileCGH$BkpInfo[,1])),
##                                       RecomputeGNL=as.integer(0),
##                                       as.integer(nbsigma),                                      
##                                       PACKAGE="GLAD")


##             profileCGH$BkpInfo[,c("BkpToDel",
##                                   "Gap",
##                                   "LogRatio",
##                                   "NextPosOrder",
##                                   "PosOrder",
##                                   "Side",
##                                   "Sigma",
##                                   "Smoothing",
##                                   "SmoothingNext",
##                                   "Weight")] <- deleteContiguousBkp[c("BkpToDel",
##                                                                       "Gap",
##                                                                       "LogRatio",
##                                                                       "NextPosOrder",
##                                                                       "PosOrder",
##                                                                       "Side",
##                                                                       "Sigma",
##                                                                       "Smoothing",
##                                                                       "SmoothingNext",
##                                                                       "Weight")]

##             RecomputeGNL <- deleteContiguousBkp$RecomputeGNL

##           }


##         ## ICI on peut encore optimiser en supprimant la boucle for
##         indexBPtoDel <- which(profileCGH$BkpInfo$BkpToDel==1)
##         if (length(indexBPtoDel)>0)
##           {

##             profile <- profileCGH$profileValues
##             profile <- profile[order(profile$PosOrder),]

##             for (ind in indexBPtoDel)
##               {

##                 profile <- profile[order(profile$PosOrder),]
##                 indextochange <- profileCGH$BkpInfo$PosOrder[ind]
##                 profile$Breakpoints[indextochange] <- 0

##                 if (profileCGH$BkpInfo$Side[ind]==0)
##                   {
##                     profile$Level[indextochange+1] <- profile$Level[indextochange]
##                     profile$Smoothing[indextochange+1] <- profile$Smoothing[indextochange]
##                   }
##                 else
##                   {
##                     profile$Level[indextochange] <- profile$Level[indextochange+1]
##                     profile$Smoothing[indextochange] <- profile$Smoothing[indextochange+1] 
##                   }

##               }


##             profileCGH$profileValues <- profile

##             profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexBPtoDel,]

##           }


##         profileCGH$BkpInfo <- profileCGH$BkpInfo[,setdiff(names(profileCGH$BkpInfo),c("Side","BkpToDel","NextPosOrder"))]

##       }


##     if (verbose) print("daglad - step recomputeGNL")
##     if (RecomputeGNL)
##       {


##         profileCGH$profileValues <- profileCGH$profileValues[,setdiff(names(profileCGH$profileValues),c("ZoneGNL","ZoneGen"))]

##         ## on prend comme référence ceux qui sont compris entre certaines valeurs
##         l <- length(profileCGH$profileValues$NormalRef)
##         NormalRange <- .C("compute_NormalRange",
##                           as.double(profileCGH$profileValues$Smoothing),
##                           as.double(profileCGH$profileValues$NormalRef),
##                           as.integer(profileCGH$profileValues$Level),
##                           NormalRange=integer(l),
##                           as.double(deltaN),
##                           as.integer(l),
##                           PACKAGE="GLAD")


##         ## le clustering est fait sur les niveaux NormalRange
##         class(profileCGH) <- "profileChr"
##         profileCGH <- findCluster(profileCGH, region="NormalRange", method=method, genome=TRUE,
##                                   lambda=lambdaclusterGen, nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)
##         class(profileCGH) <- "profileCGH"

##         ## le cluster correspondant au normal est celui qui comprend
##         ## le NormalRange 0
##         lengthDest <- length(profileCGH$profileValues$ZoneGen)
##         myZoneGNL <- .C("compute_cluster_LossNormalGain",
##                         ## variables pour la jointure
##                         as.integer(profileCGH$profileValues$ZoneGen),
##                         ZoneGNL=integer(lengthDest),
##                         as.integer(lengthDest),
##                         as.double(profileCGH$profileValues$Smoothing),
##                         as.double(profileCGH$forceGL[1]),
##                         as.double(profileCGH$forceGL[2]),
##                         as.double(profileCGH$NormalRef),
##                         as.double(profileCGH$amplicon),
##                         as.double(profileCGH$deletion),                                                                                    
##                         ## variables pour le calcul de la médiane par cluster
##                         as.double(profileCGH$profileValues$LogRatio),
##                         as.integer(profileCGH$profileValues$NormalRange),
##                         PACKAGE="GLAD")


##         profileCGH$profileValues$ZoneGNL <- myZoneGNL$ZoneGNL       

##       }

###################    
### FIN SUPPRESSION    
###################    


  ## ## ajout phupe
  ## print("ajout phupe")
  ## calc.range <- function(profileCGH = NULL)
  ## {
  ##   ind <- which(profileCGH$profileValues[["OutliersTot"]] == 0 )
  ##   LogRatio <- profileCGH$profileValues[["Smoothing"]][ind]
  ##   ZoneGNL <- profileCGH$profileValues[["ZoneGNL"]][ind]

  ##   print("max")
  ##   print(aggregate(LogRatio, by=list("GNL"=ZoneGNL), max))

  ##   print("min")
  ##   print(aggregate(LogRatio, by=list("GNL"=ZoneGNL), min))    
    
  ##   return(profileCGH)
  ## }

  ##   ## fin ajout phupe


OptimBkpFindCluster <- function(profileCGH = NULL)
  {

    ## choix de la méthode de clustering
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method <- pmatch(profileCGH$method, METHODS)


    if (is.na(method)) 
      stop("invalid clustering method")
    if (method == -1) 
      stop("ambiguous clustering method")
    
    
    startChr <- profileCGH$PosOrderRange$MinPosOrder - 1 ### car on commence à compter à 0
    sizeChr <- profileCGH$PosOrderRange$MaxPosOrder - profileCGH$PosOrderRange$MinPosOrder + 1

    
    NbChr <- length(startChr)
    l <- profileCGH$NbProbes

    resLoopChr <- .C("daglad_OptmisationBreakpoints_findCluster",
                     as.integer(profileCGH$profileValues[["Chromosome"]]),
                     Smoothing = double(l), ## valeur de sortie
                     NormalRange = integer(l),
                     as.double(profileCGH$NormalRef),
                     as.double(profileCGH$deltaN),
                     as.double(profileCGH$profileValues[["LogRatio"]]),
                     NextLogRatio = as.double(profileCGH$profileValues[["NextLogRatio"]]),   ## valeur de sortie
                     as.integer(profileCGH$profileValues[["PosOrder"]]),
                     Level = as.integer(profileCGH$profileValues[["Level"]]),                ## valeur de sortie
                     OutliersAws = as.integer(profileCGH$profileValues[["OutliersAws"]]),    ## valeur de sortie
                     OutliersMad = as.integer(profileCGH$profileValues[["OutliersMad"]]),    ## valeur de sortie
                     OutliersTot = as.integer(profileCGH$profileValues[["OutliersTot"]]),    ## valeur de sortie
                     Breakpoints = as.integer(profileCGH$profileValues[["Breakpoints"]]),    ## valeur de sortie
                     as.integer(profileCGH$msize),
                     as.double(qnorm(1-profileCGH$alpha/2)),
                     as.double(profileCGH$lambdabreak),
                     as.double(profileCGH$param["d"]),
                     as.double(profileCGH$SigmaC$Value),
                     as.integer(NbChr),   ## Nombre de chromosome a analyser
                     as.integer(sizeChr), ## taille de chaque chromosome
                     as.integer(startChr),## position pour le debut des valeurs de chaque chromosome
                     as.integer(profileCGH$BkpDetected$BkpDetected),
                     ## paramètres pour findCluster
                     as.integer(method),
                     as.double(profileCGH$findClusterSigma),
                     as.double(profileCGH$lambdaclusterGen),
                     as.integer(profileCGH$nmin),
                     as.integer(profileCGH$nmax),
                     ZoneGen = integer(l),  ## valeur de sortie
                     nbclasses = integer(1),
                     ## paramètres pour le calcul du GNL
                     ZoneGNL = integer(l),
                     as.double(profileCGH$forceGL[1]),
                     as.double(profileCGH$forceGL[2]),
                     as.double(profileCGH$NormalRef),
                     as.double(profileCGH$amplicon),
                     as.double(profileCGH$deletion),                             
                     as.integer(l), ## nombre total de sondes
                     PACKAGE = "GLAD")
    

    ## #########################################
    ## Récupération des résultats
    ## #########################################
    
    fields.replaced <- c("Smoothing", "NextLogRatio","Level", "OutliersAws", "OutliersMad", "OutliersTot", "Breakpoints", "NormalRange", "ZoneGen", "ZoneGNL")
    profileCGH$profileValues[fields.replaced] <- resLoopChr[fields.replaced]

    profileCGH$NbClusterOpt <- resLoopChr[["nbclasses"]]


    return(profileCGH)

    
  }
