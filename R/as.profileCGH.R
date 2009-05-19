as.profileCGH <- function(object, ...)
  {
    UseMethod("as.profileCGH")
  }

as.profileCGH.data.frame <- function(object, infaction=c("value","empty"), value=20, ...)
  {

    print("changer la description de l'objet profileCGH dans la doc")
    
    profileCGH <- object
    infaction <- match.arg(infaction)

    champ.obligatoire <- c("LogRatio","PosOrder","Chromosome")
    
    nomchamp <- intersect(c(champ.obligatoire, "PosBase"), names(profileCGH))

    nomchamp.supp <- setdiff(names(profileCGH), nomchamp)
    print("Les champs en plus sont:")
    print(nomchamp.supp)

    addedfields <- c("ChromosomeChar","Smoothing","Region","Level",
                     "OutliersAws","Breakpoints","OutliersMad","OutliersTot",
                     "ZoneChr","ZoneGen","ZoneGNL")
    interadded <- intersect(addedfields,names(profileCGH))


    if (length(interadded) != 0)
      {
        print(paste("Error in as.profileCGH.data.frame: the following fields already exist:"))
        print(interadded)
        stop("Please remove those fields from your data.frame.")
      }
    if (length(intersect(champ.obligatoire,names(profileCGH))) != 3)
      {
        stop(paste("Error in as.profileCGH.data.frame: the following fields are missing:", paste(setdiff(champ.obligatoire,names(profileCGH))), sep=", "))
      }

    if (!is.numeric(profileCGH$profileValues$Chromosome))
      {
        ChromosomeNum <- ChrNumeric(profileCGH$Chromosome)
        profileCGH$ChromosomeChar <- profileCGH$Chromosome        
        profileCGH$Chromosome <- ChromosomeNum

      }

    
    ## Suppression des valeurs manquantes et des LogRatio avec Inf    
    indexInf <- which(is.infinite(profileCGH$LogRatio)==TRUE)
    if (length(indexInf)>0)
      {          
        print("The LogRatio with following rows index contains infinite value:")
        print(indexInf)
        
        profileCGH$LogRatio[indexInf] <- sign(profileCGH$LogRatio[indexInf]) * switch(infaction,
                                                                                      empty = NA,
                                                                                      value = value)
        

        printinfaction <- switch(infaction,
                                 empty = "LogRatio with infinite values have been replaced by NA",
                                 value = paste("LogRatio with infinite values have been replaced by + or - ",value))

        print(printinfaction)
        

      }

    indexNA <- attr(na.omit(profileCGH[,nomchamp]),"na.action")

    if (!is.null(indexNA))
      {
        if (length(nomchamp.supp) > 0)
          {
            profileCGH <- list(profileValues = profileCGH[-indexNA, nomchamp],
                               profileValuesSUPP = profileCGH[-indexNA, nomchamp.supp],
                               profileValuesNA = profileCGH[indexNA,])            
          }
        else
          {
            profileCGH <- list(profileValues = profileCGH[-indexNA,],
                               profileValuesSUPP = profileCGH[-indexNA, nomchamp.supp],
                               profileValuesNA = profileCGH[indexNA,])
          }
        class(profileCGH) <- "profileCGH"

      }

    else
      {
        if (length(nomchamp.supp) > 0)
          {
          }
        else
          {        
            profileCGH <- list(profileValues = profileCGH)
          }
        
        class(profileCGH) <- "profileCGH"
      }


    return(profileCGH)
    
  }


as.data.frame.profileCGH <- function(x, row.names = NULL, optional = FALSE, ...)
  {

    print("Il faut modifier cette fonction")
    profileCGH <- x
    if (!is.null(profileCGH$profileValuesNA))
      {
        values <- profileCGH$profileValues
        valuesNA <- profileCGH$profileValuesNA
        nbNA <- length(valuesNA[,1])
        valuesAux <- values[1:nbNA,]
        CommonFields <- intersect(names(values),names(valuesNA))
        MissingFields <- setdiff(names(values),names(valuesNA))
        valuesAux[,CommonFields] <- valuesNA[,CommonFields]
        valuesAux[,MissingFields] <- NA
        values <- rbind(values,valuesAux)
        values <- values[order(values$Chromosome,values$PosOrder),]
        rownames(values) <- 1:length(values[,1])
      }
    else
      {
        values <- profileCGH$profileValues
        values <- values[order(values$Chromosome,values$PosOrder),]
        rownames(values) <- 1:length(values[,1])
      }

    if (length(which(names(profileCGH$profileValues) == "ChromosomeChar")))
      {
        values$Chromosome <- values$ChromosomeChar
        values <- subset(values, select=setdiff(names(values),"ChromosomeChar"))
      }


    return(values)

  }
