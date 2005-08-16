as.profileCGH <- function(object)
  {
    UseMethod("as.profileCGH")
  }

as.profileCGH.data.frame <- function(object)
  {

    profileCGH <- object
    nomchamp <- c("LogRatio","PosOrder","Chromosome")

    addedfields <- c("ChromosomeChar","Smoothing","Region","Level",
                     "OutliersAws","Breakpoints","OutliersMad","OutliersTot",
                     "ZoneChr","ZoneGen","ZoneGNL")
    interadded <- intersect(addedfields,names(profileCGH))


    if (length(interadded)!=0)
      {
        print(paste("Error in as.profileCGH.data.frame: the following fields already exist:"))
        print(interadded)
        stop("Please remove those fields from your data.frame.")
      }
    if (length(intersect(nomchamp,names(profileCGH)))!=3)
      {
        stop(paste("Error in as.profileCGH.data.frame: the following fields are missing:", paste(setdiff(nomchamp,names(profileCGH))), sep=", "))
      }

    if (!is.numeric(profileCGH$profileValues$Chromosome))
      {
        ChromosomeNum <- ChrNumeric(profileCGH$Chromosome)
        profileCGH$ChromosomeChar <- profileCGH$Chromosome        
        profileCGH$Chromosome <- ChromosomeNum

      }

    
### Suppression des valeurs manquantes et des LogRatio avec Inf    
    indexNA <- attr(na.omit(profileCGH[,nomchamp]),"na.action")
    indexInf <- which(is.finite(profileCGH$LogRatio)==FALSE)
    indexNA <- c(indexNA,indexInf)
    
    if (!is.null(indexNA))
      {
        profileCGH <- list(profileValues=profileCGH[-indexNA,], profileValuesNA=profileCGH[indexNA,])
        class(profileCGH) <- "profileCGH"

      }

    else
      {
        profileCGH <- list(profileValues=profileCGH)
        class(profileCGH) <- "profileCGH"
      }


###

    return(profileCGH)
    
  }


as.data.frame.profileCGH <- function(x, row.names = NULL, optional = FALSE)
  {
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

    if (length(which(names(profileCGH$profileValues)=="ChromosomeChar")))
      {
        values$Chromosome <- values$ChromosomeChar
        values <- subset(values, select=-ChromosomeChar)
      }


    return(values)

  }
