### Copyright (C) 2006 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2006
### Contact: glad@curie.fr
### http://bioinfo.curie.fr



CheckData <- function(...)
  {
    UseMethod("CheckData")
  }

CheckData.profileCGH <- function(profileCGH = profileCGH, bandwidth = bandwidth,
                                 smoothfunc = smoothfunc, weights.name = NULL, OnlyOptimCall = FALSE,...)
  {

    n <- dim(profileCGH$profileValues)[1]
    nb.chr <- length(unique(profileCGH$profileValues[,"Chromosome"]))
    
    if(n == 0)
      {
        stop("Error: the data contains only missing values. Check that the fields LogRatio, Chromosome or PosOrder are not empty.")
      }


    if(smoothfunc != "haarseg")
      {
        if ((n / nb.chr) > 2000) ### 2000 ~ 50K/24      
          {
            if (bandwidth > 1)
              {
                print("You can set bandwitdth to 1 to decrease computation time")
              }
          }
      }

    
    if (!is.null(weights.name))
      {
        if(length(which(names(profileCGH$profileValues) == weights.name)) != 1)
          {
            stop(paste("Variable", weights.name, "used for weights has not been found"))
          }
        else
          {
            ind <- which(profileCGH$profileValues[[weights.name]] < 0)
            if(length(ind) > 0)
              {
                stop("Weights must be positive")
              }
          }
      }

    if(OnlyOptimCall)
      {
        if(length(which(names(profileCGH$profileValues) == "Smoothing")) != 1)
          {
            stop("When OnlyOptimCall = TRUE, the fields Smoothing must exist in profileCGH$profileValues")
          }        
      }
    
  }
