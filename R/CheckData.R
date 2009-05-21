### Copyright (C) 2006 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2006
### Contact: glad@curie.fr
### http://bioinfo.curie.fr



CheckData <- function(...)
  {
    UseMethod("CheckData")
  }

CheckData.profileCGH <- function(profileCGH=profileCGH, bandwidth=bandwidth, smoothfunc=smoothfunc, ...)
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

##     if(!is.numeric(profileCGH$profileValues$Chromosome))
##       stop("Error: profileCGH$profileValues$Chromosome must be numeric")

    
  }
