### Copyright (C) 2006 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2006
### Contact: glad@curie.fr
### http://bioinfo.curie.fr



CheckData <- function(...)
  {
    UseMethod("CheckData")
  }

CheckData.profileCGH <- function(profileCGH=profileCGH, bandwidth=bandwidth, ...)
  {

    n <- dim(profileCGH$profileValues)[1]    
    if(n==0)
      {
        stop("Error: the data contains only missing values. Check that the fields LogRatio, Chromosome or PosOrder are not empty.")
      }


    if (n>200)
      {
        if (bandwidth>1)
          {
            print("You can set bandwitdth to 1 to decrease computation time")
          }
      }
    
  }
