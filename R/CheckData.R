### Copyright (C) 2006 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2006
### Contact: glad@curie.fr
### http://bioinfo.curie.fr



CheckData <- function(...)
  {
    UseMethod("CheckData")
  }

CheckData.profileCGH <- function(profileCGH=profileCGH, ...)
  {
    if(dim(profileCGH$profileValues)[1]==0)
      {
        stop("Error: the data contains only missing values. Check that the fields LogRatio, Chromosome or PosOrder are not empty.")
      }
  }
