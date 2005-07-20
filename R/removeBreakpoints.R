### This function detects chromosomal breakpoints along genome

### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr


removeBreakpoints <- function(...)
  {
    UseMethod("removeBreakpoints")
  }



removeBreakpoints.profileChr <- function(profileChr, lambda=10, type="tricubic", param=c(d=6), verbose=FALSE, msize=5, alpha=0.001,...)
  {


    if (verbose)
      {
        print("removeBreakpoints: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }
    
    

    sigma <- profileChr$findClusterSigma


### Appel de la fonction loopRemove
        profileChr <- loopRemove(profileChr, sigma, lambda=lambda,
                                  type=type, param=param, verbose=verbose, msize=msize, alpha=alpha)


        return(profileChr)


    if (verbose) print("removeBreakpoints: ending function")

  }

