### This function detects chromosomal breakpoints along genome

### Copyright (C) 2003 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2003
### Contact: glad@curie.fr


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


    subset <- profileChr
        class(subset) <- "profileChr"
        

#########################################################################################
###
###  la fonction IQR permet d'estimer l'écart-type qui sera utilisé dans la fonction kernelpen
###
##########################################################################################

    sigma <- profileChr$findClusterSigma


### Appel de la fonction loopRemove
        CGH <- loopRemove(subset, CGH, sigma, lambda=lambda, type=type, param=param, verbose=verbose, ...)


        profileChr$profileValues <- CGH
        return(profileChr)


    if (verbose) print("removeBreakpoints: ending function")

  }

