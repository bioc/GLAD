.First.lib <- function(lib, pkg){

    library.dynam("GLAD", pkg, lib)
#    if(!require(aws))stop("You must install the package named aws")
    print("Have fun with GLAD")
    print("For smoothing it is possible to use either")
    print("the AWS algorithm (Polzehl and Spokoiny, 2002)")
    print("or the HaarSeg algorithm (Ben-Yaacov and Eldar, 2008)")
    print("")
    print("If you use the package with AWS, please cite:")
    print("Hupe et al. (2004) and Polzehl and Spokoiny (2002)")
    print("")    
    print("If you use the package with HaarSeg, please cite:")
    print("Hupe et al. (2004) and (Ben-Yaacov and Eldar, 2008)")
    
}
