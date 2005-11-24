.First.lib <- function(lib, pkg){

    library.dynam("GLAD", pkg, lib)
#    if(!require(aws))stop("You must install the package named aws")
    print("Have fun with GLAD")
    
    
}
