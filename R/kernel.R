# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

#param is a named vector for the parameter of the kernel

kernelpen <- function(x, type="tricubic", param)
	{
		k <- rep(0,length(x))
		if (type=="tricubic")
		{
			if (missing(param)) stop ("set parameters for your kernel")
			if (length(which(x<0)>0)) stop("kernel function is not defined for negative numbers")
			index <- which(x<=param["d"])
			k[index] <- (1-(x[index]/param["d"])^3)^3
		}
	
    return(k)
    
}
