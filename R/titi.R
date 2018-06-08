## Hierarchical clustering, on raw input data; we will use Euclidean
## distance.  A range of criteria are supported; also there is a
## storage-economic option.
##
## We use the general routine, `hc', which caters for 7 criteria,
## using a half dissimilarity matrix; (BTW, this uses the very efficient
## nearest neighbor chain algorithm, which makes this algorithm of
## O(n^2) computational time, and differentiates it from the less
## efficient -- i.e. O(n^3) -- implementations in all commercial
## statistical packages -- as far as I am aware -- except Clustan.)
##
## Clustering Methods:
##
## 1. Ward's minimum variance or error sum of squares method.
## 2. single linkage or nearest neighbor method.
## 3. complete linkage or diameter.
## 4. average linkage, group average, or UPGMA method.
## 5. McQuitty's or WPGMA method.
## 6. median, Gower's or WPGMC method.
## 7. centroid or UPGMC method (7).
##
## Original author: F. Murtagh, May 1992
## R Modifications: Ross Ihaka, Dec 1996
##		    Friedrich Leisch, Apr 1998, Jun 2000
##                  Antoine Lucas, Oct 2002

hclustglad <-  function (d, method = "complete", members = NULL) 
  { 
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
      stop("invalid clustering method")
    if (method == -1) 
      stop("ambiguous clustering method")
    n <- as.integer(attr(d, "Size"))
    if (is.null(n)) 
      stop("invalid dissimilarities")
    if (n < 2) 
      stop("Must have n >= 2 objects to cluster")
    labels <- attr(d, "Labels")
    len <- n * (n - 1)/2
    if (is.null(members)) 
      members <- rep(1, n)
    if (length(members) != n) 
      stop("Invalid length of members")

    ## Ajout Philippe HupÃ©
    ## modification lorsque ward est utilisÃ©
    if (method == 1)
      {
      	coeff <- matrix(rep(members,n),n,n)
	coeff <- coeff*t(coeff)/(coeff + t(coeff))
	d <- as.dist(as.matrix(d*d)*coeff)
      }

##     print("method")
##     print(method)
    
    hcl <- .C("hclust",
              n = as.integer(n),
              len = as.integer(len), 
              method = as.integer(method),
              ia = integer(n),
              ib = integer(n),
              order = integer(n),
              crit = double(n),
              members = as.double(members),
              diss = as.double(d),
              res  = as.integer (1),
              PACKAGE ="GLAD" )
    
    if(hcl$res == 2)
      stop("Cannot allocate memory")
    if(hcl$res == 1)
      stop("Error")

    tree <- list(merge = cbind(hcl$ia[1:(n - 1)],
                   hcl$ib[1:(n - 1)]),
                 height = hcl$crit[1:(n - 1)],
                 order = hcl$order, 
                 labels = attr(d, "Labels"),
                 method = METHODS[method], 
                 call = match.call())
    if (!is.null(attr(d, "method"))) {
      tree$dist.method <- attr(d, "method")
    }

    class(tree) <- "hclust"
    tree
  }

