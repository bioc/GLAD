
lawsglad <- function(y,x=NULL,qlambda=NULL,eta=0.5,lkern="Triangle",model="Poisson",
                 shape=NULL,hinit=NULL,hincr=NULL,hmax=10,NN=FALSE,u=NULL,
                 graph=FALSE,demo=FALSE,symmetric=FALSE,wghts=NULL)
  {

###
###    first check arguments and initialize
###
    args <- match.call()
    eps <- 1.e-10
    if(is.null(qlambda)) if(symmetric) qlambda <- switch(model,
                                                         Gaussian=.985,
                                                         Bernoulli=.985,
                                                         Exponential=.985,
                                                         Poisson=.985,
                                                         Weibull=.985,
                                                         Volatility=.995)
    else qlambda <- switch(model,   Gaussian=.966,
                           Bernoulli=.966,
                           Exponential=.966,
                           Poisson=.966,
                           Weibull=.966,
                           Volatility=.98)
    if(qlambda>=1 || qlambda<.6) return("Inappropriate value of qlambda")
    if(eta<eps || eta>=1) return("Inappropriate value of eta")
    
    if(is.null(hinit)||hinit<1) hinit <- 1
    if(is.null(hincr)) hincr <- 1.25
    lamakt <- 5*qchisq(qlambda,1)
    if(model=="Gaussian")
      {
        lamakt <- lamakt*shape*2
      }
    
### tabulation du noyau
    getkern <- function(x,kern)
      switch(kern,Triangle=pmax(0,(1-x)),
             Quadratic=pmax(0,(1-x))^2,
             Cubic=pmax(0,(1-x))^3,
             Uniform=as.numeric(abs(x)<=1),
             Exponential=exp(-5*x),
             {
               cat("Triangle kernel is used as default\n");
               pmax(0,(1-x))
             })

    kernl <- getkern(seq(0,1.01,.01),lkern)
    kerns <- getkern(seq(0,1.01,.01),"Exponential")

    #print(length(kernl))

    n <- length(y)

##     z <- .C("iawsuni",
##             as.double(y),
##             as.integer(n),
##             as.double(hinit),
##             bi=double(n),
##             ai=double(n),
##             as.double(kernl),PACKAGE="GLAD")[c("bi","ai")]

##     bi <- z$bi
##     ai <- z$ai
##     theta <- ai/bi
    
##     Gaussian <- .C("gawsuni",
##                    as.double(y),
##                    as.integer(n),
##                    as.double(hinit),
##                    as.double(hincr),
##                    as.double(hmax),
##                    as.double(lamakt),
##                    as.double(eta),
##                    theta=as.double(theta),
##                    as.double(bi),
##                    as.double(ai),
##                    as.double(kernl),
##                    as.double(kerns),
##                    as.double(bi),
##                    PACKAGE="GLAD")$theta

    Gaussian <- .C("lawsglad",
                   as.double(y),
                   as.integer(n),
                   as.double(hinit),
                   as.double(hincr),
                   as.double(hmax),
                   as.double(lamakt),
                   as.double(eta),
                   theta=double(n),
                   double(n),
                   double(n),
                   as.double(kernl),
                   as.double(kerns),
                   double(n),
                   PACKAGE="GLAD")$theta

    
    return(Gaussian)

  }

