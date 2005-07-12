#Array plot functions

# Copyright (C) 2003 Institut Curie
# Author(s): Philippe Hupé (Institut Curie) 2003
# Contact: glad@curie.fr

arrayPlot <- function(...)
  {
    UseMethod("arrayPlot")
  }

arrayPlot.arrayCGH <- function(arrayCGH, variable, mediancenter=FALSE,
                               col=myPalette("green", "red", "yellow"),
                               contour=FALSE, nlevels=5, zlim=NULL,
                               bar=c("none", "horizontal", "vertical"),
                               layout=TRUE, ...)
  {
  bar <- match.arg(bar)
  if (bar!="none") 
    omar <- par()$mar  # store old margins
  olas <- par()$las    # store old las

  if (nrow(arrayCGH$arrayValues) != prod(arrayCGH$arrayDesign)) 
    stop("Array dimension does not match with input data")
  if (length(which(names(arrayCGH$arrayValues) == variable)) < 1) 
    stop(paste("variable", variable, "not found", sep = " "))
  arrayCGH$arrayValues <- arrayCGH$arrayValues[order(arrayCGH$arrayValues$Col, -arrayCGH$arrayValues$Row), ]
  if (mediancenter == FALSE) 
    z <- arrayCGH$arrayValues[[variable]]
  else z <- arrayCGH$arrayValues[[variable]] - median(na.omit(arrayCGH$arrayValues[[variable]]))
  
  if (is.null(zlim)) 
    zlim <- c(min(na.omit(z)), max(na.omit(z)))
  else {
    z[which(z <= zlim[1])] <- zlim[1]
    z[which(z >= zlim[2])] <- zlim[2]
  }
  ncol <- arrayCGH$arrayDesign[1] * arrayCGH$arrayDesign[3]
  nrow <- arrayCGH$arrayDesign[2] * arrayCGH$arrayDesign[4]
  z <- matrix(z, ncol, nrow, byrow=TRUE)
  
  if (bar=="horizontal") {
    par(mar = c(1, 1, 4, 1))
    if(layout)
      layout(matrix(1:2, 2, 1), height=c(8,1))
  }
  if (bar=="vertical") {
    par(mar = c(1, 1, 4, 1))
    if (layout)
      layout(matrix(1:2, 1, 2), width=c(8,1))
  }
  
  ## ajust 'col' vector to 'zlim'
  l <- length(col)
  s <- seq(zlim[1], zlim[2], length=l)
  minz <- min(z, na.rm=TRUE)
  maxz <- max(z, na.rm=TRUE)
  if (minz!=maxz)
    {
      ind1 <- floor((minz-zlim[1])/(zlim[2]-zlim[1])*l)
      ind2 <- floor((maxz-zlim[1])/(zlim[2]-zlim[1])*l)
      colimage <- col[ind1:ind2]
    }
  else
    {
      colimage <- col
    }


  image(1:ncol, 1:nrow, z, axes = FALSE, col = colimage, xlab = "", ylab = "", ...)

  if (contour) 
    contour(1:ncol, 1:nrow, z, nlevels = nlevels, add = TRUE)
  box(lwd = 1)
  abline(v = ((1:arrayCGH$arrayDesign[1] - 1) * arrayCGH$arrayDesign[3] + 0.5), lwd = 0.5)
  abline(h = ((1:arrayCGH$arrayDesign[2] - 1) * arrayCGH$arrayDesign[4] + 0.5), lwd = 0.5)

  if (bar!="none") {
    if (zlim[1] != zlim[2]) 
      x.bar <- seq(zlim[1], zlim[2], length = 41)
    else x.bar <- seq((zlim[1] - 1), (zlim[2] + 1), length = 41)

    if (bar=="horizontal")
      par(mar = c(4, 1, 0, 1), las=2)
    if (bar=="vertical")
      par(mar = c(1, 0, 4, 4))
    
    ColorBar(x.bar, horizontal = (bar=="horizontal"), col = col, main = "", k=7)
    par(mar = omar, las=olas) ## restore old margins
  }
}

arrayPlot.default <- function(Statistic, Col, Row, ArrCol, ArrRow, SpotCol, SpotRow,mediancenter=FALSE, col=myPalette("green", "red", "yellow"), contour=FALSE, nlevels=5, zlim=NULL, bar=c("none", "horizontal", "vertical"), layout=TRUE, ...)
{
  if(length(Statistic)!=ArrCol*ArrRow*SpotCol*SpotRow)
    stop("Array dimension does not match with input data")

### data must be sorted by column and descending row
  data <- data.frame(Statistic=Statistic, Col=Col, Row=Row)
  data <- data[order(data$Col,-data$Row),]

  arrayCGH <- list(arrayValues=data, arrayDesign=c(ArrCol, ArrRow, SpotCol, SpotRow))
  class(arrayCGH) <- "arrayCGH"

  arrayPlot.arrayCGH(arrayCGH, "Statistic", mediancenter=mediancenter, col=col, contour=contour, nlevels=nlevels, zlim=zlim, bar=bar, layout=layout, ...)
}

arrayPersp <- function(...)
  {
    UseMethod("arrayPersp")
  }


arrayPersp.default<-function(Statistic, Col, Row,
                             ArrCol, ArrRow, SpotCol, SpotRow,
                             mediancenter=FALSE,
                             col=myPalette("green","red","yellow"),
                             zlim=zlim, bar=TRUE, ...)
  {

  if(length(Statistic)!=ArrCol*ArrRow*SpotCol*SpotRow)
    stop("Array dimension does not match with input data")

                                        #data must be sorted by column and descending row
  data <- data.frame(Statistic=Statistic, Col=Col, Row=Row)
  data <- data[order(data$Col,-data$Row),]



  if (mediancenter==FALSE)
    z <- data$Statistic

  else
    z <- data$Statistic - median(na.omit(data$Statistic))

  if(missing(zlim))
    zlim <- c(min(na.omit(z)), max(na.omit(z)))

  else
    {
      z[which(z<=zlim[1])] <- zlim[1]
      z[which(z>=zlim[2])] <- zlim[2]
    }


  z <- matrix(z, ArrCol*SpotCol, ArrRow*SpotRow, byrow=TRUE)
  ncol <- ArrCol*SpotCol
  nrow <- ArrRow*SpotRow

  par(mar=c(4,4,4,4))

  if(bar)
    {
      layout(matrix(c(1,2),1,2),width=c(9,1))

      if(names(dev.cur())!="pdf")
        par(mar=c(4,3,5,3)) # settings for graphics device

      else
        par(mar=c(4,3,5,1)) # settings for writing pdf files
    }

  fcol <- matrix("green", nr=nrow(z), nc=ncol(z))
  zi <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z), -ncol(z)])/4
  breaks <- seq(zlim[1],zlim[2],length=length(col))
  fcol <- col[cut(zi, breaks, include.lowest = TRUE)]
  
  persp(1:ncol, 1:nrow, z, col=fcol, ...)


  if(bar)
    {
      if (zlim[1]!=zlim[2])
        x.bar <- seq(zlim[1], zlim[2], length=41)

      else
        x.bar <- seq((zlim[1]-1), (zlim[2]+1), length=41)
      

      if(names(dev.cur())!="pdf")
        par(mar=c(4,0,5,3))  # settings for graphics device

      else
        par(mar=c(4,0,5,5)) # settings for writing pdf files

      ColorBar(x.bar,horizontal=FALSE,col=col,main="")
    }

  layout(1)
  par(mar=c(5, 4, 4, 2) + 0.1)

}


arrayPersp.arrayCGH<-function(arrayCGH, variable,
			      mediancenter=FALSE,
			      col=myPalette("green","red","yellow"),
			      zlim=zlim, bar=TRUE, ...)
  {
  if(nrow(arrayCGH$arrayValues)!=prod(arrayCGH$arrayDesign))
    stop("Array dimension does not match with input data")
  if(length(which(names(arrayCGH$arrayValues)==variable))<1)stop(paste("variable",variable,"not found",sep=" "))


                                        #data must be sorted by column and descending row
  arrayCGH$arrayValues <- arrayCGH$arrayValues[order(arrayCGH$arrayValues$Col,-arrayCGH$arrayValues$Row),]
                                        #data <- data.frame(Statistic=Statistic, Col=Col, Row=Row)
                                        #data <- data[order(data$Col,-data$Row),]



  if (mediancenter==FALSE)
    z <- arrayCGH$arrayValues[[variable]]

  else
    z <-  arrayCGH$arrayValues[[variable]] - median(na.omit( arrayCGH$arrayValues[[variable]]))

  if(missing(zlim))
    zlim <- c(min(na.omit(z)), max(na.omit(z)))

  else
    {
      z[which(z<=zlim[1])] <- zlim[1]
      z[which(z>=zlim[2])] <- zlim[2]
    }

  ncol <- arrayCGH$arrayDesign[1]*arrayCGH$arrayDesign[3]
  nrow <- arrayCGH$arrayDesign[2]*arrayCGH$arrayDesign[4]
  z <- matrix(z, ncol, nrow, byrow=TRUE)

  par(mar=c(1, 1, 4, 1))

  if(bar)
    {
      layout(matrix(c(1,2),1,2),width=c(9,1))

      if(names(dev.cur())!="pdf")
        par(mar=c(4,3,5,3)) # settings for graphics device

      else
        par(mar=c(4,3,5,1)) # settings for writing pdf files
    }

  fcol <- matrix("green", nr=nrow(z), nc=ncol(z))
  zi <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z), -ncol(z)])/4
  breaks <- seq(zlim[1],zlim[2],length=length(col))
  fcol <- col[cut(zi, breaks, include.lowest = TRUE)]
  
  persp(1:ncol, 1:nrow, z, col=fcol, ...)


  if(bar)
    {
      if (zlim[1]!=zlim[2])
        x.bar <- seq(zlim[1], zlim[2], length=41)

      else
        x.bar <- seq((zlim[1]-1), (zlim[2]+1), length=41)
      

      if(names(dev.cur())!="pdf")
        par(mar=c(4,0,5,3))  # settings for graphics device

      else
        par(mar=c(4,0,5,5)) # settings for writing pdf files

      ColorBar(x.bar,horizontal=FALSE,col=col,main="")
    }

  layout(1)
  par(mar=c(5, 4, 4, 2) + 0.1)

}

myPalette <- function(low = "white",
                      high = c("green", "red"),
                      mid=NULL,
                      k =50)
  {

    low <- col2rgb(low)/255
    high <- col2rgb(high)/255
    
    if(is.null(mid))
      {
        r <- seq(low[1], high[1], len = k)
        g <- seq(low[2], high[2], len = k)
        b <- seq(low[3], high[3], len = k)
      }

    if(!is.null(mid))
      {
        k2 <- round(k/2)
        mid <- col2rgb(mid)/255
        r <- c(seq(low[1], mid[1], len = k2),
               seq(mid[1], high[1], len = k2))
        g <- c(seq(low[2], mid[2], len = k2),
               seq(mid[2], high[2], len = k2))
        b <- c(seq(low[3], mid[3], len = k2),
               seq(mid[3], high[3], len = k2))
      }
    rgb(r, g, b)
  }


ColorBar<-function(x, horizontal = TRUE, col=heat.colors(50),
                   scale=1:length(x), k=10,  ...)
  {

    if(is.numeric(x))
      {
        x <- x
        colmap <- col
      }

    else
      {
        colmap <- x
        low<-range(scale)[1]
        high<-range(scale)[2]
        x <- seq(low, high, length=length(x))
      }

    if(length(x)>k)
      x.small<-seq(x[1], x[length(x)],length=k)

    else
      x.small<-x

    if(horizontal)
      {
        image(x, 1, matrix(x,length(x),1), axes=FALSE, xlab="", ylab="", col=colmap, ...)
        axis(1, at=rev(x.small), labels=signif(rev(x.small),2), srt=270)
      }

    if(!horizontal)
      {
        image(1, x, matrix(x,1,length(x)), axes=FALSE, xlab="", ylab="", col=colmap, ...)
        par(las=1)
        axis(4, at=rev(x.small), labels=signif(rev(x.small), 2))
        par(las=0) # Back to default
      }

    box()
    
  }
