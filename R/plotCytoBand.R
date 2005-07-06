
plotCytoBand <- function(...)
  {
    UseMethod("plotCytoBand")
  }


plotCytoBand.default <- function(cytoband, y=-1, Chromosome=1, labels=TRUE, height=1,
                                 colCytoBand=c("white","darkblue"), colCentro="red", ...)
  {

    Color <- unique(cytoband$Col)
    NbColor <- length(Color)
    pal <- myPalette(low=colCytoBand[1], high=colCytoBand[2], k=NbColor)

    info <- data.frame(Color=Color, ColorName=pal)
    cytoband <- merge(cytoband,info,by="Color")
    
    
### Info sur les cytobandes
    indexChr <- which(cytoband$Chromosome==Chromosome)
    dataChr <- cytoband[indexChr,]

    
    CytoPos <- 0.5*(dataChr$Start + dataChr$End)
    CytoLength <- (dataChr$End - dataChr$Start)+1
    NbCyto <- length(dataChr[,1])
    HeightPlot <- rep(height, NbCyto)

    SizeCyto <- matrix(c(CytoLength,HeightPlot),NbCyto, 2, byrow=FALSE)

    symbols(CytoPos, y=rep(min(unique(y)),NbCyto), rectangles=SizeCyto, inches=FALSE, bg=as.character(dataChr$ColorName), add=TRUE, ...)

### Ajout de la position du centromère
    indexCentro <- which(dataChr$Centro==1)
    centroPos <- min(dataChr$End[indexCentro])

    arrows(centroPos, min(unique(y)) + height/2, centroPos, min(unique(y)) - height/2, col=colCentro, code=3, angle=120)

### Ajout du label de la bande

    dataChr$Band <- paste(dataChr$Chromosome,dataChr$Band,sep="")

    if (labels)
      {
        axis(side=3, at=CytoPos, labels=as.character(dataChr$Band), las=2)
      }

    
    
  }
