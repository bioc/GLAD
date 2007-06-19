### This function detects chromosomal breakpoints along genome

### Copyright (C) 2005 Institut Curie
### Author(s): Philippe Hupé (Institut Curie) 2005
### Contact: glad@curie.fr


plotProfile <- function(...)
  {
    UseMethod("plotProfile")
  }

plotProfile.profileCGH <- function(profileCGH, variable="LogRatio", Chromosome=NULL,
                                   Smoothing=NULL, GNL="ZoneGNL", Bkp=FALSE,
                                   labels=TRUE, plotband=TRUE, unit=0,
                                   colDAGLAD=c("black","blue","red","green","yellow"),
                                   pchSymbol=c(20,13),
                                   colCytoBand=c("white","darkblue"),
                                   colCentro="red", text=NULL,
                                   main="", ...)
  {

    
### il faudrait vérifier si le champ Length existe déjà en entrée de la table, sinon
### pb de jointure
    if (length(intersect(names(profileCGH$profileValues),"PosBase"))<1)
      {
        stop("Error in plotProfile.profileCGH: PosBase is not available")
      }

    if (!is.null(Smoothing))
      {
            if (length(intersect(names(profileCGH$profileValues),Smoothing))<1)
              {
                print(paste("Warning in plotProfile.profileCGH:", Smoothing," is not available"))
              }
      }
    
    if (Bkp)
      {
        if (length(intersect(names(profileCGH$profileValues),"Breakpoints"))<1)
          {
            print("Warning in plotProfile.profileCGH: Breakpoints is not available")
          }
      }



    profileCGH$profileValues$VarToPlot <- profileCGH$profileValues[,variable]
    

    profileCGH$profileValues$Chromosome <- ChrNumeric(profileCGH$profileValues$Chromosome)
    
    ChrNum <- TRUE


    indexna <- attr(na.omit(profileCGH$profileValues[variable]),"na.action")
    if(!is.null(indexna))
      {
        profileCGH$profileValues <- profileCGH$profileValues[-indexna,]
      }

    data("cytoband")

    if (!is.null(Chromosome))
      {
        ind <- NULL
        for (Chr in Chromosome)
          {
            indChr <- which(profileCGH$profileValues$Chromosome==Chr)
            ind <- c(ind,indChr)
          }
        profileCGH$profileValues <- profileCGH$profileValues[ind,]
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
      }
    

    
    LabelChr <- unique(na.omit(profileCGH$profileValues$Chromosome))
    NbChr <- length(LabelChr)
    LabelChr <- data.frame(Chromosome=LabelChr)

### Information dans les cytobandes
    genomeInfo <- aggregate(cytoband$End, list(Chromosome=cytoband$Chromosome, ChrNumeric=cytoband$ChrNumeric), max, na.rm=TRUE)
    names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length")
    genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome)
    genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric))

    if (ChrNum)
      {        
        LabelChr <- merge(LabelChr, genomeInfo[,c("ChrNumeric","Length")], by.x="Chromosome", by.y="ChrNumeric", all.x=TRUE)
        LabelChr <- LabelChr[order(LabelChr$Chromosome),]
      }
    else
      {
        LabelChr <- merge(LabelChr, genomeInfo, by="Chromosome", all.x=TRUE)
        LabelChr <- LabelChr[order(LabelChr$ChrNumeric),]
      }
   
  

    if (NbChr > 1)
      {

        gap <- 100000000/3
        LabelChr$Length <- LabelChr$Length + gap

        cumulLength <- cumsum(LabelChr$Length)
        LabelChr$Length <- c(0, cumulLength[1:(NbChr-1)])                                                             
        LabelChr$Length <- LabelChr$Length/(10^unit)
       
      }

    else
      {
        LabelChr$Length <- 0
      }


    if (ChrNum)
      {
        cytobandNew <- subset(cytoband, select=setdiff(names(cytoband),"Chromosome"))
        cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome", by.y="ChrNumeric")
      }
    else
      {
        cytobandNew <- subset(cytoband, select=setdiff(names(cytoband),"ChrNumeric"))
        cytobandNew <- merge(LabelChr, cytobandNew, by="Chromosome") 
      }


    cytobandNew$Start <- cytobandNew$Start/(10^unit)
    cytobandNew$End <- cytobandNew$End/(10^unit)


    cytobandNew$Start <- cytobandNew$Start + cytobandNew$Length
    cytobandNew$End <- cytobandNew$End +  cytobandNew$Length


    profileCGH$profileValues <- merge(profileCGH$profileValues, LabelChr, by="Chromosome")

    profileCGH$profileValues$NewPosBase <- profileCGH$profileValues$PosBase + profileCGH$profileValues$Length

    def.par <- par(no.readonly = TRUE)

 ##    if (plotband)
##       {
##         layout(c(1,2), heights=c(4,1))
##         par(mar=c(0,4,5,2))
##       }


    if (plotband)
      {
### Cytobandes

        layout(c(1,2), heights=c(1,4))

        par(mar=c(0,4,4,2))

        plot(0, type="n", xlim=c(0, max(cytobandNew$End)),
        ylim=c(-1.5,1.5), xaxt="n", yaxt="n", ylab="", xlab="")
        
        LabelChrCyto <- unique(cytobandNew$Chromosome)

        
        for (i in 1:NbChr)
          {
            plotCytoBand(cytobandNew, Chromosome=LabelChrCyto[i], labels=labels, y=0, height=2, colCytoBand=colCytoBand, colCentro=colCentro)
          }


        par(mar=c(4,4,0,2))

      }
    
    ### graphe des segments


    if (!is.null(Smoothing))
      {
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$Chromosome,profileCGH$profileValues$PosBase),]
        NbPos <- length(profileCGH$profileValues[,1])
        PosMax <- max(profileCGH$profileValues$NewPosBase) + 1
        Pos <- profileCGH$profileValues$NewPosBase[1:(NbPos-1)]
        PosNext <- profileCGH$profileValues$NewPosBase[2:NbPos]
        InterPos <- Pos + (PosNext-Pos)/2
        InterPos <- c(0, InterPos, PosMax)

        SmtStart <- profileCGH$profileValues[,Smoothing][1]
        SmtEnd <- profileCGH$profileValues[,Smoothing][NbPos]

        Smt1 <- profileCGH$profileValues[,Smoothing][1:(NbPos-1)]
        Smt1 <- c(SmtStart, Smt1, SmtEnd)

        Smt2 <- profileCGH$profileValues[,Smoothing][2:NbPos]
        Smt2 <- c(SmtStart, Smt2, SmtEnd)
        

        datasmt <- data.frame(PosBase=c(InterPos,InterPos),Smoothing=c(Smt1,Smt2))
        datasmt <- unique(datasmt)
        datasmt <- datasmt[order(datasmt$PosBase),]
      }


    
    if (length(intersect(names(profileCGH$profileValues),GNL))>=1)
      {

        col <- rep(colDAGLAD[5],length(profileCGH$profileValues$PosOrder))
        col[which(profileCGH$profileValues[GNL]==-1)] <- colDAGLAD[4]
        col[which(profileCGH$profileValues[GNL]==1)] <- colDAGLAD[3]
        col[which(profileCGH$profileValues[GNL]==2)] <- colDAGLAD[2]
        col[which(profileCGH$profileValues[GNL]==-10)] <- colDAGLAD[1]


        outliers <- rep(pchSymbol[1],length(profileCGH$profileValues$PosOrder))
        outliers[which(profileCGH$profileValues$OutliersTot!=0)] <- pchSymbol[2]

        if (plotband)
          {
            
            plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
                 pch=outliers, col=col, xaxt="n", xlab=main,
                 ylab=variable, ...)
          }
        else
          {
            plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
                 pch=outliers, col=col, xaxt="n", xlab="", ylab=variable, main=main)
          }


        if (!is.null(Smoothing))
          {
            lines(datasmt$Smoothing ~ datasmt$PosBase, col="black")
          }
        
        if (Bkp)
          {
            if (is.data.frame(profileCGH$BkpInfo))
              {
                profileCGH$BkpInfo <- merge(profileCGH$BkpInfo, LabelChr, by="Chromosome")
                profileCGH$BkpInfo$NewPosBase <- profileCGH$BkpInfo$PosBase + profileCGH$BkpInfo$Length
                abline(v=profileCGH$BkpInfo$NewPosBase+0.5, col="red", lty=2)
              }
            

          }
      }
    else
      {
        if (plotband)
          {
            plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
                 pch=20, xaxt="n", xlab=main, ylab=variable, ...)
          }
        else
          {
            plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
                 pch=20, xaxt="n", xlab="", ylab=variable, main=main, ...)
          }
        
        if (Bkp)          
          {
            if (is.data.frame(profileCGH$BkpInfo))
              {
                profileCGH$BkpInfo <- merge(profileCGH$BkpInfo, LabelChr, by="Chromosome")
                profileCGH$BkpInfo$NewPosBase <- profileCGH$BkpInfo$PosBase + profileCGH$BkpInfo$Length
                abline(v=profileCGH$BkpInfo$NewPosBase+0.5, col="red", lty=2)
              }

          }

        if (!is.null(Smoothing))
          {
            lines(datasmt$Smoothing ~ datasmt$PosBase, col="red")
          }
      }


    if (!is.null(text))
      {
        text(text$x, text$y, labels=text$labels, cex=text$cex)
      }




#    par(def.par)

  }
