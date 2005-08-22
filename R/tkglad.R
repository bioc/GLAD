tkglad <- function(...)
  {
    UseMethod("tkglad")
  }

tkglad.default <- function(list=character(0), ...)
  {

    require(tcltk)

    frameDaglad <- tktoplevel()
    frameOverall <- tkframe(frameDaglad)
    framePlotGlobal <- tkframe(frameDaglad)
    tkwm.title(frameDaglad,"Array CGH Analysis")
    tkpack(tklabel(frameDaglad, text="GLAD (Gain and Loss Analysis of DNA)"))
    img.path <- paste(system.file(package="GLAD"),"/doc/institut.gif",sep="")
    img.tk <- tclVar()
    tkimage.create("photo",img.tk,file=img.path)
    

    mediancenter.frm <- tkframe(frameOverall, borderwidth=2, relief="groove")
    tkpack(tklabel(mediancenter.frm,text="Median centering before processing data"), anchor="w")
    
    mediancenter <- tclVar(0)
    cb.mediancenter <- tkcheckbutton(mediancenter.frm, variable=mediancenter,text="mediancenter")
    tkpack(cb.mediancenter)
    tkpack(mediancenter.frm, fill="x")
   

    frameUpper <- tkframe(frameOverall, borderwidth=2, relief="groove")
    frameUpperLeft <- tkframe(frameUpper)
    frameUpperRight <- tkframe(frameUpper)
    smoothfunc.frm <- tkframe(frameUpperLeft, borderwidth=2, relief="groove")
    lkern.frm <- tkframe(frameUpperRight, borderwidth=2, relief="groove")
    
    
    tkpack(tklabel(frameUpper,text="Smoothing function parameters"), anchor="w")
    tkpack(tklabel(smoothfunc.frm,text="smoothfunc:"))
    tkpack(tklabel(lkern.frm,text="lkern:"))
       
    smoothfunc <- tclVar("lawsglad")
    for (i in c("lawsglad","laws","aws"))
      {
        tmp <- tkradiobutton(smoothfunc.frm,variable=smoothfunc,value=i,text=i)
        tkpack(tmp,anchor="w")
      }

    kernels <- c("Exponential","Uniform","Triangle","Quadratic","Cubic")
    lkern <- tclVar("Exponential")
    for (i in kernels)
      {
        tmp <- tkradiobutton(lkern.frm,variable=lkern,value=i,text=i)
        tkpack(tmp,anchor="w")
      }


    frameUpperLeft.left <- tkframe(frameUpperLeft)
    frameUpperLeft.right <- tkframe(frameUpperLeft)
    qlambda <- tclVar(0.999)
    l.qlambda <- tklabel(frameUpperLeft.left, text="qlambda:")
    e.qlambda <- tkentry(frameUpperLeft.right, width=10, textvariable=qlambda)

    bandwidth <- tclVar(10)
    l.bandwidth <- tklabel(frameUpperLeft.left, text="bandwidth:")
    e.bandwidth <- tkentry(frameUpperLeft.right, width=10, textvariable=bandwidth)

    round <- tclVar(1.5)
    l.round <- tklabel(frameUpperLeft.left, text="round:")
    e.round <- tkentry(frameUpperLeft.right, width=10, textvariable=round)

    tkpack(l.qlambda,anchor="e")
    tkpack(e.qlambda,anchor="w")
    tkpack(l.bandwidth,anchor="e")
    tkpack(e.bandwidth,anchor="e")
    tkpack(l.round,anchor="e")
    tkpack(e.round,anchor="e")
    tkpack(smoothfunc.frm)
    tkpack(frameUpperLeft.left,frameUpperLeft.right,side="left")
    tkpack(lkern.frm)
    tkpack(frameUpperLeft,frameUpperRight, side="left", anchor="n")
    tkpack(frameUpper, fill="x")



    lambdabreak.frm <- tkframe(frameOverall, borderwidth=2, relief="groove")
    tkpack(tklabel(lambdabreak.frm,text="Optimization of the breakpoints"), anchor="w")
 
    lambdabreak <- tclVar(8)
    l.lambdabreak <- tklabel(lambdabreak.frm, text="lambdabreak:")
    e.lambdabreak <- tkentry(lambdabreak.frm, width=10, textvariable=lambdabreak)

    param <- tclVar(6)
    l.param <- tklabel(lambdabreak.frm, text="param:")
    e.param <- tkentry(lambdabreak.frm, width=10, textvariable=param)
    tkpack(l.lambdabreak,e.lambdabreak,l.param,e.param, side="left")
    
    tkpack(lambdabreak.frm, fill="x")

    regionall.frm <- tkframe(frameOverall, borderwidth=2, relief="groove")
    tkpack(tklabel(regionall.frm, text="Region Assigment"), anchor="w")

    region.frm <- tkframe(regionall.frm)

    region.frm.left <- tkframe(region.frm)
    region.frm.right <- tkframe(region.frm)
    

    lambdacluster <- tclVar(8)
    l.lambdacluster <- tklabel(region.frm.left, text="lambdacluster:")
    e.lambdacluster <- tkentry(region.frm.right, width=10, textvariable=lambdacluster)
    tkpack(l.lambdacluster, anchor="e")
    tkpack(e.lambdacluster, anchor="w")

   
    lambdaclusterGen <- tclVar(40)
    l.lambdaclusterGen <- tklabel(region.frm.left, text="lambdaclusterGen:")
    e.lambdaclusterGen <- tkentry(region.frm.right, width=10, textvariable=lambdaclusterGen)
    tkpack(l.lambdaclusterGen, anchor="e")
    tkpack(e.lambdaclusterGen, anchor="w")

    
    tkpack(region.frm.left,region.frm.right, side="left")
    tkpack(regionall.frm, fill="x")
    tkpack(region.frm)

    
    outliers.frm <- tkframe(frameOverall, borderwidth=2, relief="groove")
    tkpack(tklabel(outliers.frm,text="Outliers Detection"), anchor="w")
    alpha <- tclVar(0.001)
    l.alpha <- tklabel(outliers.frm, text="alpha:")
    e.alpha <- tkentry(outliers.frm, width=10, textvariable=alpha)

    msize <- tclVar(5)
    l.msize <- tklabel(outliers.frm, text="msize:")
    e.msize <- tkentry(outliers.frm, width=10, textvariable=msize)

    tkpack(l.alpha,e.alpha,l.msize,e.msize, side="left")
    tkpack(outliers.frm, fill="x")


    regionall.frm <- tkframe(frameOverall, borderwidth=2, relief="groove")
    tkpack(tklabel(regionall.frm, text="Region Assigment"), anchor="w")

    region.frm <- tkframe(regionall.frm)

    region.frm.left <- tkframe(region.frm)
    region.frm.right <- tkframe(region.frm)
    

    lambdaclusterGen <- tclVar(40)
    l.lambdaclusterGen <- tklabel(region.frm.left, text="lambdaclusterGen:")
    e.lambdaclusterGen <- tkentry(region.frm.right, width=10, textvariable=lambdaclusterGen)
    tkpack(l.lambdaclusterGen, e.lambdaclusterGen)

 
  

    
    OnAnalysis <- function()
      {

        mediancenterVal <- FALSE
        if (tclvalue(mediancenter)!="0")
          mediancenterVal <- TRUE


        smoothfuncVal <- as.character(tclvalue(smoothfunc))
        lkernVal <- as.character(tclvalue(lkern))       
        qlambdaVal <- as.numeric(tclvalue(qlambda))
        bandwidthVal <- as.numeric(tclvalue(bandwidth))
        roundVal <- as.numeric(tclvalue(round))
        lambdabreakVal <- as.numeric(tclvalue(lambdabreak))       
        lambdaclusterVal <- as.numeric(tclvalue(lambdacluster))    
        lambdaclusterGenVal <- as.numeric(tclvalue(lambdaclusterGen))    
        paramVal <- as.numeric(tclvalue(param))
        alphaVal <- as.numeric(tclvalue(alpha))
        msizeVal <- as.numeric(tclvalue(msize))        

        for (array in list)
          {
            tmp <- get(array)
            
            tmp <- glad(tmp,
                        mediancenter=mediancenterVal,
                        smoothfunc=smoothfuncVal,
                        lkern=lkernVal,
                        qlambda=qlambdaVal,
                        bandwidth=bandwidthVal,
                        round=roundVal,
                        lambdabreak=lambdabreakVal,
                        lambdacluster=lambdaclusterVal,
                        lambdaclusterGen=lambdaclusterGenVal,
                        param=c(d=paramVal),
                        alpha=alphaVal,
                        msize=msizeVal)

            assign(array,tmp, envir=globalenv())

 
          }

        MakeListBox()
        tkmessageBox(title="GLAD Analysis",message="The analysis is finished",icon="info",type="ok")
      }
    
    OnQuit <- function()
      {
        tkdestroy(frameDaglad)
      }


    OnDefault <- function()
      {
        tclvalue(mediancenter) <- 0
        tclvalue(smoothfunc) <- "lawsglad"
        tclvalue(lkern) <- "Exponential"
        tclvalue(qlambda) <- 0.999
        tclvalue(bandwidth) <- 10
        tclvalue(lambdabreak) <- 8
        tclvalue(param) <- 6
        tclvalue(alpha) <- 0.001
        tclvalue(msize) <- 5
        tclvalue(lambdaclusterGen) <- 40
      }


    OnPlot <- function()
      {
        
        BkpVal <- FALSE
        if (tclvalue(Bkp)!="0")
          BkpVal <- TRUE


        LogRatioVal <- as.character(tkget(lb.LogRatio,tkcurselection(lb.LogRatio)))

        SmoothingVal <- as.character(tkget(lb.Smoothing,tkcurselection(lb.Smoothing)))

        if (SmoothingVal=="<None>")
          SmoothingVal <- NULL

        unitVal <- as.numeric(tclvalue(unit))

        plotbandVal <- FALSE
        if (tclvalue(plotband)!="0")
          plotbandVal <- TRUE

        ChrVal <- as.numeric(tclvalue(Chr))
        if (ChrVal==0)
          ChrVal <- NULL

        for (array in list)
          {
            tmp <- get(array)
            X11()
            plotProfile(tmp,
                        unit=unitVal,
                        variable=LogRatioVal,
                        Smoothing=SmoothingVal,
                        Bkp=BkpVal,
                        plotband=plotbandVal,
                        Chromosome=ChrVal)
          }


      }

    
    OnLogo <- function()
      {

        tkmessageBox(title="GLAD package",message="http://bioinfo.curie.fr",icon="info",type="ok")

      }

    control.frm <- tkframe(frameOverall)
    but.analysis <- tkbutton(control.frm,text="Analyze", command=OnAnalysis)
    but.quit <- tkbutton(control.frm,text="Quit", command=OnQuit)
    but.default <- tkbutton(control.frm,text="Default", command=OnDefault)

    ### Plot options
        

    framePlot <- tkframe(framePlotGlobal, borderwidth=2, relief="groove")
    tkpack(tklabel(framePlot, text="Plot options"))
    frame1 <- tkframe(framePlot)
    tkpack(frame1)
    plot.frm.left <- tkframe(framePlot)
    plot.frm.right <- tkframe(framePlot)
    
    Bkp <- tclVar(1)
    cb.Bkp <- tkcheckbutton(framePlot, variable=Bkp, text="Breakpoints")
    tkpack(cb.Bkp, anchor="w")


    unit <- tclVar(0)
    l.unit <- tklabel(plot.frm.left, text="unit:")
    e.unit <- tkentry(plot.frm.right, width=10, textvariable=unit)
    tkpack(l.unit, anchor="e")
    tkpack(e.unit, anchor="w")


    
    plotband <- tclVar(1)
    cb.plotband <- tkcheckbutton(framePlot, variable=plotband, text="ploband")
    tkpack(cb.plotband, anchor="w")

    Chr <- tclVar(0)
    l.Chr <- tklabel(plot.frm.left, text="Chr:")
    e.Chr <- tkentry(plot.frm.right, width=10, textvariable=Chr)
    tkpack(l.Chr, anchor="e")
    tkpack(e.Chr, anchor="w")

    tkpack(plot.frm.left, plot.frm.right, side="left")

    tkpack(framePlot)

    frame2 <- 
    l.LogRatio <- tklabel(frame1, text="variable:")
    variable.frm <- tkframe(frame1)  
    sb.LogRatio <- tkscrollbar(variable.frm, command=function(...)tkyview(lb.LogRatio,...))
    lb.LogRatio<-tklistbox(variable.frm,height=4,selectmode="single",background="white",
                           yscrollcommand=function(...)tkset(sb.LogRatio,...), exportselection=FALSE)

    l.Smoothing <- tklabel(frame1, text="Smoothing:")
    smoothing.frm <- tkframe(frame1)
    sb.Smoothing <- tkscrollbar(smoothing.frm, command=function(...)tkyview(lb.Smoothing,...))
    lb.Smoothing<-tklistbox(smoothing.frm,height=4,selectmode="single",background="white",
                           yscrollcommand=function(...)tkset(sb.Smoothing,...), exportselection=FALSE)
  




    MakeListBox <- function()
      {

        tkdelete(lb.LogRatio,0,tksize(lb.LogRatio))
        tkdelete(lb.Smoothing,0,tksize(lb.Smoothing))
        Champ <- names(get(list[1])$profileValues)
        tkinsert(lb.Smoothing,"end","<None>")
        for (i in Champ)
          {
            tkinsert(lb.LogRatio,"end",i)
            tkinsert(lb.Smoothing,"end",i)
          }

        ind.LogRatio <- which(Champ=="LogRatio")-1
        ind.Smoothing <- which(Champ=="Smoothing")
        
        if (length(ind.LogRatio)>0)
          tkselection.set(lb.LogRatio,ind.LogRatio)

        if (length(ind.Smoothing)>0)
          tkselection.set(lb.Smoothing,ind.Smoothing)
        else
          tkselection.set(lb.Smoothing,0)
      }

    MakeListBox()
       
    tkpack(l.LogRatio, anchor="w")
    tkpack(lb.LogRatio,sb.LogRatio, side="left", fill="y")
    tkpack(variable.frm)
    tkpack(l.Smoothing, anchor="w")
    tkpack(lb.Smoothing,sb.Smoothing, side="left", fill="y")
    tkpack(smoothing.frm)

   
    but.plot <- tkbutton(framePlotGlobal,text="Plot", command=OnPlot)
    tkpack(but.plot)

    tkpack(tkbutton(framePlotGlobal, image=img.tk, command=OnLogo), side="bottom")

    

    tkpack(but.analysis, but.default, but.quit, side="left")
    tkpack(control.frm)

    tkpack(frameOverall, framePlotGlobal, side="left", anchor="n")    


  }


