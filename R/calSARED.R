#####
calSARED<-
function(Data, model="gok", origin=FALSE, 
         ErrorMethod="mc", nsim=1000, weight=TRUE, 
         trace=TRUE, outpdf=NULL, outfile=NULL) {
    UseMethod("calSARED")
} #
### 2016.07.16.
calSARED.default<-
function(Data, model="gok", origin=FALSE,   
         ErrorMethod="mc", nsim=1000, weight=TRUE, 
         trace=TRUE, outpdf=NULL, outfile=NULL) {
    ### Stop if not.
    stopifnot(ncol(Data)==5L, nrow(Data)>=5L,
              is.numeric(Data[,1L,drop=TRUE]), all(abs(Data[,1L]-round(Data[,1L]))<.Machine$double.eps^0.5),
              is.numeric(Data[,3L,drop=TRUE]), is.numeric(Data[,4L,drop=TRUE]), is.numeric(Data[,5L,drop=TRUE]),
              all(Data[,3L,drop=TRUE]>=0), all(Data[,5L,drop=TRUE]>0),
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(ErrorMethod)==1L, is.character(ErrorMethod), ErrorMethod %in% c("sp","mc"),
              is.numeric(nsim), length(nsim)==1L, nsim>=50L, nsim<=3000L,
              length(weight)==1L, is.logical(weight),
              length(trace)==1L, is.logical(trace),
              is.null(outpdf) || (is.character(outpdf) && length(outpdf)==1L),
              is.null(outfile) || (is.character(outfile) && length(outfile)==1L))
    ###
    colnames(Data) <- c("Grain.NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    ###
    ### Check Grain.NO and SAR.Cycle for Data.
    GrainNumber <- as.numeric(levels(factor(Data[,"Grain.NO",drop=TRUE])))
    n <- length(GrainNumber)
    nag <- n
    ###
    for (i in seq(n)) {
        GrainIndex <- which(Data[,"Grain.NO",drop=TRUE]==GrainNumber[i])
        SarCyclei <- substr(Data[GrainIndex,"SAR.Cycle",drop=TRUE], start=1L, stop=1L)   
        ###
        if (!all(diff(GrainIndex)==1L)) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data appears in discontinuous locations!", sep=""))
        } ### end if. 
        ###  
        if (!all(SarCyclei %in% c("N","R"))) {
             stop(paste("Error: Grain.NO", GrainNumber[i], 
                  " of Data contains incorrect SAR.Cycle!", sep=""))
        } # end if.
        ###
        if (all(SarCyclei=="N")) {
           stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data should contain SAR.Cycle of R!", sep=""))
        } # end if. 
        ###
        if (all(SarCyclei=="R")) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data should contain SAR.Cycle of N!", sep=""))
        } # end if.
        ###
        if (sum(SarCyclei=="N")>1L) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data should contain only one SAR.Cycle of N!", sep=""))
        } # end if.
    } # end for.
    ###
    ###
    DataList <- vector(mode="list", length=n)
    for (i in seq(n)) {
        GrainIndex <- which(Data[,"Grain.NO",drop=TRUE]==GrainNumber[i])  
        DataList[[i]] <- Data[GrainIndex,,drop=FALSE]
    } # end for.
    ###
    ###
    if (trace==TRUE) {
        print("Analyzing growth curves for all Grain.NO!")
    } # end if.
    ###
    failFitGrainNumber <- c()
    saturateGrainNumber <- c()
    for (i in seq(n)) {
        Data4L <- DataList[[i]][,c("SAR.Cycle","Dose","Signal","Signal.Err"),drop=FALSE]
        DataN <- as.numeric(Data4L[Data4L[,"SAR.Cycle",drop=TRUE]=="N",-1L])
        DataR <- Data4L[Data4L[,"SAR.Cycle",drop=TRUE]!="N",-1L,drop=FALSE]
        ###
        res1 <- try(fitGrowth(DataR,model=model,origin=origin,
                    weight=weight,plot=FALSE), silent=TRUE)
        ###
        if (class(res1)=="try-error") {
            failFitGrainNumber <- c(failFitGrainNumber, GrainNumber[i])
        } else {
            ###
            ### Test saturated natural signals.
            pars <- res1$LMpars[,1L,drop=TRUE]
            cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
            if (model %in% c("exp","lexp","dexp","gok")) {
                if (model=="exp") {
                    aa <- pars[1L]
                    bb <- pars[2L]
                    cc <- cst
                    Xm <- -log(1.0e-4/aa/bb)/bb
                    Ym <- aa*(1.0-exp(-bb*Xm)) + cc
                } else if (model=="lexp") {
                    aa <- pars[1L]
                    bb <- pars[2L]
                    cc <- pars[3L]
                    dd <- cst
                    Xm <- ifelse(cc>=1.0e-4, 1e5, -log((1.0e-4-cc)/aa/bb)/bb)
                    Ym <- aa*(1.0-exp(-bb*Xm)) + cc*Xm + dd  
                } else if (model=="dexp") {
                    aa <- ifelse(pars[2L]>pars[4L], pars[3L], pars[1L])
                    bb <- ifelse(pars[2L]>pars[4L], pars[4L], pars[2L])
                    cc <- ifelse(pars[2L]>pars[4L], pars[1L], pars[3L])
                    dd <- ifelse(pars[2L]>pars[4L], pars[2L], pars[4L])
                    ee <- cst
                    Xm <- -log(1.0e-4/aa/bb)/bb
                    Ym <- aa*(1.0-exp(-bb*Xm)) + cc*(1.0-exp(-dd*Xm)) + ee
                } else if (model=="gok") {
                    aa <- pars[1L]
                    bb <- pars[2L]
                    cc <- pars[3L]
                    dd <- cst
                    Xm <- (1.0e-4)^(-cc/(1.0+cc))*(1.0-(1.0e-4)^(cc/(1.0+cc))*
                          (1.0/aa/bb)^(cc/(1.0+cc)))*(1.0/aa/bb)^(-cc/(1.0+cc))/bb/cc
                    Ym <- aa*(1.0-(1.0+bb*cc*Xm)^(-1.0/cc)) + dd
                } # end if.
                ###
                maxSig <- DataN[2L] + DataN[3L]
                if (maxSig>=Ym) {
                    saturateGrainNumber <- c(saturateGrainNumber, GrainNumber[i])
                } # end if.
            } # end if.
            ###
        } # end if.
    } # end for.
    ###
    ###
    nfailFitGrainNumber <- length(failFitGrainNumber)
    if (nfailFitGrainNumber>=1L) {
        if (nfailFitGrainNumber==n) stop("Error: fail in fitting all growth curves!")
        Data <- Data[!Data[,"Grain.NO",drop=TRUE] %in% failFitGrainNumber,,drop=FALSE]
        GrainNumber <- as.numeric(levels(factor(Data[,"Grain.NO",drop=TRUE])))
        n <- length(GrainNumber)
    } # end if. 
    ###
    ###
    nsaturateGrainNumber <- length(saturateGrainNumber)
    if (nsaturateGrainNumber>=1L) {
        if (nsaturateGrainNumber==n) stop("Error: all natural OSL signals are saturated!")
        Data <- Data[!Data[,"Grain.NO",drop=TRUE] %in% saturateGrainNumber,,drop=FALSE]
        GrainNumber <- as.numeric(levels(factor(Data[,"Grain.NO",drop=TRUE])))
        n <- length(GrainNumber) 
    } # end if.
    ###
    ###
    DataList <- vector(mode="list", length=n)
    for (i in seq(n)) {
        GrainIndex <- which(Data[,"Grain.NO",drop=TRUE]==GrainNumber[i])  
        DataList[[i]] <- Data[GrainIndex,,drop=FALSE]
    } # end for.
    ###
    ###   
    Plot1 <- function(xvalue,yvalue,simED,mainName,dose,
                      doseltx,sdoseltx,pars,model,origin) {
        ###
        lowerX <- min(min(dose,xvalue),0)*1.2
        upperX <- max(dose,xvalue)*1.2
        lowerY <- min(min(doseltx,yvalue),0)*1.2
        upperY <- max(doseltx,yvalue)*1.2
        ###
        par(mar=c(5,5,4,1)+0.1)
        plot(NA, NA, main=mainName, xlab="Dose (Gy)", ylab="Standardised OSL",
             las=0, cex.lab=1.5, cex.main=1.5, xlim=c(lowerX,upperX),
             ylim=c(lowerY,upperY), xaxs="i", yaxs="i", lab=c(7,7,9))
        if (!is.null(simED) && all(yvalue>0)) {
            dmcED <- density(simED)
            dxy <- cbind(dmcED$x,dmcED$y)
            dxy[,2L] <- (dxy[,2L,drop=TRUE]-min(dxy[,2L,drop=TRUE]))/
                        (max(dxy[,2L,drop=TRUE])-min(dxy[,2L,drop=TRUE]))*yvalue[1L]*0.8
            polygon(dxy, col="grey")
            #rug(simED, quiet=TRUE)
        } # end if
        ###
        points(dose, doseltx, pch=21, cex=1.5, bg="black")
        ###
        x <- NULL
        cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
        ###
        if(model=="line") {
            curve(pars[1L]*x+cst, type="l", add=TRUE, 
                  lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="exp") {
           curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, type="l", add=TRUE, 
                 lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="lexp")  {
           curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, type="l", 
                 add=TRUE, lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="dexp") {
           curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                 type="l", add=TRUE, lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="gok") {
            curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                 type="l", add=TRUE, lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } # end if.
        ###
        ###
        points(x=xvalue[1L], y=yvalue[1L], pch=23, cex=1.5, bg="grey")
        ###
        arrowIndex <- which(sdoseltx>0.05 & sdoseltx/doseltx>0.01)
        if (length(arrowIndex)>=1L) {
            arrows(x0=dose[arrowIndex], y0=doseltx[arrowIndex]-sdoseltx[arrowIndex]/2.0, 
                   x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sdoseltx[arrowIndex]/2.0,
                   code=3, lwd=1, angle=90, length=0.05, col="black")
        } # end if.
        ###
        if (xvalue[2L]>0.05 && xvalue[2L]/xvalue[1L]>0.01) {
            arrows(x0=xvalue[1L]-xvalue[2L]/2.0, y0=yvalue[1L],
                   x1=xvalue[1L]+xvalue[2L]/2.0, y1=yvalue[1L],
                   code=3, lwd=1, angle=90, length=0.05, col="black")
        } # end if.
        ###
        lines(x=c(0,xvalue[1L],xvalue[1L]), y=c(yvalue[1L],yvalue[1L],0), lty="dashed", lwd=1)    
        legend("topleft", legend=paste("ED=",round(xvalue[1L],2L), "+/-",round(xvalue[2L],2L)," (Gy)",sep=""), 
               yjust=2, ncol=1L, cex=1.5, bty="n")
        grid()
        box(lwd=1)
        ###
        par(mar=c(5,4,4,2)+0.1)
    } # end function Plot1.
    ###
    ###
    tab <- data.frame(matrix(nrow=n, ncol=12L))
    rownames(tab) <- paste("Grain.NO",GrainNumber,sep="")
    colnames(tab) <- c("ED","Std.ED", "rsdED",
                       "RecyclingRatio","Std.RecyclingRatio",
                       "Recuperation1","Std.Recuperation1",
                       "Recuperation2","Std.Recuperation2",
                       "Method","FOM","RCS")
    ###
    if (!is.null(outpdf)) pdf(file=paste(outpdf,".pdf",sep=""))
    ###
    ###
    failEDIndex <- c()
    LMpars <- list()
    for (i in seq(n)) {
        Data4L <- DataList[[i]][,c("SAR.Cycle","Dose","Signal","Signal.Err"),drop=FALSE]
        DataN <- as.numeric(Data4L[Data4L[,"SAR.Cycle",drop=TRUE]=="N",-1L])
        DataR <- Data4L[Data4L[,"SAR.Cycle",drop=TRUE]!="N",-1L,drop=FALSE]
        ###
        if (trace==TRUE) {
            print(paste("Calculating SAR ED for Grain.NO ",GrainNumber[i],sep=""))
        } # end if.
        ###
        res2 <- try(calED(Curvedata=DataR, Ltx=DataN[-1L], model=model, origin=origin, 
                    ErrorMethod=ErrorMethod, nsim=nsim, weight=weight, plot=FALSE), silent=TRUE)
        ###
        ### Check error.
        if (class(res2)=="try-error") {
            tab[i,] <- NA
            failEDIndex <- c(failEDIndex, i)
            print(attr(res2,"condition"))
        } else {
            characterGrain.NO <- paste("Grain.NO",GrainNumber[i],sep="")
            LMpars[[characterGrain.NO]] <- res2$LMpars
            ###
            ### Calculated ED values.
            tab[i,c(1L,2L)] <- res2$ED 
            ###
            ###
            ### Calculate relative standard error of ED values.
            tab[i,3L] <- (res2$ED[2L]/res2$ED[1L])*100.0
            ###
            ###
            ### Calculate recycling ratio.
            dose <- DataR[,"Dose",drop=TRUE]
            doseltx <- DataR[,"Signal",drop=TRUE]
            sdoseltx <- DataR[,"Signal.Err",drop=TRUE]
            ###
            lvl.dose <- as.numeric(levels(factor(dose)))
            existrpd <- length(dose)>length(lvl.dose)
            if (existrpd==TRUE) {
                RepeatIndex <- apply(as.matrix(lvl.dose), MARGIN=1L, function(x,y)
                                     which(abs(x-y)<=.Machine$double.eps^0.5), dose)
                RepeatIndex <- unlist(RepeatIndex[sapply(RepeatIndex,length)>=2L])
                RecycleRatio <- doseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[1L]]
                sRecycleRatio <- abs(RecycleRatio)*sqrt((sdoseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[2L]])^2L+
                    (sdoseltx[RepeatIndex[1L]]/doseltx[RepeatIndex[1L]])^2L)
                tab[i,c(4L,5L)] <- c(RecycleRatio, sRecycleRatio)
            } else {
                tab[i,c(4L,5L)] <- NA
            } # end if.
            ###
            ###
            ### Calculate recuperation.
            exist0d <- which(abs(dose)<=.Machine$double.eps^0.5)
            if (length(exist0d)>=1L) {
                Recuperation1 <- (doseltx[exist0d[1L]]/DataN[2L])*100.0
                sRecuperation1 <- abs(Recuperation1)*sqrt((sdoseltx[exist0d[1L]]/doseltx[exist0d[1L]])^2L+
                    (DataN[3L]/DataN[2L])^2L)
                Recuperation2 <- (doseltx[exist0d[1L]]/max(doseltx))*100.0
                sRecuperation2 <- abs(Recuperation2)*sqrt((sdoseltx[exist0d[1L]]/doseltx[exist0d[1L]])^2L+
                    (sdoseltx[which.max(doseltx)]/max(doseltx))^2L)
                tab[i,c(6L,7L)] <- c(Recuperation1, sRecuperation1)
                tab[i,c(8L,9L)] <- c(Recuperation2, sRecuperation2)
            } else {
                tab[i,c(6L,7L)] <- NA
                tab[i,c(8L,9L)] <- NA
            } # end if.  
            ###
            ###
            ### Test whether equivalent dose is calculated 
            ### by Interpolating or Extrapolating.
            pars <- res2$LMpars[,1L,drop=TRUE]
            maxDose <- max(dose)
            cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
            if (model=="line") {
                maxLtx <- pars[1L]*maxDose+cst
            } else if (model=="exp") {
                maxLtx <- pars[1L]*(1.0-exp(-pars[2L]*maxDose))+cst
            } else if (model=="lexp") {
                maxLtx <- pars[1L]*(1.0-exp(-pars[2L]*maxDose))+pars[3L]*maxDose+cst
            } else if (model=="dexp") {
                maxLtx <- pars[1L]*(1.0-exp(-pars[2L]*maxDose))+pars[3L]*(1.0-exp(-pars[4L]*maxDose))+cst
            } else if (model=="gok") {
                maxLtx <- pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*maxDose)^(-1.0/pars[3L]))+cst
            } # end if.
            tab[i,10L] <- ifelse(DataN[2L]>maxLtx, "Extra", "Inter")
            ###
            ###
            ### Calculate FOM and RCS.
            tab[i,11L] <- 100*sum(abs(res2$fit.value[,2L,drop=TRUE]-
                res2$fit.value[,3L,drop=TRUE]))/sum(res2$fit.value[,3L,drop=TRUE])
            tab[i,12L] <- res2$value/(nrow(DataR)-nrow(res2$LMpars))
            ###
            ###
            ### Plot the results of ED calculation.
            if (!is.null(outpdf)) {
                xvalue <- res2$ED
                yvalue <- DataN[-1L]
                simED <- if (is.null(res2$mcED)) NULL else res2$mcED
                mainName <- paste("Grain.NO",GrainNumber[i],sep="")
                sdoseltx <- DataR[,"Signal.Err",drop=TRUE]
                ###
                Plot1(xvalue, yvalue, simED, mainName, dose,
                      doseltx, sdoseltx, pars, model, origin)
            } # end if.
        } # end if.    
    } # end for. 
    ###
    ###
    if (!is.null(outpdf)) dev.off()
    ###
    if (!is.null(outfile)) {
        write.csv(tab, file=paste(outfile,".csv",sep=""))
    } # end if.
    ###
    ###
    failEDGrainNumber <- NULL
    if (length(failEDIndex)>=1L) {
        if (length(failEDIndex)==n) stop("Error: fail in equivalent dose calculation for all Grain.NO!")
        failEDGrainNumber <- GrainNumber[failEDIndex]
        GrainNumber <- GrainNumber[-failEDIndex]
        tab <- tab[-failEDIndex,,drop=FALSE]
    } # end if.
    ###
    extrapolateGrainNumber <- GrainNumber[tab[,"Method",drop=TRUE]=="Extra"]
    if (length(extrapolateGrainNumber)==0L) extrapolateGrainNumber <- NULL
    ###
    largeRecyclingRatioGrainNumber <- GrainNumber[complete.cases(tab) &
        (tab[,"RecyclingRatio",drop=TRUE]>1.1 | tab[,"RecyclingRatio",drop=TRUE]<0.9)]
    if (length(largeRecyclingRatioGrainNumber)==0L) largeRecyclingRatioGrainNumber <- NULL
    ###
    largeRecuperatuionGrainNumber1 <- GrainNumber[complete.cases(tab) & tab[,"Recuperation1",drop=TRUE]>5.0]
    if (length(largeRecuperatuionGrainNumber1)==0L) largeRecuperatuionGrainNumber1 <- NULL
    ###
    largeRecuperatuionGrainNumber2 <- GrainNumber[complete.cases(tab) & tab[,"Recuperation1",drop=TRUE]>10.0]
    if (length(largeRecuperatuionGrainNumber2)==0L) largeRecuperatuionGrainNumber2 <- NULL
    ###
    ###
    output <- list("LMpars"=LMpars,
                   "N"=nag,
                   "failFit.NO"=failFitGrainNumber,
                   "saturate.NO"=saturateGrainNumber,
                   "failED.NO"=failEDGrainNumber,
                   "extrapolate.NO"=extrapolateGrainNumber,
                   "largeRcy.NO"=largeRecyclingRatioGrainNumber,
                   "largeRcp5.NO"=largeRecuperatuionGrainNumber1,
                   "largeRcp10.NO"=largeRecuperatuionGrainNumber2,
                   "tab"=tab)
    class(output) <- "calSARED"   
    ###            
    return(invisible(output))
} # end function calSARED.default. 
###
###
summary.calSARED <- function(object,...) {
    UseMethod("summary.calSARED")
} #
### 2016.06.28.
summary.calSARED.default <- function(object,...) {
    stopifnot(class(object)=="calSARED", is.list(object), length(object)==10L,
              all(names(object) %in% c("LMpars","N","failFit.NO","saturate.NO","failED.NO",
              "extrapolate.NO","largeRcy.NO","largeRcp5.NO","largeRcp10.NO","tab")))
    ###
    cat("\n")
    cat("=====================================Summary=====================================\n")
    cat(paste("Number of analyzed grains: N=", object$N, "\n\n", sep=""))
    cat(paste("Number of rejected grains: N=", length(object$failFit.NO)+length(object$saturate.NO)+length(object$failED.NO), "\n\n", sep=""))
    ###
    cat(paste("Grain.NO failed in DRC fitting (rejected N=",length(object$failFit.NO),")",":\n", sep=""))
    print(object$failFit.NO)
    cat("\n")
    ###
    cat(paste("Grain.NO saturated (rejected N=",length(object$saturate.NO),")",":\n", sep=""))
    print(object$saturate.NO)
    cat("\n")
    ###
    cat(paste("Grain.NO failed in ED (error) calculation (rejected N=",length(object$failED.NO),")",":\n", sep=""))
    print(object$failED.NO)
    cat("\n")
    ###
    cat(paste("Grain.NO calculated by extrapolation (N=",length(object$extrapolate.NO),")",":\n", sep=""))
    print(object$extrapolate.NO)
    cat("\n")
    ### 
    cat(paste("Grain.NO with large recylcing ratio [outside the range 0.9-1.1] (N=",length(object$largeRcy.NO),")",":\n", sep=""))
    print(object$largeRcy.NO)
    cat("\n")
    ###
    cat(paste("Grain.NO with recuperation1 > 5% (N=",length(object$largeRcp5.NO),")",":\n", sep=""))
    print(object$largeRcp5.NO)
    cat("\n")
    ###   
    cat(paste("Grain.NO with recuperation1 > 10% (N=",length(object$largeRcp10.NO),")",":\n", sep=""))
    print(object$largeRcp10.NO)
    cat("=======================================End=======================================\n\n")
    ###
} # end function summary.calSARED.default.
##### 
