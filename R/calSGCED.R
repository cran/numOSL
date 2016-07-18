#####
calSGCED<-
function(Data1,pars,model,origin,
         ErrorMethod="sp",outpdf=NULL) {
    UseMethod("calSGCED")
} ###
### 2016.07.17.
calSGCED.default<-
function(Data1,pars,model,origin,
         ErrorMethod="sp",outpdf=NULL) {
    ### Stop if not.
    stopifnot(ncol(Data1)==5L, nrow(Data1)>=5L,
              is.numeric(Data1[,1L,drop=TRUE]), all(abs(Data1[,1L]-round(Data1[,1L]))<.Machine$double.eps^0.5),
              is.numeric(Data1[,3L,drop=TRUE]), is.numeric(Data1[,4L,drop=TRUE]), is.numeric(Data1[,5L,drop=TRUE]),
              all(Data1[,3L,drop=TRUE]>=0), all(Data1[,5L,drop=TRUE]>0),
              is.vector(pars), is.numeric(pars), 
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(ErrorMethod)==1L, is.character(ErrorMethod), ErrorMethod=="sp",
              is.null(outpdf) || (is.character(outpdf) && length(outpdf)==1L))
    ###
    ### Check argument pars.
    if (model=="line") {
        if (origin==TRUE && length(pars)!=1L) stop("Error: need provide one parameter!")
        if (origin==FALSE && length(pars)!=2L) stop("Error: need provide two parameter!")
        if (pars[1L]<=0) stop("Error: improper parameters!")
    } else if (model=="exp") {
        if (origin==TRUE && length(pars)!=2L) stop("Error: need provide two parameter!")
        if (origin==FALSE && length(pars)!=3L) stop("Error: need provide three parameter!")
        if (any(pars[seq(2L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="lexp") {
        if (origin==TRUE && length(pars)!=3L) stop("Error: need provide three parameter!")
        if (origin==FALSE && length(pars)!=4L) stop("Error: need provide four parameter!")
        if (any(pars[seq(3L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="dexp") {
        if (origin==TRUE && length(pars)!=4L) stop("Error: need provide four parameter!")
        if (origin==FALSE && length(pars)!=5L) stop("Error: need provide five parameter!")
        if (any(pars[seq(4L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="gok") {
        if (origin==TRUE && length(pars)!=3L) stop("Error: need provide three parameter!")
        if (origin==FALSE && length(pars)!=4L) stop("Error: need provide four parameter!")
        if (any(pars[seq(3L)]<=0)) stop("Error: improper parameters!")
    } # end if.       
    ###
    ###
    ### Check Grain.NO and SAR.Cycle for Data1.
    colnames(Data1) <- c("Grain.NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    GrainNumber <- as.numeric(levels(factor(Data1[,"Grain.NO",drop=TRUE])))
    n <- length(GrainNumber)
    ###
    for (i in seq(n)) {
        GrainIndex <- which(Data1[,"Grain.NO",drop=TRUE]==GrainNumber[i])
        SarCyclei <- substr(Data1[GrainIndex,"SAR.Cycle",drop=TRUE], start=1L, stop=1L)
        ###
        if (!all(diff(GrainIndex)==1L)) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data1 appears in discontinuous locations!", sep=""))
        } ### end if.
        ###
        if (length(SarCyclei)!=2L) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                 " of Data1 should contain two SAR.Cycle!", sep=""))
        } # end if.
        ###
        if (!("N" %in% SarCyclei && "R" %in% SarCyclei)) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                 " of Data1 contains incorrect SAR.Cycle!", sep=""))
        } # end if.
    } # end for.
    ###
    ###
    cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
    fcn <- function(x) {
        if(model=="line") {
           pars[1L]*x+cst
        } else if(model=="exp") {
           pars[1L]*(1.0-exp(-pars[2L]*x))+cst
        } else if(model=="lexp")  {
           pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst
        } else if(model=="dexp") {
           pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst
        } else if(model=="gok") {
           pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst
        } # end if.
    } # end function fcn.
    ###
    ###
    Mat <- matrix(nrow=n, ncol=4L)
    ScaledNaturalSignal <- matrix(nrow=n, ncol=2L)
    ###
    for(i in seq(n)) {
        Mat[i,c(1L,2L)] <- as.numeric(Data1[Data1[,"Grain.NO",drop=TRUE]==GrainNumber[i] &
                           substr(Data1[,"SAR.Cycle",drop=TRUE],start=1L,stop=1L)=="R",c("Dose","Signal")])
        Mat[i,c(3L,4L)] <- as.numeric(Data1[Data1[,"Grain.NO",drop=TRUE]==GrainNumber[i] &
                           substr(Data1[,"SAR.Cycle",drop=TRUE],start=1L,stop=1L)=="N",c("Signal","Signal.Err")])
        ###
        ScaledNaturalSignal[i,] <- (fcn(Mat[i,1L])/Mat[i,2L])*Mat[i,c(3L,4L),drop=TRUE]
    } # end for.
    ###
    ###
    ### Test satureated natural signals.
    saturateGrainNumber <- NULL
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
        saturateIndex <- which(ScaledNaturalSignal[,1L,drop=TRUE]+ScaledNaturalSignal[,2L,drop=TRUE]>=Ym)
        if (length(saturateIndex)>=1L) {
            if (length(saturateIndex)==n) stop("Error: all natural OSL signals are saturated!")
            saturateGrainNumber <- GrainNumber[saturateIndex]
            ScaledNaturalSignal <- ScaledNaturalSignal[-saturateIndex,,drop=FALSE]
            GrainNumber <- GrainNumber[-saturateIndex]
            n <- length(GrainNumber)
        } # end if.     
    } # end if.
    ###
    ###
    if (model=="line") {
        n2 <- 1L+!origin
        mdl <- 0L
    } else if (model=="exp") {
        n2 <- 2L+!origin
        mdl <- 1L
    } else if (model=="lexp") {
        n2 <- 3L+!origin
        mdl <- 2L
    } else if (model=="dexp") {
        n2 <- 4L+!origin
        mdl <- 3L
    } else if (model=="gok") {
        n2 <- 3L+!origin
        mdl <- 7L
    } # end if.
    outDose <- matrix(0, nrow=n, ncol=2L)
    message <- whichErr <- 0
    ###
    res <- .Fortran("calSGCED",as.integer(n),as.integer(n2),as.double(ScaledNaturalSignal),
                    as.double(pars),outDose=as.double(outDose),as.integer(mdl),
                    message=as.integer(message),whichErr=as.integer(whichErr),
                    PACKAGE="numOSL")
    ###
    if (res$message==1) {
        stop(paste("Error: natural OSL signal saturated for Grain.NO", GrainNumber[res$whichErr], sep=""))
    } else if (res$message==2) {
        stop(paste("Error: fail in equivalent dose calculation for Grain.NO", GrainNumber[res$whichErr], sep=""))
    } else if (res$message==3L) {
        stop(paste("Error: fail in equivalent dose error estimation by Simple Transformation for Grain.NO", GrainNumber[res$whichErr], sep=""))
    } # end if.
    ###
    ###
    rownames(ScaledNaturalSignal) <- paste("Grain.NO",GrainNumber,sep="")
    colnames(ScaledNaturalSignal) <- c("Ltx", "Std.Ltx")
    ###
    sgcED <- matrix(res$outDose, ncol=2L)
    rownames(sgcED) <- paste("Grain.NO",GrainNumber,sep="")
    colnames(sgcED) <- c("ED", "Std.ED")
    ###
    output <- list("scaleLtx"=ScaledNaturalSignal,
                   "sgcED"=sgcED,
                   "saturate.NO"=saturateGrainNumber)
    
    Plot1 <- function(xvalue,yvalue,mainName,dose,pars,model,origin) {
        ###
        lowerX <- min(min(dose,xvalue),0)*1.2
        upperX <- max(dose,xvalue)*1.2
        ###
        par(mar=c(5,5,4,1)+0.1)
        x<-NULL
        cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
        ###
        if(model=="line") {
            curve(pars[1L]*x+cst, type="l", main=mainName, xlab="Dose", ylab="Normalised standardised OSL", 
                  las=0, cex.lab=1.5, cex.main=1.5, xaxs="i", yaxs="i", lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="exp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, type="l", main=mainName, xlab="Dose", ylab="Normalised standardised OSL", 
                  las=0, cex.lab=1.5, cex.main=1.5, xaxs="i", yaxs="i", lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="lexp")  {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, type="l", main=mainName, xlab="Dose", ylab="Normalised standardised OSL", 
                  las=0, cex.lab=1.5, cex.main=1.5, xaxs="i", yaxs="i", lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="dexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, type="l", main=mainName, xlab="Dose", ylab="Normalised standardised OSL", 
                  las=0, cex.lab=1.5, cex.main=1.5, xaxs="i", yaxs="i", lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } else if(model=="gok") {
            curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, type="l", main=mainName, xlab="Dose", ylab="Normalised standardised OSL", 
                  las=0, cex.lab=1.5, cex.main=1.5, xaxs="i", yaxs="i", lwd=1.5, from=lowerX, to=upperX, col="skyblue")
        } # end if.
        ###
        points(x=xvalue[1L], y=yvalue[1L], pch=23, cex=1.5, bg="grey")
        ###
        if (xvalue[2L]>0.05 && xvalue[2L]/xvalue[1L]>0.01) {
            arrows(x0=xvalue[1L]-xvalue[2L]/2.0, y0=yvalue[1L],
                   x1=xvalue[1L]+xvalue[2L]/2.0, y1=yvalue[1L],
                   code=3, lwd=1, angle=90, length=0.05, col="black")
        } # end if.
        ###
        lines(x=c(0, xvalue[1L], xvalue[1L]), y=c(yvalue[1L],yvalue[1L],0), lty="dashed", lwd=1)    
        legend("topleft", legend=paste("ED=",round(xvalue[1L],2L), "+/-",round(xvalue[2L],2L)," (Gy)",sep=""), 
               yjust=2, ncol=1L, cex=1.5, bty="n")
        grid()
        box(lwd=1)
        ###
        par(mar=c(5,4,4,2)+0.1)
    } # end function Plot1.
    ###
    ###
    ###
    if (!is.null(outpdf)) {
        pdf(file=paste(outpdf,".pdf",sep=""))
        ###
        dose <- Data1[Data1[,"Grain.NO",drop=TRUE] %in% GrainNumber,3L,drop=TRUE]
        ###
        for (i in seq(length(GrainNumber))) {
            xvalue <- sgcED[i,,drop=TRUE] 
            yvalue <- ScaledNaturalSignal[i,,drop=TRUE]
            mainName <- paste("Grain.NO", GrainNumber[i], sep="")
            Plot1(xvalue,yvalue,mainName,dose,pars,model,origin)
        } # end for.
        ###
        dev.off()
    } # end if.
    ### 
    return(output)
} # end function calSGCED.default.
#####
