#####
calED <-
function(Curvedata, Ltx, model="exp", origin=FALSE, 
         ErrorMethod="mc", nsim=1000, weight=TRUE, 
         plot=TRUE) {
    UseMethod("calED")
} ###
### 2016.07.16.
calED.default <-
function(Curvedata, Ltx, model="exp", origin=FALSE, 
         ErrorMethod="mc", nsim=1000, weight=TRUE, 
         plot=TRUE) {
    ### Stop if not.
    stopifnot(ncol(Curvedata)==3L, all(Curvedata[,1L,drop=TRUE]>=0), all(Curvedata[,3L,drop=TRUE]>0),
              is.vector(Ltx), length(Ltx)==2L,
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(ErrorMethod)==1L, is.character(ErrorMethod), ErrorMethod %in% c("mc","sp"),
              length(nsim)==1L, is.numeric(nsim), nsim>=50L, nsim<=3000L,
              length(weight)==1L, is.logical(weight), 
              length(plot)==1L, is.logical(plot))
    ###
    dose <- as.numeric(Curvedata[,1L,drop=TRUE])
    doseltx <- as.numeric(Curvedata[,2L,drop=TRUE])
    sdoseltx <- as.numeric(Curvedata[,3L,drop=TRUE])
    ndat <- nrow(Curvedata)
    n2 <- if (model=="line") {
            1L+!origin
        } else if (model=="exp") {
            2L+!origin
        } else if (model=="lexp") {
            3L+!origin
        } else if (model=="dexp") {
            4L+!origin
        } else if (model=="gok") {
            3L+!origin
        } # end if. 
    ### 
    ### Check if data points is enough for fitting the model.
    if (ndat<n2) {
        stop("Error: data points is not enough for model optimization!")
    } # end if.
    ###
    inltx <- matrix(Ltx, nrow=1L) 
    ###
    outDose <- matrix(0, nrow=1L, ncol=2L)
    mcED <- vector(length=nsim)
    pars <- stdp <- vector(length=n2)
    model1 <- if (model=="line") {
        0L } else if (model=="exp") {
        1L } else if (model=="lexp") {
        2L } else if (model=="dexp") {
        3L } # end if.
    uw <- ifelse(weight==FALSE, 0L, 1L)
    method <- ifelse(ErrorMethod=="sp", 0L, 1L)
    fvec1 <- vector(length=ndat)
    fmin <- 0
    message <- 0
    ###
    if (model=="gok") {
        res <- .Fortran("calED1",as.double(dose),as.double(doseltx),as.double(sdoseltx),
                        as.integer(ndat),as.integer(n2),as.double(inltx),outDose=as.double(outDose),
                        mcED=as.double(mcED),pars=as.double(pars),stdp=as.double(stdp),as.integer(uw),
                        as.integer(method),as.integer(nsim),fvec1=as.double(fvec1),fmin=as.double(fmin),
                        message=as.integer(message),PACKAGE="numOSL")
    } else {
        res <- .Fortran("calED",as.double(dose),as.double(doseltx),as.double(sdoseltx),
                        as.integer(ndat),as.integer(n2),as.double(inltx),outDose=as.double(outDose),
                        mcED=as.double(mcED),pars=as.double(pars),stdp=as.double(stdp),as.integer(model1),
                        as.integer(uw),as.integer(method),as.integer(nsim),fvec1=as.double(fvec1),
                        fmin=as.double(fmin),message=as.integer(message),PACKAGE="numOSL")
    } # end if.
    ###
    if (res$message==1) {
        stop("Error: fail in growth curve fitting!")
    } else if (res$message==2) {
        stop("Error: natural OSL signal saturated!")
    } else if (res$message==3) {
        stop("Error: fail in equivalent dose calculation!")
    } else if (res$message==4) {
        stop("Error: fail in equivalent dose error estimation by Simple Transformation!")
    } else if (res$message==5) {
        stop("Error: fail in equivalent dose error estimation by Monte Carlo Simulation!")
    } # end if.
    ###
    LMpars <- cbind(res$pars,res$stdp)
    colnames(LMpars) <- c("Pars","Std.Pars")
    rownames(LMpars) <- (c("a","b","c","d","e"))[seq(n2)]
    ###
    fit.value <- cbind(dose, doseltx, res$fvec1)
    colnames(fit.value) <- c("Redose", "Lx/Tx", "Fit.Lx/Tx")
    rownames(fit.value) <- paste("Redose", seq(ndat), sep="")
    ###
    ED <- matrix(res$outDose, ncol=2L)
    colnames(ED) <- c("ED", "Std.ED")
    ###
    if (ErrorMethod=="mc") {
        mcED <- res$mcED
    } else {
        mcED <- NULL
    } # end if
    ###
    output<-list("mcED"=mcED,
                 "LMpars"=LMpars, 
                 "value"=res$fmin,
                 "fit.value"=fit.value, 
                 "ED"=ED)
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
        ###
        plot(NA, NA, main=mainName, xlab="Dose (Gy)", ylab="Standardised OSL",
             las=0, cex.lab=1.5, cex.main=1.5, xlim=c(lowerX,upperX),
             ylim=c(lowerY,upperY), xaxs="i", yaxs="i", lab=c(7,7,9))
        ###
        if (!is.null(simED) && all(yvalue>0)) {
            dmcED <- density(simED)
            dxy <- cbind(dmcED$x,dmcED$y)
            dxy[,2L] <- (dxy[,2L,drop=TRUE]-min(dxy[,2L,drop=TRUE]))/
                        (max(dxy[,2L,drop=TRUE])-min(dxy[,2L,drop=TRUE]))*
                        yvalue[1L]*0.8
            polygon(dxy, col="grey")
            #rug(simED, quiet=TRUE)
        } # end if
        ###
        points(dose, doseltx, pch=21, cex=1.5, bg="black")
        ###
        x<-NULL
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
        ### 
        legend("topleft", legend=paste("ED=",round(xvalue[1L],2L), "+/-", round(xvalue[2L],2L), " (Gy)", sep=""), 
               yjust=2, ncol=1L, cex=1.5, bty="n")
        ###
        grid()
        box(lwd=1)
        ###
        par(mar=c(5,4,4,2)+0.1)
    } # end function Plot1.
    ###
    ###
    if (plot==TRUE) {
        pars <- LMpars[,1L,drop=TRUE]
        xvalue <- ED[1L,,drop=TRUE]
        yvalue <- inltx[1L,,drop=TRUE]
        simED <- mcED
        mainName <- "Growth Curve"
        Plot1(xvalue, yvalue, simED, mainName, dose,
              doseltx, sdoseltx, pars, model, origin)
    } # end if.
    ###
    return(output)
} # end function calED.default.
#####
