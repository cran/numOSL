#####
fitGrowth<-
function(Curvedata, model="exp", origin=FALSE, 
         weight=TRUE, plot=TRUE) {
    UseMethod("fitGrowth")
} #
### 2016.07.06.
fitGrowth.default<-
function(Curvedata, model="exp", origin=FALSE, 
         weight=TRUE, plot=TRUE) {
    ### Stop if not.
    stopifnot(ncol(Curvedata)==3L, all(Curvedata[,1L,drop=TRUE]>=0), all(Curvedata[,3L,drop=TRUE]>0),
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
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
    if (ndat<n2) {
        stop("Error: data points is not enough for fitting the model!")
    } # end if.
    ###
    pars <- stdp <- vector(length=n2)
    model1 <- if (model=="line") {
        0L } else if (model=="exp") {
        1L } else if (model=="lexp") {
        2L } else if (model=="dexp") {
        3L } # end if.
    uw <- ifelse(weight==FALSE, 0L, 1L)
    fvec1 <- vector(length=ndat)
    fmin <- 0
    message <- 0
    ###
    if (model=="gok") {
        res<-.Fortran("fitGOK",as.double(dose),as.double(doseltx),as.double(sdoseltx),
                      as.integer(ndat),as.integer(n2),pars=as.double(pars),stdp=as.double(stdp),
                      as.integer(uw),fvec1=as.double(fvec1),fmin=as.double(fmin),
                      message=as.integer(message),PACKAGE="numOSL")
    } else {
        res<-.Fortran("fitGrowth",as.double(dose),as.double(doseltx),as.double(sdoseltx),
                      as.integer(ndat),as.integer(n2),pars=as.double(pars),stdp=as.double(stdp),
                      as.integer(model1),as.integer(uw),fvec1=as.double(fvec1),
                      fmin=as.double(fmin),message=as.integer(message),
                      PACKAGE="numOSL")
    } # end if.
    ###
    if (res$message!=0) {
        stop("Error: fail in growth curve fitting!")
    } # end if
    ###
    LMpars <- cbind(res$pars,res$stdp)
    colnames(LMpars) <- c("Pars","Std.Pars")
    rownames(LMpars) <- (c("a","b","c","d","e"))[seq(n2)]
    ###
    fit.value <- cbind(dose, doseltx, res$fvec1)
    colnames(fit.value) <- c("Redose", "Lx/Tx", "Fit.Lx/Tx")
    rownames(fit.value) <- paste("Redose", seq(ndat), sep="")
    ###
    output<-list("LMpars"=LMpars,
                 "value"=res$fmin,
                 "fit.value"=fit.value)
    ###
    if (plot==TRUE) {
        par(mar=c(5,5,4,1)+0.1)
        plot(dose, doseltx, main="Growth Curve", xlab="Dose (Gy)", 
             ylab="Standardised OSL", pch=21, bg="black", cex=1.5, 
             cex.lab=1.5, cex.main=1.5)
        ###
        pars <- LMpars[,1L,drop=TRUE]
        x <- NULL
        cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
        ###
        if(model=="line") {
            curve(pars[1L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if(model=="exp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if(model=="lexp")  {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if(model=="dexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if(model=="gok") {
            curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } # end if.
        ###
        arrowIndex <- which(sdoseltx>0.05 & sdoseltx/doseltx>0.01)
        if (length(arrowIndex)>=1L) {
            arrows(x0=dose[arrowIndex], y0=doseltx[arrowIndex]-sdoseltx[arrowIndex]/2.0, 
                   x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sdoseltx[arrowIndex]/2.0,
                   code=3, lwd=1, angle=90, length=0.05, col="black")
        } # end if.
        ###
        grid()
        box(lwd=2)
        ###
        par(mar=c(5,4,4,2)+0.1)
    } # end if.
    ###
    return(output)
} # end function fitGrowth.default.
#####
