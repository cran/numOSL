#####
fastED <-
function(Sigdata, Redose, ncomp=2, constant=TRUE,
         control.args=list(), typ="cw", model="exp",
         origin=FALSE, ErrorMethod="mc", 
         nsim=1000, weight=TRUE) {
    UseMethod("fastED")
} ###
### 2016.07.15.
fastED.default <-
function(Sigdata, Redose, ncomp=2, constant=TRUE,
         control.args=list(), typ="cw", model="exp", 
         origin=FALSE, ErrorMethod="mc", 
         nsim=1000, weight=TRUE) {
    ### Stop if not.
    stopifnot(ncol(Sigdata)>=5L, ncol(Sigdata)%%2L==1L, all(Sigdata[,-1L,drop=TRUE]>0),
              is.vector(Redose), length(Redose)==(ncol(Sigdata)-3L)/2L, is.numeric(Redose), all(Redose>=0),
              length(ncomp)==1L, ncomp %in% seq(4L),
              length(constant)==1L, is.logical(constant),
              class(control.args)=="list", all(names(control.args) %in% list("factor","f","cr","maxiter","tol")),
              length(typ)==1L, is.character(typ), typ=="cw",
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(ErrorMethod)==1L, is.character(ErrorMethod), ErrorMethod %in% c("mc","sp"),
              length(nsim)==1L, is.numeric(nsim), nsim>=50L, nsim<=3000L,
              length(weight)>=1L, is.logical(weight))
    ###
    ###
    Plot2<-function(tim,sig,pars,cvalue,samplename) {
        ###
        par(mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
        ###
        plot(tim, sig, main=samplename, log="x", las=0, cex.main=1.25,
             lab=c(7,7,9), ylim=c(-max(sig)*0.01,max(sig)*1.01), cex.lab=1,
             xlab="Stimulated Time (s)", ylab="Photon Counts", xaxs="r", 
             yaxs="i", type="p", pch=21, cex=1.5, bg="white", col="black")
        ### 
        colors<-c("blue", "green", "red", "deepskyblue")
        ###
        nrp <- nrow(pars)
        x <- seq(min(tim), max(tim), by=(max(tim)-min(tim))/length(tim)/100L)
        lines(x, eval(parse(text=paste("pars[",seq(nrp),",1L,drop=TRUE]*pars[",
              seq(nrp),",3L,drop=TRUE]*exp(-pars[",seq(nrp),",3L,drop=TRUE]*x)",
              collapse="+",sep="")))+cvalue,lwd=3,col="black",lty="solid")
        ###
        for (i in seq(nrp)) {
            curve(pars[i,1L,drop=TRUE]*pars[i,3L,drop=TRUE]*
                  exp(-pars[i,3L,drop=TRUE]*x),lwd=3,col=colors[i],
                  lty="solid",add=TRUE)
        } # end for.
        ###
        if (cvalue>0) {
            abline(h=cvalue, lty="dashed", lwd=3)
            legend("topright",legend=c("Fitted.Curve",paste("Comp.",
                   seq(nrp),sep=""),"Constant"),col=c("black", 
                   colors[seq(nrp)],"black"),pch=c(21,rep(NA,nrp+1L)), 
                   lty=c(rep("solid",nrp+1L),"dashed"),yjust=2,ncol=1L, 
                   cex=par("cex"),bty="o",lwd=3,pt.bg="white")
        } else {
            legend("topright",legend=c("Fitted.Curve",paste("Comp.",seq(nrp),sep="")), 
                   col=c("black", colors[seq(nrp)]),pch=c(21,rep(NA,nrp)),lty="solid", 
                   yjust=2,ncol=1L,cex=par("cex"),bty="o",lwd=3,pt.bg="white")
        } # end if.
        ###
        grid(equilogs=FALSE)
        box(lwd=2L)
        par(mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
        ###
    } # end function Plot2.
    ###
    ###
    ncs1 <- ncol(Sigdata)-1L
    fLtx <- matrix(nrow=ncs1, ncol=2L)
    pars <- vector("list", length=ncs1)    
    ###
    args <- list(factor=10L,f=0.5,cr=0.99,maxiter=500L,tol=0.1)
    args[names(control.args)] <- control.args
    factor <- args$factor
    f <- args$f
    cr <- args$cr
    maxiter <- args$maxiter
    tol <- args$tol
    stopifnot(length(factor)==1L, is.numeric(factor), factor>=5L, factor<=50L,
              length(f)==1L, is.numeric(f), f>0.0, f<=1.2,
              length(cr)==1L, is.numeric(cr), cr>0.0, cr<=1.0,
              length(maxiter)==1L, is.numeric(maxiter), maxiter>=10L, maxiter<=5000L,
              length(tol)==1L, is.numeric(tol), tol>0.0)
    ###
    res <- try(decomp(Sigdata[,c(1L,2L),drop=FALSE], 
               ncomp=ncomp, constant=constant, typ=typ, 
               control.args=args, transf=TRUE, weight=FALSE,
               plot=FALSE), silent=TRUE)
    if (class(res)=="try-error") {
        stop("Error: fail in decomposing the natural decay curve!")
    } # end if.
    ###
    ###
    par(mfrow=c(ifelse(ncs1%%4L==0L, ncs1/4L+1L, ncs1%/%4L+1L), 4L))
    pars[[1L]] <- res$LMpars
    fLtx[1L,] <- res$LMpars[which.max(res$LMpars[,3L,drop=TRUE]),seq(2L),drop=TRUE]
    ###
    Plot2(tim=Sigdata[,1L,drop=TRUE], sig=Sigdata[,2L,drop=TRUE], pars=res$LMpars, 
          cvalue=ifelse(constant==FALSE,0,res$constant[1L]), samplename="Natural")   
    ###
    ###
    for (i in 3L:ncol(Sigdata)) {
        res <- try(decomp(Sigdata[,c(1L,i),drop=FALSE], ncomp=ncomp, constant=constant, 
                   typ=typ, control.args=args, transf=TRUE, weight=FALSE, plot=FALSE), silent=TRUE)
        ###
        samplename <- ifelse(i==3L,"Test[Natural]",
                      ifelse(i%%2L==0L,paste("Redose.",i/2L-1L,sep=""), 
                      paste("Test[Redose.",(i-1L)/2L-1L,"]",sep="")))
        ###
        if (class(res)=="try-error") {
            cat(paste("Note: fail in decomposing the ",i-1L,"th decay curve!\n",sep=""))
            par(mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
            ###
            plot(x=Sigdata[,1L,drop=TRUE], y=Sigdata[,i,drop=TRUE], 
                 main=samplename, log="x", las=0, cex.main=1.25,
                 lab=c(7,7,9), ylim=c(-max(Sigdata[,i,drop=TRUE])*0.01, 
                 max(Sigdata[,i,drop=TRUE])*1.2), cex.lab=1, 
                 xlab="Stimulated Time (s)", ylab="Photon Counts",
                 xaxs="r", yaxs="i", typ="p", pch=21, cex=1.5, 
                 bg="white", col="black")
            grid(equilogs=FALSE)
            box(lwd=2L)
            par(mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
        } else {
            pars[[i-1L]] <- res$LMpars
            fLtx[i-1L,] <- res$LMpars[which.max(res$LMpars[,3L,drop=TRUE]),seq(2L),drop=TRUE]
            ###
            Plot2(tim=Sigdata[,1L,drop=TRUE], sig=Sigdata[,i,drop=TRUE], pars=res$LMpars, 
                  cvalue=ifelse(constant==FALSE,0,res$constant[1L]), samplename=samplename)   
        } # end if.
    } # end for.
    ###
    ###
    fLtx[,2L] <- fLtx[,2L,drop=TRUE]/fLtx[,1L,drop=TRUE]
    fLxTx <- fLtx[seq(nrow(fLtx))%%2L==1L,1L,drop=TRUE]/fLtx[seq(nrow(fLtx))%%2L==0L,1L,drop=TRUE]
    sfLxTx <- fLxTx*sqrt((fLtx[seq(nrow(fLtx))%%2L==1L,2L,drop=TRUE])^2L+
        (fLtx[seq(nrow(fLtx))%%2L==0L,2L,drop=TRUE])^2L)
    ###
    Curvedata <- data.frame("Redose"=Redose, "OSL"=fLxTx[-1L], "Std.OSL"=sfLxTx[-1L])
    Curvedata <- Curvedata[complete.cases(Curvedata),,drop=FALSE]
    ###
    NatureLxTx <- c(fLxTx[1L],sfLxTx[1L])
    if(any(!is.finite(NatureLxTx))) { 
        stop("Error: fail in calculating the standardised natural OSL!")
    } # end if.
    ###
    ###
    dose <- Curvedata[,1L,drop=TRUE]
    doseltx <- Curvedata[,2L,drop=TRUE]
    sdoseltx <- Curvedata[,3L,drop=TRUE]
    ###
    ### Calculate recycling ratio.
    lvl.dose <- as.numeric(levels(factor(dose)))
    existrpd <- length(Curvedata[,1L,drop=TRUE])>length(lvl.dose)
    if (existrpd==TRUE) {
        RepeatIndex <- apply(as.matrix(lvl.dose), MARGIN=1L, function(x,y) 
                             which(abs(x-y)<=.Machine$double.eps^0.5), dose)
        RepeatIndex <- unlist(RepeatIndex[sapply(RepeatIndex,length)>=2L])
        RecycleRatio <- doseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[1L]]
        sRecycleRatio <- abs(RecycleRatio)*sqrt((sdoseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[2L]])^2L+
            (sdoseltx[RepeatIndex[1L]]/doseltx[RepeatIndex[1L]])^2L)
    } else {
        RecycleRatio <- sRecycleRatio <- NA
    } # end if
    ###
    ### Calculate recuperation.
    exist0d <- which(abs(dose)<=.Machine$double.eps^0.5)
    if (length(exist0d)>=1L) {
        Recuperation <- (doseltx[exist0d[1L]]/NatureLxTx[1L])*100.0
        sRecuperation <- abs(Recuperation)*sqrt((sdoseltx[exist0d[1L]]/doseltx[exist0d[1L]])^2L+
            (NatureLxTx[2L]/NatureLxTx[1L])^2L)
    } else {
        Recuperation <- sRecuperation <- NA
    } # end if
    ###
    ###
    res <- calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=model, origin=origin, 
                 ErrorMethod=ErrorMethod, nsim=nsim, weight=weight, plot=TRUE)
    ###
    par(mfrow=c(1L,1L))
    ###
    output <- list("decayPars"=pars,
                   "Curvedata"=Curvedata, 
                   "Ltx"=NatureLxTx, 
                   "LMpars"=res$LMpars, 
                   "value"=res$value, 
                   "fastED"=res$ED, 
                   "RecyclingRatio"=c(RecycleRatio,sRecycleRatio), 
                   "Recuperation"=c(Recuperation,sRecuperation))
    ### 
    return(output)
    ###
} # end function fastED.default.
#####                    
