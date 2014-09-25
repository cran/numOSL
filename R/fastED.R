#####
fastED<-
function(Sigdata, Redose, ncomp=2, constant=TRUE,
         control.args=list(), typ="cw", nstart=100, upb=0.5,
         ErrorMethod=c("mc","sp"), model=NULL, origin=NULL) {
    UseMethod("fastED")
} ###
### 2014.09.14.
fastED.default<-
function(Sigdata, Redose, ncomp=2, constant=TRUE,
         control.args=list(), typ="cw", nstart=100, upb=0.5,
         ErrorMethod=c("mc","sp"), model=NULL, origin=NULL) {
    ### Stop if not.
    stopifnot(is.data.frame(Sigdata), 
              ncol(Sigdata)>=5L, ncol(Sigdata)%%2L==1L,
              is.vector(Redose), length(Redose)==(ncol(Sigdata)-3L)/2L,
              length(ncomp)==1L, ncomp %in% (1L:4L),
              length(constant)==1L, is.logical(constant),
              class(control.args)=="list",
              all(names(control.args) %in% 
              list("factor","f","cr","maxiter","tol")),
              length(typ)==1L, typ=="cw",
              is.numeric(nstart), nstart>=10L, nstart<=5000L,
              length(upb)==1L, is.numeric(upb), upb>0, upb<=100,
              all(ErrorMethod %in% c("mc","sp")),
              all(model %in% c("line","exp","lexp","dexp")),
              is.null(origin) || is.logical(origin))
    ###
    ###
    Plot2<-function(tim,sig,pars,cvalue,samplename) {
        ###
        par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
        ###
        plot(tim, sig, main=samplename, log="x", las=0, cex.main=1.25,
             lab=c(7,7,9), ylim=c(-max(sig)*0.01,max(sig)*1.01), cex.lab=1,
             xlab="Stimulated Time (s)", ylab="Photon Counts", xaxs="r", 
             yaxs="i", type="p", pch=21, cex=1.5, bg="white", col="black")
        ### 
        colors<-c("blue", "green", "red", "deepskyblue")
        ###
        nrp<-nrow(pars)
        x<-seq(min(tim), max(tim), by=(max(tim)-min(tim))/length(tim)/100L)
        lines(x, eval(parse(text=paste("pars[",seq(nrp),",1L,drop=TRUE]*pars[",
              seq(nrp),",3L,drop=TRUE]*exp(-pars[",seq(nrp),",3L,drop=TRUE]*x)",
              collapse="+",sep="")))+cvalue,lwd=3,col="black",lty="solid")
        ###
        for (i in seq(nrp)) {
            curve(pars[i,1L,drop=TRUE]*pars[i,3L,drop=TRUE]*
                  exp(-pars[i,3L,drop=TRUE]*x),lwd=3,col=colors[i],
                  lty="solid",add=TRUE)
        } # end for
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
        } # end if
        ###
        grid(equilogs=FALSE)
        box(lwd=2L)
        par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
        ###
    } # end function Plot2.
    ###
    ###
    ncs1<-ncol(Sigdata)-1L
    fLtx<-matrix(nrow=ncs1, ncol=2L)
    pars<-vector("list", length=ncs1)    
    ###
    args<-list(factor=10L,f=0.5,cr=0.99,maxiter=500L,tol=0.1)
    args[names(control.args)]<-control.args
    factor<-args$factor
    f<-args$f
    cr<-args$cr
    maxiter<-args$maxiter
    tol<-args$tol
    stopifnot(is.numeric(factor), factor>=3L, factor<=50L,
              is.numeric(f), f>0.0, f<=1.2,
              is.numeric(cr), cr>0.0, cr<=1.0,
              is.numeric(maxiter), maxiter>=10L, maxiter<=5000L,
              is.numeric(tol), tol>0.0)
    ###
    res<-try(decomp(Sigdata[,c(1L,2L),drop=FALSE], 
             ncomp=ncomp, constant=constant, typ=typ, 
             control.args=args, transf=TRUE, 
             plot=FALSE), silent=TRUE)
    if (class(res)=="try-error") {
        stop("Error: fail in decomposing the natural decay curve!")
    } # end if
    ###
    par(mfrow=c(ifelse(ncs1%%4L==0L, ncs1/4L+1L, ncs1%/%4L+1L), 4L))
    pars[[1L]]<-res$pars[,-5L,drop=FALSE]
    fLtx[1L,]<-res$pars[which.max(res$pars[,3L,drop=TRUE]),seq(2L),drop=TRUE]
    ###
    Plot2(tim=Sigdata[,1L,drop=TRUE], sig=Sigdata[,2L,drop=TRUE],
          pars=res$pars, cvalue=ifelse(constant==FALSE,0,
          res$constant[1L]), samplename="Natural")   
    ###
    ###
    for (i in 3L:ncol(Sigdata)) {
        res<-try(decomp(Sigdata[,c(1L,i),drop=FALSE], 
                 ncomp=ncomp, constant=constant, typ=typ, 
                 control.args=args, transf=TRUE, 
                 plot=FALSE), silent=TRUE)
        ###
        samplename<-ifelse(i==3L,"Test[Natural]",
                    ifelse(i%%2L==0L,paste("Redose.",i/2L-1L,sep=""), 
                    paste("Test[Redose.",(i-1L)/2L-1L,"]",sep="")))
        ###
        if (class(res)=="try-error") {
            cat(paste("Note: fail in decomposing the ",i-1L,"th decay curve!\n",sep=""))
            par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
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
            par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
        } else {
            pars[[i-1L]]<-res$pars[,-5L,drop=FALSE]
            fLtx[i-1L,]<-res$pars[which.max(res$pars[,3L,drop=TRUE]),seq(2L),drop=TRUE]
            ###
            Plot2(tim=Sigdata[,1L,drop=TRUE], sig=Sigdata[,i,drop=TRUE],
                  pars=res$pars, cvalue=ifelse(constant==FALSE,0,
                  res$constant[1L]), samplename=samplename)   
        } # end if
    } # end for
    ###
    ###
    decayRateMat<-matrix(unlist(sapply(pars,function(x) 
                  x[,3L,drop=TRUE])), ncol=ncomp, byrow=TRUE)
    ###
    fLtx[,2L]<-fLtx[,2L,drop=TRUE]/fLtx[,1L,drop=TRUE]
    fLxTx<-fLtx[seq(nrow(fLtx))%%2L==1L,1L,drop=TRUE]/
           fLtx[seq(nrow(fLtx))%%2L==0L,1L,drop=TRUE]
    sfLxTx<-fLxTx*sqrt((fLtx[seq(nrow(fLtx))%%2L==1L,2L,drop=TRUE])^2L+ 
                       (fLtx[seq(nrow(fLtx))%%2L==0L,2L,drop=TRUE])^2L)
    Curvedata<-data.frame("Redose"=Redose,
                          "OSL"=fLxTx[-1L],
                          "Std.OSL"=sfLxTx[-1L])
    Curvedata<-Curvedata[complete.cases(Curvedata),,drop=FALSE]
    NatureLxTx<-c(fLxTx[1L],sfLxTx[1L])
    if(any(!is.finite(NatureLxTx))) { 
        stop("Error: fail in calculating the standardised natural OSL!")
    } # end if
    ###
    ### Calculate recycling ratio.
    lvl.dose<-as.numeric(levels(factor(Curvedata[,1L,drop=TRUE])))
    existrpd<-length(Curvedata[,1L,drop=TRUE])==length(lvl.dose)+1L
    if (existrpd==TRUE) {
        RepeatIndex<-apply(as.matrix(lvl.dose), MARGIN=1L, function(x,y) 
                     which(abs(x-y)<=.Machine$double.eps^0.5), Curvedata[,1L,drop=TRUE])
        RepeatIndex<-unlist(RepeatIndex[sapply(RepeatIndex,length)==2L])
        RecycleRatio<-Curvedata[,2L,drop=TRUE][RepeatIndex[2L]]/
                      Curvedata[,2L,drop=TRUE][RepeatIndex[1L]]
        if (RecycleRatio>1.3 || RecycleRatio<0.7) {
            cat("Note: recycling ratio is large (small)!\n")
        } # end if
    } else {
        cat("Note: recycling ratio is not available!\n")
        RecycleRatio<-NULL
    } # end if
    ###
    ### Calculate recuperation.
    exist0d<-which(abs(Curvedata[,1L,drop=TRUE])<=.Machine$double.eps^0.5)
    if (length(exist0d)>0L) {
        Recuperation<-Curvedata[exist0d[1L],2L,drop=TRUE]/NatureLxTx[1L]
        if (Recuperation>0.3) {
            cat("Note: recuperation is large!\n")
        } # end if
    } else {
        cat("Note: recuperation is not available!\n")
        Recuperation<-NULL
    } # end if
    ###
    ###
    Models<-if(is.null(model)) c("line","exp","lexp") else model
    Origins<-if(is.null(origin)) c(TRUE,FALSE) else origin
    min.chi.square<-1e20
    ###
    for (i in Models) {
        for (j in Origins) {
            res<-try(calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=i,
                           origin=j, nstart=nstart, upb=upb, 
                           ErrorMethod=ErrorMethod, nsim=100L, 
                           plot=FALSE), silent=TRUE)
            ###
            if (class(res)!="try-error") {
                if (res$value<min.chi.square) {
                    model<-i
                    origin<-j
                    min.chi.square<-res$value
                    OK<-1L
                } # end if
            } # end if
        } # end for
    } # end for
    ###
    if (exists("OK")) {
        res<-calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=model, origin=origin, 
                   nstart=nstart, upb=upb, ErrorMethod=ErrorMethod, 
                   nsim=1000L, plot=TRUE)
    } else {
        par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
        plot(Curvedata[,c(1L,2L)], main="Growth Curve", pch=21, cex=3, 
             bg="white", xlab="Dose (Gy)", ylab="Standardised OSL", las=0, 
             cex.main=1.25, cex.lab=1, xlim=c(0,max(Curvedata[,1L,drop=TRUE])*1.2), 
             ylim=c(0,max(Curvedata[,2L,drop=TRUE])*1.2), 
             xaxs="i", yaxs="i", lab=c(7,7,9))
        grid(equilogs=FALSE)
        box(lwd=2L)
        par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
        stop("Error: fail in calculating the fast-component ED!")
    } # end if
    ###
    par(mfrow=c(1L,1L))
    output<-list("pars"=pars, 
                 "decayRateMat"=decayRateMat,
                 "Curvedata"=Curvedata, 
                 "Ltx"=NatureLxTx, 
                 "model"=paste(model,"(origin=",origin,")",sep=""), 
                 "LMpars"=res$LMpars, 
                 "value"=res$value, 
                 "ED"=res$ED, 
                 "RecyclingRatio"=RecycleRatio, 
                 "Recuperation"=Recuperation)
    ### 
    return(output)
    ###
} # end function fastED.
#####                    
