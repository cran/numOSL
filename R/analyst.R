#####
analyst<-
function(Sigdata, Redose, sig.channel=NULL, back.channel=NULL, 
         mr=0.01, disb=c("p","op"), typ="cw", nstart=100,  
         upb=0.5, ErrorMethod=c("mc","sp"), nsim=1000,
         plot=TRUE, model=NULL, origin=NULL) {
    UseMethod("analyst")
} #
### 2014.09.23.
analyst.default<-
function(Sigdata, Redose, sig.channel=NULL, back.channel=NULL, 
         mr=0.01, disb=c("p","op"), typ="cw", nstart=100,  
         upb=0.5, ErrorMethod=c("mc","sp"), nsim=1000,
         plot=TRUE, model=NULL, origin=NULL) {
    ### Stop if not.
    stopifnot(is.data.frame(Sigdata), nrow(Sigdata)>=25L,
              ncol(Sigdata)>=5L, ncol(Sigdata)%%2L==1L,
              is.vector(Redose), length(Redose)==(ncol(Sigdata)-3L)/2L,
              is.null(sig.channel) || is.vector(sig.channel),
              all(sig.channel %in% seq(nrow(Sigdata))),
              is.null(back.channel) || is.vector(back.channel),
              all(back.channel %in% seq(nrow(Sigdata))),
              length(mr)==1L, is.numeric(mr),
              all(disb %in% c("p","op")),
              length(typ)==1L, typ=="cw",
              is.numeric(nstart), nstart>=10L, nstart<=5000L,
              length(upb)==1L, is.numeric(upb), upb>0, upb<=100,
              all(ErrorMethod %in% c("mc","sp")),
              is.numeric(nsim), nsim>=100L, nsim<=3000L,
              length(plot)==1L, is.logical(plot), 
              all(model %in% c("line","exp","lexp","dexp")),
              is.null(origin) || is.logical(origin))
    ###
    ndat<-nrow(Sigdata)
    if (is.null(sig.channel)) {
        sig.channel<-seq(4L)
    } # end if
    if (is.null(back.channel)) {
        back.channel<-(ndat-19L):ndat
    } # end if
    ###
    n<-length(sig.channel)
    m<-length(back.channel)
    k<-m/n
    if (!k %in% seq(10000L)) {
        stop(paste("Error: n.back.channel/n.sig.channel=",round(k,3L),sep=""))
    } # end if
    ### 
    ### R function for calculating Lx/Tx.
    calLxTx<-function(sig1,sig2) {
        ###
        ### Total signal.
        Lx<-sum(sig1[sig.channel])
        Tx<-sum(sig2[sig.channel])
        ###
        ### Background signal.
        Yl<-Yt<-vector(length=k)
        backMat<-matrix(back.channel,nrow=k,byrow=TRUE)
        for (i in seq(k)) {
            Yl[i]<-sum(sig1[backMat[i,,drop=TRUE]])
            Yt[i]<-sum(sig2[backMat[i,,drop=TRUE]])
        } # end if
        bLx<-mean(Yl)
        bTx<-mean(Yt)
        ###
        ### Net signal.
        netLx<-Lx-bLx      
        netTx<-Tx-bTx
        ###
        sigmal<-abs(var(Yl)-mean(Yl))
        sigmat<-abs(var(Yt)-mean(Yt))
        ### 
        if (disb[1L]=="p") {
            ### Eqn.3 of Galbraith (2002).
            rsnetLx<-sqrt(Lx+bLx/k)/(Lx-bLx)
            rsnetTx<-sqrt(Tx+bTx/k)/(Tx-bTx)
        } else if (disb[1L]=="op") {
            ### Eqn.6 of Galbraith (2002).
            rsnetLx<-sqrt(Lx+bLx/k+sigmal*(1.0+1.0/k))/
                     (Lx-bLx)
            rsnetTx<-sqrt(Tx+bTx/k+sigmat*(1.0+1.0/k))/
                     (Tx-bTx)
        } # 
        ###
        rsnetLx<-sqrt(rsnetLx^2L+mr^2L)
        rsnetTx<-sqrt(rsnetTx^2L+mr^2L)
        ###
        LxTx<-netLx/netTx
        sLxTx<-abs(LxTx)*sqrt(rsnetLx^2L+rsnetTx^2L)
        ###
        return(c(LxTx,sLxTx))
    } # end function calLxTx.
    ###
    ###
    nLxTx<-(ncol(Sigdata)-1L)/2L
    matLxTx<-matrix(nrow=nLxTx,ncol=2L)
    for (i in seq(nLxTx)) {
        matLxTx[i,]<-calLxTx(Sigdata[,2*i,drop=TRUE],
                             Sigdata[,2*i+1L,drop=TRUE])
    } # end for
    ###
    if (any(!is.finite(matLxTx))) {
        stop("Fail in calculating Lx/Tx!")
    } # end if
    ###
    Curvedata<-data.frame("Redose"=Redose,
                          "OSL"=matLxTx[-1L,1L,drop=TRUE],
                          "Std.OSL"=matLxTx[-1L,2L,drop=TRUE])
    NatureLxTx<-matLxTx[1L,,drop=TRUE]
    ###
    ### Calculate recycling ratio.
    lvl.dose<-as.numeric(levels(factor(Redose)))
    existrpd<-length(Redose)==length(lvl.dose)+1L
    if (existrpd==TRUE) {
        RepeatIndex<-apply(as.matrix(lvl.dose), MARGIN=1L, function(x,y) 
                     which(abs(x-y)<=.Machine$double.eps^0.5), Redose)
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
    exist0d<-which(abs(Redose)<=.Machine$double.eps^0.5)
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
        res<-calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=model, 
                   origin=origin, nstart=nstart, upb=upb, 
                   ErrorMethod=ErrorMethod, 
                   nsim=nsim, plot=plot)
    } else {
        stop("Error: fail in calculating ED!")
    } # end if   
    ###
    output<-list("Curvedata"=Curvedata, 
                 "Ltx"=NatureLxTx, 
                 "model"=paste(model,"(origin=",origin,")",sep=""), 
                 "LMpars"=res$LMpars, 
                 "value"=res$value, 
                 "ED"=res$ED, 
                 "RecyclingRatio"=RecycleRatio, 
                 "Recuperation"=Recuperation)
    ###
    return(output)
} # end function analyst.default.  
