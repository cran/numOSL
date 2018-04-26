#####
fastED <-
function(Sigdata, Redose, delay.off=c(0,0), ncomp=2, constant=TRUE, control.args=list(), 
         typ="cw", model="gok", origin=FALSE, errMethod="sp",nsim=500, weight.decomp=FALSE, 
         weight.fitGrowth=TRUE, trial=TRUE, nofit.rgd=NULL, outpdf=NULL, log="x", lwd=2, 
         test.dose=NULL, agID=NULL) {
    UseMethod("fastED")
} ###
### 2018.04.23.
fastED.default <-
function(Sigdata, Redose, delay.off=c(0,0), ncomp=2, constant=TRUE, control.args=list(), 
         typ="cw", model="gok", origin=FALSE, errMethod="sp",nsim=500, weight.decomp=FALSE, 
         weight.fitGrowth=TRUE, trial=TRUE, nofit.rgd=NULL, outpdf=NULL, log="x", lwd=2, 
         test.dose=NULL, agID=NULL) {
    ### Stop if not.
    stopifnot(ncol(Sigdata)>=5L, ncol(Sigdata)%%2L==1L,
              is.numeric(Redose), all(Redose>=0),
              length(delay.off)==2L, all(delay.off>=0L),
              length(ncomp)==1L, ncomp %in% seq(5L),
              length(constant)==1L, is.logical(constant),
              class(control.args)=="list", 
              all(names(control.args) %in% list("factor","f","cr","maxiter","tol")),
              length(typ)==1L, typ=="cw",
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(errMethod)==1L, errMethod %in% c("sp","mc"),
              length(nsim)==1L, is.numeric(nsim), nsim>=50L, nsim<=3000L,
              length(weight.decomp)==1L, is.logical(weight.decomp),
              length(weight.fitGrowth)==1L, is.logical(weight.fitGrowth),
              length(trial)==1L, is.logical(trial),
              is.null(nofit.rgd) || is.numeric(nofit.rgd),
              is.null(outpdf) || (length(outpdf)==1L && is.character(outpdf)),
              length(log)==1L, log %in% c("", "x", "y", "xy", "yx"),
              length(lwd)==1L, is.numeric(lwd),
              is.null(test.dose) || (length(test.dose)==1L && is.numeric(test.dose)),
              is.null(agID) || (length(agID)==3L && is.numeric(agID)))
    ###
    if (any(Sigdata[,,drop=TRUE]<=0)) stop("Error: values in 'Sigdata' should larger than zero!")
    ###
    if (length(Redose)!=(ncol(Sigdata)-3L)/2L) stop("Error: incorrect length of 'Redose'!")
    ###
    length_curve <- ncol(Sigdata)-1L
    fast_Ltx <- matrix(nrow=length_curve, ncol=2L)
    decomp_pars <- list()
    ###
    args <- list(factor=20L,f=0.5,cr=0.99,maxiter=500L,tol=0.1)
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
    if (!is.null(outpdf)) {
        pdf(paste(outpdf,".pdf",sep=""))
        if_plot <- TRUE
    } else {
        if_plot <- FALSE
    } # end if.
    ###
    LnTn.curve <- list()
    ###
    ###-------------------------------------------------------------------------------
    res <- try(decomp(Sigdata[,c(1L,2L),drop=FALSE], delay.off=delay.off, ncomp=ncomp, 
                      constant=constant, typ=typ, control.args=args, weight=weight.decomp, 
                      plot=if_plot, log=log, lwd=lwd, curve.no=1L, SAR.Cycle="Natural", 
                      irr.dose=0, outfile=NULL, transf=TRUE), silent=TRUE)
    ###
    if (class(res)!="try-error" && res$message==0L) {
        decomp_pars[["N"]] <- res$LMpars
        fast_Ltx[1L,] <- res$LMpars[which.max(res$LMpars[,3L,drop=TRUE]),seq(2L),drop=TRUE]    
    } else {
        if (!is.null(outpdf)) dev.off()
        if (class(res)=="try-error")  print(attr(res,"condition"))
        stop("Error: fail in natural (1st) decay curve fitting!")
    } # end if.
    LnTn.curve[["Ln.x"]] <- res$comp.sig[,"Time",drop=TRUE]
    LnTn.curve[["Ln.y"]] <- res$comp.sig[,"Comp.1",drop=TRUE]
    ###--------------------------------------------------------------------------------
    res <- try(decomp(Sigdata[,c(1L,3L),drop=FALSE], delay.off=delay.off, ncomp=ncomp, 
                      constant=constant, typ=typ, control.args=args, weight=weight.decomp, 
                      plot=if_plot, log=log, lwd=lwd, curve.no=2L, SAR.Cycle="Test [Natural]", 
                      irr.dose=test.dose, outfile=NULL, transf=TRUE), silent=TRUE)
    ###
    if (class(res)!="try-error" && res$message==0L) {
        decomp_pars[["TN"]] <- res$LMpars
        fast_Ltx[2L,] <- res$LMpars[which.max(res$LMpars[,3L,drop=TRUE]),seq(2L),drop=TRUE]    
    } else {
        if (!is.null(outpdf)) dev.off()
        if (class(res)=="try-error")  print(attr(res,"condition"))
        stop("Error: fail in natural test dose (2st) decay curve fitting!")
    } # end if.
    LnTn.curve[["Tn.x"]] <- res$comp.sig[,"Time",drop=TRUE]
    LnTn.curve[["Tn.y"]] <- res$comp.sig[,"Comp.1",drop=TRUE]
    ###--------------------------------------------------------------------------------
    ###
    for (i in 4L:ncol(Sigdata)) {
        iSAR.Cycle <- ifelse(i%%2L==0L,paste("Redose.",i/2L-1L,sep=""), 
                             paste("Test [Redose.",(i-1L)/2L-1L,"]",sep=""))
        if (i%%2L==0L) i_irr.dose <- Redose[i/2L-1L] else i_irr.dose <- test.dose
        ###
        res <- try(decomp(Sigdata[,c(1L,i),drop=FALSE], delay.off=delay.off, ncomp=ncomp, 
                          constant=constant, typ=typ, control.args=args, weight=weight.decomp, 
                          plot=if_plot, log=log, lwd=lwd, curve.no=i-1L, SAR.Cycle=iSAR.Cycle, 
                          irr.dose=i_irr.dose, outfile=NULL, transf=TRUE), silent=TRUE)
        ###
        if (class(res)!="try-error" && res$message==0L) {
            characterNO <- ifelse(i%%2L==0L,paste("R",i/2L-1L,sep=""), 
                                  paste("TR",(i-1L)/2L-1L,sep=""))
            decomp_pars[[characterNO]] <- res$LMpars
            fast_Ltx[i-1L,] <- res$LMpars[which.max(res$LMpars[,3L,drop=TRUE]),seq(2L),drop=TRUE]
            ###
        } else {
            if (class(res)=="try-error")  print(attr(res,"condition"))
            cat(paste("Note: fail in fit the ",i-1L,"th decay curve!\n",sep=""))
            fast_Ltx[i-1L,] <- NA
        } # end if.
    } # end for.
    ###--------------------------------------------------------------------------------------
    ###
    odd_index <- which(seq(length_curve)%%2L==1L)
    even_index <- which(seq(length_curve)%%2L==0L)
    ###
    Lx_vec <- fast_Ltx[odd_index]
    Tx_vec <- fast_Ltx[even_index]
    rseLx_vec <- (fast_Ltx[,2L,drop=TRUE]/fast_Ltx[,1L,drop=TRUE])[odd_index]
    rseTx_vec <- (fast_Ltx[,2L,drop=TRUE]/fast_Ltx[,1L,drop=TRUE])[even_index]
    ###
    TxTn_vec <- Tx_vec[is.finite(Tx_vec)]
    TxTn_vec <- TxTn_vec/TxTn_vec[1L]
    ###
    LxTx_vec <- Lx_vec/Tx_vec
    seLxTx_vec <- abs(Lx_vec/Tx_vec)*sqrt(rseLx_vec^2L+rseTx_vec^2L)
    ###
    Curvedata <- data.frame("Redose"=Redose, "LxTx"=LxTx_vec[-1L], "seLxTx"=seLxTx_vec[-1L])
    Curvedata <- Curvedata[complete.cases(Curvedata),,drop=FALSE]
    if (nrow(Curvedata)==0L) {
        if (!is.null(outpdf)) dev.off() 
        stop("Error: no data can be used to build the fast-component growth curve!")
    } # end if.
    ###
    Nature_LxTx <- c(LxTx_vec[1L], seLxTx_vec[1L])
    ###
    res <- try(calED(Curvedata=Curvedata, Ltx=Nature_LxTx, model=model, origin=origin, 
                     errMethod=errMethod, nsim=nsim, weight=weight.fitGrowth, trial=trial, 
                     plot=if_plot, nofit.rgd=nofit.rgd, agID=agID, Tn=NULL, Tn3BG=NULL, 
                     TnBG.ratio=NULL, rseTn=NULL, FR=NULL, LnTn.curve=LnTn.curve, TxTn=TxTn_vec), 
                     silent=TRUE)
    ###
    if (!is.null(outpdf)) dev.off() 
    ### 
    if (class(res)=="try-error") {
        print(attr(res, "condition"))
        stop("Error: fast-component equivalent dose calculation fails!")
    } # end if.
    ###
    output <- list("decomp.pars"=decomp_pars,
                   "Curvedata"=Curvedata[res$fitIDX,], 
                   "Ltx"=Nature_LxTx, 
                   "LMpars"=res$LMpars, 
                   "value"=res$value, 
                   "avg.error"=res$avg.error,
                   "RCS"=res$RCS,
                   "FOM"=res$FOM, 
                   "calED.method"=res$calED.method,
                   "mcED"=res$mcED,
                   "ED"=res$ED,
                   "ConfInt"=res$ConfInt,
                   "RecyclingRatio1"=res$RecyclingRatio1,
                   "RecyclingRatio2"=res$RecyclingRatio2,
                   "RecyclingRatio3"=res$RecyclingRatio3,
                   "Recuperation1"=res$Recuperation1,
                   "Recuperation2"=res$Recuperation2)
    ###
    invisible(output)
    ###
} # end function fastED.default.
#####
