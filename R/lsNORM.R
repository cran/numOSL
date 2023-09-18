#####
lsNORM <- 
function(SARdata, model="gok", origin=FALSE, weight=FALSE, 
         natural.rm=TRUE, norm.dose=NULL, maxiter=10, plot=TRUE) {
    UseMethod("lsNORM")
} ###	
### 2023.08.31.
lsNORM.default <- 
function(SARdata, model="gok", origin=FALSE, weight=FALSE, 
         natural.rm=TRUE, norm.dose=NULL, maxiter=10, plot=TRUE) {
    ### Stop if not.
    stopifnot(ncol(SARdata)==5L, nrow(SARdata)>=5L,
              is.numeric(SARdata[,1L,drop=TRUE]), 
              is.numeric(SARdata[,3L,drop=TRUE]), all(SARdata[,3L,drop=TRUE]>=0),
              is.numeric(SARdata[,4L,drop=TRUE]), 
              is.numeric(SARdata[,5L,drop=TRUE]), all(SARdata[,5L,drop=TRUE]>0),
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(maxiter)==1L, is.numeric(maxiter), maxiter>0, maxiter<=1e3,
              length(weight)==1L, is.logical(weight),
              length(natural.rm)==1L, is.logical(natural.rm),
              is.null(norm.dose) || (length(norm.dose)==1L && is.numeric(norm.dose)),
              length(plot)==1L, is.logical(plot))
    ###
    colnames(SARdata) <- c("NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    NO <- sort(as.numeric(levels(factor(SARdata[,"NO",drop=TRUE]))))
    n <- length(NO)
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
    ###
    ### Check Grain.NO and SAR.Cycle for SARdata.
    for (i in seq(n)) {  
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i])  
        ###
        if (!all(diff(iIndex)==1L)) {
            stop(paste("[NO=", NO[i], 
                       "]: 'NO' appears in discontinuous locations!", sep=""))
        } ### end if. 
        ###     
        iSAR.Cycle <- substr(SARdata[iIndex,"SAR.Cycle",drop=TRUE], start=1L, stop=1L) 
        ###
        ###
        if (!all(iSAR.Cycle %in% c("N","R"))) {
            stop(paste("[NO=", NO[i], 
                       "]: incorrect SAR.Cycle!", sep=""))
        } # end if.
        ###
        if (sum(iSAR.Cycle=="N")>1L) {
            stop(paste("[NO=", NO[i],
                       "]: should contain not more than one SAR.Cycle of 'N'!", sep=""))
        } # end if.
        ###
    } # end for.
    ###
    ### Remove data with SAR.Cycle=N if possible.
    if (natural.rm==TRUE) {
        SARdata <- SARdata[SARdata[,"SAR.Cycle",drop=TRUE]!="N",,drop=FALSE]
    } # end if. 
    ###
    origin_SARdata <- SARdata 
    origin_Curvedata <- SARdata[SARdata[,"SAR.Cycle",drop=TRUE]!="N",
                                c("Dose","Signal","Signal.Err"),drop=FALSE]
    ###
    vary_SARdata <- SARdata 
    vary_Curvedata <- origin_Curvedata
    ###
    iter <- 0
    tol <- 1.0e-06
    ax <- 1.0e-05
    bx <- 1.0e+05
    ###
    ###----------------------------------------------------------------
    cat("LS-normalisation is in progress, please wait, ...\n")
    ###
    repeat {  
        old_rsd <- sd(vary_Curvedata[,"Signal",drop=TRUE])/
                 mean(vary_Curvedata[,"Signal",drop=TRUE])
        ###
        res1 <- try(fitGrowth(Curvedata=vary_Curvedata, model=model, origin=origin, weight=weight, 
                              trial=FALSE, plot=FALSE, nofit.rgd=NULL, agID=NULL, Tn=NULL, Tn3BG=NULL, 
                              TnBG.ratio=NULL, rseTn=NULL, FR=NULL, RecyclingRatio1=NULL, 
                              RecyclingRatio2=NULL, RecyclingRatio3=NULL, Recuperation1=NULL, 
                              Recuperation2=NULL, LnTn.curve=NULL, TxTn=NULL), silent=TRUE)
        ###
        if (inherits(res1,what="try-error")==TRUE || res1$message==1L) {
            if (inherits(res1,what="try-error")==TRUE) print(attr(res1,"condition"))
            stop("Error: growth curve fitting failed!")
        } # end if.
        ###
        pars <- res1$LMpars[,1L,drop=TRUE]
        ###
        for (i in seq(n)) {
            iIndex <- which(vary_SARdata[,"NO",drop=TRUE]==NO[i] &
                            vary_SARdata[,"SAR.Cycle",drop=TRUE]!="N") 
            xx <- vary_SARdata[iIndex,"Dose",drop=TRUE]
            yy <- vary_SARdata[iIndex,"Signal",drop=TRUE]
            nd <- length(xx)
            SF <- fmin <- 0.0
            ###    
            res2 <- .Fortran("calcSF",as.double(ax),as.double(bx),as.double(xx),as.double(yy),
                             as.double(pars),as.integer(nd),as.integer(n2),as.integer(mdl),
                             SF=as.double(SF),fmin=as.double(fmin),PACKAGE="numOSL")
            ###
            iIndex <- which(vary_SARdata[,"NO",drop=TRUE]==NO[i])
            vary_SARdata[iIndex, c("Signal","Signal.Err")] <- 
            vary_SARdata[iIndex, c("Signal","Signal.Err")]*res2$SF
        } # end for. 
        ###
        ###
        vary_Curvedata <- vary_SARdata[vary_SARdata[,"SAR.Cycle",drop=TRUE]!="N",
                                       c("Dose","Signal","Signal.Err"),drop=FALSE]
        ###
        new_rsd <- sd(vary_Curvedata[,"Signal",drop=TRUE])/
                 mean(vary_Curvedata[,"Signal",drop=TRUE])
        ###
        iter <- iter + 1L
        ###
        cat(paste("Iteration=", iter,  ": RSD of SARdata=", old_rsd, "\n",sep=""))
        ###
        if (iter==1L) {
            LMpars1 <- res1$LMpars
            value1 <- res1$value
            avg.error1 <- res1$avg.error
            RCS1 <- res1$RCS
            FOM1 <- res1$FOM
        } # end if
        ###
        if (iter==maxiter) break 
        if (abs(new_rsd-old_rsd)<=tol) break
    } # end repeat.
    ###----------------------------------------------------------------
    ###
    LMpars2 <- res1$LMpars
    value2 <- res1$value
    avg.error2 <- res1$avg.error
    RCS2 <- res1$RCS
    FOM2 <- res1$FOM
    ###
    if (!is.null(norm.dose)) {
        ###
        pars <- res1$LMpars[,1L,drop=TRUE]
        se_pars <- res1$LMpars[,2L,drop=TRUE]
        ###
        if (model=="line") {
            ###
            if (origin==TRUE) {
                norm.signal <- pars[1L]*norm.dose
                pars[1L] <- pars[1L]/norm.signal
                se_pars[1L] <- se_pars[1L]/norm.signal
            } else {
                norm.signal <- pars[1L]*norm.dose+pars[2L]
                pars[c(1L,2L)] <- pars[c(1L,2L)]/norm.signal
                se_pars[c(1L,2L)] <- se_pars[c(1L,2L)]/norm.signal
            } # end if.
            ###
        } else if (model=="exp") {
            ###
            if (origin==TRUE) {
                norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))
                pars[1L] <- pars[1L]/norm.signal
                se_pars[1L] <- se_pars[1L]/norm.signal
            } else {
                norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+pars[3L]
                pars[c(1L,3L)] <- pars[c(1L,3L)]/norm.signal
                se_pars[c(1L,3L)] <- se_pars[c(1L,3L)]/norm.signal
            } # end if.
            ###
        } else if (model=="lexp")  {
            ###
            if (origin==TRUE) {
                norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+
                               pars[3L]*norm.dose
                pars[c(1L,3L)] <- pars[c(1L,3L)]/norm.signal
                se_pars[c(1L,3L)] <- se_pars[c(1L,3L)]/norm.signal
            } else {
                norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+
                               pars[3L]*norm.dose+pars[4L]
                pars[c(1L,3L,4L)] <- pars[c(1L,3L,4L)]/norm.signal
                se_pars[c(1L,3L,4L)] <- se_pars[c(1L,3L,4L)]/norm.signal
            } # end if.
            ###
        } else if(model=="dexp") {
            ###
            if (origin==TRUE) {
                norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+
                               pars[3L]*(1.0-exp(-pars[4L]*norm.dose))
                pars[c(1L,3L)] <- pars[c(1L,3L)]/norm.signal
                se_pars[c(1L,3L)] <- se_pars[c(1L,3L)]/norm.signal
            } else {
                norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+
                               pars[3L]*(1.0-exp(-pars[4L]*norm.dose))+pars[5L]
                pars[c(1L,3L,5L)] <- pars[c(1L,3L,5L)]/norm.signal
                se_pars[c(1L,3L,5L)] <- se_pars[c(1L,3L,5L)]/norm.signal
            } # end if.
            ###
        } else if(model=="gok") {
            ###
            if (origin==TRUE) {
                norm.signal <- pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*norm.dose)^(-1.0/pars[3L]))
                pars[1L] <- pars[1L]/norm.signal
                se_pars[1L] <- se_pars[1L]/norm.signal
            } else {
                norm.signal <- pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*norm.dose)^(-1.0/pars[3L]))+pars[4L]
                pars[c(1L,4L)] <- pars[c(1L,4L)]/norm.signal
                se_pars[c(1L,4L)] <- se_pars[c(1L,4L)]/norm.signal
            } # end if.
            ###
        } # end if.
        ###
        vary_SARdata[, c("Signal","Signal.Err")] <- 
        vary_SARdata[, c("Signal","Signal.Err")]/norm.signal
        ###
        vary_Curvedata <- vary_SARdata[vary_SARdata[,"SAR.Cycle",drop=TRUE]!="N",
                                       c("Dose","Signal","Signal.Err"),drop=FALSE]
        ###
        LMpars2 <- cbind(pars, se_pars)
        colnames(LMpars2) <- c("Pars", "sePars")
        ###
        avg.error2 <- res1$avg.error/norm.signal
        attr(avg.error2, "names") <- NULL
    } # end if.
    ###
    new_SARdata <- vary_SARdata
    new_Curvedata <- vary_Curvedata
    ###
    ###
    sf_vec <- vector(length=n)
    for (i in seq(n)) {
        iIndex <- which(new_SARdata[,"NO",drop=TRUE]==NO[i]) 
        sf_vec[i] <- new_SARdata[iIndex,"Signal"][1L]/
                     origin_SARdata[iIndex,"Signal"][1L]
    } # end if.
    ###   
    ###
    output <- list("norm.SARdata"=new_SARdata,
                   "sf"=sf_vec, "iter"=iter,
                   "LMpars1"=LMpars1, "value1"=value1,
                   "avg.error1"=avg.error1, "RCS1"=RCS1, "FOM1"=FOM1,
                   "LMpars2"=LMpars2, "value2"=value2,
                   "avg.error2"=avg.error2, "RCS2"=RCS2, "FOM2"=FOM2)
    ###
    ###
    if (plot==TRUE) {
        opar <- par("mfrow", "mgp", "mar")
        on.exit(par(opar))
        ###
        layout(matrix(c(1L,1L,1L,2L,2L,2L,1L,1L,1L,2L,2L,2L,
               3L,3L,3L,4L,4L,4L),nrow=6L), respect=FALSE)
        par(mgp=c(2.5,1,0)) 
        ### 
        ### The first plot.
        ###---------------------------------------------------------------------------
        par(mar=c(4,4,2,0.5)+0.1)
        ###
        ylim <- c(min(origin_Curvedata[,"Signal",drop=TRUE],0)*1.1, 
                  max(origin_Curvedata[,"Signal",drop=TRUE])*1.1)
        ###
        plot(origin_Curvedata[,c("Dose","Signal")], type="p", pch=23, bg="red", cex=1.5,
             ylim=ylim, xlab="Regenerative dose (<Gy>|<s>)", ylab="Standardised OSL", 
             las=0, xaxs="r", yaxs="i", cex.lab=1.5, cex.axis=1.5) 
        legend("topleft", legend="A", yjust=2, ncol=1, cex=1.5, bty="o")
        ###
        dose <- origin_Curvedata[,"Dose",drop=TRUE]
        doseltx <- origin_Curvedata[,"Signal",drop=TRUE]
        sdoseltx <- origin_Curvedata[,"Signal.Err",drop=TRUE]
        ###
        arrowIndex <- which(sdoseltx/doseltx>0.001)
        if (length(arrowIndex)>=1L) {
            suppressWarnings(arrows(x0=dose[arrowIndex], 
                y0=doseltx[arrowIndex]-sdoseltx[arrowIndex]/2.0, 
                x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sdoseltx[arrowIndex]/2.0,
                code=3, lwd=1, angle=90, length=0.05, col="gray30"))
        } # end if.
        ###
        x <- NULL
        pars <- LMpars1[,1L,drop=TRUE]
        cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
        ###
        if (model=="line") {
            curve(pars[1L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="exp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="lexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="dexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="gok") {
            curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } # end if.
        ###
        ###------------------------------------------------------------------------------
        ### The second plot.
        par(mar=c(4,4,0.5,0.5)+0.1)
        ylim <- c(min(new_Curvedata[,"Signal",drop=TRUE],0)*1.1, 
                  max(new_Curvedata[,"Signal",drop=TRUE])*1.1)
        plot(new_Curvedata[,c("Dose","Signal")], type="p", pch=23, bg="blue", cex=1.5,
             ylim=ylim, xlab="Regenerative dose (<Gy>|<s>)", ylab="Normalised standardised OSL", 
             las=0, xaxs="r", yaxs="i", cex.lab=1.5, cex.axis=1.5) 
        legend("topleft", legend="B", yjust=2, ncol=1, cex=1.5, bty="o")
        ###
        dose <- new_Curvedata[,"Dose",drop=TRUE]
        doseltx <- new_Curvedata[,"Signal",drop=TRUE]
        sdoseltx <- new_Curvedata[,"Signal.Err",drop=TRUE]
        ###
        arrowIndex <- which(sdoseltx/doseltx>0.001)
        if (length(arrowIndex)>=1L) {
            suppressWarnings(arrows(x0=dose[arrowIndex], 
                y0=doseltx[arrowIndex]-sdoseltx[arrowIndex]/2.0, 
                x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sdoseltx[arrowIndex]/2.0,
                code=3, lwd=1, angle=90, length=0.05, col="gray30"))
        } # end if.
        ###
        x <- NULL
        pars <- LMpars2[,1L,drop=TRUE]
        cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
        ###
        if (model=="line") {
            curve(pars[1L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="exp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="lexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="dexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="gok") {
            curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } # end if.
        ###
        ###--------------------------------------------------------------------------
        ### The thrid plot.
        par(mar=c(4,0.5,2,0.5)+0.1)
        par(mgp=c(1,1,0))
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary A", ylab="", cex.lab=1.5)
        ###
        pars <- LMpars1[,1L,drop=TRUE]
        se_pars <-  LMpars1[,2L,drop=TRUE]
        ###
        legend("center", 
               legend=c("Before LS-normalisation", 
               "======================",
               paste("Fit model: ", model, sep=""),
               paste("Pass origin: ", origin, sep=""),
               paste("Weighted fit: ", weight, sep=""),
               paste("Minimized value: ", round(value1,2L), sep=""),
               paste("Average error in fit: ", round(avg.error1,2L), sep=""),
               paste("Reduced Chi-Square: ", round(RCS1,2L), sep=""),
               paste("Figure Of Merit: ", round(FOM1,2L)," (%)", sep=""),
               "======================",
               paste("a=", signif(pars[1L],2L), " +/- ", signif(se_pars[1L],2L), sep=""),
               paste("b=", ifelse(length(pars)>=2L, paste(signif(pars[2L],2L),
                     " +/- ",signif(se_pars[2L],2L),sep=""), "NULL"), sep=""),
               paste("c=", ifelse(length(pars)>=3L, paste(signif(pars[3L],2L),
                     " +/- ",signif(se_pars[3L],2L),sep=""), "NULL"), sep=""),
               paste("d=", ifelse(length(pars)>=4L, paste(signif(pars[4L],2L),
                     " +/- ",signif(se_pars[4L],2L),sep=""), "NULL"), sep=""),
               paste("e=", ifelse(length(pars)>=5L, paste(signif(pars[5L],2L),
                     " +/- ",signif(se_pars[5L],2L),sep=""), "NULL"), sep="")),  
               yjust=2, ncol=1, cex=1.2, bty="n")
        box(lwd=1L)
        ###
        ###---------------------------------------------------------------------------
        ### The fourth plot.
        par(mar=c(4,0.5,0.5,0.5)+0.1)
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary B", ylab="", cex.lab=1.5)
        ###
        pars <- LMpars2[,1L,drop=TRUE]
        se_pars <-  LMpars2[,2L,drop=TRUE]
        ###
        legend("center", 
               legend=c("After LS-normalisation", 
               "======================",
               paste("Fit model: ", model, sep=""),
               paste("Pass origin: ", origin, sep=""),
               paste("Weighted fit: ", weight, sep=""),
               paste("Minimized value: ", round(value2,2L), sep=""),
               paste("Average error in fit: ", round(avg.error2,2L), sep=""),
               paste("Reduced Chi-Square: ", round(RCS2,2L), sep=""),
               paste("Figure Of Merit: ", round(FOM2,2L)," (%)", sep=""),
               "======================",
               paste("a=", signif(pars[1L],2L), " +/- ", signif(se_pars[1L],2L), sep=""),
               paste("b=", ifelse(length(pars)>=2L, paste(signif(pars[2L],2L),
                     " +/- ",signif(se_pars[2L],2L),sep=""), "NULL"), sep=""),
               paste("c=", ifelse(length(pars)>=3L, paste(signif(pars[3L],2L),
                     " +/- ",signif(se_pars[3L],2L),sep=""), "NULL"), sep=""),
               paste("d=", ifelse(length(pars)>=4L, paste(signif(pars[4L],2L),
                     " +/- ",signif(se_pars[4L],2L),sep=""), "NULL"), sep=""),
               paste("e=", ifelse(length(pars)>=5L, paste(signif(pars[5L],2L),
                     " +/- ",signif(se_pars[5L],2L),sep=""), "NULL"), sep="")),
               yjust=2, ncol=1, cex=1.2, bty="n")
        box(lwd=1L)
        ###
    } # end if. 
    ### 
    invisible(output)
} # end function lsNORM.default.
#####
