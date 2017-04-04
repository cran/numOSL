#####
decomp <-
function(Sigdata, delay.off=c(0,0), ncomp=2, constant=TRUE, 
         typ="cw", control.args=list(), weight=FALSE, plot=TRUE, 
         log="x", lwd=2, curve.no=NULL, SAR.Cycle=NULL, 
         irr.dose=NULL, outfile=NULL, transf=TRUE) {
    UseMethod("decomp")
} ###
### 2017.03.30.
decomp.default <-
function(Sigdata, delay.off=c(0,0), ncomp=2, constant=TRUE, 
         typ="cw", control.args=list(), weight=FALSE, plot=TRUE, 
         log="x", lwd=2, curve.no=NULL, SAR.Cycle=NULL, 
         irr.dose=NULL, outfile=NULL, transf=TRUE) {
    #### Stop if not.
    stopifnot(ncol(Sigdata)==2L,
              length(delay.off)==2L, all(delay.off>=0L),
              length(ncomp)==1L, ncomp %in% seq(7L),
              length(constant)==1L, is.logical(constant),
              length(typ)==1L, typ %in% c("cw","lm"),
              class(control.args)=="list",
              all(names(control.args) %in% list("factor","f","cr","maxiter","tol")),
              length(weight)==1L, is.logical(weight),
              length(plot)==1L, is.logical(plot),
              length(log)==1L, log %in% c("", "x", "y", "xy", "yx"),
              length(lwd)==1L, is.numeric(lwd),
              is.null(curve.no) || (length(curve.no)==1L && is.numeric(curve.no)),
              is.null(SAR.Cycle) || (length(SAR.Cycle)==1L && is.character(SAR.Cycle)),
              is.null(irr.dose) || (length(irr.dose)==1L && is.numeric(irr.dose)),
              is.null(outfile) || (length(outfile)==1L && is.character(outfile)),
              length(transf)==1L, is.logical(transf))
    ###
    tim <- as.numeric(Sigdata[,1L,drop=TRUE])
    if (any(tim<=0)) stop("Error: time value in [Sigdata] should larger than zero!")
    ###
    sig <- as.numeric(Sigdata[,2L,drop=TRUE])
    if (any(sig<=0)) stop("Error: signal value in [Sigdata] should larger than zero!")
    ###
    ntim <- length(tim)
    if (sum(delay.off)>=ntim) stop("Error: sum of 'delay.off' exceeds number of data points!")
    ###
    if (delay.off[1L]>0L || delay.off[2L]>0L) {
        ntim <- ntim - sum(delay.off)
        tim <- seq(from=tim[1L], to=ntim*tim[1L], by=tim[1L])
        sig <- sig[(delay.off[1L]+1L):(delay.off[1L]+ntim)]
    } # end if.
    ###
    addc <- ifelse(constant==TRUE, 1L, 0L)
    n2 <- 2L*ncomp + addc
    uw <- ifelse(weight==FALSE, 0L, 1L)
    pars <- stdp <- vector(length=n2)
    typ1 <- ifelse(typ=="cw",1L, 2L)
    ###
    ### 
    if (ntim<n2) stop("Error: data points is not enough for model fitting!")
    ###
    ###
    ### Default parameters for Deffirential Evolution algorithm.
    args <- list(factor=20L, f=0.5, cr=0.99, maxiter=500L, tol=0.1)
    args[names(control.args)] <- control.args
    ###
    factor <- args$factor
    f <- args$f
    cr <- args$cr
    maxiter <- args$maxiter
    tol <- args$tol
    ###
    stopifnot(length(factor)==1L, is.numeric(factor), factor>=5L, factor<=50L,
              length(f)==1L, is.numeric(f), f>0.0, f<=1.2,
              length(cr)==1L, is.numeric(cr), cr>0.0, cr<=1.0,
              length(maxiter)==1L, is.numeric(maxiter),  maxiter>=10L, maxiter<=5000L,
              length(tol)==1L, is.numeric(tol), tol>0.0)
    ####
    fvec1 <- vector(length=ntim)
    fmin <- 0.0
    message <- 0
    ###
    res<-.Fortran("decomp_fort",as.double(tim),as.double(sig),as.integer(ntim), 
                  pars=as.double(pars),stdp=as.double(stdp),as.integer(n2),
                  as.integer(uw),as.integer(addc),as.integer(typ1),as.integer(factor), 
                  as.double(f),as.double(cr),as.integer(maxiter),as.double(tol), 
                  fvec1=as.double(fvec1),fmin=as.double(fmin),message=as.integer(message),
                  PACKAGE="numOSL")
    message <- res$message
    ###
    if (message==0L) {
        pars <- cbind(res$pars[seq(ncomp)], 
                      res$stdp[seq(ncomp)],
                      res$pars[(ncomp+1L):(2L*ncomp)], 
                      res$stdp[(ncomp+1L):(2L*ncomp)])
        pars <- pars[order(pars[,3L,drop=TRUE],decreasing=TRUE), , drop=FALSE]
        colnames(pars) <- c("Ithn", "seIthn", "Lamda", "seLamda")
        rownames(pars) <- paste("Comp.", seq(ncomp), sep="")
        ###
        if (constant==TRUE) {
            constant_d <- c("Constant"=res$pars[2L*ncomp+1L],
                            "seConstant"=res$stdp[2L*ncomp+1L])
        } else {
            constant_d <- 0.0
        } # end if.
        ###
        comp.sig <- apply(cbind(pars[,1L,drop=TRUE], pars[,3L,drop=TRUE]), MARGIN=1L,
                          function(x) if(typ=="cw") x[1L]*x[2L]*exp(-x[2L]*tim) else 
                          x[1L]*x[2L]*(tim/max(tim))*exp(-x[2L]*tim^2L/2L/max(tim)))
        ###
        if (constant==TRUE) {
            comp.sig <- cbind(tim, sig, res$fvec1, comp.sig, if (typ=="cw") constant_d[1L] else constant_d[1L]*tim/max(tim))
            colnames(comp.sig)<-c("Time","Signal","Comp.Sum", paste("Comp.",seq(ncomp),sep=""),"Background")
        } else {
            comp.sig <- cbind(tim, sig, res$fvec1, comp.sig)
            colnames(comp.sig)<-c("Time","Signal","Comp.Sum", paste("Comp.",seq(ncomp),sep=""))
        } # end if.
        ###
        if (!is.null(outfile)) write.csv(comp.sig, file=paste(outfile, ".csv", sep=""))
        ###
        value <- res$fmin
        FOM <- 100.0*sum(abs(sig-res$fvec1))/sum(res$fvec1)
    } else {
        pars <- NULL
        constant_d <- NULL
        comp.sig <- NULL
        value <- NULL
        FOM <- NULL
    } # end if.   
    ### 
    ###
    if (plot==TRUE) {
        layout(matrix(c(1L,1L,1L,1L,2L,2L,
                        1L,1L,1L,1L,2L,2L,
                        3L,3L,3L,4L,4L,4L),nrow=6L), 
                        respect=FALSE)
        par(mgp=c(2.5,1,0))
        ###--------------------------------------------------------------------
        par(mar=c(4,4,2,0.5)+0.1)
        plot(x=tim, y=sig, log=log, las=0, lab=c(7,7,9), ylim=c(1.0, max(sig)*1.2),
             xlab="Stimulation time (s)", ylab="Photon counts", cex.lab=1.5, 
             xaxs="r", yaxs="r", type="p", pch=21, cex=ifelse(typ=="cw",1.5,1.25),
             bg="white", col="black")
        grid(col="pink3")
        box(lwd=1L)
        ###
        if (message==0L) {
            colors<-c("blue", "green", "red", "deepskyblue","purple", "orange", "brown")
            x <- seq(min(tim),max(tim),by=(max(tim)-min(tim))/ntim/100L)
            ###
            if (typ=="cw") {
                lines(x, eval(parse(text=paste("pars[",seq(ncomp),
                      ",1L,drop=TRUE]*pars[",seq(ncomp),
                      ",3L,drop=TRUE]*exp(-pars[",
                      seq(ncomp),",3L,drop=TRUE]*x)", 
                      collapse="+",sep="")))+constant_d[1L], 
                      lwd=lwd, col="black", lty="solid")
                ###
                for (i in seq(ncomp)) {
                    curve(pars[i,1L,drop=TRUE]*pars[i,3L,drop=TRUE]*
                          exp(-pars[i,3L,drop=TRUE]*x)+1e-7, 
                          lwd=lwd, col=colors[i], lty="solid", add=TRUE)
                } # end for.
                ###
                if (constant_d[1L]>0.0) abline(h=constant_d[1L], lty="dashed", lwd=lwd)
                ###
            } else if (typ=="lm") {
                lines(x, eval(parse(text=paste("pars[",seq(ncomp),
                      ",1L,drop=TRUE]*pars[", seq(ncomp),
                      ",3L,drop=TRUE]*(x/max(tim))*exp(-pars[",
                      seq(ncomp),",3L,drop=TRUE]*x^2/2L/max(tim))",
                      collapse="+",sep="")))+constant_d[1L]*x/max(tim), 
                      lwd=lwd, col="black", lty="solid")
                ###
                for (i in seq(ncomp)) {
                    curve(pars[i,1L,drop=TRUE]*pars[i,3L,drop=TRUE]*
                          (x/max(tim))*exp(-pars[i,3L,drop=TRUE]*x^2L/2L/max(tim))+1e-7, 
                           lwd=lwd, col=colors[i], lty="solid", add=TRUE)
                } # end for.
                ###
                if (constant_d[1L]>0.0) {
                    points(tim, constant_d[1L]*tim/max(tim), type="l", lty="dashed", lwd=lwd)
                } # end if.
            } # end if.
        } # end if. 
        ###
        ###--------------------------------------------------------------
        par(mar=c(4,4,0.5,0.5)+0.1)
        if (log=="y") {
           t_log <-  ""
        } else if (log %in% c("xy","yx")) {
            t_log <- "x"
        } else {
            t_log <- log
        } # end if.
        ###
        if (message==0L) {
            plot(x=tim, y=sig-comp.sig[,"Comp.Sum",drop=TRUE], log=t_log, las=0, lab=c(7,7,9), 
                 xlab="Stimulation time (s)", ylab="Residuals", cex.lab=1.5, 
                 xaxs="r", yaxs="r", type="o", pch=21, cex=ifelse(typ=="cw",0.75,0.5), 
                 bg="black", col="gray")
        } else {
            plot(x=tim, y=sig, log=t_log, las=0, lab=c(7,7,9),
                 xlab="Stimulation time (s)", ylab="Residuals", cex.lab=1.5, 
                 xaxs="r", yaxs="r", type="n", pch=21, cex=ifelse(typ=="cw",0.75,0.5), 
                 bg="black", col="gray")
        } # end if.
        grid(col="pink3")
        box(lwd=1L)
        ###
        ###---------------------------------------------------------------
        par(mar=c(3,0.5,7,0.5)+0.1)
        par(mgp=c(1,1,0))
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Legend", ylab="",cex.lab=1.5)
        box(lwd=1L)
        if (message==0L) {
            ###
            if (constant[1L]>0.0) {
                legend("center",
                       legend=c("Comp.Sum", paste("Comp.",seq(ncomp),sep=""), "Background"),
                       col=c("black", colors[seq(ncomp)], "black"),
                       pch=c(21, rep(NA,ncomp+1L)),
                       lty=c(rep("solid",ncomp+1L),"dashed"),
                       lwd=2.5, yjust=2, ncol=1L, cex=1.5, bty="n")
            } else {
                legend("center",
                       legend=c("Comp.Sum",paste("Comp.",seq(ncomp),sep="")),
                       col=c("black", colors[seq(ncomp)]),
                       pch=c(21, rep(NA,ncomp)),
                       lty=rep("solid", ncomp+1L),
                       lwd=2.5, yjust=2, ncol=1L, cex=1.5, bty="n")
            } # end if. 
        } else {
            legend("center",
                   legend="Data Points",
                   col="black",
                   pch=21, 
                   yjust=2, ncol=1L, cex=1.5, bty="n")
        } # end if.
        ###
        ###-----------------------------------------------------------   
        par(mar=c(6,0.5,4,0.5)+0.1)
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary", ylab="",cex.lab=1.5)
        box(lwd=1L)
        ###
        if (message==0L) {
            legend("top",
                   legend=c(paste("Curve NO: ", ifelse(is.null(curve.no),"NULL",curve.no), sep=""),
                   paste("SAR Cycle: ", ifelse(is.null(SAR.Cycle),"NULL",SAR.Cycle), sep=""),
                   paste("Irradiation dose: ", ifelse(is.null(irr.dose),"NULL",irr.dose), " (Gy|s)",sep=""),
                   "=========================",
                   "Status: OK",
                   "=========================",
                   paste("Type: ", ifelse(typ=="cw","CW","LM"), sep=""),
                   paste("ncomp: ", ncomp, sep=""),
                   paste("Add constant: ", constant, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   "=========================",
                   paste("Minimized value: ", round(value,2L), sep=""),
                   paste("Figure Of Merit: ", round(FOM,2L), " (%)", sep="")),
                   yjust=2, ncol=1L, cex=1.1, bty="n")          
        } else {
            legend("top",
                   legend=c(paste("Curve NO: ", ifelse(is.null(curve.no),"NULL",curve.no), sep=""),
                   paste("SAR Cycle: ", ifelse(is.null(SAR.Cycle),"NULL",SAR.Cycle), sep=""),
                   paste("Irradiation dose: ", ifelse(is.null(irr.dose),"NULL",irr.dose), " (Gy|s)",sep=""),
                   "=========================",
                   "Status: Model fit failed",
                   "=========================",
                   paste("Type: ", ifelse(typ=="cw","CW","LM"), sep=""),
                   paste("ncomp: ", ncomp, sep=""),
                   paste("Add constant: ", constant, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   "=========================",
                   "Minimized value: NULL",
                   "Figure Of Merit: NULL"),
                   yjust=2, ncol=1L, cex=1.1, bty="n")
        } # end if. 
        ###
        on.exit(par(mar=c(5,4,4,2)+0.1,
                    mgp=c(3,1,0),
                    mfrow=c(1L,1L)))
    } # end if. 
    ###  
    ###
    output <- list("message"=message,
                   "comp.sig"=comp.sig,
                   "LMpars"=pars, 
                   "constant"=constant_d, 
                   "value"=value, 
                   "FOM"=FOM)
    invisible(output)
} # end function decomp.default.
#####
