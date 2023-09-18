#####
calED <-
function(Curvedata, Ltx, model="gok", origin=FALSE, errMethod="sp", 
         nsim=500, weight=TRUE, trial=FALSE, plot=TRUE, nofit.rgd=NULL, 
         agID=NULL, Tn=NULL, Tn3BG=NULL, TnBG.ratio=NULL, rseTn=NULL,  
         FR=NULL, LnTn.curve=NULL, TxTn=NULL) {
    UseMethod("calED")
} ###
### 2023.09.01. 
calED.default <-
function(Curvedata, Ltx, model="gok", origin=FALSE, errMethod="sp", 
         nsim=500, weight=TRUE, trial=FALSE, plot=TRUE, nofit.rgd=NULL, 
         agID=NULL, Tn=NULL, Tn3BG=NULL, TnBG.ratio=NULL, rseTn=NULL, 
         FR=NULL, LnTn.curve=NULL, TxTn=NULL) {
    ### Stop if not.
    stopifnot(ncol(Curvedata)==3L, nrow(Curvedata)>=1L,
              length(Ltx)==2L, is.numeric(Ltx),
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(errMethod)==1L, errMethod %in% c("sp","mc"),
              length(nsim)==1L, is.numeric(nsim), nsim>=50L, nsim<=3000L,
              length(weight)==1L, is.logical(weight), 
              length(trial)==1L, is.logical(trial),
              length(plot)==1L, is.logical(plot),
              is.null(nofit.rgd) || is.numeric(nofit.rgd),
              is.null(agID) || (length(agID)==3L && is.numeric(agID)),
              is.null(Tn) || (length(Tn)==2L && is.numeric(Tn)),
              is.null(Tn3BG) || (length(Tn3BG)==1L && Tn3BG %in% c(0L,1L)),
              is.null(TnBG.ratio) || (length(TnBG.ratio)==2L && is.numeric(TnBG.ratio)),
              is.null(rseTn) || (length(rseTn)==1L && is.numeric(rseTn)),
              is.null(FR) || (length(FR)==2L && is.numeric(FR)),
              is.null(LnTn.curve) || is.list(LnTn.curve),
              is.null(TxTn) || is.numeric(TxTn))
    ###
    dose <- as.numeric(Curvedata[,1L,drop=TRUE])
    if (any(dose<0.0)) {
        stop("Error: dose value in [Curvedata] should not less than zero!")
    } # end if.
    ###
    doseltx <- as.numeric(Curvedata[,2L,drop=TRUE])
    ###
    sedoseltx <- as.numeric(Curvedata[,3L,drop=TRUE])
    if (any(sedoseltx<=0.0)) {
        stop("Error: standard error of signal in [Curvedata] should larger than zero!")
    } # end if.
    ###

    ###
    NCD <- seq(nrow(Curvedata))
    ###
    if (is.null(nofit.rgd)) {
        Curvedata1 <- Curvedata
        fitIDX <- NCD
    } else {
        Curvedata1 <- Curvedata[-nofit.rgd,,drop=FALSE]
        fitIDX <- NCD[-nofit.rgd]
    } # end if.
    ###    

    ###
    ndat <- nrow(Curvedata1)
    dose1 <- as.numeric(Curvedata1[,1L,drop=TRUE])
    doseltx1 <- as.numeric(Curvedata1[,2L,drop=TRUE])
    sedoseltx1 <- as.numeric(Curvedata1[,3L,drop=TRUE])

    ###
    if (model=="line") {
        require_npars <- 1L+!origin
    } else if (model=="exp") {
        require_npars <- 2L+!origin
    } else if (model=="lexp") {
        require_npars <- 3L+!origin
    } else if (model=="dexp") {
        require_npars <- 4L+!origin
    } else if (model=="gok") {
        require_npars <- 3L+!origin
    } # end if. 
    ### 
    n_lost <- require_npars - ndat
    no_enough_data <- FALSE
    ###
    if (model=="line" && origin==FALSE) {
        if (n_lost<=0L) { 
            model_vec <- 0L
            npars_vec <- 2L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="line" && origin==TRUE) {
        if (n_lost<=0L) {
            model_vec <- 0L
            npars_vec <- 1L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="exp" && origin==FALSE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(3L,2L) else 3L
        } else if (n_lost==1L) {
            model_vec  <- 0L 
            npars_vec  <- 2L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="exp" && origin==TRUE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(2L,1L) else 2L
        } else if (n_lost==1L) {
            model_vec <- 0L
            npars_vec <- 1L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="lexp" && origin==FALSE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(2L,7L,1L,0L) else 2L
            npars_vec <- if (trial==TRUE) c(4L,4L,3L,2L) else 4L
        } else if (n_lost==1L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(3L,2L) else 3L
        } else if (n_lost==2L) {
            model_vec <- 0L
            npars_vec <- 2L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="lexp" && origin==TRUE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(2L,7L,1L,0L) else 2L
            npars_vec <- if (trial==TRUE) c(3L,3L,2L,1L) else 3L
        } else if (n_lost==1L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(2L,1L) else 2L
        } else if (n_lost==2L) {
            model_vec <- 0L
            npars_vec <- 1L
        } else {
            no_enough_data <- TRUE 
        } # end if.
    } else if (model=="dexp" && origin==FALSE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(3L,7L,1L,0L) else 3L
            npars_vec <- if (trial==TRUE) c(5L,4L,3L,2L) else 5L
        } else if (n_lost==1L) {
            model_vec <- if (trial==TRUE) c(7L,1L,0L) else 7L
            npars_vec <- if (trial==TRUE) c(4L,3L,2L) else 4L
        } else if (n_lost==2L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(3L,2L) else 3L
        } else if (n_lost==3L) {
            model_vec <- 0L
            npars_vec <- 2L
        } else {
            no_enough_data <- TRUE
        } # end if. 
    } else if (model=="dexp" && origin==TRUE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(3L,7L,1L,0L) else 3L
            npars_vec <- if (trial==TRUE) c(4L,3L,2L,1L) else 4L
        } else if (n_lost==1L) {
            model_vec <- if (trial==TRUE) c(7L,1L,0L) else 7L
            npars_vec <- if (trial==TRUE) c(3L,2L,1L) else 3L
        } else if (n_lost==2L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(2L,1L) else 2L
        } else if (n_lost==3L) {
            model_vec <- 0L
            npars_vec <- 1L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="gok" && origin==FALSE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(7L,1L,0L) else 7L
            npars_vec <- if (trial==TRUE) c(4L,3L,2L) else 4L
        } else if (n_lost==1L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(3L,2L) else 3L
        } else if (n_lost==2L) {
            model_vec <- 0L
            npars_vec <- 2L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } else if (model=="gok" && origin==TRUE) {
        if (n_lost<=0L) {
            model_vec <- if (trial==TRUE) c(7L,1L,0L) else 7L
            npars_vec <- if (trial==TRUE) c(3L,2L,1L) else 3L
        } else if (n_lost==1L) {
            model_vec <- if (trial==TRUE) c(1L,0L) else 1L
            npars_vec <- if (trial==TRUE) c(2L,1L) else 2L
        } else if (n_lost==2L) {
            model_vec <- 0L
            npars_vec <- 1L
        } else {
            no_enough_data <- TRUE
        } # end if.
    } # end if.
    ###
    if (no_enough_data==TRUE)  {
        stop("Error: data points is not enough for model optimization!")
    } # end if.
    ###
    if (n_lost>0L) {
        if (!is.null(agID)) {
            cat(paste("[NO=", agID[1L], ",Position=",agID[2L], ",Grain=",agID[3L],"]: ", ndat, 
                      " data points are not enough for fitting the ", 
                      model, " (origin=",origin,") model!\n",sep=""))
        } else {
            cat(paste("NOTE: ", ndat, " data points are not enough for fitting the ", 
                      model, " (origin=",origin,") model!\n",sep=""))
        } # end if.
    } # end if.
    ###
    ###
    inltx <- matrix(Ltx, nrow=1L, ncol=2L) 
    outDose <- matrix(0, nrow=1L, ncol=2L)
    eemm <- ifelse(errMethod=="sp", 0L, 1L)
    mcED <- vector(length=nsim)
    uw <- ifelse(weight==FALSE, 0L, 1L)
    fvec1 <- vector(length=ndat)
    fmin <- 0.0
    saturateDose <- -99.0
    acceptRate <- 0.0
    message <- 0
    ###
    n_loop <- 0L
    max_loop <- length(model_vec)
    ###
    repeat {
        n_loop <- n_loop + 1L
        mdl <- model_vec[n_loop]
        n2 <- npars_vec[n_loop]
        pars <- stdp <- vector(length=n2)
        ###
        if (mdl==7L) {
            res <- .Fortran("calED1",as.double(dose1),as.double(doseltx1),as.double(sedoseltx1),
                       as.integer(ndat),as.integer(n2),as.double(inltx),outDose=as.double(outDose),
                       as.integer(eemm),mcED=as.double(mcED),pars=as.double(pars),stdp=as.double(stdp),
                       as.integer(uw),as.integer(nsim),fvec1=as.double(fvec1),fmin=as.double(fmin),
                       saturateDose=as.double(saturateDose),acceptRate=as.double(acceptRate),
                       message=as.integer(message),PACKAGE="numOSL")
        } else if (mdl %in% c(0L,1L,2L,3L)) {
            res <- .Fortran("calED_fort",as.double(dose1),as.double(doseltx1),as.double(sedoseltx1),
                       as.integer(ndat),as.integer(n2),as.double(inltx),outDose=as.double(outDose),
                       as.integer(eemm),mcED=as.double(mcED),pars=as.double(pars),stdp=as.double(stdp),
                       as.integer(mdl),as.integer(uw),as.integer(nsim),fvec1=as.double(fvec1),
                       fmin=as.double(fmin),saturateDose=as.double(saturateDose),
                       acceptRate=as.double(acceptRate),message=as.integer(message),
                       PACKAGE="numOSL")
        } # end if.
        ###
        if (res$message!=1L) break
        if (n_loop==max_loop) break
    } # end repeat.
    ###
    message <- res$message
    if (message %in% c(0L,4L)) ED <- res$outDose[1L] else ED <- NA
    if (message==0L) seED <- res$outDose[2L] else  seED <- NA
    if (message==0L && errMethod=="mc") mcED <- res$mcED else  mcED <- NULL 
    if (message==0L) {
        if (errMethod=="mc") {
            ConfInt <- quantile(mcED, probs=c(0.16,0.84,0.025,0.975)) 
        } # end if.
        if (errMethod=="sp") {
            ConfInt <- vector(length=4L)
            ConfInt[1L] <- ED - 0.9944579*seED
            ConfInt[2L] <- ED + 0.9944579*seED
            ConfInt[3L] <- ED - 1.959964*seED
            ConfInt[4L] <- ED + 1.959964*seED
        } # end if.
        names(ConfInt) <- c("lower68", "upper68", "lower95", "upper95")
    } else {
        ConfInt <- NA
    } # end if.
    ###
    ###
    if (message!=1L) {
        min_obj <- res$fmin
        avg_fit_error <- sqrt(sum((doseltx1-res$fvec1)^2L))/ndat
        FOM <- 100.0*sum(abs(doseltx1-res$fvec1))/sum(res$fvec1)
        if ((ndat-n2)>0L) RCS <- res$fmin/(ndat-n2) else RCS <- NA
        ###
        LMpars <- cbind(res$pars,res$stdp)
        colnames(LMpars) <- c("Pars","sePars")
        rownames(LMpars) <- (c("a","b","c","d","e"))[seq(n2)]
    } else {
        min_obj <- avg_fit_error <- FOM <- RCS <- NA
        LMpars <- NA
    } # end if.
    ###
    res_calRcyRcp <- calRcyRcp(Curvedata, Ltx)
    ###
    RecyclingRatio1 <- res_calRcyRcp$RecyclingRatio1[1L]
    seRecyclingRatio1 <- res_calRcyRcp$RecyclingRatio1[2L]
    ###
    RecyclingRatio2 <- res_calRcyRcp$RecyclingRatio2[1L]
    seRecyclingRatio2 <- res_calRcyRcp$RecyclingRatio2[2L]
    ###
    RecyclingRatio3 <- res_calRcyRcp$RecyclingRatio3[1L]
    seRecyclingRatio3 <- res_calRcyRcp$RecyclingRatio3[2L]
    ###
    Recuperation1 <- res_calRcyRcp$Recuperation1[1L]
    seRecuperation1 <- res_calRcyRcp$Recuperation1[2L]
    ###
    Recuperation2 <- res_calRcyRcp$Recuperation2[1L]
    seRecuperation2 <- res_calRcyRcp$Recuperation2[2L]
    ### 
    if (mdl==0L) {
        model <- "line"
    } else if (mdl==1L) {
        model <- "exp"
    } else if (mdl==2L) {
        model <- "lexp"
    } else if (mdl==3L) {
        model <- "dexp"
    } else if (mdl==7L) {
        model <- "gok"
    } # end if.
    ###
    maxDose <- max(dose)
    pars <- res$pars
    cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
    ###
    if (message!=1L) {
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
        calED_method <- ifelse(Ltx[1L]>=maxLtx, "Extrapolation", "Interpolation")
    } else {
        calED_method <- "NULL"
    } # end if.
    ###
    ###
    saturateDose <- res$saturateDose
    acceptRate <- res$acceptRate
    ###
    ###-----------------------------------------------------------------------------
    if (plot==TRUE) {
        opar <- par("mfrow", "mgp", "mar")
        on.exit(par(opar))
        ###
        layout(matrix(c(1L,1L,2L,1L,1L,3L,4L,4L,4L),nrow=3L), respect=TRUE)
        par(mgp=c(2.5,1,0))
        ###
        if (message %in% c(0L,4L)) {
            # message=0L: ED and Error calculation succeeded.
            # message=4L: ED Error calculation failed.
            lowerX <- min(dose,ED,0)*1.2
            upperX <- max(dose,ED)*1.2
            lowerY <- min(doseltx,Ltx[1L],0)*1.2
            upperY <- max(doseltx,Ltx[1L])*1.2
        } else if (message %in% c(1L,3L)) {
            # message=1L: growth curve fitting failed.
            # message=3L: ED calculation failed.
            lowerX <- min(dose,0)*1.2
            upperX <- max(dose)*1.2
            lowerY <- min(doseltx,Ltx[1L],0)*1.2
            upperY <- max(doseltx,Ltx[1L])*1.2
        } else if (message==2L) {
            ### message=2L: natural signal saturated. 
            lowerX <- 0.0
            upperX <- max(dose,saturateDose)*1.5
            lowerY <- min(doseltx,Ltx[1L],0)*1.2
            upperY <- max(doseltx,Ltx[1L])*1.2
        } # end if.
        ###
        par(mar=c(4,4,2,0.5)+0.1)
        plot(NA, NA, main=NULL, xlab="Regenerative dose (<Gy>|<s>)", 
             ylab="Sensitivity-corrected OSL",
             las=0, cex.lab=1.5, cex.main=1.5, xlim=c(lowerX,upperX),
             ylim=c(lowerY,upperY), xaxs="i", yaxs="i", lab=c(7,7,9))
        ###
        points(dose, doseltx, pch=21, cex=1.5, bg="black")
        ###
        arrowIndex <- which(sedoseltx/doseltx>0.001)
        if (length(arrowIndex)>=1L) {
            suppressWarnings(arrows(x0=dose[arrowIndex], y0=doseltx[arrowIndex]-sedoseltx[arrowIndex]/2.0, 
                x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sedoseltx[arrowIndex]/2.0,
                code=3, lwd=1, angle=90, length=0.05, col="black"))
        } # end if.
        ###
        ###
        if (message!=1L) {
            ### message!=1L: growth curve fitting succeeded.
            x<-NULL
            if (model=="line") {
                curve(pars[1L]*x+cst, type="l", add=TRUE, 
                      lwd=2.0, from=lowerX, to=upperX, col="skyblue")
            } else if (model=="exp") {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, type="l", add=TRUE, 
                      lwd=2.0, from=lowerX, to=upperX, col="skyblue")
            } else if (model=="lexp")  {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, type="l", 
                      add=TRUE, lwd=1.5, from=lowerX, to=upperX, col="skyblue")
            } else if (model=="dexp") {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                      type="l", add=TRUE, lwd=2.0, from=lowerX, to=upperX, col="skyblue")
            } else if (model=="gok") {
                curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                     type="l", add=TRUE, lwd=2.0, from=lowerX, to=upperX, col="skyblue")
            } # end if.
        } # end if.
        ###
        if (message==0L) {
            ### message=0L: ED and Error calculation succeeded.
            points(x=ED[1L], y=Ltx[1L], pch=23, cex=1.5, bg="grey")
            ###
            if (Ltx[1L]>0 && errMethod=="mc") {
                dmcED <- density(mcED)
                dxy <- cbind(dmcED$x,dmcED$y)
                dxy[,2L] <- (dxy[,2L,drop=TRUE]-min(dxy[,2L,drop=TRUE]))/
                            (max(dxy[,2L,drop=TRUE])-min(dxy[,2L,drop=TRUE]))*Ltx[1L]*0.8
                polygon(dxy, col="grey")  
            } # end if
            ###  
            lines(x=c(0,ED,ED), y=c(Ltx[1L],Ltx[1L],0), lty="dashed", lwd=1.5) 
            ###     
            if (seED/ED>0.001) {
                suppressWarnings(arrows(x0=ED-seED/2.0, y0=Ltx[1L],
                    x1=ED+seED/2.0, y1=Ltx[1L], code=3, lwd=1, angle=90, length=0.05, col="black"))
            } # end if.
            ###           
        } else if (message==2L) {
            ### message=2L: natural signal saturated.
            abline(h=Ltx[1L], lty="dashed", col="red", lwd=1.5)
            ###
        } else if (message==4L) {
            ### message=4L: ED Error calculation failed.
            points(x=ED[1L], y=Ltx[1L], pch=23, cex=1.5, bg="grey")
            lines(x=c(0,ED,ED), y=c(Ltx[1L],Ltx[1L],0), lty="dashed", lwd=1.5) 
        } # end if.
        ###
        grid(col="pink3")
        box(lwd=1)
        ###
        ###-------------------------------------------------------------------------------------- 
        par(mar=c(4,4,0.5,0.5)+0.1)
        if (!is.null(LnTn.curve) && length(LnTn.curve[["Ln.x"]])>1L) {
            x_max <- max(max(LnTn.curve[["Ln.x"]]), max(LnTn.curve[["Tn.x"]]), na.rm=TRUE)
            y_max <- max(max(LnTn.curve[["Ln.y"]]), max(LnTn.curve[["Tn.y"]]), na.rm=TRUE)
            ###
            plot(x=LnTn.curve[["Ln.x"]], y=LnTn.curve[["Ln.y"]], type="l", lwd=1.5, col="blue",
                 main=NULL, xlim=c(0, x_max), ylim=c(0, y_max), xlab="Stimulation time (s)", 
                 ylab="Photon counts", las=0, xaxt="n", cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
            ###
            if (length(LnTn.curve[["Tn.x"]])>1L) {
                # Both Tn.x and Tn.y are not equal to NA of length 1.
                points(x=LnTn.curve[["Tn.x"]], y=LnTn.curve[["Tn.y"]], type="l", lwd=1.5, col="red")
            } # end if.
            ###
            x_axis_location <- axTicks(side=1L)
            axis(side=1L, at=x_axis_location, labels=as.character(x_axis_location))
            legend("topright", legend=c("Ln decay curve", "Tn decay curve"), col=c("blue","red"), 
                   lwd=1.5, yjust=2, ncol=1L, cex=1.0, bty="n")
        } else {
            plot(x=1L, y=1.0, type="n", lwd=1.5,
                 main=NULL, xlab="Stimulation time (s)", ylab="Photon counts", las=0, 
                 xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
            axis(side=1L, at=c(0.7, 1, 1.3), labels=c("x1", "x2", "x3"))
            axis(side=2L, at=c(0.7, 1, 1.3), labels=c("y1", "y2", "y3"))
        } # end if.
        ###
        ###--------------------------------------------------------------------------------------
        par(mar=c(4,4,0.5,0.5)+0.1)
        if (!is.null(TxTn) && all(is.finite(TxTn))) {
            plot(x=seq(length(TxTn)), y=TxTn, type="p", pch=21, bg="blue", cex=2.0, 
                 main=NULL, xlab="SAR Cycle", ylab="Tx/Tn", las=0, xaxt="n", 
                 ylim=c(0, max(TxTn)*1.5), cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
            axis(side=1L, at=seq(length(TxTn)), labels=as.character(seq(length(TxTn))))
            abline(h=1.0, lty="dashed", col="red", lwd=1.5)
        } else {
            plot(x=1L, y=1.0, type="n", pch=21, bg="blue", cex=2.0, 
                 main=NULL, xlab="SAR Cycle", ylab="Tx/Tn", las=0, 
                 xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
            axis(side=1L, at=c(0.7, 1, 1.3), labels=c("x1", "x2", "x3"))
            axis(side=2L, at=c(0.7, 1, 1.3), labels=c("y1", "y2", "y3"))
        } # end if.
        ###
        ###--------------------------------------------------------------------------------------
        par(mar=c(11,0.5,11,0.5)+0.1)
        par(mgp=c(1,1,0))
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary", ylab="",cex.lab=1.5)
        ###
        if (is.null(agID)) {
            NO_Position_Grain <- "NULL"
        } else {
            NO_Position_Grain <- paste("[NO=",agID[1L],", Position=",agID[2L],", Grain=",agID[3L],"]",sep="")
        } # end if
        ###
        character_model_vec <- vector(length=n_loop)
        for (i in seq(n_loop)) {
            if (model_vec[i]==0L) character_model_vec[i] <- "line"
            if (model_vec[i]==1L) character_model_vec[i] <- "exp"
            if (model_vec[i]==2L) character_model_vec[i] <- "lexp"
            if (model_vec[i]==3L) character_model_vec[i] <- "dexp"
            if (model_vec[i]==7L) character_model_vec[i] <- "gok"
        } # end for.
        ###
        if (message==0L) {
            ### message=0L: ED and Error calculation succeeded.
            legend("center", 
                   legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                   "=========================",
                   "Status: OK",
                   "=========================",
                   paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                         " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                   paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""),
                   paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                         " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                   paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""),
                   paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                         " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                   "=========================",
                   paste("Recycling ratio 1: ", round(RecyclingRatio1,2L), 
                         " +/- ", round(seRecyclingRatio1,2L), sep=""),
                   paste("Recycling ratio 2: ", round(RecyclingRatio2,2L), 
                         " +/- ", round(seRecyclingRatio2,2L), sep=""),  
                   paste("Recycling ratio 3: ", round(RecyclingRatio3,2L), 
                         " +/- ", round(seRecyclingRatio3,2L), sep=""),
                   "=========================",
                   paste("Recuperation 1: ", round(Recuperation1,2L), 
                         " +/- ", round(seRecuperation1,2L), " (%)", sep=""),
                   paste("Recuperation 2: ", round(Recuperation2,2L), 
                         " +/- ", round(seRecuperation2,2L), " (%)",sep=""),
                   "=========================",
                   paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                   paste("Pass origin: ", origin, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   paste("Minimized value: ", round(min_obj,2L),sep=""),
                   paste("Average error in fit: ", round(avg_fit_error,2L),sep=""), 
                   paste("Reduced Chi-Square: ", round(RCS,2L),sep=""),
                   paste("Figure Of Merit: ", round(FOM,2L)," (%)",sep=""),
                   "=========================",
                   paste("calED method: ",calED_method, sep=""),
                   paste("MC acceptance-rate: ", ifelse(errMethod=="sp", "NULL",
                         round(acceptRate,2L)), " (%)", sep=""),
                   paste("ED: ",round(ED,2L), " +/- ",round(seED,2L)," (<Gy>|<s>)",sep=""),
                   paste("RSE of ED: ",round(seED/abs(ED)*100.0,2L)," (%)",sep=""),
                   paste("95% interval: [",round(ConfInt[3L],2L),", ",round(ConfInt[4L],2L),"]", sep=""), 
                   paste("68% interval: [",round(ConfInt[1L],2L),", ",round(ConfInt[2L],2L),"]", sep="")),                    
                   yjust=2, ncol=1L, cex=0.9, bty="n")
        } else if (message==1L) {
            ### message=1L: Growth curve fitting failed.
            legend("center", 
                   legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                   "=========================",
                   "Status: model fit failed",
                   "=========================",
                   paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                         " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                   paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""), 
                   paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                         " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                   paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""), 
                   paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                         " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                   "=========================",
                   paste("Recycling ratio 1: ", round(RecyclingRatio1,2L), 
                         " +/- ", round(seRecyclingRatio1,2L), sep=""),
                   paste("Recycling ratio 2: ", round(RecyclingRatio2,2L), 
                         " +/- ", round(seRecyclingRatio2,2L), sep=""),  
                   paste("Recycling ratio 3: ", round(RecyclingRatio3,2L), 
                         " +/- ", round(seRecyclingRatio3,2L), sep=""),
                   "=========================",
                   paste("Recuperation 1: ", round(Recuperation1,2L), 
                         " +/- ", round(seRecuperation1,2L), " (%)", sep=""),
                   paste("Recuperation 2: ", round(Recuperation2,2L), 
                         " +/- ", round(seRecuperation2,2L), " (%)",sep=""),  
                   "=========================",
                   paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                   paste("Pass origin: ", origin, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   "Minimized value: NULL",
                   "Average error in fit: NULL", 
                   "Reduced Chi-Square: NULL",
                   "Figure Of Merit: NULL (%)",
                   "=========================",
                   paste("calED method: ",calED_method, sep=""),
                   "MC acceptance-rate: NULL (%)",
                   paste("ED: ",NA, " +/- ",NA," (<Gy>|<s>)",sep=""),
                   "RSE of ED: NA (%)",
                   "95% interval: [NA, NA]", 
                   "68% interval: [NA, NA]"),                      
                   yjust=2, ncol=1L, cex=0.9, bty="n")
        } else if (message==2L) {
            ### message=2L: natural signal saturated.
            legend("center", 
                   legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                   "=========================",
                   "Status: Ln/Tn saturated",
                   "=========================",
                   paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                         " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                   paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""),
                   paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                         " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                   paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""),
                   paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                         " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                   "=========================",
                   paste("Recycling ratio 1: ", round(RecyclingRatio1,2L), 
                         " +/- ", round(seRecyclingRatio1,2L), sep=""),
                   paste("Recycling ratio 2: ", round(RecyclingRatio2,2L), 
                         " +/- ", round(seRecyclingRatio2,2L), sep=""),  
                   paste("Recycling ratio 3: ", round(RecyclingRatio3,2L), 
                         " +/- ", round(seRecyclingRatio3,2L), sep=""),
                   "=========================",
                   paste("Recuperation 1: ", round(Recuperation1,2L), 
                         " +/- ", round(seRecuperation1,2L), " (%)", sep=""),
                   paste("Recuperation 2: ", round(Recuperation2,2L), 
                         " +/- ", round(seRecuperation2,2L), " (%)",sep=""), 
                   "=========================",
                   paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                   paste("Pass origin: ", origin, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   paste("Minimized value: ", round(min_obj,2L),sep=""),
                   paste("Average error in fit: ", round(avg_fit_error,2L),sep=""), 
                   paste("Reduced Chi-Square: ", round(RCS,2L),sep=""),
                   paste("Figure Of Merit: ", round(FOM,2L)," (%)",sep=""),
                   "=========================",
                   paste("calED method: ",calED_method, sep=""),
                   "MC acceptance-rate: NULL (%)",
                   paste("ED: ",Inf, " +/- ",NA," (<Gy>|<s>)",sep=""),
                   "RSE of ED: NA (%)",
                   "95% interval: [-Inf, +Inf]", 
                   "68% interval: [-Inf, +Inf]"),                                
                   yjust=2, ncol=1L, cex=0.9, bty="n")
        } else if (message==3L) {
            ### message=3L: ED calculation failed.
            legend("center", 
                  legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                  "=========================",
                  "Status: ED failed (ED < -50 Gy|s)",
                  "=========================",
                  paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                        " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                  paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""),
                  paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                        " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                  paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""), 
                  paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                        " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                  "=========================",
                  paste("Recycling ratio 1: ", round(RecyclingRatio1,2L), 
                        " +/- ", round(seRecyclingRatio1,2L), sep=""),
                  paste("Recycling ratio 2: ", round(RecyclingRatio2,2L), 
                        " +/- ", round(seRecyclingRatio2,2L), sep=""),  
                  paste("Recycling ratio 3: ", round(RecyclingRatio3,2L), 
                        " +/- ", round(seRecyclingRatio3,2L), sep=""),
                  "=========================",
                  paste("Recuperation 1: ", round(Recuperation1,2L), 
                        " +/- ", round(seRecuperation1,2L), " (%)", sep=""),
                  paste("Recuperation 2: ", round(Recuperation2,2L), 
                        " +/- ", round(seRecuperation2,2L), " (%)",sep=""), 
                  "=========================",
                  paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                  paste("Pass origin: ", origin, sep=""),
                  paste("Weighted fit: ", weight, sep=""),
                  paste("Minimized value: ", round(min_obj,2L),sep=""),
                  paste("Average error in fit: ", round(avg_fit_error,2L),sep=""), 
                  paste("Reduced Chi-Square: ", round(RCS,2L),sep=""),
                  paste("Figure Of Merit: ", round(FOM,2L)," (%)",sep=""),
                  "=========================",
                  paste("calED method: ",calED_method, sep=""),
                  "MC acceptance-rate: NULL (%)",
                  paste("ED: ",NA, " +/- ",NA," (<Gy>|<s>)",sep=""),
                  "RSE of ED: NA (%)",
                  "95% interval: [NA, NA]", 
                  "68% interval: [NA, NA]"),                               
                  yjust=2, ncol=1L, cex=0.9, bty="n")
        } else if (message==4L) {
            ### message=4L: ED Error calculation failed.
            legend("center", 
                  legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                  "=========================",
                  ifelse(errMethod=="mc",
                  "Status: ED Error failed (MC AR < 1%)",
                  "Status: ED Error failed (infinite upper ED)"),
                  "=========================",
                  paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                        " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                  paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""),
                  paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                        " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                  paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""),
                  paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                        " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                  "=========================",
                  paste("Recycling ratio 1: ", round(RecyclingRatio1,2L), 
                        " +/- ", round(seRecyclingRatio1,2L), sep=""),
                  paste("Recycling ratio 2: ", round(RecyclingRatio2,2L), 
                        " +/- ", round(seRecyclingRatio2,2L), sep=""),  
                  paste("Recycling ratio 3: ", round(RecyclingRatio3,2L), 
                        " +/- ", round(seRecyclingRatio3,2L), sep=""),
                  "=========================",
                  paste("Recuperation 1: ", round(Recuperation1,2L), 
                        " +/- ", round(seRecuperation1,2L), " (%)", sep=""),
                  paste("Recuperation 2: ", round(Recuperation2,2L), 
                        " +/- ", round(seRecuperation2,2L), " (%)",sep=""), 
                  "=========================",
                  paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                  paste("Pass origin: ", origin, sep=""),
                  paste("Weighted fit: ", weight, sep=""),
                  paste("Minimized value: ", round(min_obj,2L),sep=""),
                  paste("Average error in fit: ", round(avg_fit_error,2L),sep=""),
                  paste("Reduced Chi-Square: ", round(RCS,2L),sep=""),
                  paste("Figure Of Merit: ", round(FOM,2L)," (%)",sep=""),
                  "=========================", 
                  paste("calED method: ",calED_method, sep=""),
                  "MC acceptance-rate: <1 (%)",         
                  paste("ED: ",round(ED,2L), " +/- ",NA," (<Gy>|<s>)",sep=""),
                  "RSE of ED: NA (%)",
                  "95% interval: [NA, NA]", 
                  "68% interval: [NA, NA]"),                          
                  yjust=2, ncol=1L, cex=0.9, bty="n")
        } # end if.
        ###
    } # end if.
    ###--------------------------------------------------------------------------------------
    ###
    output<-list("message"=message,
                 "fitIDX"=fitIDX,
                 "LMpars"=LMpars, 
                 "value"=min_obj,
                 "avg.error"=avg_fit_error,
                 "RCS"=RCS,
                 "FOM"=FOM, 
                 "calED.method"=calED_method,
                 "mcED"=mcED,
                 "ED"=c("ED"=ED, "seED"=seED),
                 "ConfInt"=ConfInt,
                 "RecyclingRatio1"=res_calRcyRcp$RecyclingRatio1,
                 "RecyclingRatio2"=res_calRcyRcp$RecyclingRatio2,
                 "RecyclingRatio3"=res_calRcyRcp$RecyclingRatio3,
                 "Recuperation1"=res_calRcyRcp$Recuperation1,
                 "Recuperation2"=res_calRcyRcp$Recuperation2)
    ###
    invisible(output)
} # end function calED.default.
#####
