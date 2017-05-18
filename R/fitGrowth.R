#####
fitGrowth <-
function(Curvedata, model="gok", origin=FALSE, weight=TRUE,
         trial=FALSE, plot=TRUE, agID=NULL, Tn=NULL, Tn3BG=NULL, 
         TnBG.ratio=NULL, rseTn=NULL, FR=NULL, RecyclingRatio1=NULL, 
         RecyclingRatio2=NULL, RecyclingRatio3=NULL, Recuperation1=NULL, 
         Recuperation2=NULL, LnTn.curve=NULL, TxTn=NULL) {
    UseMethod("fitGrowth")
} #
### 2017.05.18.
fitGrowth.default <-
function(Curvedata, model="gok", origin=FALSE, weight=TRUE,
         trial=FALSE, plot=TRUE, agID=NULL, Tn=NULL, Tn3BG=NULL, 
         TnBG.ratio=NULL, rseTn=NULL, FR=NULL, RecyclingRatio1=NULL, 
         RecyclingRatio2=NULL, RecyclingRatio3=NULL, Recuperation1=NULL, 
         Recuperation2=NULL, LnTn.curve=NULL, TxTn=NULL) {
    ### Stop if not.
    stopifnot(ncol(Curvedata)==3L, nrow(Curvedata)>=1L,
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(weight)==1L, is.logical(weight),
              length(trial)==1L, is.logical(trial),
              length(plot)==1L, is.logical(plot),
              is.null(agID) || (length(agID)==3L && is.numeric(agID)),
              is.null(Tn) || (length(Tn)==2L && is.numeric(Tn)),
              is.null(Tn3BG) || (length(Tn3BG)==1L && Tn3BG %in% c(0L,1L)),
              is.null(TnBG.ratio) || (length(TnBG.ratio)==2L && is.numeric(TnBG.ratio)),
              is.null(rseTn) || (length(rseTn)==1L && is.numeric(rseTn)),
              is.null(FR) || (length(FR)==2L && is.numeric(FR)),
              is.null(RecyclingRatio1) || (length(RecyclingRatio1)==2L && is.numeric(RecyclingRatio1)),
              is.null(RecyclingRatio2) || (length(RecyclingRatio2)==2L && is.numeric(RecyclingRatio2)),
              is.null(RecyclingRatio3) || (length(RecyclingRatio3)==2L && is.numeric(RecyclingRatio3)),
              is.null(Recuperation1) || (length(Recuperation1)==2L && is.numeric(Recuperation1)),
              is.null(Recuperation2) || (length(Recuperation2)==2L && is.numeric(Recuperation2)),
              is.null(LnTn.curve) || is.list(LnTn.curve), 
              is.null(TxTn) || is.numeric(TxTn))
    ###
    dose <- as.numeric(Curvedata[,1L,drop=TRUE])
    if (any(dose<0.0)) stop("Error: dose value in [Curvedata] should not less than zero!")
    ###
    doseltx <- as.numeric(Curvedata[,2L,drop=TRUE])
    ###
    sdoseltx <- as.numeric(Curvedata[,3L,drop=TRUE])
    if (any(sdoseltx<=0.0)) stop("Error: standard error of signal in [Curvedata] should larger than zero!")
    ###
    ndat <- nrow(Curvedata)
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
    if (no_enough_data==TRUE)  stop("Error: data points is not enough for model optimization!")
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
    uw <- ifelse(weight==FALSE, 0L, 1L)
    fvec1 <- vector(length=ndat)
    fmin <- 0.0
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
            res<-.Fortran("fitGOK",as.double(dose),as.double(doseltx),as.double(sdoseltx),
                          as.integer(ndat),as.integer(n2),pars=as.double(pars),stdp=as.double(stdp),
                          as.integer(uw),fvec1=as.double(fvec1),fmin=as.double(fmin),
                          message=as.integer(message),PACKAGE="numOSL")
            ###
        } else if (mdl %in% c(0L,1L,2L,3L)) {
            res<-.Fortran("fitGrowth_fort",as.double(dose),as.double(doseltx),as.double(sdoseltx),
                          as.integer(ndat),as.integer(n2),pars=as.double(pars),stdp=as.double(stdp),
                          as.integer(mdl),as.integer(uw),fvec1=as.double(fvec1),fmin=as.double(fmin),
                          message=as.integer(message),PACKAGE="numOSL")
        } # end if.
        if (res$message==0L) break
        if (n_loop==max_loop) break
    } # end repeat. 
    ###
    message <- res$message
    ###
    if (message==0L) {
        min_obj <- res$fmin
        avg_fit_error <- sqrt(sum((doseltx-res$fvec1)^2L))/ndat
        FOM <- 100.0*sum(abs(doseltx-res$fvec1))/sum(res$fvec1)
        if ((ndat-n2-1L)>0L) RCS <- res$fmin/(ndat-n2-1L) else RCS <- NA
        ###
        LMpars <- cbind(res$pars,res$stdp)
        colnames(LMpars) <- c("Pars","sePars")
        rownames(LMpars) <- (c("a","b","c","d","e"))[seq(n2)]
    } else {
        min_obj <- avg_fit_error <- FOM <- RCS <- NA
        LMpars <- NULL
    } # end if.
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
    ###-----------------------------------------------------------------
    if (plot==TRUE) {
        layout(matrix(c(1L,1L,2L,1L,1L,3L,4L,4L,4L),nrow=3L), respect=TRUE)
        par(mgp=c(2.5,1,0))
        ###      
        lowerX <- min(dose,0)*1.2
        upperX <- max(dose)*1.2
        lowerY <- min(doseltx,0)*1.2
        upperY <- max(doseltx)*1.2
        ###
        par(mar=c(4,4,2,0.5)+0.1)
        plot(NA, NA, main=NULL, xlab="Regenerative dose (Gy|s)", 
             ylab="Sensitivity-corrected OSL",
             las=0, cex.lab=1.5, cex.main=1.5, xlim=c(lowerX,upperX),
             ylim=c(lowerY,upperY), xaxs="i", yaxs="i", lab=c(7,7,9))
        ### 
        points(dose, doseltx, pch=21, cex=1.5, bg="black")
        ###
        arrowIndex <- which(sdoseltx/doseltx>0.001)
        if (length(arrowIndex)>=1L) {
            suppressWarnings(arrows(x0=dose[arrowIndex], y0=doseltx[arrowIndex]-sdoseltx[arrowIndex]/2.0, 
                x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sdoseltx[arrowIndex]/2.0,
                code=3, lwd=1, angle=90, length=0.05, col="black"))
        } # end if.
        ###
        if (message==0L) { 
            pars <- res$pars
            se_pars <- res$stdp
            x <- NULL
            cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
            ###
            if(model=="line") {
                curve(pars[1L]*x+cst, from=0,
                      type="l", add=TRUE, lwd=2, col="skyblue")
            } else if(model=="exp") {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, from=0,
                      type="l", add=TRUE, lwd=2, col="skyblue")
            } else if(model=="lexp")  {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, 
                      from=0, type="l", add=TRUE, lwd=2, col="skyblue")
            } else if(model=="dexp") {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                      from=0, type="l", add=TRUE, lwd=2, col="skyblue")
            } else if(model=="gok") {
                curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                      from=0, type="l", add=TRUE, lwd=2, col="skyblue")
            } # end if.
            ###
        } # end if.
        ###
        grid(col="pink3")
        box(lwd=1)
        ###
        ###--------------------------------------------------------------
        par(mar=c(4,4,0.5,0.5)+0.1)
        if (!is.null(LnTn.curve)) {
            x_max <- max(max(LnTn.curve[["Ln.x"]]), max(LnTn.curve[["Tn.x"]]), na.rm=TRUE)
            y_max <- max(max(LnTn.curve[["Ln.y"]]), max(LnTn.curve[["Tn.y"]]), na.rm=TRUE)
            ###
            plot(x=LnTn.curve[["Ln.x"]], y=LnTn.curve[["Ln.y"]], type="l", lwd=1.5, col="blue",
                 main=NULL, xlim=c(0, x_max), ylim=c(0, y_max), xlab="Stimulation time (s)", 
                 ylab="Photon counts", las=0, xaxt="n", cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
            ###
            if (length(LnTn.curve[["Tn.x"]])>1L) {
                ### Both Tn.x and Tn.y are not equal to NA of length 1.
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
        ###--------------------------------------------------------------
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
        ###--------------------------------------------------------------
        par(mar=c(10,0.5,10,0.5)+0.1)
        par(mgp=c(1,1,0))
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary", ylab="", cex.lab=1.5)
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
            ### message=0L: Growth curve fitting succeeded.
            legend("center", 
                   legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                   "========================",
                   "Status: OK",
                   "========================",
                   paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                         " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                   paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""),
                   paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                         " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                   paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""), 
                   paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                         " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                   "========================",
                   paste("Recycling ratio 1: ", if(!is.null(RecyclingRatio1)) paste(round(RecyclingRatio1[1L],2L), 
                         " +/- ", round(RecyclingRatio1[2L],2L), sep="") else "NULL", sep=""),
                   paste("Recycling ratio 2: ", if(!is.null(RecyclingRatio2)) paste(round(RecyclingRatio2[1L],2L), 
                         " +/- ", round(RecyclingRatio2[2L],2L), sep="") else "NULL", sep=""), 
                   paste("Recycling ratio 3: ", if(!is.null(RecyclingRatio3)) paste(round(RecyclingRatio3[1L],2L), 
                         " +/- ", round(RecyclingRatio3[2L],2L), sep="") else "NULL", sep=""),
                   "========================",
                   paste("Recuperation 1: ", if(!is.null(Recuperation1)) paste(round(Recuperation1[1L],2L), 
                         " +/- ", round(Recuperation1[2L],2L), " (%)", sep="") else "NULL (%)", sep=""),
                   paste("Recuperation 2: ", if(!is.null(Recuperation2)) paste(round(Recuperation2[1L],2L), 
                         " +/- ", round(Recuperation2[2L],2L), " (%)", sep="") else "NULL (%)", sep=""), 
                   "========================",
                   paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                   paste("Pass origin: ", origin, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   paste("Minimized value: ", round(min_obj,2L),sep=""),
                   paste("Average error in fit: ", round(avg_fit_error,2L),sep=""),
                   paste("Reduced Chi-Square: ", round(RCS,2L), sep=""),
                   paste("Figure Of Merit: ", round(FOM,2L)," (%)", sep=""),
                   "========================",
                   paste("a=", signif(pars[1L],2L), " +/- ", signif(se_pars[1L],2L), sep=""),
                   paste("b=", ifelse(length(pars)>=2L, paste(signif(pars[2L],2L),
                         " +/- ",signif(se_pars[2L],2L),sep=""), "NULL"), sep=""),
                   paste("c=", ifelse(length(pars)>=3L, paste(signif(pars[3L],2L),
                         " +/- ",signif(se_pars[3L],2L),sep=""), "NULL"), sep=""),
                   paste("d=", ifelse(length(pars)>=4L, paste(signif(pars[4L],2L),
                         " +/- ",signif(se_pars[4L],2L),sep=""), "NULL"), sep=""),
                   paste("e=", ifelse(length(pars)>=5L, paste(signif(pars[5L],2L),
                         " +/- ",signif(se_pars[5L],2L),sep=""), "NULL"), sep="")),                  
                   yjust=2, ncol=1L, cex=1.0, bty="n")
            ###
        } else if (message==1L) {
            ### message=1L: Growth curve fitting failed.
            legend("center", 
                   legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                   "========================",
                   "Status: model fit failed",
                   "========================",
                   paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[1L],2L), 
                         " +/- ", round(Tn[2L],2L), sep="") else "NULL", sep=""),
                   paste("Tn above 3 sigma BG: ", ifelse(!is.null(Tn3BG), as.logical(Tn3BG), "NULL"), sep=""),
                   paste("Ratio of Tn to BG: ", if(!is.null(TnBG.ratio)) paste(round(TnBG.ratio[1L],2L), 
                         " +/- ", round(TnBG.ratio[2L],2L), sep="") else "NULL", sep=""),
                   paste("RSE of Tn: ", ifelse(!is.null(rseTn), round(rseTn,2L), "NULL"), " (%)",sep=""), 
                   paste("Fast ratio of Tn: ", if(!is.null(FR)) paste(round(FR[1L],2L), 
                         " +/- ", round(FR[2L],2L), sep="") else "NULL", sep=""),
                   "========================",
                   paste("Recycling ratio 1: ", if(!is.null(RecyclingRatio1)) paste(round(RecyclingRatio1[1L],2L), 
                         " +/- ", round(RecyclingRatio1[2L],2L), sep="") else "NULL", sep=""),
                   paste("Recycling ratio 2: ", if(!is.null(RecyclingRatio2)) paste(round(RecyclingRatio2[1L],2L), 
                         " +/- ", round(RecyclingRatio2[2L],2L), sep="") else "NULL", sep=""), 
                   paste("Recycling ratio 3: ", if(!is.null(RecyclingRatio3)) paste(round(RecyclingRatio3[1L],2L), 
                         " +/- ", round(RecyclingRatio3[2L],2L), sep="") else "NULL", sep=""),
                   "========================",
                   paste("Recuperation 1: ", if(!is.null(Recuperation1)) paste(round(Recuperation1[1L],2L), 
                         " +/- ", round(Recuperation1[2L],2L), " (%)", sep="") else "NULL (%)", sep=""),
                   paste("Recuperation 2: ", if(!is.null(Recuperation2)) paste(round(Recuperation2[1L],2L), 
                         " +/- ", round(Recuperation2[2L],2L), " (%)", sep="") else "NULL (%)", sep=""), 
                   "========================",
                   paste("Fit model: ", paste(character_model_vec, collapse="-"), sep=""),
                   paste("Pass origin: ", origin, sep=""),
                   paste("Weighted fit: ", weight, sep=""),
                   "Minimized value: NULL",
                   "Average error in fit: NULL",
                   "Reduced Chi-Square: NULL",
                   "Figure Of Merit: NULL (%)",
                   "========================",
                   "a=NULL",
                   "b=NULL",
                   "c=NULL",
                   "d=NULL",
                   "e=NULL"),           
                   yjust=2, ncol=1L, cex=1.0, bty="n")
            ###
        } # end if.
        ###
        on.exit(par(mar=c(5,4,4,2)+0.1,
                    mgp=c(3,1,0),
                    mfrow=c(1L,1L)))
    } # end if.
    ###---------------------------------------------------------------
    ###
    output<-list("message"=message,
                 "LMpars"=LMpars, 
                 "value"=min_obj,
                 "avg.error"=avg_fit_error,
                 "RCS"=RCS,
                 "FOM"=FOM)
    ###
    return(output)
} # end function fitGrowth.default.
#####
