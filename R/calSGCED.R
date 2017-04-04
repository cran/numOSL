#####
calSGCED <-
function(SGCdata, pars, model, origin, avgDev, method="gSGC", 
         errMethod="sp", nsim=500, outpdf=NULL) {
    UseMethod("calSGCED")
} ###
### 2017.04.04. 
calSGCED.default <-
function(SGCdata, pars, model, origin, avgDev, method="gSGC", 
         errMethod="sp", nsim=500, outpdf=NULL) {
    ### Stop if not.
    stopifnot(ncol(SGCdata)==5L, nrow(SGCdata)>=1L,
              is.numeric(SGCdata[,1L,drop=TRUE]),
              is.numeric(SGCdata[,3L,drop=TRUE]), all(SGCdata[,3L,drop=TRUE]>=0),
              is.numeric(SGCdata[,4L,drop=TRUE]), 
              is.numeric(SGCdata[,5L,drop=TRUE]), all(SGCdata[,5L,drop=TRUE]>0),
              is.numeric(pars), 
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(avgDev)==1L, is.numeric(avgDev), avgDev>0.0,
              length(method)==1L, method %in% c("SGC", "gSGC"),
              length(errMethod)==1L, errMethod %in% c("sp","mc"),
              length(nsim)==1L, is.numeric(nsim), nsim>=50L, nsim<=3000L,
              is.null(outpdf) || (length(outpdf)==1L && is.character(outpdf)))
    ###
    ### Check argument pars.
    if (model=="line") {
        if (origin==TRUE && length(pars)!=1L) stop("Error: need provide one parameter!")
        if (origin==FALSE && length(pars)!=2L) stop("Error: need provide two parameter!")
        if (pars[1L]<=0) stop("Error: improper parameters!")
    } else if (model=="exp") {
        if (origin==TRUE && length(pars)!=2L) stop("Error: need provide two parameter!")
        if (origin==FALSE && length(pars)!=3L) stop("Error: need provide three parameter!")
        if (any(pars[seq(2L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="lexp") {
        if (origin==TRUE && length(pars)!=3L) stop("Error: need provide three parameter!")
        if (origin==FALSE && length(pars)!=4L) stop("Error: need provide four parameter!")
        if (any(pars[seq(3L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="dexp") {
        if (origin==TRUE && length(pars)!=4L) stop("Error: need provide four parameter!")
        if (origin==FALSE && length(pars)!=5L) stop("Error: need provide five parameter!")
        if (any(pars[seq(4L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="gok") {
        if (origin==TRUE && length(pars)!=3L) stop("Error: need provide three parameter!")
        if (origin==FALSE && length(pars)!=4L) stop("Error: need provide four parameter!")
        if (any(pars[seq(3L)]<=0)) stop("Error: improper parameters!")
    } # end if.       
    ###
    ###
    ### Check NO and SAR.Cycle for SGCdata.
    colnames(SGCdata) <- c("NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    NO <- sort(as.numeric(levels(factor(SGCdata[,"NO",drop=TRUE]))))
    n <- length(NO)
    ###
    for (i in seq(n)) {
        iIndex <- which(SGCdata[,"NO",drop=TRUE]==NO[i])
        iSAR.Cycle <- substr(SGCdata[iIndex,"SAR.Cycle",drop=TRUE], start=1L, stop=1L)
        ###
        if (!all(diff(iIndex)==1L)) {
            stop(paste("[NO=", NO[i], 
                       "]: 'NO' appears in discontinuous locations!", sep=""))
        } ### end if.
        ###
        if (method=="SGC") {
            if (length(iSAR.Cycle)!=1L) {
                stop(paste("[NO=", NO[i], 
                           "]: should contain only one 'SAR.Cycle'!", sep=""))
            } # end if.
            if (!("N" %in% iSAR.Cycle)) {
                stop(paste("[NO=", NO[i], 
                           "]: should contain 'SAR.Cycle' of 'N'!", sep=""))
            } # end if.
            ###
        } else if (method=="gSGC") {
            if (length(iSAR.Cycle)!=2L) {
                stop(paste("[NO=", NO[i], 
                           "]: should contain two 'SAR.Cycle'!", sep=""))
            } # end if.
            if (!("N" %in% iSAR.Cycle && "R" %in% iSAR.Cycle)) {
                stop(paste("[NO=", NO[i], 
                           "]: should contain 'SAR.Cycle' of 'N' and 'R'!", sep=""))
            } # end if.
        } # end if.
        ###
    } # end for.
    ###
    ###
    cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
    fcn <- function(x) {
        if(model=="line") {
           pars[1L]*x+cst
        } else if(model=="exp") {
           pars[1L]*(1.0-exp(-pars[2L]*x))+cst
        } else if(model=="lexp")  {
           pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst
        } else if(model=="dexp") {
           pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst
        } else if(model=="gok") {
           pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst
        } # end if.
    } # end function fcn.
    ###
    ###
    ScaledNaturalSignal <- matrix(nrow=n, ncol=2L)
    ###
    for(i in seq(n)) {
        iIndex <- which(SGCdata[,"NO",drop=TRUE]==NO[i])
        iSAR.Cycle <- substr(SGCdata[iIndex,"SAR.Cycle",drop=TRUE],start=1L,stop=1L)
        ###
        regenerativeDose   <- SGCdata[iIndex,"Dose",drop=TRUE][iSAR.Cycle=="R"]
        regenerativeSignal <- SGCdata[iIndex,"Signal",drop=TRUE][iSAR.Cycle=="R"]
        ###
        naturalSignal      <- SGCdata[iIndex,"Signal",drop=TRUE][iSAR.Cycle=="N"]
        naturalSignalError <- SGCdata[iIndex,"Signal.Err",drop=TRUE][iSAR.Cycle=="N"]
        ###
        if (method=="SGC") {
            scalingFactor <- 1.0 
        } else if (method=="gSGC") {
            scalingFactor <- fcn(regenerativeDose)/regenerativeSignal
        } # end if.
        ###
        ScaledNaturalSignal[i,] <- scalingFactor*c(naturalSignal,naturalSignalError)       
    } # end for.
    ###
    ###
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
    ### 
    scaleLtx <- sgcED <- ConfInt <- c()
    accept_ID <- saturate_ID <- failED_ID <- failEDError_ID <- c() 
    ###
    if (!is.null(outpdf)) pdf(file=paste(outpdf,".pdf",sep=""))
    ###
    ###
    for (i in seq(n)) {     
        LxTx_seLxTx <- ScaledNaturalSignal[i,,drop=FALSE]
        outDose <- matrix(0, nrow=1L, ncol=2L)
        eemm <- ifelse(errMethod=="sp", 0L, 1L)
        mcED <- vector(length=nsim)
        saturateDose <- -99.0
        acceptRate <- 0.0
        message <- 100
        ###
        res <- .Fortran("calSGCED_fort",as.integer(n2),as.double(LxTx_seLxTx),
                    as.double(pars),outDose=as.double(outDose),as.integer(eemm),
                    as.double(avgDev),mcED=as.double(mcED),saturateDose=as.double(saturateDose),
                    as.integer(mdl),as.integer(nsim),acceptRate=as.double(acceptRate),
                    message=as.integer(message),PACKAGE="numOSL")
        ###
        message <- res$message
        ED <- res$outDose[1L]
        seED <- res$outDose[2L]   
        mcED <- res$mcED
        saturateDose <- res$saturateDose
        acceptRate <- res$acceptRate
        ###
        if (message==0L) {
            scaleLtx <- rbind(scaleLtx, LxTx_seLxTx)
            sgcED <- rbind(sgcED, c(ED, seED))
            if (errMethod=="mc") {
                i_ConfInt <- quantile(mcED,probs=c(0.16,0.84,0.025,0.975))
            } # end if.
            if (errMethod=="sp") {
                i_ConfInt <- vector(length=4L)
                i_ConfInt[1L] <- ED - 0.9944579*seED
                i_ConfInt[2L] <- ED + 0.9944579*seED
                i_ConfInt[3L] <- ED - 1.959964*seED
                i_ConfInt[4L] <- ED + 1.959964*seED
            } # end if.
            ConfInt <- rbind(ConfInt, i_ConfInt)
            accept_ID <- c(accept_ID, NO[i])
        } else if (message==1L) {
            saturate_ID <- c(saturate_ID, NO[i])
        } else if (message==2L) {
            failED_ID <- c(failED_ID, NO[i])
        } else if (message==3L) {
            failEDError_ID <- c(failEDError_ID, NO[i])
        } # end if.
        ###
        ###
        if (!is.null(outpdf)) {
            layout(matrix(c(1L,1L,1L,1L,2L,2L),nrow=2L), respect=TRUE)
            par(mgp=c(2.5,1,0))
            ###
            if (message %in% c(0L, 3L)) {
                lowerX <- min(ED,0.0)*2.0
                upperX <- abs(ED)*2.0
                lowerY <- min(LxTx_seLxTx[1L],0.0)*2.0
                upperY <- abs(LxTx_seLxTx[1L])*1.3
            } else {
                lowerX <- 0.0
                upperX <- saturateDose*1.5
                lowerY <- 0.0
                upperY <- abs(LxTx_seLxTx[1L])*1.3              
            } # end if.
            ###
            par(mar=c(4,4,2,0.5)+0.1)
            plot(NA, NA, main=NULL, xlim=c(lowerX, upperX),  
                 ylim=c(lowerY, upperY), xlab="Regenerative dose (Gy|s)", 
                 ylab=ifelse(method=="gSGC","Normalised standardised OSL","Standardised OSL"),
                 las=0, cex.lab=1.25, cex.main=1.05, xaxs="i", yaxs="i")
            x <- NULL
            ###
            if(model=="line") {
                curve(pars[1L]*x+cst, type="l", add=TRUE, 
                      lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="exp") {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, type="l",   
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="lexp")  {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, type="l",                     
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="dexp") {
                curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, type="l",                      
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="gok") {
                curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, type="l",   
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } # end if.
            ###
            ###
            if (message==0L) {
                ### message=0L: ED and Error calculation succeeded.
                points(x=ED, y=LxTx_seLxTx[1L], pch=23, cex=1.5, bg="grey")
                ###
                if (LxTx_seLxTx[1L]>0 && errMethod=="mc") {
                    dmcED <- density(mcED)
                    dxy <- cbind(dmcED$x,dmcED$y)
                    dxy[,2L] <- (dxy[,2L,drop=TRUE]-min(dxy[,2L,drop=TRUE]))/
                                (max(dxy[,2L,drop=TRUE])-min(dxy[,2L,drop=TRUE]))*
                                 LxTx_seLxTx[1L]*0.8
                    polygon(dxy, col="grey")
                } # end if.
                ###
                if (seED/ED>0.001) {
                    suppressWarnings(arrows(x0=ED-seED/2.0, y0=LxTx_seLxTx[1L],
                        x1=ED+seED/2.0, y1=LxTx_seLxTx[1L],code=3, lwd=1, angle=90, length=0.05, col="black"))
                } # end if.
                lines(x=c(0, ED, ED), y=c(LxTx_seLxTx[1L], LxTx_seLxTx[1L], 0), lty="dashed", lwd=1.5) 
                ###                
                ###
            } else if (message %in% c(1L,2L)) {
                ### message=1L: natural signal saturated.
                ### message=2L: ED calculation failed.
                abline(h=LxTx_seLxTx[1L], lty="dashed", col="red", lwd=1.5)
                ###
            } else if (message==3L) {
                ### message==3L: ED error calculation failed. 
                points(x=ED, y=LxTx_seLxTx[1L], pch=23, cex=1.5, bg="grey")
                lines(x=c(0, ED, ED), y=c(LxTx_seLxTx[1L], LxTx_seLxTx[1L], 0), lty="dashed", lwd=1.5)
            } # end if.
            ###
            grid()
            box(lwd=1)
            ###
            ###
            par(mar=c(9,0.5,9,1)+0.1)
            par(mgp=c(1,1,0))
            plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary", ylab="", cex.lab=1.5)
            if (message==0L) {
                legend("center", 
                       legend=c(paste("ID: [NO=", i, "]", sep=""),
                       "======================",
                       "Status: OK",
                       "======================",
                       paste("Method: ", method, sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""),
                       paste("MC accept-rate: ", ifelse(errMethod=="sp","NULL", round(acceptRate,2L)), " (%)", sep=""),
                       "======================",
                       paste("ED: ",round(ED,2L), " +/- ",round(seED,2L)," (Gy|s)",sep=""), 
                       paste("RSE of ED: ",round(seED/abs(ED)*100.0,2L), " (%)",sep=""),
                       paste("95% interval: [",round(i_ConfInt[3L],2L),", ",round(i_ConfInt[4L],2L),"]", sep=""), 
                       paste("68% interval: [",round(i_ConfInt[1L],2L),", ",round(i_ConfInt[2L],2L),"]", sep=""),
                       "======================"),           
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } else if (message==1L) {
                legend("center",
                       legend=c(paste("Aliquot (grain) ID: [NO=", i, "]", sep=""),
                       "======================",
                       "Status: Ln/Tn saturated",
                       "======================",
                       paste("Method: ", method, sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""),
                       "======================",
                       paste("ED: ",Inf, " +/- ",NA," (Gy|s)",sep=""), 
                       "RSE of ED: NA (%)",  
                       "95% interval: [-Inf, +Inf]", 
                       "68% interval: [-Inf, +Inf]",  
                       "======================"),               
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } else if (message==2L) {
                legend("center",
                       legend=c(paste("ID: [NO=", i, "]", sep=""),
                       "======================",
                       "Status: ED failed (ED < -50 Gy)",
                       "======================",
                       paste("Method: ", method, sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""), 
                       "======================",               
                       paste("ED: ",NA, " +/- ",NA," (Gy|s)",sep=""), 
                       "RSE of ED: NA (%)", 
                       "95% interval: [NA, NA]", 
                       "68% interval: [NA, NA]",   
                       "======================"),             
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } else if (message==3L) {
                legend("center", 
                       legend=c(paste("ID: [NO=", i, "]", sep=""),
                       "======================",
                       ifelse(errMethod=="mc",
                       "Status: ED Error failed (MC accept-rate < 1%)",
                       "Status: ED Error failed (infinite upper ED)"),
                       "======================",
                       paste("Method: ", method, sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""),
                       "======================",
                       paste("ED: ",round(ED,2L), "+/-",NA," (Gy|s)",sep=""), 
                       "RSE of ED: NA (%)", 
                       "95% interval: [NA, NA]", 
                       "68% interval: [NA, NA]",   
                       "======================"),            
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } # end if.
            ###
        } # end if.
        ###
        ###
    } # end for.
    ###
    ###
    if (!is.null(outpdf)) dev.off()
    ###
    if (is.null(sgcED)) stop("Error: fail in SGC ED calculation!")
    ###
    ###
    rownames(scaleLtx) <- paste("[NO=",accept_ID,"]",sep="")
    colnames(scaleLtx) <- c("Ltx", "seLtx")
    ###
    rownames(sgcED) <- paste("[NO=",accept_ID,"]",sep="")
    colnames(sgcED) <- c("ED", "seED")
    ###
    rownames(ConfInt) <- paste("[NO=",accept_ID,"]",sep="")
    colnames(ConfInt) <- c("lower68", "upper68", "lower95", "upper95")
    ###       
    output <- list("saturate.NO"=saturate_ID,
                   "failED.NO"=failED_ID,
                   "failEDError.NO"=failEDError_ID,
                   "scale.Ltx"=scaleLtx,
                   "sgcED"=sgcED,
                   "ConfInt"=ConfInt)
    invisible(output)
} # end function calSGCED.default.
#####
