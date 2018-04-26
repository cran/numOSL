#####
analyseBINdata <- 
function(obj_pickBIN, nfchn, nlchn, bg="late", me=2.0, 
         distp="p", kph=NULL, kdc=NULL, dcr=NULL, 
         FR.fchn=NULL, FR.mchn=NULL, FR.lchn=NULL, 
         signal.type="LxTx", outfile=NULL) {
    UseMethod("analyseBINdata")
} #
### 2018.04.21.
analyseBINdata.default <- 
function(obj_pickBIN, nfchn, nlchn, bg="late", me=2.0, 
         distp="p", kph=NULL, kdc=NULL, dcr=NULL, 
         FR.fchn=NULL, FR.mchn=NULL, FR.lchn=NULL, 
         signal.type="LxTx", outfile=NULL) {
    ### Stop if not.
    stopifnot(class(obj_pickBIN)=="pickBIN", 
              names(obj_pickBIN)==c("BINdata","agID"),
              length(nfchn)==1L, is.numeric(nfchn), 
              length(nlchn)==1L, is.numeric(nlchn),
              length(bg)==1L, bg %in% c("early", "late"),
              length(me)==1L, is.numeric(me),
              length(distp)==1L, distp %in% c("p","op"),
              is.null(kph) || (length(kph)==1L && is.numeric(kph)),
              is.null(kdc) || (length(kdc)==1L && is.numeric(kdc)),
              is.null(dcr) || (length(dcr)==1L && is.numeric(dcr)),
              is.null(FR.fchn) || is.numeric(FR.fchn),
              is.null(FR.mchn) || is.numeric(FR.mchn),
              is.null(FR.lchn) || is.numeric(FR.lchn),
              length(signal.type)==1L, signal.type %in% c("LxTx", "Lx", "Tx"),
              is.null(outfile) || (length(outfile)==1L && is.character(outfile)))
    ###
    if (is.null(obj_pickBIN$agID)) {
        stop("Error: set force.matrix=FALSE in function pickBINdata()!") 
    } # end if.
    ###
    if (distp=="op") {
        if (is.null(kph) || is.null(kdc) || is.null(dcr)) {
            stop("Error: parameter [kph,kdc,dcr] should be provided!")
        } # end if.
    } # end if.
    ###
    recordList <- obj_pickBIN$BINdata
    nList <- length(recordList)
    ###
    ###-----------------------------------
    ### R function for calculating 
    ### signal and associated error. 
    analyseSig <- function(sigData, n_fchn, n_lchn) {
        NPoints <- attr(sigData, "NPoints")
        Delay <- attr(sigData, "Delay")
        Off <- attr(sigData, "Off")
        On <- NPoints - Delay - Off
        if (On<=1L) return(rep(NA,7L))
        ###
        Low <- attr(sigData, "Low")
        High <- attr(sigData, "High")       
        vstep <- (High-Low)/NPoints
        ###
        IRRTime <- attr(sigData, "IRRTime")
        ###   
        ###
        ### Calculate Fast Ratio and its standard error. 
        #-------------------------------------------------------------------------
        if (!is.null(FR.fchn) && !is.null(FR.mchn) && !is.null(FR.lchn)) {
            range_FR.fchn <- range(FR.fchn)
            if (range_FR.fchn[1L]<1L) range_FR.fchn[1L] <- 1L
            if (range_FR.fchn[2L]>On) range_FR.fchn[2L] <- On
            ###
            range_FR.mchn <- range(FR.mchn)
            if (range_FR.mchn[1L]<1L) range_FR.mchn[1L] <- 1L
            if (range_FR.mchn[2L]>On) range_FR.mchn[2L] <- On
            ###
            range_FR.lchn <- range(FR.lchn)
            if (range_FR.lchn[1L]<1L) range_FR.lchn[1L] <- 1L
            if (range_FR.lchn[2L]>On) range_FR.lchn[2L] <- On
            ###
            n_FR.fchn <- range_FR.fchn[2L] - range_FR.fchn[1L] + 1L
            n_FR.mchn <- range_FR.mchn[2L] - range_FR.mchn[1L] + 1L
            n_FR.lchn <- range_FR.lchn[2L] - range_FR.lchn[1L] + 1L
            ###
            kf <- n_FR.lchn/n_FR.fchn
            km <- n_FR.lchn/n_FR.mchn
            ###
            all_On_channels <- (Delay+1L):(Delay+On)
            photon_FR.fchn <- sum((sigData[all_On_channels])[range_FR.fchn[1L]:range_FR.fchn[2L]])
            photon_FR.mchn <- sum((sigData[all_On_channels])[range_FR.mchn[1L]:range_FR.mchn[2L]])
            ###
            photon_FR.lchn <- sum((sigData[all_On_channels])[range_FR.lchn[1L]:range_FR.lchn[2L]])
            photon_FR.lchn_f <- n_FR.fchn*mean((sigData[all_On_channels])[range_FR.lchn[1L]:range_FR.lchn[2L]])
            photon_FR.lchn_m <- n_FR.mchn*mean((sigData[all_On_channels])[range_FR.lchn[1L]:range_FR.lchn[2L]])
            ###       
            net_photon_FR.fchn <- photon_FR.fchn - photon_FR.lchn_f
            net_photon_FR.mchn <- photon_FR.mchn - photon_FR.lchn_m 
            if (abs(net_photon_FR.fchn)<=.Machine$double.eps^0.5) net_photon_FR.fchn <- 1.0e-5  
            if (abs(net_photon_FR.mchn)<=.Machine$double.eps^0.5) net_photon_FR.mchn <- 1.0e-5 
            ###
            FR <- net_photon_FR.fchn / net_photon_FR.mchn
            ###
            if (distp=="p") {              
                ### Eqn.3 of Galbraith (2002).
                rse_net_photon_FR.fchn <- sqrt(photon_FR.fchn+photon_FR.lchn_f/kf)/abs(net_photon_FR.fchn)
                rse_net_photon_FR.mchn <- sqrt(photon_FR.mchn+photon_FR.lchn_m/km)/abs(net_photon_FR.mchn)             
            } else if (distp=="op") {
                ### Eqn.10 of Bluszcz et al.(2015).
                vdif <- kdc^2L-kph^2L
                rse_net_photon_FR.fchn <- sqrt(kph^2L*photon_FR.fchn+vdif*dcr*n_FR.fchn*vstep+
                                         (kph^2L*photon_FR.lchn+vdif*dcr*n_FR.lchn*vstep)/kf^2L)/net_photon_FR.fchn
                ###
                rse_net_photon_FR.mchn <- sqrt(kph^2L*photon_FR.mchn+vdif*dcr*n_FR.mchn*vstep+
                                         (kph^2L*photon_FR.lchn+vdif*dcr*n_FR.lchn*vstep)/km^2L)/net_photon_FR.mchn
            } # end if.
            ###
            rse_net_photon_FR.fchn <- sqrt(rse_net_photon_FR.fchn^2L+(me/100.0)^2L)
            rse_net_photon_FR.mchn <- sqrt(rse_net_photon_FR.mchn^2L+(me/100.0)^2L)
            ###
            seFR <- abs(FR)*sqrt(rse_net_photon_FR.fchn^2L+rse_net_photon_FR.mchn^2L)
        } else {
            FR <- seFR <- NA
        } # end if.
        ###--------------------------------------------------------------
        ###
        ###
        ### Calculate net OSL signal and its standard error.
        #------------------------------------------------------------------
        ###  
        if (n_fchn<On) {
            Lx <- sum(sigData[(Delay+1L):(Delay+n_fchn)])
            n_fchn <- n_fchn
        } else {
            Lx <- sum(sigData[(Delay+1L):(Delay+On-1L)])
            n_fchn <- On-1L
        } # end if.
        ###
        if (bg=="early" && n_lchn<=(On-n_fchn)) {
            bLx <- n_fchn*mean(sigData[(Delay+n_fchn+1L):(Delay+n_fchn+n_lchn)])
            sum_bLx <- sum(sigData[(Delay+n_fchn+1L):(Delay+n_fchn+n_lchn)])
            n_lchn <- n_lchn
        } else if (bg=="late" && n_lchn<=(On-n_fchn)) {
            bLx <- n_fchn*mean(sigData[(Delay+On-n_lchn+1L):(Delay+On)])
            sum_bLx <- sum(sigData[(Delay+On-n_lchn+1L):(Delay+On)])
            n_lchn <- n_lchn
        } else {
            bLx <- n_fchn*mean(sigData[(Delay+n_fchn+1L):(Delay+On)])
            sum_bLx <- sum(sigData[(Delay+n_fchn+1L):(Delay+On)])
            n_lchn <- On-n_fchn 
        } # end if.
        ###
        k <- n_lchn/n_fchn
        ###
        netLx <- Lx-bLx 
        if (abs(netLx)<=.Machine$double.eps^0.5) netLx <- 1.0e-5    
        ###
        if (distp=="p") {              
            ### Eqn.3 of Galbraith (2002).
            rse_netLx <- sqrt(Lx+bLx/k)/abs(netLx)
        } else if (distp=="op") {
            ### Eqn.10 of Bluszcz et al.(2015).
            vdif <- kdc^2L-kph^2L
            rse_netLx <- sqrt(kph^2L*Lx+vdif*dcr*n_fchn*vstep+
                        (kph^2L*sum_bLx+vdif*dcr*n_lchn*vstep)/k^2L)/abs(netLx)
        } # end if.
        rse_netLx <- sqrt(rse_netLx^2L+(me/100.0)^2L)
        ###---------------------------------------------------------------------
        ###
        values <- c(IRRTime, Lx, bLx, netLx, rse_netLx, FR, seFR)
        ###
        return(values)
    } # end function analyseSig.
    ###---------------------------------------------
    ###
    ###
    NO <- obj_pickBIN$agID[,"NO",drop=TRUE]
    LxTx_list <- vector(length=length(NO), mode="list")
    ###
    LnTn.curve <- list()
    ###
    ###
    for (i in seq(NO)) {
        iRecord <- recordList[[i]]
        length_iRecord <- length(iRecord)
        ###
        ###
        Ln_NPoints <- attr(iRecord[[1L]], "NPoints")
        ###
        if (Ln_NPoints>0L) {
            Ln_Low <- attr(iRecord[[1L]], "Low")
            Ln_High <- attr(iRecord[[1L]], "High")
            Ln_vstep <- (Ln_High-Ln_Low)/Ln_NPoints
            ###
            iLn_curve_x <- seq(from=Ln_vstep, to=Ln_High, by=Ln_vstep)                  
            iLn_curve_y <- as.numeric(iRecord[[1L]])
        } else {
            iLn_curve_x <- iLn_curve_y <- NA
        } # end if.
        ###
        ###
        if (length_iRecord>=2L) {
            Tn_NPoints <- attr(iRecord[[2L]], "NPoints")
            ###
            if (Tn_NPoints>0L) {
                Tn_Low <- attr(iRecord[[2L]], "Low")
                Tn_High <- attr(iRecord[[2L]], "High")
                Tn_vstep <- (Tn_High-Tn_Low)/Tn_NPoints
                ###
                iTn_curve_x <- seq(from=Tn_vstep, to=Tn_High, by=Tn_vstep) 
                iTn_curve_y <- as.numeric(iRecord[[2L]])
            } else {
                iTn_curve_x <- iTn_curve_y <- NA
            } # end if.
        } else {
            iTn_curve_x <- iTn_curve_y <- NA
        } # end if.
        ###
        characterNO <- paste("NO", i, sep="")
        LnTn.curve[[characterNO]] <- list("Ln.x"=iLn_curve_x, "Ln.y"=iLn_curve_y,
                                          "Tn.x"=iTn_curve_x, "Tn.y"=iTn_curve_y)
        ###
        ###
        i_matrix <- c()
        for (j in seq(length_iRecord)) {
            i_matrix <- rbind(i_matrix, analyseSig(iRecord[[j]], nfchn, nlchn))
        } # end for.
        colnames(i_matrix) <- c("IRRTime","OSL","BG","netOSL","rse.netOSL","FR","seFR")
        LxTx_list[[i]] <- i_matrix
    } # end for.
    ###
    ###
    Position <- obj_pickBIN$agID[,"Position",drop=TRUE]
    Grain <- obj_pickBIN$agID[,"Grain",drop=TRUE]
    ###
    ###
    SARdata <- ALLdata <- data.frame(stringsAsFactors=FALSE)
    ###
    Tn3BG_vec <- rseTn_vec <- 
    TnBG.ratio_vec <- seTnBG.ratio_vec <-
    FR_vec <- seFR_vec <- c()
    ###
    Tn_vec <- seTn_vec <- c()
    ###
    TxTn <- list()
    ###
    agID <- c()
    ###
    accept_characterNO <- c()
    ###
    for (i in seq(NO)) {
        iLxTx <- LxTx_list[[i]]
        selectedIndex <- which(is.finite(iLxTx[,"IRRTime",drop=TRUE]) &
                               is.finite(iLxTx[,"netOSL",drop=TRUE]) &
                               is.finite(iLxTx[,"rse.netOSL",drop=TRUE]))
        iLxTx <- iLxTx[selectedIndex,,drop=FALSE]
        nr <- nrow(iLxTx)
        nc <- ncol(iLxTx)
        ###
        ###
        if (nr==0L) {
            cat(paste("[NO=",NO[i],",Position=",Position[i],",Grain=",Grain[i],
                      "]: analysis failed!\n",sep=""))
            ###
        } else if (nr==1L) {
            cat(paste("[NO=",NO[i],",Position=",Position[i],",Grain=",Grain[i],
                      "]: only one Lx[Tx] (omitted)!\n",sep=""))
            ###
        } else if (!(selectedIndex[1L]==1L && selectedIndex[2L]==2L)) {
            cat(paste("[NO=",NO[i],",Position=",Position[i],",Grain=",Grain[i],
                      "]: analysis failed!\n",sep=""))
        } else {
            ###
            accept_characterNO <- c(accept_characterNO, paste("NO", i, sep=""))
            ###
            if (nr%%2L!=0L) {
                cat(paste("[NO=",NO[i],",Position=",Position[i],",Grain=",Grain[i],
                          "]: unpaired Lx-Tx (omit the last Lx[Tx])!\n",sep="")) 
                iLxTx <- iLxTx[-nr,,drop=FALSE]
            } # end if.
            ###
            nPairedLxTx <- nrow(iLxTx)/2L
            ###
            iNO <- rep(NO[i], nPairedLxTx)
            iPosition <- rep(Position[i], nPairedLxTx)
            iGrain <- rep(Grain[i], nPairedLxTx)
            ###
            if (nPairedLxTx==1L) {
                iSAR.Cycle <- "N"
            } else {
                iSAR.Cycle <- c("N", paste("R",seq(nPairedLxTx-1L),sep=""))
            } # end if.
            ###
            ###
            Dose <- 
            Init <-   BG <- netLx <- rse_netLx <- 
            TInit <- TBG <- netTx <- rse_netTx <- 
            TFR <- seTFR <- LxTx <- seLxTx <- vector(length=nPairedLxTx) 
            ###
            for (j in seq(nPairedLxTx)) {
                Dose[j] <- iLxTx[2L*j-1L,"IRRTime"]
                ###
                Init[j] <- iLxTx[2L*j-1L,"OSL"]
                BG[j] <- iLxTx[2L*j-1L,"BG"]
                ###
                netLx[j] <- iLxTx[2L*j-1L,"netOSL"]
                rse_netLx[j] <- iLxTx[2L*j-1L,"rse.netOSL"]
                ###
                TInit[j] <- iLxTx[2L*j,"OSL"]
                TBG[j] <- iLxTx[2L*j,"BG"]
                ###
                netTx[j] <- iLxTx[2L*j,"netOSL"]
                rse_netTx[j] <- iLxTx[2L*j,"rse.netOSL"]
                ###
                TFR[j] <- iLxTx[2L*j,"FR"]
                seTFR[j] <- iLxTx[2L*j,"seFR"]
                ###
                LxTx[j] <- netLx[j]/netTx[j]
                seLxTx[j] <- abs(LxTx[j])*sqrt((rse_netLx[j])^2L+
                                               (rse_netTx[j])^2L)    
            } # end for.
            ### 
            ###
            if (distp=="p") {
                Tn3BG_vec <- c(Tn3BG_vec, netTx[1L]>0.0 && netTx[1L]>3.0*sqrt(TBG[1L]))
            } else {
                if (bg=="early") {
                    Tn3BG_vec <- c(Tn3BG_vec, netTx[1L]>0.0 && netTx[1L]>3.0*kph*sqrt(TBG[1L]))
                } else if (bg=="late") {
                    Tn3BG_vec <- c(Tn3BG_vec, netTx[1L]>0.0 && netTx[1L]>3.0*kdc*sqrt(TBG[1L]))
                } # end if.
            } # end if.
            ###
            rseTn_vec <- c(rseTn_vec, rse_netTx[1L]*100.0)
            ###
            TnBG.ratio_vec <- c(TnBG.ratio_vec, TInit[1L]/TBG[1L])
            if (distp=="p") {
                useitpass <- abs(TInit[1L]/TBG[1L])*sqrt(1.0/abs(TInit[1L]) + 1.0/abs(TBG[1L]))
            } else if (distp=="op") {
                useitpass <- abs(TInit[1L]/TBG[1L])*sqrt(kph^2L/abs(TInit[1L]) + kdc^2L/abs(TBG[1L]))
            } # end if.
            seTnBG.ratio_vec <- c(seTnBG.ratio_vec, useitpass)
            ###
            FR_vec <- c(FR_vec, TFR[1L])
            seFR_vec <- c(seFR_vec, seTFR[1L])
            ###
            Tn_vec <- c(Tn_vec, netTx[1L])
            seTn_vec <- c(seTn_vec, abs(netTx[1L])*rse_netTx[1L])
            ###
            characterNO <- paste("NO", i, sep="")
            TxTn[[characterNO]] <- netTx/netTx[1L]
            ###
            agID <- rbind(agID, c(NO[i], Position[i], Grain[i]))
            ###
            if (signal.type=="LxTx") {
                DataFrame1 <- data.frame(iNO, iSAR.Cycle, Dose, LxTx, seLxTx,
                                         stringsAsFactors=FALSE) 
            } else if (signal.type=="Lx") {
                DataFrame1 <- data.frame(iNO, iSAR.Cycle, Dose, netLx, abs(netLx)*rse_netLx,
                                         stringsAsFactors=FALSE) 
            } else if (signal.type=="Tx") {
                DataFrame1 <- data.frame(iNO, iSAR.Cycle, Dose, netTx, abs(netTx)*rse_netTx,
                                         stringsAsFactors=FALSE) 
            } # end if.
            SARdata <- rbind(SARdata, DataFrame1) 
            ###    
            if (!is.null(outfile)) {
                DataFrame2 <- data.frame(iNO, iPosition, iGrain, 
                                         iSAR.Cycle, Dose, 
                                         Init, BG, netLx, abs(netLx)*rse_netLx,
                                         TInit, TBG, netTx, abs(netTx)*rse_netTx, 
                                         LxTx, seLxTx,
                                         stringsAsFactors=FALSE)
                ALLdata <- rbind(ALLdata, DataFrame2)
            } # end if.
            ###  
        } # end if.
        ###
    } # end for.
    ###
    ###
    if (nrow(SARdata)==0L) stop("Error: no SAR data can be returned!")
    ###
    LnTn.curve <- LnTn.curve[accept_characterNO]
    ###
    if (!is.null(outfile)) {
        rownames(ALLdata) <- NULL
        colnames(ALLdata) <- c("NO","Position","Grain",
                               "SAR.Cycle","Dose",
                               "Init","BG","Lx","seLx",
                               "TInit","TBG","Tx","seTx",
                               "LxTx","seLxTx")
        write.csv(ALLdata, file=paste(outfile, ".csv", sep=""))
    } # end if.
    ###
    rownames(SARdata) <- NULL
    if (signal.type=="LxTx") {
        colnames(SARdata) <- c("NO","SAR.Cycle","Dose","LxTx","seLxTx")
    } else if (signal.type=="Lx") {
        colnames(SARdata) <- c("NO","SAR.Cycle","Dose","Lx","seLx")
    } else if (signal.type=="Tx") {
        colnames(SARdata) <- c("NO","SAR.Cycle","Dose","Tx","seTx")
    } # end if.
    ###
    criteria <- cbind(Tn3BG_vec, TnBG.ratio_vec, seTnBG.ratio_vec, rseTn_vec, FR_vec, seFR_vec)
    rownames(criteria) <- NULL
    colnames(criteria) <- c("Tn3BG", "TnBG.ratio", "seTnBG.ratio", "rseTn", "FR", "seFR")
    ###
    Tn <- cbind(Tn_vec, seTn_vec)
    rownames(Tn) <- NULL
    colnames(Tn) <- c("Tn", "seTn")
    ###
    rownames(agID) <- NULL
    colnames(agID) <- c("NO","Position","Grain")
    ###
    output <- list("SARdata"=SARdata,
                   "criteria"=criteria,
                   "Tn"=Tn,
                   "LnTn.curve"=LnTn.curve,
                   "TxTn"=TxTn,
                   "agID"=agID) 
    ###
    class(output) <- "analyseBIN"
    ###
    invisible(output)
} # end function analyseBINdata.default.
#####
