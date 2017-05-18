#####
calSGCED <-
function(obj_analyseBIN, SGCpars, model, origin, avgDev, 
         method="SGC", SAR.Cycle="N", errMethod="sp", Tn.above.3BG=TRUE, 
         TnBG.ratio.low=NULL, rseTn.up=NULL, FR.low=NULL, rseED.up=NULL, 
         use.se=TRUE, outpdf=NULL, outfile=NULL) {
    UseMethod("calSGCED")
} ###
### 2017.05.18. 
calSGCED.default <-
function(obj_analyseBIN, SGCpars, model, origin, avgDev, 
         method="SGC", SAR.Cycle="N", errMethod="sp", Tn.above.3BG=TRUE, 
         TnBG.ratio.low=NULL, rseTn.up=NULL, FR.low=NULL, rseED.up=NULL, 
         use.se=TRUE, outpdf=NULL, outfile=NULL) {
    ### Stop if not.
    stopifnot(class(obj_analyseBIN)=="analyseBIN", 
              names(obj_analyseBIN)==c("SARdata","criteria","Tn","LnTn.curve","TxTn","agID"),
              is.numeric(SGCpars), 
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(avgDev)==1L, is.numeric(avgDev), avgDev>=0.0,
              length(method)==1L, method %in% c("SGC", "gSGC"),
              length(SAR.Cycle) %in% c(1L, 2L), all(substr(SAR.Cycle,1L,1L) %in% c("N","R")),
              length(errMethod)==1L, errMethod =="sp",
              length(Tn.above.3BG)==1L, is.logical(Tn.above.3BG),
              is.null(TnBG.ratio.low) || (length(TnBG.ratio.low)==1L && is.numeric(TnBG.ratio.low)),
              is.null(rseTn.up) || (length(rseTn.up)==1L && is.numeric(rseTn.up)),
              is.null(FR.low) || (length(FR.low)==1L && is.numeric(FR.low)),
              is.null(rseED.up) || (length(rseED.up)==1L && is.numeric(rseED.up)),
              length(use.se)==1L, is.logical(use.se),
              is.null(outpdf) || (length(outpdf)==1L && is.character(outpdf)),
              is.null(outfile) || (length(outfile)==1L && is.character(outfile)))
    ###
    ### Check argument SGCpars.
    if (model=="line") {
        if (origin==TRUE && length(SGCpars)!=1L) stop("Error: need provide one parameter!")
        if (origin==FALSE && length(SGCpars)!=2L) stop("Error: need provide two parameter!")
        if (SGCpars[1L]<=0) stop("Error: improper parameters!")
    } else if (model=="exp") {
        if (origin==TRUE && length(SGCpars)!=2L) stop("Error: need provide two parameter!")
        if (origin==FALSE && length(SGCpars)!=3L) stop("Error: need provide three parameter!")
        if (any(SGCpars[seq(2L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="lexp") {
        if (origin==TRUE && length(SGCpars)!=3L) stop("Error: need provide three parameter!")
        if (origin==FALSE && length(SGCpars)!=4L) stop("Error: need provide four parameter!")
        if (any(SGCpars[seq(3L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="dexp") {
        if (origin==TRUE && length(SGCpars)!=4L) stop("Error: need provide four parameter!")
        if (origin==FALSE && length(SGCpars)!=5L) stop("Error: need provide five parameter!")
        if (any(SGCpars[seq(4L)]<=0)) stop("Error: improper parameters!")
    } else if (model=="gok") {
        if (origin==TRUE && length(SGCpars)!=3L) stop("Error: need provide three parameter!")
        if (origin==FALSE && length(SGCpars)!=4L) stop("Error: need provide four parameter!")
        if (any(SGCpars[seq(3L)]<=0)) stop("Error: improper parameters!")
    } # end if.       
    ###
    if (method=="SGC" && length(SAR.Cycle)==2L) {
        stop("Error: the original SGC method needs only SAR cycle of 'N'!")
    } # end if.
    ###
    if (method=="gSGC" && length(SAR.Cycle)==1L) {
        stop("Error: the improved SGC method needs both SAR cycle of 'N' and 'R'!")
    } # end if.
    ###
    SARdata <- obj_analyseBIN$SARdata
    criteria <- obj_analyseBIN$criteria
    Tn <- obj_analyseBIN$Tn
    LnTn.curve <- obj_analyseBIN$LnTn.curve
    TxTn <- obj_analyseBIN$TxTn
    agID <- obj_analyseBIN$agID
    ###
    ###
    colnames(SARdata) <- c("NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    ###
    ### Check NO and SAR.Cycle for SARdata.
    NO <- agID[,"NO",drop=TRUE]
    Position <- agID[,"Position",drop=TRUE]
    Grain <- agID[,"Grain",drop=TRUE]
    n <- length(NO)
    ###
    nag <- n
    ###
    for (i in seq(n)) {
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i])
        iSAR.Cycle <- substr(SARdata[iIndex,"SAR.Cycle",drop=TRUE], start=1L, stop=1L)   
        ###
        if (!all(diff(iIndex)==1L)) {
            stop(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                       "]: 'NO' appears in discontinuous locations!", sep=""))
        } ### end if. 
        ###  
        if (!all(iSAR.Cycle %in% c("N","R"))) {
             stop(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                        "]: incorrect 'SAR.Cycle'!", sep=""))
        } # end if.
        ###
        if (all(iSAR.Cycle=="R")) {
            stop(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                       "]: should contain 'SAR.Cycle' of 'N'!", sep=""))
        } # end if.
        ###
        if (sum(iSAR.Cycle=="N")>1L) {
            stop(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                       "]: should contain only one 'SAR.Cycle' of 'N'!", sep=""))
        } # end if.
    } # end for.
    ###
    ###
    action_character <- "Total number of analyzed aliquots (grains)"
    step_reject_N <- nag
    ###
    sigma <- 2.0
    ###
    NPG <- function(x) paste("[NO=",x[1L],",Position=",x[2L],",Grain=",x[3L],"]",sep="")
    ###
    is_forced_object <- (is.null(criteria)) && (is.null(Tn)) && 
                        (is.null(LnTn.curve)) && (is.null(TxTn))
    ###
    ### Apply signal related rejection criteria. 
    ###------------------------------------------------------------------
    if (is_forced_object==FALSE) { 
        ### No se consideration.
        if (Tn.above.3BG==TRUE) {
            all_value <- criteria[,"Tn3BG",drop=TRUE]
            ###
            select_index <- which(all_value==1L)
            ###
            if (length(select_index)==0L) {
                stop("Error: no acceptable SGC ED if [Tn.above.3BG] is applied!")
            } # end if.
            ###
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                "Rejection criterion: Tn below 3 sigma BG")
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            ###
            Tn3BG_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            Tn3BG_reject <- NULL
        } # end if.
        ###

        ### Have se consideration.
        if (!is.null(TnBG.ratio.low)) {
            all_value <- criteria[,"TnBG.ratio",drop=TRUE]
            ###
            if (use.se==FALSE) {
                select_index <- which(all_value>TnBG.ratio.low)
            } else {
                all_se_value <- criteria[,"seTnBG.ratio",drop=TRUE]
                select_index <- which(all_value+sigma*all_se_value>TnBG.ratio.low)
            } # end if.
            ###
            if (length(select_index)==0L) {
                stop("Error: no acceptable SGC ED if [TnBG.ratio.low] is applied!")
            } # end if.
            ###
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: ratio of Tn to BG below ", TnBG.ratio.low, sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            ###
            TnBG.ratio_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            TnBG.ratio_reject <- NULL
        } # end if.
        ###

        ### No se consideration.
        if (!is.null(rseTn.up)) {
            all_value <- criteria[,"rseTn",drop=TRUE]
            ###
            select_index <- which(abs(all_value)<rseTn.up)
            ###
            if (length(select_index)==0L) {
                stop("Error: no acceptable SGC ED if [rseTn.up] is applied!")
            } # end if.
            ###
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: RSE of Tn exceeds ", rseTn.up, "%", sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            ###
            rseTn_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            rseTn_reject <- NULL
        } # end if.
        ###

        ### Have se consideration.
        if (!is.null(FR.low)) {
            all_value <- criteria[,"FR",drop=TRUE]
            ###
            if (use.se==FALSE) {
                select_index <- which(all_value>FR.low)
            } else {
                all_se_value <- criteria[,"seFR",drop=TRUE]
                select_index <- which(all_value+sigma*all_se_value>FR.low)
            } # end if.
            ###
            if (length(select_index)==0L) {
                stop("Error: no acceptable SGC ED if [FR.low] is applied!")
            } # end if.
            ###
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: fast ratio of Tn below ", FR.low, sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ### 
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            ###
            FR_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            FR_reject <- NULL
        } # end if.
        ###
    } else {
        if ((Tn.above.3BG==TRUE) || (!is.null(TnBG.ratio.low)) ||
            (!is.null(rseTn.up)) || (!is.null(FR.low))) {
            cat("Note: signal-related rejection criteria cannot be applied!\n")
        } # end if.
    } # end if.
    ###---------------------------------------------------------------------

    ###
    cst <- ifelse(origin==TRUE, 0, SGCpars[length(SGCpars)])
    fcn <- function(x) {
        if(model=="line") {
           SGCpars[1L]*x+cst
        } else if(model=="exp") {
           SGCpars[1L]*(1.0-exp(-SGCpars[2L]*x))+cst
        } else if(model=="lexp")  {
           SGCpars[1L]*(1.0-exp(-SGCpars[2L]*x))+SGCpars[3L]*x+cst
        } else if(model=="dexp") {
           SGCpars[1L]*(1.0-exp(-SGCpars[2L]*x))+SGCpars[3L]*(1.0-exp(-SGCpars[4L]*x))+cst
        } else if(model=="gok") {
           SGCpars[1L]*(1.0-(1.0+SGCpars[2L]*SGCpars[3L]*x)^(-1.0/SGCpars[3L]))+cst
        } # end if.
    } # end function fcn.
    ###
    
    ### Re-write NO, Position, and Grain.
    NO <- agID[,"NO",drop=TRUE]
    Position <- agID[,"Position",drop=TRUE]
    Grain <- agID[,"Grain",drop=TRUE]
    n <- length(NO)
    ###
    
    ###
    ScaledNaturalSignal <- c()
    ###
    extractLOOP <- c()
    ###
    for (i in seq(n)) {
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i])  
        ith_SARdata <- SARdata[iIndex,,drop=FALSE]
        ###
        select_index <- which(ith_SARdata[,"SAR.Cycle",drop=TRUE] %in% SAR.Cycle)
        ###
        if (length(select_index)==length(SAR.Cycle)) {
            ith_SGCdata <- ith_SARdata[select_index,,drop=FALSE]
            iSAR.Cycle <- substr(ith_SGCdata[,"SAR.Cycle",drop=TRUE],start=1L,stop=1L)
            ###
            naturalSignal      <- ith_SGCdata[iSAR.Cycle=="N","Signal",drop=TRUE]
            naturalSignalError <- ith_SGCdata[iSAR.Cycle=="N","Signal.Err",drop=TRUE]
            ###
            if (method=="SGC") {
                scalingFactor <- 1.0 
            } else if (method=="gSGC") {
                regenerativeDose   <- ith_SGCdata[iSAR.Cycle=="R","Dose",drop=TRUE]
                regenerativeSignal <- ith_SGCdata[iSAR.Cycle=="R","Signal",drop=TRUE]
                ###
                scalingFactor <- fcn(regenerativeDose)/regenerativeSignal
            } # end if.
            ###
            ScaledNaturalSignal <- rbind(ScaledNaturalSignal, 
                scalingFactor*c(naturalSignal,naturalSignalError))
            ###
            extractLOOP <- c(extractLOOP, i)
        } else {
             cat(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                       "]: desired 'SAR.cycle' cannot be extracted (omitted)!\n", sep=""))
        } # end if.
    } # end for.
    ###
    if (is.null(extractLOOP)) stop("Error: data used for SGC ED calculation is not available!")
    ### 
    ### Re-write agID.
    agID <- agID[extractLOOP,,drop=FALSE]
    ###
    ### Re-write NO, Position, and Grain.
    NO <- agID[,"NO",drop=TRUE]
    Position <- agID[,"Position",drop=TRUE]
    Grain <- agID[,"Grain",drop=TRUE]
    n <- length(NO)
    ###
    ### Re-write criteria, Tn, and LnTn.curve.
    criteria <- criteria[extractLOOP,,drop=FALSE]
    Tn <- Tn[extractLOOP,,drop=FALSE]
    LnTn.curve <- LnTn.curve[extractLOOP]
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
    acceptNO <- acceptPosition <- acceptGrain <-  
    Tn_vec <- seTn_vec <- Tn3BG_vec <- TnBG.ratio_vec <- 
    seTnBG.ratio_vec <- rseTn_vec <- FR_vec <- seFR_vec <- 
    ScaleLtx_vec <- seScaleLtx_vec <- ED_vec <- seED_vec <- rseED_vec <- 
    lower68_vec <- upper68_vec <- lower95_vec <- upper95_vec <- 
    saturate_ID <- failED_ID <- failEDError_ID <- c()
    ###
    if (!is.null(outpdf)) pdf(file=paste(outpdf,".pdf",sep=""))
    ###
    ###
    nsim <- 600L
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
                    as.double(SGCpars),outDose=as.double(outDose),as.integer(eemm),
                    as.double(avgDev),mcED=as.double(mcED),saturateDose=as.double(saturateDose),
                    as.integer(mdl),as.integer(nsim),acceptRate=as.double(acceptRate),
                    message=as.integer(message),PACKAGE="numOSL")
        ###
        message <- res$message
        saturateDose <- res$saturateDose
        ###
        if (message==0L) {
            ###
            ### Both SGC ED calculation and error assessment succeed.
            acceptNO <- c(acceptNO, NO[i])
            acceptPosition <- c(acceptPosition, Position[i])
            acceptGrain <- c(acceptGrain, Grain[i])
            ###
            Tn_vec <- c(Tn_vec, if(!is.null(Tn)) Tn[i,"Tn",drop=TRUE] else NA)
            seTn_vec <- c(seTn_vec, if(!is.null(Tn)) Tn[i,"seTn",drop=TRUE] else NA)
            ###
            Tn3BG_vec <- c(Tn3BG_vec, if(!is.null(criteria)) criteria[i,"Tn3BG",drop=TRUE] else NA)
            TnBG.ratio_vec <- c(TnBG.ratio_vec, if(!is.null(criteria)) criteria[i,"TnBG.ratio",drop=TRUE] else NA)
            seTnBG.ratio_vec <- c(seTnBG.ratio_vec, if(!is.null(criteria)) criteria[i,"seTnBG.ratio",drop=TRUE] else NA)
            rseTn_vec <- c(rseTn_vec, if(!is.null(criteria)) criteria[i,"rseTn",drop=TRUE] else NA)           
            FR_vec <- c(FR_vec, if(!is.null(criteria)) criteria[i,"FR",drop=TRUE] else NA)
            seFR_vec <- c(seFR_vec,  if(!is.null(criteria)) criteria[i,"seFR",drop=TRUE] else NA)
            ### 
            ScaleLtx_vec <- c(ScaleLtx_vec, LxTx_seLxTx[1L])
            seScaleLtx_vec <- c(seScaleLtx_vec, LxTx_seLxTx[2L])        
            ###
            ED_vec <- c(ED_vec, res$outDose[1L])
            seED_vec <- c(seED_vec, res$outDose[2L])
            ###
            rseED_vec <- c(rseED_vec, res$outDose[2L]/abs(res$outDose[1L])*100.0)
            ###
            i_ConfInt <- vector(length=4L)
            i_ConfInt[1L] <- res$outDose[1L]-0.9944579*res$outDose[2L]
            i_ConfInt[2L] <- res$outDose[1L]+0.9944579*res$outDose[2L]
            i_ConfInt[3L] <- res$outDose[1L]-1.959964*res$outDose[2L]
            i_ConfInt[4L] <- res$outDose[1L]+1.959964*res$outDose[2L]
            ###
            lower68_vec <- c(lower68_vec, i_ConfInt[1L])
            upper68_vec <- c(upper68_vec, i_ConfInt[2L])
            lower95_vec <- c(lower95_vec, i_ConfInt[3L])
            upper95_vec <- c(upper95_vec, i_ConfInt[4L])
            ###
        } else if (message==1L) {
            ###
            ### Natural signal saturates.
            saturate_ID <- rbind(saturate_ID, agID[i,,drop=TRUE])
            ###
        } else if (message==2L) {
            ###
            ### SGC ED calculation fails.
            failED_ID <- rbind(failED_ID, agID[i,,drop=TRUE])
            ###
        } else if (message==3L) {
            ###
            ### SGC ED error assessment fails.
            failEDError_ID <- rbind(failEDError_ID, agID[i,,drop=TRUE])
            ###
        } # end if.
        ###
        ###===================================================================================
        if (!is.null(outpdf)) {
            layout(matrix(c(1L,1L,2L,1L,1L,2L,3L,3L,3L),nrow=3L), respect=TRUE)
            par(mgp=c(2.5,1,0))
            ###
            if (message %in% c(0L, 3L)) {
                lowerX <- min(res$outDose[1L],0.0)*2.0
                upperX <- abs(res$outDose[1L])*2.0
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
                curve(SGCpars[1L]*x+cst, type="l", add=TRUE, 
                      lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="exp") {
                curve(SGCpars[1L]*(1.0-exp(-SGCpars[2L]*x))+cst, type="l",   
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="lexp")  {
                curve(SGCpars[1L]*(1.0-exp(-SGCpars[2L]*x))+SGCpars[3L]*x+cst, type="l",                     
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="dexp") {
                curve(SGCpars[1L]*(1.0-exp(-SGCpars[2L]*x))+SGCpars[3L]*(1.0-exp(-SGCpars[4L]*x))+cst, type="l",                      
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } else if(model=="gok") {
                curve(SGCpars[1L]*(1.0-(1.0+SGCpars[2L]*SGCpars[3L]*x)^(-1.0/SGCpars[3L]))+cst, type="l",   
                      add=TRUE, lwd=2.5, from=lowerX, to=upperX, col="skyblue")
            } # end if.
            ###
            ###
            if (message==0L) {
                ### message=0L: ED and Error calculation succeeded.
                points(x=res$outDose[1L], y=LxTx_seLxTx[1L], pch=23, cex=1.5, bg="grey")
                ###
                if (res$outDose[2L]/res$outDose[1L]>0.001) {
                    suppressWarnings(arrows(x0=res$outDose[1L]-res$outDose[2L]/2.0, y0=LxTx_seLxTx[1L],
                        x1=res$outDose[1L]+res$outDose[2L]/2.0, y1=LxTx_seLxTx[1L], code=3, lwd=1, 
                        angle=90, length=0.05, col="black"))
                } # end if.
                lines(x=c(0, res$outDose[1L], res$outDose[1L]), 
                      y=c(LxTx_seLxTx[1L], LxTx_seLxTx[1L], 0), lty="dashed", lwd=1.5) 
                ###                
                ###
            } else if (message %in% c(1L,2L)) {
                ### message=1L: natural signal saturated.
                ### message=2L: ED calculation failed.
                abline(h=LxTx_seLxTx[1L], lty="dashed", col="red", lwd=1.5)
                ###
            } else if (message==3L) {
                ### message==3L: ED error calculation failed. 
                points(x=res$outDose[1L], y=LxTx_seLxTx[1L], pch=23, cex=1.5, bg="grey")
                lines(x=c(0, res$outDose[1L], res$outDose[1L]), 
                      y=c(LxTx_seLxTx[1L], LxTx_seLxTx[1L], 0), lty="dashed", lwd=1.5)
            } # end if.
            ###
            grid()
            box(lwd=1)
            ###-------------------------------------------------------------------------------------- 
            par(mar=c(4,4,0.5,0.5)+0.1)
            if (!is.null(LnTn.curve)) {
                x_max <- max(max(LnTn.curve[[i]][["Ln.x"]]), max(LnTn.curve[[i]][["Tn.x"]]), na.rm=TRUE)
                y_max <- max(max(LnTn.curve[[i]][["Ln.y"]]), max(LnTn.curve[[i]][["Tn.y"]]), na.rm=TRUE)
                ###
                plot(x=LnTn.curve[[i]][["Ln.x"]], y=LnTn.curve[[i]][["Ln.y"]], type="l", 
                     lwd=1.5, col="blue", main=NULL, xlim=c(0, x_max), ylim=c(0, y_max),
                     xlab="Stimulation time (s)", ylab="Photon counts", las=0, 
                     xaxt="n", cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
                ###
                if (length(LnTn.curve[[i]][["Tn.x"]])>1L) {
                    # Both Tn.x and Tn.y are not equal to NA of length 1.
                    points(x=LnTn.curve[[i]][["Tn.x"]], y=LnTn.curve[[i]][["Tn.y"]], 
                           type="l", lwd=1.5, col="red")
                } # end if.
                ###
                x_axis_location <- axTicks(side=1L)
                axis(side=1L, at=x_axis_location, labels=as.character(x_axis_location))
                legend("topright", legend=c("Ln decay curve", "Tn decay curve"), col=c("blue","red"), 
                       lwd=1.5, yjust=2, ncol=1L, cex=1.2, bty="n")
            } else {
                plot(x=1L, y=1.0, type="n", lwd=1.5,
                     main=NULL, xlab="Stimulation time (s)", ylab="Photon counts", las=0, 
                     xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5, lab=c(7,7,9))
                axis(side=1L, at=c(0.7, 1, 1.3), labels=c("x1", "x2", "x3"))
                axis(side=2L, at=c(0.7, 1, 1.3), labels=c("y1", "y2", "y3"))
            } # end if.
            ###-------------------------------------------------------------------------------------- 
            ###
            par(mar=c(15,0.5,15,0.5)+0.1)
            par(mgp=c(1,1,0))
            plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="Summary", ylab="", cex.lab=1.5)
            ###
            NO_Position_Grain <- paste("[NO=",NO[i],", Position=",Position[i],", Grain=",Grain[i],"]",sep="")
            ###
            if (message==0L) {
                legend("center", 
                       legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                       "========================",
                       "Status: OK",
                       "========================",
                       paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[i,"Tn",drop=TRUE],2L), " +/- ", 
                             round(Tn[i,"seTn",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("Tn above 3 sigma BG: ", ifelse(!is.null(criteria), 
                             as.logical(criteria[i,"Tn3BG",drop=TRUE]), "NULL"), sep=""),
                       paste("Ratio of Tn to BG: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"TnBG.ratio",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seTnBG.ratio",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("RSE of Tn: ", ifelse(!is.null(criteria), 
                             round(criteria[i,"rseTn",drop=TRUE],2L), "NULL"), " (%)",sep=""),
                       paste("Fast ratio of Tn: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"FR",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seFR",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       "========================",
                       paste("Method: ", ifelse(method=="SGC","original SGC","improved SGC"), sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""),
                       "========================",
                       paste("ED: ",round(res$outDose[1L],2L), " +/- ",round(res$outDose[2L],2L)," (Gy|s)",sep=""), 
                       paste("RSE of ED: ",round(res$outDose[2L]/abs(res$outDose[1L])*100.0,2L), " (%)",sep=""),
                       paste("95% interval: [",round(i_ConfInt[3L],2L),", ",round(i_ConfInt[4L],2L),"]", sep=""), 
                       paste("68% interval: [",round(i_ConfInt[1L],2L),", ",round(i_ConfInt[2L],2L),"]", sep=""),
                       "========================"),           
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } else if (message==1L) {
                legend("center",
                       legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                       "========================",
                       "Status: Ln/Tn saturated",
                       "========================",
                       paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[i,"Tn",drop=TRUE],2L), " +/- ", 
                             round(Tn[i,"seTn",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("Tn above 3 sigma BG: ", ifelse(!is.null(criteria), 
                             as.logical(criteria[i,"Tn3BG",drop=TRUE]), "NULL"), sep=""),
                       paste("Ratio of Tn to BG: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"TnBG.ratio",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seTnBG.ratio",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("RSE of Tn: ", ifelse(!is.null(criteria), 
                             round(criteria[i,"rseTn",drop=TRUE],2L), "NULL"), " (%)",sep=""),
                       paste("Fast ratio of Tn: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"FR",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seFR",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       "========================",
                       paste("Method: ", ifelse(method=="SGC","original SGC","improved SGC"), sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""),
                       "========================",
                       paste("ED: ",Inf, " +/- ",NA," (Gy|s)",sep=""), 
                       "RSE of ED: NA (%)",  
                       "95% interval: [-Inf, +Inf]", 
                       "68% interval: [-Inf, +Inf]",  
                       "========================"),               
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } else if (message==2L) {
                legend("center",
                       legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                       "========================",
                       "Status: ED failed (ED < -50 Gy)",
                       "========================",
                       paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[i,"Tn",drop=TRUE],2L), " +/- ", 
                             round(Tn[i,"seTn",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("Tn above 3 sigma BG: ", ifelse(!is.null(criteria), 
                             as.logical(criteria[i,"Tn3BG",drop=TRUE]), "NULL"), sep=""),
                       paste("Ratio of Tn to BG: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"TnBG.ratio",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seTnBG.ratio",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("RSE of Tn: ", ifelse(!is.null(criteria), 
                             round(criteria[i,"rseTn",drop=TRUE],2L), "NULL"), " (%)",sep=""),
                       paste("Fast ratio of Tn: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"FR",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seFR",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       "========================",
                       paste("Method: ", ifelse(method=="SGC","original SGC","improved SGC"), sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""), 
                       "========================",               
                       paste("ED: ",NA, " +/- ",NA," (Gy|s)",sep=""), 
                       "RSE of ED: NA (%)", 
                       "95% interval: [NA, NA]", 
                       "68% interval: [NA, NA]",   
                       "========================"),             
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } else if (message==3L) {
                legend("center", 
                       legend=c(paste("ID: ", NO_Position_Grain, sep=""),
                       "========================",
                       "Status: ED Error failed (infinite upper ED)",
                       "========================",
                       paste("Tn: ", if(!is.null(Tn)) paste(round(Tn[i,"Tn",drop=TRUE],2L), " +/- ", 
                             round(Tn[i,"seTn",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("Tn above 3 sigma BG: ", ifelse(!is.null(criteria), 
                             as.logical(criteria[i,"Tn3BG",drop=TRUE]), "NULL"), sep=""),
                       paste("Ratio of Tn to BG: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"TnBG.ratio",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seTnBG.ratio",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       paste("RSE of Tn: ", ifelse(!is.null(criteria), 
                             round(criteria[i,"rseTn",drop=TRUE],2L), "NULL"), " (%)",sep=""),
                       paste("Fast ratio of Tn: ", if(!is.null(criteria)) 
                             paste(round(criteria[i,"FR",drop=TRUE],2L), " +/- ", 
                             round(criteria[i,"seFR",drop=TRUE],2L), sep="") else "NULL", sep=""),
                       "========================",
                       paste("Method: ", ifelse(method=="SGC","original SGC","improved SGC"), sep=""),
                       paste("Fit model: ", model, sep=""),
                       paste("Pass origin: ", origin, sep=""),
                       "========================",
                       paste("ED: ",round(res$outDose[1L],2L), "+/-",NA," (Gy|s)",sep=""), 
                       "RSE of ED: NA (%)", 
                       "95% interval: [NA, NA]", 
                       "68% interval: [NA, NA]",   
                       "========================"),            
                       yjust=2, ncol=1L, cex=1.0, bty="n")
            } # end if.
            ###
        } # end if.
        ###===================================================================================
        ###
    } # end for.
    ###
    ###
    if (!is.null(outpdf)) dev.off()
    ###
    if (!is.null(acceptNO)) {
        SGCED.table <- data.frame("NO"=acceptNO, "Position"=acceptPosition, "Grain"=acceptGrain,
                       "Tn"=Tn_vec, "seTn"=seTn_vec, "Tn3BG"=Tn3BG_vec, "TnBG.ratio"=TnBG.ratio_vec, 
                       "seTnBG.ratio"=seTnBG.ratio_vec, "rseTn"=rseTn_vec, "FR"=FR_vec, "seFR"=seFR_vec, 
                       "ScaleLtx"=ScaleLtx_vec, "seScaleLtx"=seScaleLtx_vec,
                       "rseED"=rseED_vec,"ED"=ED_vec, "seED"=seED_vec, 
                       "lower68"=lower68_vec, "upper68"=upper68_vec, 
                       "lower95"=lower95_vec, "upper95"=upper95_vec,
                       stringsAsFactors=FALSE)
        ###
        agID <- cbind("NO"=acceptNO, "Position"=acceptPosition, "Grain"=acceptGrain)
        ###

        ###
        ### No se consideration.
        if (!is.null(rseED.up)) {
            all_value <- SGCED.table[,"rseED",drop=TRUE]
            ###
            select_index <- which(abs(all_value)<rseED.up)
            ###
            if (length(select_index)==0L) {
                stop("Error: no acceptable SGC ED if [rseED.up] is applied!")
            } # end if.
            ###
            reject_N <- nrow(SGCED.table) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: RSE of ED exceeds ", rseED.up, "%",sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            SGCED.table <- SGCED.table[select_index,,drop=FALSE]
            rseED_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            rseED_reject <- NULL
        } # end if.
        ###
        ###
        if (!is.null(outfile)) write.csv(SGCED.table, file=paste(outfile,".csv",sep=""))
        ###
        scaleLtx <- as.matrix(SGCED.table[,c("ScaleLtx","seScaleLtx"),drop=FALSE])
        rownames(scaleLtx) <- paste("NO",SGCED.table[,"NO",drop=TRUE],sep="")
        ###
        sgcED <- as.matrix(SGCED.table[,c("ED","seED"),drop=FALSE])
        rownames(sgcED) <- paste("NO",SGCED.table[,"NO",drop=TRUE],sep="")
        ###
        ConfInt <- as.matrix(SGCED.table[,c("lower68","upper68","lower95","upper95"),drop=FALSE])
        rownames(ConfInt) <- paste("NO",SGCED.table[,"NO",drop=TRUE],sep="")
        ### 
        output <- list("scale.Ltx"=scaleLtx,
                       "sgcED"=sgcED,
                       "ConfInt"=ConfInt,
                       "agID"=agID)
    } # end if.
    ###
    ###------------------------------------------------------
    if (is_forced_object==FALSE) { 
        if (length(Tn3BG_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [Tn.above.3BG]:\n")
            print(Tn3BG_reject)
            cat("\n")
        } # end if.
        ###
        if (length(TnBG.ratio_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [TnBG.ratio]:\n")
            print(TnBG.ratio_reject)
            cat("\n")
        } # end if.
        ###
        if (length(rseTn_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rseTn]:\n")
            print(rseTn_reject)
            cat("\n")
        } # end if.
        ### 
        if (length(FR_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [FR]:\n")
            print(FR_reject)
            cat("\n")
        } # end if.
        ###
    } # end if.
    ###------------------------------------------------------
    ###
    ###------------------------------------------------------
    if (!is.null(acceptNO)) {
        if (length(rseED_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rseED]:\n")
            print(rseED_reject)
            cat("\n")
        } # end if.
    } # end if.
    ###--------------------------------------------------------
    ###
    if (!is.null(saturate_ID)) {
        cat("\n")
        cat("Aliquot (grain) ID saturated in Ln/Tn:\n")
        print(apply(saturate_ID, MARGIN=1L, NPG))
        cat("\n")
    } # end if.
    ###
    if (!is.null(failED_ID)) {
        cat("\n")
        cat("Aliquot (grain) ID failed in ED calculation:\n")
        print(apply(failED_ID, MARGIN=1L, NPG))
        cat("\n")
    } # end if.
    ###
    ###
    if (!is.null(failEDError_ID)) {
        cat("\n")
        cat("Aliquot (grain) ID failed in ED error estimation:\n")
        print(apply(failEDError_ID, MARGIN=1L, NPG))
        cat("\n")
    } # end if.
    ###
    action_character <- c(action_character,
                          "Saturated in Ln/Tn",
                          "Failed in ED calculation",
                          "Failed in ED error estimation",
                          "Total number of rejected aliquots (grains)",
                          "Total number of accepted aliquots (grains)")
    ###
    step_reject_N <- c(step_reject_N,
                       ifelse(is.null(saturate_ID), 0L, nrow(saturate_ID)),
                       ifelse(is.null(failED_ID), 0L, nrow(failED_ID)),
                       ifelse(is.null(failEDError_ID), 0L, nrow(failEDError_ID)),
                       ifelse(is.null(acceptNO), nag, nag-nrow(agID)), 
                       ifelse(is.null(acceptNO), 0L, nrow(agID)))
    ###
    summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
    print(summary_info)
    ###
    if (!is.null(acceptNO)) {
        return(invisible(output))
    } else {
        return(invisible(NULL))
    } # end if.
    ###
} # end function calSGCED.default.
#####
