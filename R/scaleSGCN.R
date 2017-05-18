#####
scaleSGCN <-
function(obj_analyseBIN, SGCpars, model, origin, SAR.Cycle, 
         Tn.above.3BG=TRUE, TnBG.ratio.low=NULL, rseTn.up=NULL, 
         FR.low=NULL, use.se=TRUE, outfile=NULL) {
    UseMethod("scaleSGCN")
} ###
### 2017.05.18.
scaleSGCN.default <-
function(obj_analyseBIN, SGCpars, model, origin, SAR.Cycle, 
         Tn.above.3BG=TRUE, TnBG.ratio.low=NULL, rseTn.up=NULL, 
         FR.low=NULL, use.se=TRUE, outfile=NULL) {
    ### Stop if not.
    stopifnot(class(obj_analyseBIN)=="analyseBIN", 
              names(obj_analyseBIN)==c("SARdata","criteria","Tn","LnTn.curve","TxTn","agID"),
              is.numeric(SGCpars), 
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(SAR.Cycle)==2L, all(substr(SAR.Cycle,1L,1L) %in% c("N","R")),
              length(Tn.above.3BG)==1L, is.logical(Tn.above.3BG),
              is.null(TnBG.ratio.low) || (length(TnBG.ratio.low)==1L && is.numeric(TnBG.ratio.low)),
              is.null(rseTn.up) || (length(rseTn.up)==1L && is.numeric(rseTn.up)),
              is.null(FR.low) || (length(FR.low)==1L && is.numeric(FR.low)),
              length(use.se)==1L, is.logical(use.se),
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
    ###
    SARdata <- obj_analyseBIN$SARdata
    criteria <- obj_analyseBIN$criteria
    Tn <- obj_analyseBIN$Tn
    LnTn.curve <- obj_analyseBIN$LnTn.curve
    TxTn <- obj_analyseBIN$TxTn
    agID <- obj_analyseBIN$agID
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
        if (all(iSAR.Cycle=="N")) {
            stop(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                       "]: should contain 'SAR.Cycle' of 'R'!", sep=""))
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
                stop("Error: no acceptable scaled natrual signal if [Tn.above.3BG] is applied!")
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
                stop("Error: no acceptable scaled natrual signal if [TnBG.ratio.low] is applied!")
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
                stop("Error: no acceptable scaled natrual signal if [rseTn.up] is applied!")
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
                stop("Error: no acceptable scaled natrual signal if [FR.low] is applied!")
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
    acceptNO <- acceptPosition <- acceptGrain <-  
    Tn_vec <- seTn_vec <- Tn3BG_vec <- TnBG.ratio_vec <- 
    seTnBG.ratio_vec <- rseTn_vec <- FR_vec <- seFR_vec <- 
    ScaleLtx_vec <- seScaleLtx_vec <- c()
    ###
    for (i in seq(n)) {
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i])  
        ith_SARdata <- SARdata[iIndex,,drop=FALSE]
        ###
        select_index <- which(ith_SARdata[,"SAR.Cycle",drop=TRUE] %in% SAR.Cycle)
        ###
        if (length(select_index)==2L) {
            ###
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
            ith_SGCdata <- ith_SARdata[select_index,,drop=FALSE]
            iSAR.Cycle <- substr(ith_SGCdata[,"SAR.Cycle",drop=TRUE],start=1L,stop=1L)
            ###
            naturalSignal      <- ith_SGCdata[iSAR.Cycle=="N","Signal",drop=TRUE]
            naturalSignalError <- ith_SGCdata[iSAR.Cycle=="N","Signal.Err",drop=TRUE]
            ###
            regenerativeDose   <- ith_SGCdata[iSAR.Cycle=="R","Dose",drop=TRUE]
            regenerativeSignal <- ith_SGCdata[iSAR.Cycle=="R","Signal",drop=TRUE]
            ###
            scalingFactor <- fcn(regenerativeDose)/regenerativeSignal
            ###
            ScaleLtx_vec <- c(ScaleLtx_vec, scalingFactor*naturalSignal)
            seScaleLtx_vec <- c(seScaleLtx_vec, scalingFactor*naturalSignalError)
            ###
        } else {
             cat(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                       "]: desired 'SAR.cycle' cannot be extracted (omitted)!\n", sep=""))
        } # end if.
    } # end for.
    ###
    ###
    if (!is.null(acceptNO)) {
        ###
        scaleSGCN.table <- data.frame("NO"=acceptNO, "Position"=acceptPosition, "Grain"=acceptGrain,
                           "Tn"=Tn_vec, "seTn"=seTn_vec, "Tn3BG"=Tn3BG_vec, "TnBG.ratio"=TnBG.ratio_vec, 
                           "seTnBG.ratio"=seTnBG.ratio_vec, "rseTn"=rseTn_vec, "FR"=FR_vec, "seFR"=seFR_vec, 
                           "ScaleLtx"=ScaleLtx_vec, "seScaleLtx"=seScaleLtx_vec,
                           stringsAsFactors=FALSE)
        ### Re-write agID.
        agID <- cbind("NO"=acceptNO, "Position"=acceptPosition, "Grain"=acceptGrain)
        ###
        ###
        if (!is.null(outfile)) write.csv(scaleSGCN.table, file=paste(outfile,".csv",sep=""))
        ###
        scaleLtx <- as.matrix(scaleSGCN.table[,c("ScaleLtx","seScaleLtx"),drop=FALSE])
        rownames(scaleLtx) <- paste("NO",scaleSGCN.table[,"NO",drop=TRUE],sep="")
        ###
        output <- list("scale.Ltx"=scaleLtx,
                       "agID"=agID)
    } # end if.
    ###
    ###
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
    ###
    action_character <- c(action_character,
                          "Total number of rejected aliquots (grains)",
                          "Total number of accepted aliquots (grains)")
     ###
    step_reject_N <- c(step_reject_N,
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
} # end function scaleSGCN.default.
#####
