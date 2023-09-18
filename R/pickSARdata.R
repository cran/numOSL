#####
pickSARdata <-
function(obj_analyseBIN, model="gok", origin=FALSE, weight=TRUE, trial=TRUE, 
         nofit.rgd=NULL, Tn.above.3BG=TRUE, TnBG.ratio.low=NULL, rseTn.up=NULL,
         FR.low=NULL, rcy1.range=NULL, rcy2.range=NULL, rcy3.range=NULL, 
         rcp1.up=NULL, rcp2.up=NULL, fom.up=NULL, rcs.up=NULL, 
         use.se=TRUE, norm.dose=NULL, outpdf=NULL, outfile=NULL) {
    UseMethod("pickSARdata")
} #
### 2023.09.09.
pickSARdata.default <-
function(obj_analyseBIN, model="gok", origin=FALSE, weight=TRUE, trial=TRUE, 
         nofit.rgd=NULL, Tn.above.3BG=TRUE, TnBG.ratio.low=NULL, rseTn.up=NULL,
         FR.low=NULL, rcy1.range=NULL, rcy2.range=NULL, rcy3.range=NULL, 
         rcp1.up=NULL, rcp2.up=NULL, fom.up=NULL, rcs.up=NULL, 
         use.se=TRUE, norm.dose=NULL, outpdf=NULL, outfile=NULL) {
    ### Stop if not.
    stopifnot(class(obj_analyseBIN)=="analyseBIN", 
              names(obj_analyseBIN)==c("SARdata","criteria","Tn","LnTn.curve","TxTn","agID"),
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(weight)==1L, is.logical(weight),
              length(trial)==1L, is.logical(trial),
              is.null(nofit.rgd) || is.numeric(nofit.rgd),
              length(Tn.above.3BG)==1L, is.logical(Tn.above.3BG),
              is.null(TnBG.ratio.low) || (length(TnBG.ratio.low)==1L && is.numeric(TnBG.ratio.low)),
              is.null(rseTn.up) || (length(rseTn.up)==1L && is.numeric(rseTn.up)),
              is.null(FR.low) || (length(FR.low)==1L && is.numeric(FR.low)), 
              is.null(rcy1.range) || (length(rcy1.range)==2L && is.numeric(rcy1.range)),
              is.null(rcy2.range) || (length(rcy2.range)==2L && is.numeric(rcy2.range)),
              is.null(rcy3.range) || (length(rcy3.range)==2L && is.numeric(rcy3.range)),
              is.null(rcp1.up) || (length(rcp1.up)==1L && is.numeric(rcp1.up)),
              is.null(rcp2.up) || (length(rcp2.up)==1L && is.numeric(rcp2.up)),
              is.null(fom.up) || (length(fom.up)==1L && is.numeric(fom.up)),
              is.null(rcs.up) || (length(rcs.up)==1L && is.numeric(rcs.up)), 
              length(use.se)==1L, is.logical(use.se),
              is.null(norm.dose) || (length(norm.dose)==1L && is.numeric(norm.dose)),  
              is.null(outpdf) || (length(outpdf)==1L && is.character(outpdf)),
              is.null(outfile) || (length(outfile)==1L && is.character(outfile)))
    ###
    SARdata <- obj_analyseBIN$SARdata
    criteria <- obj_analyseBIN$criteria
    Tn <- obj_analyseBIN$Tn
    LnTn.curve <- obj_analyseBIN$LnTn.curve
    TxTn <- obj_analyseBIN$TxTn
    agID <- obj_analyseBIN$agID
    ###
    if (nrow(SARdata)<5L) stop("Error: need more data points!")
    ###
    colnames(SARdata) <- c("NO","SAR.Cycle","Dose","Signal","Signal.Err")
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
        if (sum(iSAR.Cycle=="N")>=2L) {
           stop(paste("[NO=", NO[i], ",Position=",Position[i], ",Grain=",Grain[i],
                      "]: should not contain more than one 'SAR.Cycle' of 'N'!", sep=""))
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

    
    ### Apply signal related rejection criteria. 
    ###*************************************************************************************************
    if (is_forced_object==FALSE) { 

        ###
        ###-----------------------------------------------------------
        ### Tn above 3sigma BG: no se consideration.
        if (Tn.above.3BG==TRUE) {
            all_value <- criteria[,"Tn3BG",drop=TRUE]
            ###
            select_index <- which(all_value==1L)
            ###
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                "Rejection criterion: Tn below 3 sigma BG")
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            if (length(select_index)==0L) {
                ###
                cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
                ###
                action_character <- c(action_character, 
                    "Total number of rejected aliquots (grains)",
                    "Total number of accepted aliquots (grains)")
                ###
                step_reject_N <- c(step_reject_N, nag, 0L)
                ###
                summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
                ###
                print(summary_info)
                ###
                return(invisible(summary_info))
            } # end if.
            ###
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            TxTn <- TxTn[select_index]
            ###
            Tn3BG_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            Tn3BG_reject <- NULL
        } # end if.
        ###-----------------------------------------------------------
        ###

        ###
        ###-----------------------------------------------------------
        ### Ratio of Tn to BG: have se consideration.
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
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: ratio of Tn to BG below ", TnBG.ratio.low, sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            if (length(select_index)==0L) {
                ###
                cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
                ###
                action_character <- c(action_character, 
                    "Total number of rejected aliquots (grains)",
                    "Total number of accepted aliquots (grains)")
                ###
                step_reject_N <- c(step_reject_N, nag, 0L)
                ###
                summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
                ###
                print(summary_info)
                ###
                return(invisible(summary_info))
            } # end if.
            ###
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            TxTn <- TxTn[select_index]
            ###
            TnBG.ratio_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            TnBG.ratio_reject <- NULL
        } # end if.
        ###-----------------------------------------------------------
        ###

        ###
        ###-----------------------------------------------------------
        ### Relative standard error of Tn: no se consideration.
        if (!is.null(rseTn.up)) {
            all_value <- criteria[,"rseTn",drop=TRUE]
            ###
            select_index <- which(abs(all_value)<rseTn.up)
            ###
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: RSE of Tn exceeds ", rseTn.up, "%", sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            if (length(select_index)==0L) {
                ###
                cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
                ###
                action_character <- c(action_character, 
                    "Total number of rejected aliquots (grains)",
                    "Total number of accepted aliquots (grains)")
                ###
                step_reject_N <- c(step_reject_N, nag, 0L)
                ###
                summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
                ###
                print(summary_info)
                ###
                return(invisible(summary_info))
            } # end if.
            ### 
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            TxTn <- TxTn[select_index]
            ###
            rseTn_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            rseTn_reject <- NULL
        } # end if.
        ###-----------------------------------------------------------
        ###

        ###
        ###-----------------------------------------------------------
        ### Fast ratio: have se consideration.
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
            reject_N <- nrow(criteria) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: fast ratio of Tn below ", FR.low, sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            if (length(select_index)==0L) {
                ###
                cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
                ###
                action_character <- c(action_character, 
                    "Total number of rejected aliquots (grains)",
                    "Total number of accepted aliquots (grains)")
                ###
                step_reject_N <- c(step_reject_N, nag, 0L)
                ###
                summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
                ###
                print(summary_info)
                ###
                return(invisible(summary_info))
            } # end if.
            ### 
            criteria <- criteria[select_index,,drop=FALSE]
            Tn <- Tn[select_index,,drop=FALSE]
            LnTn.curve <- LnTn.curve[select_index]
            TxTn <- TxTn[select_index]
            ###
            FR_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
        } else {
            FR_reject <- NULL
        } # end if.
        ###-----------------------------------------------------------
        ###
    } else {
        if ((Tn.above.3BG==TRUE) || (!is.null(TnBG.ratio.low)) ||
            (!is.null(rseTn.up)) || (!is.null(FR.low))) {
            cat("Note: signal-related rejection criteria cannot be applied!\n")
        } # end if.
    } # end if.
    ###*************************************************************************************************
    ###

    ###
    ### Re-write NO, Position, and Grain.
    NO <- agID[,"NO",drop=TRUE]
    Position <- agID[,"Position",drop=TRUE]
    Grain <- agID[,"Grain",drop=TRUE]
    n <- length(NO)
    ###
    RcyRcp_mat <- matrix(nrow=n, ncol=10L)
    colnames(RcyRcp_mat) <- c("RecyclingRatio1", "seRecyclingRatio1",
        "RecyclingRatio2", "seRecyclingRatio2", "RecyclingRatio3", "seRecyclingRatio3",
        "Recuperation1", "seRecuperation1", "Recuperation2", "seRecuperation2")
    ###
    for (i in seq(n)) {
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i]) 
        ### 
        ith_SARdata <- SARdata[iIndex,,drop=FALSE]
        ###
        Data4L <- ith_SARdata[,c("SAR.Cycle","Dose","Signal","Signal.Err"),drop=FALSE]
        DataR <- Data4L[Data4L[,"SAR.Cycle",drop=TRUE]!="N",c("Dose","Signal","Signal.Err"),drop=FALSE]
        ###
        index_SAR.Cycle_N <- which(Data4L[,"SAR.Cycle",drop=TRUE]=="N")
        if (length(index_SAR.Cycle_N)==0L) {
            DataN <- rep(NA, 2L)
        } else {
            DataN <- as.numeric(Data4L[index_SAR.Cycle_N,c("Signal","Signal.Err"),drop=FALSE])
        } # end if.
        ###
        res_calRcyRcp <- calRcyRcp(Curvedata=DataR, Ltx=DataN)
        RcyRcp_mat[i,1L:2L] <- res_calRcyRcp$RecyclingRatio1
        RcyRcp_mat[i,3L:4L] <- res_calRcyRcp$RecyclingRatio2
        RcyRcp_mat[i,5L:6L] <- res_calRcyRcp$RecyclingRatio3
        RcyRcp_mat[i,7L:8L] <- res_calRcyRcp$Recuperation1
        RcyRcp_mat[i,9L:10L] <- res_calRcyRcp$Recuperation2  
    } # end for.
    ###

    ### Apply growth curve related rejection criteria (1).
    ###**********************************************************************************************************************

    ###
    ###---------------------------------------------------------------------------------------------
    ### The first recycling ratio: have se consideration. 
    if (!is.null(rcy1.range)) {
        all_value <- RcyRcp_mat[,"RecyclingRatio1",drop=TRUE]
        ### 
        if (use.se==FALSE) {
            select_index <- which((all_value>rcy1.range[1L] & all_value<rcy1.range[2L]) | is.na(all_value))
        } else {
            all_se_value <- RcyRcp_mat[,"seRecyclingRatio1",drop=TRUE]
            select_index <- which((all_value-sigma*all_se_value<rcy1.range[1L] & all_value+sigma*all_se_value>rcy1.range[2L]) |
                                  (all_value-sigma*all_se_value>rcy1.range[1L] & all_value+sigma*all_se_value<rcy1.range[2L]) |
                                  (all_value-sigma*all_se_value>rcy1.range[1L] & all_value-sigma*all_se_value<rcy1.range[2L]) |
                                  (all_value+sigma*all_se_value>rcy1.range[1L] & all_value+sigma*all_se_value<rcy1.range[2L]) |
                                  is.na(all_value))
        } # end if.
        ###
        reject_N <- nrow(RcyRcp_mat) - length(select_index)
        ###
        action_character <- c(action_character, 
            paste("Rejection criterion: recycling ratio 1 outsides [",
            rcy1.range[1L], ",",rcy1.range[2L],"]", sep=""))
        step_reject_N <- c(step_reject_N, reject_N)
        ###
        if (length(select_index)==0L) {
            ###
            cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
            ###
            action_character <- c(action_character, 
                "Total number of rejected aliquots (grains)",
                "Total number of accepted aliquots (grains)")
            ###
            step_reject_N <- c(step_reject_N, nag, 0L)
            ###
            summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
            ###
            print(summary_info)
            ###
            return(invisible(summary_info))
        } # end if.
        ###
        RcyRcp_mat <- RcyRcp_mat[select_index,,drop=FALSE]
        ###
        if (!is.null(criteria)) criteria <- criteria[select_index,,drop=FALSE]
        if (!is.null(Tn)) Tn <- Tn[select_index,,drop=FALSE]
        if (!is.null(LnTn.curve)) LnTn.curve <- LnTn.curve[select_index]
        if (!is.null(TxTn)) TxTn <- TxTn[select_index]
        ###
        rcy1_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
        agID <- agID[select_index,,drop=FALSE]
    } else {
        rcy1_reject <- NULL
    } # end if.
    ###---------------------------------------------------------------------------------------------
    ###

    ###
    ###---------------------------------------------------------------------------------------------
    ### The second recycling ratio: have se consideration.
    if (!is.null(rcy2.range)) {
        all_value <- RcyRcp_mat[,"RecyclingRatio2",drop=TRUE]
        ###
        if (use.se==FALSE) {
            select_index <- which((all_value>rcy2.range[1L] & all_value<rcy2.range[2L]) | is.na(all_value))
        } else {
            all_se_value <- RcyRcp_mat[,"seRecyclingRatio2",drop=TRUE]
            select_index <- which((all_value-sigma*all_se_value<rcy2.range[1L] & all_value+sigma*all_se_value>rcy2.range[2L]) |
                                  (all_value-sigma*all_se_value>rcy2.range[1L] & all_value+sigma*all_se_value<rcy2.range[2L]) |
                                  (all_value-sigma*all_se_value>rcy2.range[1L] & all_value-sigma*all_se_value<rcy2.range[2L]) |
                                  (all_value+sigma*all_se_value>rcy2.range[1L] & all_value+sigma*all_se_value<rcy2.range[2L]) |
                                  is.na(all_value))
        } # end if.
        ###
        reject_N <- nrow(RcyRcp_mat) - length(select_index)
        ###
        action_character <- c(action_character, 
            paste("Rejection criterion: recycling ratio 2 outsides [",
            rcy2.range[1L], ",",rcy2.range[2L],"]", sep=""))
        step_reject_N <- c(step_reject_N, reject_N)
        ###
        if (length(select_index)==0L) {
            ###
            cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
            ###
            action_character <- c(action_character, 
                "Total number of rejected aliquots (grains)",
                "Total number of accepted aliquots (grains)")
            ###
            step_reject_N <- c(step_reject_N, nag, 0L)
            ###
            summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
            ###
            print(summary_info)
            ###
            return(invisible(summary_info))
        } # end if.
        ###
        RcyRcp_mat <- RcyRcp_mat[select_index,,drop=FALSE]
        ###
        if (!is.null(criteria)) criteria <- criteria[select_index,,drop=FALSE]
        if (!is.null(Tn)) Tn <- Tn[select_index,,drop=FALSE]
        if (!is.null(LnTn.curve)) LnTn.curve <- LnTn.curve[select_index]
        if (!is.null(TxTn)) TxTn <- TxTn[select_index]
        ###
        rcy2_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
        agID <- agID[select_index,,drop=FALSE]
    } else {
        rcy2_reject <- NULL
    } # end if.
    ###---------------------------------------------------------------------------------------------
    ###

    ###
    ###---------------------------------------------------------------------------------------------
    ### The third recycling ratio: have se consideration.
    if (!is.null(rcy3.range)) {
        all_value <- RcyRcp_mat[,"RecyclingRatio3",drop=TRUE]
        ###
        if (use.se==FALSE) {
            select_index <- which((all_value>rcy3.range[1L] & all_value<rcy3.range[2L]) | is.na(all_value))
        } else {
            all_se_value <- RcyRcp_mat[,"seRecyclingRatio3",drop=TRUE]
            select_index <- which((all_value-sigma*all_se_value<rcy3.range[1L] & all_value+sigma*all_se_value>rcy3.range[2L]) |
                                  (all_value-sigma*all_se_value>rcy3.range[1L] & all_value+sigma*all_se_value<rcy3.range[2L]) |
                                  (all_value-sigma*all_se_value>rcy3.range[1L] & all_value-sigma*all_se_value<rcy3.range[2L]) |
                                  (all_value+sigma*all_se_value>rcy3.range[1L] & all_value+sigma*all_se_value<rcy3.range[2L]) |
                                  is.na(all_value))
        } # end if.
        ###
        reject_N <- nrow(RcyRcp_mat) - length(select_index)
        ###
        action_character <- c(action_character, 
            paste("Rejection criterion: recycling ratio 3 outsides [",
            rcy3.range[1L], ",",rcy3.range[2L],"]", sep=""))
        step_reject_N <- c(step_reject_N, reject_N)
        ###
        if (length(select_index)==0L) {
            ###
            cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
            ###
            action_character <- c(action_character, 
                "Total number of rejected aliquots (grains)",
                "Total number of accepted aliquots (grains)")
            ###
            step_reject_N <- c(step_reject_N, nag, 0L)
            ###
            summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
            ###
            print(summary_info)
            ###
            return(invisible(summary_info))
        } # end if.
        ###
        RcyRcp_mat <- RcyRcp_mat[select_index,,drop=FALSE]
        ###
        if (!is.null(criteria)) criteria <- criteria[select_index,,drop=FALSE]
        if (!is.null(Tn)) Tn <- Tn[select_index,,drop=FALSE]
        if (!is.null(LnTn.curve)) LnTn.curve <- LnTn.curve[select_index]
        if (!is.null(TxTn)) TxTn <- TxTn[select_index]
        ###
        rcy3_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
        agID <- agID[select_index,,drop=FALSE]
    } else {
        rcy3_reject <- NULL
    } # end if.
    ###---------------------------------------------------------------------------------------------
    ###

    ###
    ###---------------------------------------------------------------------------------------------
    ### The first recuperation: Have se consideration.
    if (!is.null(rcp1.up)) {
        all_value <- RcyRcp_mat[,"Recuperation1",drop=TRUE]
        ###
        if (use.se==FALSE) {
            select_index <- which((all_value> -rcp1.up & all_value<rcp1.up) | is.na(all_value))
        } else {
            all_se_value <- RcyRcp_mat[,"seRecuperation1",drop=TRUE]
            select_index <- which((all_value-sigma*all_se_value< -rcp1.up & all_value+sigma*all_se_value>rcp1.up) |
                                  (all_value-sigma*all_se_value> -rcp1.up & all_value+sigma*all_se_value<rcp1.up) |
                                  (all_value-sigma*all_se_value> -rcp1.up & all_value-sigma*all_se_value<rcp1.up) |
                                  (all_value+sigma*all_se_value> -rcp1.up & all_value+sigma*all_se_value<rcp1.up) |
                                  is.na(all_value))
        } # end if.
        ###
        reject_N <- nrow(RcyRcp_mat) - length(select_index)
        ###
        action_character <- c(action_character, 
            paste("Rejection criterion: recuperation 1 exceeds ", rcp1.up, "%",sep=""))
        step_reject_N <- c(step_reject_N, reject_N)
        ###
        if (length(select_index)==0L) {
            ###
            cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
            ###
            action_character <- c(action_character, 
                "Total number of rejected aliquots (grains)",
                "Total number of accepted aliquots (grains)")
            ###
            step_reject_N <- c(step_reject_N, nag, 0L)
            ###
            summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
            ###
            print(summary_info)
            ###
            return(invisible(summary_info))
        } # end if.
        ###
        RcyRcp_mat <- RcyRcp_mat[select_index,,drop=FALSE]
        ###
        if (!is.null(criteria)) criteria <- criteria[select_index,,drop=FALSE]
        if (!is.null(Tn)) Tn <- Tn[select_index,,drop=FALSE]
        if (!is.null(LnTn.curve)) LnTn.curve <- LnTn.curve[select_index]
        if (!is.null(TxTn)) TxTn <- TxTn[select_index]
        ###
        rcp1_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
        agID <- agID[select_index,,drop=FALSE]
    } else {
        rcp1_reject <- NULL
    } # end if.
    ###---------------------------------------------------------------------------------------------

    ###
    ###---------------------------------------------------------------------------------------------
    ### The second recuperation: have se consideration.
    if (!is.null(rcp2.up)) {
        all_value <- RcyRcp_mat[,"Recuperation2",drop=TRUE]
        ###
        if (use.se==FALSE) {
            select_index <- which((all_value> -rcp2.up & all_value<rcp2.up) | is.na(all_value))
        } else {
            all_se_value <- RcyRcp_mat[,"seRecuperation2",drop=TRUE]
            select_index <- which((all_value-sigma*all_se_value< -rcp2.up & all_value+sigma*all_se_value>rcp2.up) |
                                  (all_value-sigma*all_se_value> -rcp2.up & all_value+sigma*all_se_value<rcp2.up) |
                                  (all_value-sigma*all_se_value> -rcp2.up & all_value-sigma*all_se_value<rcp2.up) |
                                  (all_value+sigma*all_se_value> -rcp2.up & all_value+sigma*all_se_value<rcp2.up) |
                                  is.na(all_value))
        } # end if.
        ###
        reject_N <- nrow(RcyRcp_mat) - length(select_index)
        ###
        action_character <- c(action_character, 
            paste("Rejection criterion: recuperation 2 exceeds ", rcp2.up, "%",sep=""))
        step_reject_N <- c(step_reject_N, reject_N)
        ###
        if (length(select_index)==0L) {
            ###
            cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
            ###
            action_character <- c(action_character, 
                "Total number of rejected aliquots (grains)",
                "Total number of accepted aliquots (grains)")
            ###
            step_reject_N <- c(step_reject_N, nag, 0L)
            ###
            summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
            ###
            print(summary_info)
            ###
            return(invisible(summary_info))
        } # end if.
        ###
        RcyRcp_mat <- RcyRcp_mat[select_index,,drop=FALSE]
        ###
        if (!is.null(criteria)) criteria <- criteria[select_index,,drop=FALSE]
        if (!is.null(Tn)) Tn <- Tn[select_index,,drop=FALSE]
        if (!is.null(LnTn.curve)) LnTn.curve <- LnTn.curve[select_index]
        if (!is.null(TxTn)) TxTn <- TxTn[select_index]
        ###
        rcp2_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
        agID <- agID[select_index,,drop=FALSE]
    } else {
        rcp2_reject <- NULL
    } # end if.
    ###------------------------------------------------------------------------------------
    ###
    ###**********************************************************************************************************************
    ###

    ### Re-write NO, Position, and Grain.
    NO <- agID[,"NO",drop=TRUE]
    Position <- agID[,"Position",drop=TRUE]
    Grain <- agID[,"Grain",drop=TRUE]
    n <- length(NO)
    ###
    DataList <- vector(mode="list", length=n)
    ###
    for (i in seq(n)) {
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i])  
        DataList[[i]] <- SARdata[iIndex,,drop=FALSE]
    } # end for.
    ###
    ###
    passNO <- passPosition <- passGrain <-  
    Tn_vec <- seTn_vec <- Tn3BG_vec <- TnBG.ratio_vec <- 
    seTnBG.ratio_vec <- rseTn_vec <- FR_vec <- seFR_vec <- 
    RecyclingRatio1_vec <- seRecyclingRatio1_vec <- 
    RecyclingRatio2_vec <- seRecyclingRatio2_vec <-
    RecyclingRatio3_vec <- seRecyclingRatio3_vec <-
    Recuperation1_vec <- seRecuperation1_vec <-
    Recuperation2_vec <- seRecuperation2_vec <- 
    tryError_ID <- failFit_ID <- FOM_vec <- RCS_vec <- c()
    ###
    LMpars <- list()
    ###
    norm.SARdata <- c()
    SARdata1 <- c()
    ###
    if (!is.null(outpdf)) {
        pdf(paste(outpdf, ".pdf", sep=""))
        if_plot <- TRUE
    } else {
        if_plot <- FALSE
    } # end if.
    ###
    if (n>=5L) {
        pb <- txtProgressBar(min=1L, max=n, initial=1L, char="=") 
        cat("Growth curve fitting is in progress, please wait, ...\n")
    } # end if.
    ###-------------------------------------------------------------------
    ###
    for (i in seq(n)) {
        if (n>=5L) setTxtProgressBar(pb, i)
        ###
        Data4L <- DataList[[i]][,c("SAR.Cycle","Dose","Signal","Signal.Err"),drop=FALSE]
        DataR <- Data4L[Data4L[,"SAR.Cycle",drop=TRUE]!="N",c("Dose","Signal","Signal.Err"),drop=FALSE]
        ###
        index_SAR.Cycle_N <- which(Data4L[,"SAR.Cycle",drop=TRUE]=="N")
        if (length(index_SAR.Cycle_N)==0L) {
            DataN <- rep(NA, 2L)
        } else {
            DataN <- as.numeric(Data4L[index_SAR.Cycle_N,c("Signal","Signal.Err"),drop=FALSE])
        } # end if.
        ###
        fitGrowth_Tn <- if (!is.null(Tn)) Tn[i,c("Tn","seTn"),drop=TRUE] else NULL
        fitGrowth_Tn3BG <- if (!is.null(criteria)) criteria[i,"Tn3BG",drop=TRUE] else NULL
        fitGrowth_TnBG.ratio <- if (!is.null(criteria)) criteria[i,c("TnBG.ratio","seTnBG.ratio"),drop=TRUE] else NULL
        fitGrowth_rseTn <- if (!is.null(criteria)) criteria[i,"rseTn",drop=TRUE] else NULL
        fitGrowth_FR <- if (!is.null(criteria)) criteria[i,c("FR","seFR"),drop=TRUE] else NULL
        ###
        fitGrowth_RecyclingRatio1 <- RcyRcp_mat[i,c("RecyclingRatio1","seRecyclingRatio1"),drop=TRUE]
        fitGrowth_RecyclingRatio2 <- RcyRcp_mat[i,c("RecyclingRatio2","seRecyclingRatio2"),drop=TRUE]
        fitGrowth_RecyclingRatio3 <- RcyRcp_mat[i,c("RecyclingRatio3","seRecyclingRatio3"),drop=TRUE]
        ###
        fitGrowth_Recuperation1 <- RcyRcp_mat[i,c("Recuperation1","seRecuperation1"),drop=TRUE]
        fitGrowth_Recuperation2 <- RcyRcp_mat[i,c("Recuperation2","seRecuperation2"),drop=TRUE]
        ###
        fitGrowth_LnTn.curve <- if (!is.null(LnTn.curve)) LnTn.curve[[i]] else NULL
        fitGrowth_TxTn <- if (!is.null(TxTn)) TxTn[[i]] else NULL
        ###
        ### Force the fitting to pass the origin if zero regeneraive-dose is not available.
        if (is.null(nofit.rgd)) {
            rgdddo <- DataR[,"Dose",drop=TRUE]
        } else {
            rgdddo <- DataR[-nofit.rgd,"Dose",drop=TRUE]
        } # end if.
        ###
        zeroDose_index <- which(abs(rgdddo)<=.Machine$double.eps^0.5)
        ###
        if (length(zeroDose_index)==0L && origin==FALSE) {
            origin1 <- TRUE 
            cat(paste("[NO=",NO[i],",Position=",Position[i],",Grain=",Grain[i],
                "]: no zero regenerative-dose can be used and the growth curve is forced to pass the origin!\n", sep=""))
        } else {
            origin1 <- origin
        } # end if.
        ###
        res <- try(fitGrowth(Curvedata=DataR, model=model, origin=origin1, 
                             weight=weight, trial=trial, plot=if_plot, 
                             nofit.rgd=nofit.rgd, agID=agID[i,,drop=TRUE], Tn=fitGrowth_Tn, 
                             Tn3BG=fitGrowth_Tn3BG, TnBG.ratio=fitGrowth_TnBG.ratio, 
                             rseTn=fitGrowth_rseTn, FR=fitGrowth_FR, 
                             RecyclingRatio1=fitGrowth_RecyclingRatio1, 
                             RecyclingRatio2=fitGrowth_RecyclingRatio2, 
                             RecyclingRatio3=fitGrowth_RecyclingRatio3,
                             Recuperation1=fitGrowth_Recuperation1,
                             Recuperation2=fitGrowth_Recuperation2,
                             LnTn.curve=fitGrowth_LnTn.curve,
                             TxTn=fitGrowth_TxTn), silent=TRUE)
        ###
        if (inherits(res,what="try-error")==TRUE) {
             tryError_ID <- rbind(tryError_ID, agID[i,,drop=TRUE])
             cat(paste("[NO=",NO[i],",Position=",Position[i],",Grain=",Grain[i],"]:\n", sep=""))
             print(attr(res, "condition"))
        } else {
            if (res$message==0L) {
                passNO <- c(passNO, NO[i])
                passPosition <- c(passPosition, Position[i])
                passGrain <- c(passGrain, Grain[i])
                ###
                Tn_vec <- c(Tn_vec, if(!is.null(Tn)) Tn[i,"Tn",drop=TRUE] else NA)
                seTn_vec <- c(seTn_vec, if(!is.null(Tn)) Tn[i,"seTn",drop=TRUE] else NA)
                ###
                Tn3BG_vec <- c(Tn3BG_vec, if(!is.null(criteria)) criteria[i,"Tn3BG",drop=TRUE] else NA)
                TnBG.ratio_vec <- c(TnBG.ratio_vec, if(!is.null(criteria)) criteria[i,"TnBG.ratio",drop=TRUE] else NA)
                seTnBG.ratio_vec <- c(seTnBG.ratio_vec, if(!is.null(criteria)) criteria[i,"seTnBG.ratio",drop=TRUE] else NA)
                rseTn_vec <- c(rseTn_vec, if(!is.null(criteria)) criteria[i,"rseTn",drop=TRUE] else NA)
                FR_vec <- c(FR_vec, if(!is.null(criteria)) criteria[i,"FR",drop=TRUE] else NA)
                seFR_vec <- c(seFR_vec, if(!is.null(criteria)) criteria[i,"seFR",drop=TRUE] else NA)
                ###
                res_calRcyRcp <- calRcyRcp(Curvedata=DataR, Ltx=DataN)
                ###
                RecyclingRatio1_vec   <- c(RecyclingRatio1_vec,   res_calRcyRcp$RecyclingRatio1[1L])
                seRecyclingRatio1_vec <- c(seRecyclingRatio1_vec, res_calRcyRcp$RecyclingRatio1[2L])
                ###
                RecyclingRatio2_vec   <- c(RecyclingRatio2_vec,   res_calRcyRcp$RecyclingRatio2[1L])
                seRecyclingRatio2_vec <- c(seRecyclingRatio2_vec, res_calRcyRcp$RecyclingRatio2[2L])
                ###
                RecyclingRatio3_vec   <- c(RecyclingRatio3_vec,   res_calRcyRcp$RecyclingRatio3[1L])
                seRecyclingRatio3_vec <- c(seRecyclingRatio3_vec, res_calRcyRcp$RecyclingRatio3[2L])
                ###
                Recuperation1_vec   <- c(Recuperation1_vec,   res_calRcyRcp$Recuperation1[1L])
                seRecuperation1_vec <- c(seRecuperation1_vec, res_calRcyRcp$Recuperation1[2L])
                ###
                Recuperation2_vec   <- c(Recuperation2_vec,   res_calRcyRcp$Recuperation2[1L])
                seRecuperation2_vec <- c(seRecuperation2_vec, res_calRcyRcp$Recuperation2[2L])
                ###
                FOM_vec <- c(FOM_vec, res$FOM)
                RCS_vec <- c(RCS_vec, res$RCS)
                ###
                ### Growth curve fitting succeeded.
                characterNO <- paste("NO", NO[i], sep="")
                LMpars[[characterNO]] <- res$LMpars
                ###
                ###
                dlsyc <- DataList[[i]][,"SAR.Cycle",drop=TRUE]
                dddxxx1 <- DataList[[i]][dlsyc=="N",,drop=FALSE]
                dddxxx2 <- DataList[[i]][dlsyc!="N",,drop=FALSE]
                ###
                ith_SARdata1 <- rbind(dddxxx1, dddxxx2[res$fitIDX,,drop=FALSE])
                ###
                SARdata1 <- rbind(SARdata1, ith_SARdata1)
                ###
                if (!is.null(norm.dose)) {
                    pars <-  res$LMpars[,1L,drop=TRUE]
                    cst <- ifelse(origin1==TRUE, 0, pars[length(pars)])
                    if (model=="line") {
                        norm.signal <- pars[1L]*norm.dose+cst
                    } else if(model=="exp") {
                        norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+cst
                    } else if(model=="lexp")  {
                        norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+
                                       pars[3L]*norm.dose+cst
                    } else if(model=="dexp") {
                        norm.signal <- pars[1L]*(1.0-exp(-pars[2L]*norm.dose))+
                                       pars[3L]*(1.0-exp(-pars[4L]*norm.dose))+cst
                    } else if(model=="gok") {
                        norm.signal <- pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*norm.dose)^(-1.0/pars[3L]))+cst
                    } # end if.
                    ###
                    ith_norm.SARdata <- ith_SARdata1
                    ###
                    ith_norm.SARdata[,c("Signal","Signal.Err")] <- 
                        ith_norm.SARdata[,c("Signal","Signal.Err"),drop=FALSE]/norm.signal
                    ###
                    norm.SARdata <- rbind(norm.SARdata, ith_norm.SARdata)
                } # end if.
                ###
            } else if (res$message==1L) {
                failFit_ID <- rbind(failFit_ID, agID[i,,drop=TRUE])
            } # end if. 
            ###
        } # end if.
    } # end for.
    ###------------------------------------------------------------------------
    if (n>=5L) close(pb)
    ###
    if (!is.null(outpdf)) dev.off()
    ###

    ###
    ###-------------------------------------------------------------------------------
    ### Fliter on improper arguments.
    action_character <- c(action_character, 
        "Function fitGrowth(): improper input argument")
    ###
    if (!is.null(tryError_ID)) {
        reject_N <- nrow(tryError_ID)
        step_reject_N <- c(step_reject_N, reject_N)
        tryError_reject <- apply(tryError_ID, MARGIN=1L, NPG)
    } else {
        step_reject_N <- c(step_reject_N, 0L)
        tryError_reject <- NULL
    } # end if.
    ###-------------------------------------------------------------------------------
    ###

    ###
    ###-------------------------------------------------------------------------------
    ### Fliter on growth curve fitting fails.
    action_character <- c(action_character, 
        "Function fitGrowth(): failed in growth curve fitting")
    ###
    if (!is.null(failFit_ID)) {
        reject_N <- nrow(failFit_ID)
        step_reject_N <- c(step_reject_N, reject_N)
        failFit_reject <- apply(failFit_ID, MARGIN=1L, NPG)
    } else {
        step_reject_N <- c(step_reject_N, 0L)
        failFit_reject <- NULL
    } # end if.
    ###-------------------------------------------------------------------------------
    ###

    ###
    ###**************************************************************************************************************************
    ###--------------------------------------------------------------------------------------------------------------------------
    if (!is.null(passNO)) {
        ###
        SARdata.table <- data.frame("NO"=passNO, "Position"=passPosition, "Grain"=passGrain,
                         "Tn"=Tn_vec, "seTn"=seTn_vec, "Tn3BG"=Tn3BG_vec, "TnBG.ratio"=TnBG.ratio_vec, 
                         "seTnBG.ratio"=seTnBG.ratio_vec, "rseTn"=rseTn_vec, "FR"=FR_vec, "seFR"=seFR_vec,
                         "RecyclingRatio1"=RecyclingRatio1_vec, "seRecyclingRatio1"=seRecyclingRatio1_vec, 
                         "RecyclingRatio2"=RecyclingRatio2_vec, "seRecyclingRatio2"=seRecyclingRatio2_vec,
                         "RecyclingRatio3"=RecyclingRatio3_vec, "seRecyclingRatio3"=seRecyclingRatio3_vec,
                         "Recuperation1"=Recuperation1_vec, "seRecuperation1"=seRecuperation1_vec,
                         "Recuperation2"=Recuperation2_vec, "seRecuperation2"=seRecuperation2_vec,
                         "FOM"=FOM_vec, "RCS"=RCS_vec, stringsAsFactors=FALSE)
        ###
        agID <- cbind("NO"=passNO, "Position"=passPosition, "Grain"=passGrain)
        ###

        ###
        ### Apply growth curve related rejection criteria (2).
        ###*****************************************************************************************
        ###
        ###---------------------------------------------------------------------------------
        ### Figure-of-merit: no se consideration.
        if (!is.null(fom.up)) {
            all_value <- SARdata.table[,"FOM",drop=TRUE]
            ###
            select_index <- which(abs(all_value)<fom.up)
            ###
            reject_N <- nrow(SARdata.table) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: FOM of growth curve exceeds ", fom.up, "%",sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            if (length(select_index)==0L) {
                ###
                cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
                ###
                action_character <- c(action_character, 
                    "Total number of rejected aliquots (grains)",
                    "Total number of accepted aliquots (grains)")
                ###
                step_reject_N <- c(step_reject_N, nag, 0L)
                ###
                summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
                ###
                print(summary_info)
                ###
                return(invisible(summary_info))
            } # end if.
            ###
            SARdata.table <- SARdata.table[select_index,,drop=FALSE]
            fom_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
            LMpars <- LMpars[select_index]
        } else {
            fom_reject <- NULL
        } # end if.
        ###---------------------------------------------------------------------------------
        ###

        ###
        ###---------------------------------------------------------------------------------
        ### Reduced chi-square: no se consideration.
        if (!is.null(rcs.up)) {
            all_value <- SARdata.table[,"RCS",drop=TRUE]
            ###
            select_index <- which(abs(all_value)<rcs.up | is.na(all_value))
            ###
            reject_N <- nrow(SARdata.table) - length(select_index)
            ###
            action_character <- c(action_character, 
                paste("Rejection criterion: RCS of growth curve exceeds ", rcs.up, sep=""))
            step_reject_N <- c(step_reject_N, reject_N)
            ###
            if (length(select_index)==0L) {
                cat("NOTE: no acceptable SAR data if the specified rejection criteria are applied!\n")
                ###
                action_character <- c(action_character, 
                    "Total number of rejected aliquots (grains)",
                    "Total number of accepted aliquots (grains)")
                ###
                step_reject_N <- c(step_reject_N, nag, 0L)
                ###
                summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
                ###
                print(summary_info)
                ###
                return(invisible(summary_info))
            } # end if.
            ###
            SARdata.table <- SARdata.table[select_index,,drop=FALSE]
            rcs_reject <- apply(agID[-select_index,,drop=FALSE], MARGIN=1L, NPG)
            agID <- agID[select_index,,drop=FALSE]
            LMpars <- LMpars[select_index]
        } else {
            rcs_reject <- NULL
        } # end if.
        ###------------------------------------------------------------------------------------
        ###
        ###*****************************************************************************************
        ###

        ###
        if (!is.null(outfile)) write.csv(SARdata.table, file=paste(outfile,".csv",sep=""))
        ###
        select_NO <- agID[,"NO",drop=TRUE]
        ###
        SARdata1 <- SARdata1[SARdata1[,"NO",drop=TRUE] %in% select_NO,,drop=FALSE]
        ###
        if (!is.null(norm.SARdata)) {
            norm.SARdata <- norm.SARdata[norm.SARdata[,"NO",drop=TRUE] %in% select_NO,,drop=FALSE]
        } ### end if. 
        ###
        output <- list("LMpars"=LMpars,
                       "SARdata"=SARdata1,
                       "norm.SARdata"=norm.SARdata,
                       "agID"=agID)

        ###
        ###-----------------------------------------------------------------------------
        if (is_forced_object==FALSE) {
            ###
            if (length(Tn3BG_reject)>0L) {
                cat("\n")
                cat("Rejection criterion: aliquot (grain) ID rejected use [Tn.above.3BG]:\n")
                print(Tn3BG_reject)
                cat("\n")
            } # end if.
            ###

            ###
            if (length(TnBG.ratio_reject)>0L) {
                cat("\n")
                cat("Rejection criterion: aliquot (grain) ID rejected use [TnBG.ratio]:\n")
                print(TnBG.ratio_reject)
                cat("\n")
            } # end if.
            ###

            ###
            if (length(rseTn_reject)>0L) {
                cat("\n")
                cat("Rejection criterion: aliquot (grain) ID rejected use [rseTn]:\n")
                print(rseTn_reject)
                cat("\n")
            } # end if.
            ### 

            ###
            if (length(FR_reject)>0L) {
                cat("\n")
                cat("Rejection criterion: aliquot (grain) ID rejected use [FR]:\n")
                print(FR_reject)
                cat("\n")
            } # end if.
            ###
        } # end if.
        ###-----------------------------------------------------------------------------

        ###-----------------------------------------------------------------------------
        ###
        if (length(rcy1_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rcy1]:\n")
            print(rcy1_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (length(rcy2_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rcy2]:\n")
            print(rcy2_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (length(rcy3_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rcy3]:\n")
            print(rcy3_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (length(rcp1_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rcp1]:\n")
            print(rcp1_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (length(rcp2_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rcp2]:\n")
            print(rcp2_reject)
            cat("\n")
        } # end if.
        ###
        ###-----------------------------------------------------------------------------
        ###

        ###
        ###-----------------------------------------------------------------------------
        ###
        if (!is.null(tryError_ID)) {
            cat("\n")
            cat("Function fitGrowth(): aliquot (grain) ID with improper input argument:\n")
            print(tryError_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (!is.null(failFit_ID)) {
            cat("\n")
            cat("Function fitGrowth(): aliquot (grain) ID failed in growth curve fitting:\n")
            print(failFit_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (length(fom_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [fom]:\n")
            print(fom_reject)
            cat("\n")
        } # end if.
        ###

        ###
        if (length(rcs_reject)>0L) {
            cat("\n")
            cat("Rejection criterion: aliquot (grain) ID rejected use [rcs]:\n")
            print(rcs_reject)
            cat("\n")
        } # end if.
        ###

    } # end if.
    ###**************************************************************************************************************************
    ###--------------------------------------------------------------------------------------------------------------------------
   
    ###
    action_character <- c(action_character, 
                          "Total number of rejected aliquots (grains)",
                          "Total number of accepted aliquots (grains)")
    ###
    step_reject_N <- c(step_reject_N, 
                       ifelse(is.null(passNO), nag, nag-nrow(agID)),
                       ifelse(is.null(passNO), 0L, nrow(agID)))
    ###
    summary_info <- data.frame("Description"=action_character, "N"=step_reject_N)
    print(summary_info)
    ###
  
    ###
    if (!is.null(passNO)) {
        ###
        output$summary.info <- summary_info
        ###
        return(invisible(output))
    } else {
        return(invisible(summary_info))
    } # end if.
    ###
} # end function pickSARdata.default.
#####
