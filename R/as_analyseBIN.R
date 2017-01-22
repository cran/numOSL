#####
as_analyseBIN <- 
function(SARdata) {
    UseMethod("as_analyseBIN")
} #
### 2017.01.22. 
as_analyseBIN.default <- 
function(SARdata) {
    stopifnot(ncol(SARdata)==5L,
              nrow(SARdata)>=1L)
    ###
    colnames(SARdata) <- c("NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    NO <- sort(as.numeric(levels(factor(SARdata[,"NO",drop=TRUE]))))
    n <- length(NO)
    ###
    criteria <- cbind(rep(1L,n), rep(99999.0,n), rep(0.0,n), rep(0.0,n), rep(99999.0,n), rep(0,n))
    rownames(criteria) <- NULL
    colnames(criteria) <- c("Tn3BG", "TnBG.ratio", "seTnBG.ratio", "rseTn", "FR", "seFR")
    ###
    Tn <- cbind(rep(99999.0,n), rep(0.0,n))
    rownames(Tn) <- NULL
    colnames(Tn) <- c("Tn","seTn")
    ###
    TxTn <- NULL
    ###
    agID <- matrix(nrow=n, ncol=3L)
    colnames(agID) <- c("NO", "Position", "Grain")
    ###
    for (i in seq(n)) {
        iIndex <- which(SARdata[,"NO",drop=TRUE]==NO[i])
        agID[i,] <-  c(SARdata[iIndex,"NO",drop=TRUE][1L],0L,0L)
    } # end for.
    ###
    output <- list("SARdata"=SARdata,
                   "criteria"=criteria,
                   "Tn"=Tn,
                   "TxTn"=TxTn,
                   "agID"=agID)
    class(output) <- "analyseBIN"
    ###
    invisible(output)
} # end function as_analyseBIN.default.
#####
                   
