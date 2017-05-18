#####
as_analyseBIN <- 
function(SARdata) {
    UseMethod("as_analyseBIN")
} #
### 2017.05.16. 
as_analyseBIN.default <- 
function(SARdata) {
    stopifnot(ncol(SARdata)==5L, nrow(SARdata)>=1L,
              is.numeric(SARdata[,1L,drop=TRUE]), 
              is.numeric(SARdata[,3L,drop=TRUE]), all(SARdata[,3L,drop=TRUE]>=0),
              is.numeric(SARdata[,4L,drop=TRUE]), 
              is.numeric(SARdata[,5L,drop=TRUE]), all(SARdata[,5L,drop=TRUE]>0))
    ###
    NO <- sort(as.numeric(levels(factor(SARdata[,1L,drop=TRUE]))))
    n <- length(NO)
    ###
    order_SARdata <- c()
    for (i in seq(n)) {
        iIndex <- which(SARdata[,1L,drop=TRUE]==NO[i])  
        ith_SARdata <- SARdata[iIndex,,drop=FALSE]
        order_SARdata <- rbind(order_SARdata, ith_SARdata)
    } # end for.
    colnames(order_SARdata) <- c("NO","SAR.Cycle","Dose","Signal","Signal.Err")
    SARdata <- order_SARdata
    ###
    criteria <- NULL
    ###
    Tn <- NULL
    ###
    TxTn <- NULL
    ###
    LnTn.curve <- NULL
    ###
    agID <- cbind(NO, 0L, 0L)
    colnames(agID) <- c("NO", "Position", "Grain")
    ###
    output <- list("SARdata"=SARdata,
                   "criteria"=criteria,
                   "Tn"=Tn,
                   "LnTn.curve"=LnTn.curve,
                   "TxTn"=TxTn,
                   "agID"=agID)
    class(output) <- "analyseBIN"
    ###
    invisible(output)
} # end function as_analyseBIN.default.
#####
