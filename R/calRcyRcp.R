#####
calRcyRcp <-
function(Curvedata, Ltx) {
    UseMethod("calRcyRcp")
} #
### 2017.01.22.
calRcyRcp.default <- 
function(Curvedata, Ltx) {
    stopifnot(ncol(Curvedata)==3L, nrow(Curvedata)>=1L,
              length(Ltx)==2L)
    ###
    dose <- as.numeric(Curvedata[,1L,drop=TRUE])
    doseltx <- as.numeric(Curvedata[,2L,drop=TRUE])
    sdoseltx <- as.numeric(Curvedata[,3L,drop=TRUE])
    ###
    ### Calculate recycling ratio.
    dose_rmZERO <- dose[abs(dose)>.Machine$double.eps^0.5]
    dose_level <- sort(as.numeric(levels(factor(dose_rmZERO))))
    exist_repeatDose <- length(dose_rmZERO)>length(dose_level)
    if (exist_repeatDose==TRUE) {
        RepeatIndex <- apply(as.matrix(dose_level), MARGIN=1L, function(x,y) 
                             which(abs(x-y)<=.Machine$double.eps^0.5), dose)
        ###
        if (class(RepeatIndex)=="list") {
            RepeatIndex <- RepeatIndex[sapply(RepeatIndex,length)>=2L][[1L]]
        } else {
            RepeatIndex <- RepeatIndex[,1L,drop=TRUE]
        } # end if.
        ###
        ###
        RecyclingRatio1 <- doseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[1L]]
        rseV1 <- sdoseltx[RepeatIndex[1L]]/doseltx[RepeatIndex[1L]]
        rseV2 <- sdoseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[2L]]
        seRecyclingRatio1 <- abs(RecyclingRatio1)*sqrt(rseV1^2L+rseV2^2L)
        ###
        if (length(RepeatIndex)>=3L) {
            RecyclingRatio2 <- doseltx[RepeatIndex[3L]]/doseltx[RepeatIndex[1L]]
            rseV1 <- sdoseltx[RepeatIndex[1L]]/doseltx[RepeatIndex[1L]]
            rseV2 <- sdoseltx[RepeatIndex[3L]]/doseltx[RepeatIndex[3L]]
            seRecyclingRatio2 <- abs(RecyclingRatio2)*sqrt(rseV1^2L+rseV2^2L)
            ###
            RecyclingRatio3 <- doseltx[RepeatIndex[3L]]/doseltx[RepeatIndex[2L]]
            rseV1 <- sdoseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[2L]]
            rseV2 <- sdoseltx[RepeatIndex[3L]]/doseltx[RepeatIndex[3L]]
            seRecyclingRatio3 <- abs(RecyclingRatio3)*sqrt(rseV1^2L+rseV2^2L)
        } else {
            RecyclingRatio2 <- seRecyclingRatio2 <- NA
            RecyclingRatio3 <- seRecyclingRatio3 <- NA
        } # end if.      
    } else {
        RecyclingRatio1 <- seRecyclingRatio1 <- NA
        RecyclingRatio2 <- seRecyclingRatio2 <- NA
        RecyclingRatio3 <- seRecyclingRatio3 <- NA
    } # end if
    ###
    ### Calculate recuperation.
    zeroDose_index <- which(abs(dose)<=.Machine$double.eps^0.5)
    if (length(zeroDose_index)>=1L) {
        Recuperation1 <- (doseltx[zeroDose_index[1L]]/Ltx[1L])*100.0
        rseV1 <- sdoseltx[zeroDose_index[1L]]/doseltx[zeroDose_index[1L]]
        rseV2 <- Ltx[2L]/Ltx[1L]
        seRecuperation1 <- abs(Recuperation1)*sqrt(rseV1^2L+rseV2^2L)
        ###
        Recuperation2 <- (doseltx[zeroDose_index[1L]]/max(doseltx))*100.0
        rseV1 <- sdoseltx[zeroDose_index[1L]]/doseltx[zeroDose_index[1L]]
        rseV2 <- sdoseltx[which.max(doseltx)]/max(doseltx)
        seRecuperation2 <- abs(Recuperation2)*sqrt(rseV1^2L+rseV2^2L)
    } else {
        Recuperation1 <- seRecuperation1 <- NA
        Recuperation2 <- seRecuperation2 <- NA     
    } # end if.  
    ###
    RecyclingRatio1_seRecyclingRatio1 <- c("RecyclingRatio1"=RecyclingRatio1, "seRecyclingRatio1"=seRecyclingRatio1)
    RecyclingRatio2_seRecyclingRatio2 <- c("RecyclingRatio2"=RecyclingRatio2, "seRecyclingRatio2"=seRecyclingRatio2)
    RecyclingRatio3_seRecyclingRatio3 <- c("RecyclingRatio3"=RecyclingRatio3, "seRecyclingRatio3"=seRecyclingRatio3)
    Recuperation1_seRecuperation1 <- c("Recuperation1"=Recuperation1, "seRecuperation1"=seRecuperation1)
    Recuperation2_seRecuperation2 <- c("Recuperation2"=Recuperation2, "seRecuperation2"=seRecuperation2)
    ###
    output<-list("RecyclingRatio1"=RecyclingRatio1_seRecyclingRatio1,
                 "RecyclingRatio2"=RecyclingRatio2_seRecyclingRatio2,
                 "RecyclingRatio3"=RecyclingRatio3_seRecyclingRatio3,
                 "Recuperation1"=Recuperation1_seRecuperation1,
                 "Recuperation2"=Recuperation2_seRecuperation2)
    ###
    return(output)
} # end function calRcyRcp.default.
#####
