#####
pickSARdata<-
function(Data, model="gok", origin=FALSE, weight=TRUE, 
         rcy.interval=NULL, rcp1.limit=NULL, rcp2.limit=NULL, 
         fom.limit=NULL, rcs.limit=NULL, outfile=NULL) {
    UseMethod("pickSARdata")
} #
### 2016.07.14.
pickSARdata.default<-
function(Data, model="gok", origin=FALSE, weight=TRUE, 
         rcy.interval=NULL, rcp1.limit=NULL, rcp2.limit=NULL, 
         fom.limit=NULL, rcs.limit=NULL, outfile=NULL) {
    ### Stop if not.
    stopifnot(ncol(Data)==5L, nrow(Data)>=5L,
              is.numeric(Data[,1L,drop=TRUE]), all(abs(Data[,1L]-round(Data[,1L]))<.Machine$double.eps^0.5),
              is.numeric(Data[,3L,drop=TRUE]), is.numeric(Data[,4L,drop=TRUE]), is.numeric(Data[,5L,drop=TRUE]),
              all(Data[,3L,drop=TRUE]>=0), all(Data[,5L,drop=TRUE]>0),
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(weight)==1L, is.logical(weight),
              is.null(rcy.interval) || (is.numeric(rcy.interval) && length(rcy.interval)==2L),
              is.null(rcp1.limit) || (is.numeric(rcp1.limit) && length(rcp1.limit)==1L),
              is.null(rcp2.limit) || (is.numeric(rcp2.limit) && length(rcp2.limit)==1L),
              is.null(fom.limit) || (is.numeric(fom.limit) && length(fom.limit)==1L),
              is.null(rcs.limit) || (is.numeric(rcs.limit) && length(rcs.limit)==1L),
              is.null(outfile) || (is.character(outfile) && length(outfile)==1L))
    ###
    colnames(Data) <- c("Grain.NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    GrainNumber <- as.numeric(levels(factor(Data[,"Grain.NO",drop=TRUE])))
    n <- length(GrainNumber)
    originGrainNumber <- GrainNumber
    ###
    ### Check Grain.NO and SAR.Cycle for Data.
    for (i in seq(n)) {
        GrainIndex <- which(Data[,"Grain.NO",drop=TRUE]==GrainNumber[i])
        SarCyclei <- substr(Data[GrainIndex,"SAR.Cycle",drop=TRUE], start=1L, stop=1L)   
        ###
        if (!all(diff(GrainIndex)==1L)) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data appears in discontinuous locations!", sep=""))
        } ### end if. 
        ###  
        if (!all(SarCyclei %in% c("N","R"))) {
             stop(paste("Error: Grain.NO", GrainNumber[i], 
                  " of Data contains incorrect SAR.Cycle!", sep=""))
        } # end if.
        ###
        if (all(SarCyclei=="N")) {
           stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data should contain SAR.Cycle of R!", sep=""))
        } # end if. 
        ###
        if (all(SarCyclei=="R")) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data should contain SAR.Cycle of N!", sep=""))
        } # end if.
        ###
        if (sum(SarCyclei=="N")>1L) {
            stop(paste("Error: Grain.NO", GrainNumber[i], 
                " of Data should contain only one SAR.Cycle of N!", sep=""))
        } # end if.
    } # end for.
    ###
    ###
    DataList <- vector(mode="list", length=n)
    for (i in seq(n)) {
        GrainIndex <- which(Data[,"Grain.NO",drop=TRUE]==GrainNumber[i])  
        DataList[[i]] <- Data[GrainIndex,,drop=FALSE]
    } # end for.
    ###
    ###
    tab <- data.frame(matrix(nrow=n, ncol=8L))
    rownames(tab) <- paste("Grain.NO",GrainNumber,sep="")
    colnames(tab) <- c("RecyclingRatio", "Std.RecyclingRatio",
                       "Recuperation1", "Std.Recuperation1",
                       "Recuperation2", "Std.Recuperation2",
                       "FOM", "RCS")
    ###
    failFitGrainIndex <- c()
    failFitGrainNumber <- NULL
    for (i in seq(n)) {
        Data4L <- DataList[[i]][,c("SAR.Cycle","Dose","Signal","Signal.Err"),drop=FALSE]
        DataN <- as.numeric(Data4L[Data4L[,"SAR.Cycle",drop=TRUE]=="N",-1L])
        DataR <- Data4L[Data4L[,"SAR.Cycle",drop=TRUE]!="N",-1L,drop=FALSE]
        ###
        res <- try(fitGrowth(DataR, model=model, origin=origin,
                   weight=weight, plot=FALSE), silent=TRUE)
        ###
        if (class(res)=="try-error") {
            failFitGrainIndex <- c(failFitGrainIndex, i)
            tab[i,] <- NA
        } else {
            ###
            ### Calculate recycling ratio.
            dose <- DataR[,"Dose",drop=TRUE]
            doseltx <- DataR[,"Signal",drop=TRUE]
            sdoseltx <- DataR[,"Signal.Err",drop=TRUE]
            ###
            lvl.dose <- as.numeric(levels(factor(dose)))
            existrpd <- length(dose)>length(lvl.dose)
            if (existrpd==TRUE) {
                RepeatIndex <- apply(as.matrix(lvl.dose), MARGIN=1L, function(x,y)
                                     which(abs(x-y)<=.Machine$double.eps^0.5), dose)
                RepeatIndex <- unlist(RepeatIndex[sapply(RepeatIndex,length)>=2L])
                RecycleRatio <- doseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[1L]]
                sRecycleRatio <- abs(RecycleRatio)*sqrt((sdoseltx[RepeatIndex[2L]]/doseltx[RepeatIndex[2L]])^2L+
                    (sdoseltx[RepeatIndex[1L]]/doseltx[RepeatIndex[1L]])^2L)
                tab[i,c(1L,2L)] <- c(RecycleRatio, sRecycleRatio)
            } else {
                tab[i,c(1L,2L)] <- NA
            } # end if.
            ###
            ###
            ### Calculate recuperation.
            exist0d <- which(abs(dose)<=.Machine$double.eps^0.5)
            if (length(exist0d)>=1L) {
                Recuperation1 <- (doseltx[exist0d[1L]]/DataN[2L])*100.0
                sRecuperation1 <- abs(Recuperation1)*sqrt((sdoseltx[exist0d[1L]]/doseltx[exist0d[1L]])^2L+
                    (DataN[3L]/DataN[2L])^2L)
                Recuperation2 <- (doseltx[exist0d[1L]]/max(doseltx))*100.0
                sRecuperation2 <- abs(Recuperation2)*sqrt((sdoseltx[exist0d[1L]]/doseltx[exist0d[1L]])^2L+
                    (sdoseltx[which.max(doseltx)]/max(doseltx))^2L)
                tab[i,c(3L,4L)] <- c(Recuperation1, sRecuperation1)
                tab[i,c(5L,6L)] <- c(Recuperation2, sRecuperation2)
            } else {
                tab[i,c(3L,4L)] <- NA
                tab[i,c(5L,6L)] <- NA
            } # end if.  
            ###
            ###
            tab[i,7L] <- 100*sum(abs(res$fit.value[,2L,drop=TRUE]-
                res$fit.value[,3L,drop=TRUE]))/sum(res$fit.value[,3L,drop=TRUE])
            tab[i,8L] <- res$value/(nrow(DataR)-nrow(res$LMpars))
        } # end if.
    } # end for.
    ###
    ###
    nfailFitGrainIndex <- length(failFitGrainIndex)
    if (nfailFitGrainIndex>=1L) {
        if (nfailFitGrainIndex==n) stop("Error: fail in fitting al growth curves!")
        tab <- tab[-failFitGrainIndex,,drop=FALSE]
        failFitGrainNumber <- GrainNumber[failFitGrainIndex] 
        GrainNumber <- GrainNumber[-failFitGrainIndex]
    } # end if.
    ###
    ###
    ###
    if (!is.null(rcy.interval)) {
        indexValue <- tab[,"RecyclingRatio",drop=TRUE]
        selected <- which(indexValue>rcy.interval[1L] & indexValue<rcy.interval[2L])
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
        GrainNumber <- GrainNumber[selected]   
    } # end if.
    ###
    if (!is.null(rcp1.limit)) {
        indexValue <- tab[,"Recuperation1",drop=TRUE]
        selected <- which(indexValue<rcp1.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
        GrainNumber <- GrainNumber[selected]
    } # end if.
    ###
    if (!is.null(rcp2.limit)) {
        indexValue <- tab[,"Recuperation2",drop=TRUE]
        selected <- which(indexValue<rcp2.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
        GrainNumber <- GrainNumber[selected]
    } # end if.
    ###
    if (!is.null(fom.limit)) {
        indexValue <- tab[,"FOM",drop=TRUE]
        selected <- which(indexValue<fom.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
        GrainNumber <- GrainNumber[selected]
    } # end if.
    ###
    if (!is.null(rcs.limit)) {
        indexValue <- tab[,"RCS",drop=TRUE]
        selected <- which(indexValue<rcs.limit)
        tab <- tab[selected,,drop=FALSE]
        if (nrow(tab)<=0L) stop("Error: no Grain.NO satisfies the given conditions!")
        GrainNumber <- GrainNumber[selected]
    } # end if.
    ###
    Data <- Data[Data[,"Grain.NO",drop=TRUE] %in% GrainNumber,,drop=FALSE]
    ###
    if (!is.null(outfile)) {
        write.csv(tab, file=paste(outfile,".csv",sep=""))
    } # end if.  
    ###
    rejectGrainNumber <- originGrainNumber[which(!originGrainNumber %in% 
        as.numeric(substr(rownames(tab), start=9L, stop=10000L)))]
    if (length(rejectGrainNumber)==0L) rejectGrainNumber <- NULL
    ###
    output <- list("Data"=Data,
                   "tab"=tab,
                   "failFit.NO"=failFitGrainNumber,
                   "reject.NO"=rejectGrainNumber)
    ###
    return(output)
} # end function pickSARdata.default.
#####
