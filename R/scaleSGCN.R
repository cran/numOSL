#####
scaleSGCN<-
function(SGCdata, pars, model, origin) {
    UseMethod("scaleSGCN")
} ###
### 2017.01.12.
scaleSGCN.default<-
function(SGCdata, pars, model, origin) {
    ### Stop if not.
    stopifnot(ncol(SGCdata)==5L, nrow(SGCdata)>=1L,
              is.numeric(SGCdata[,1L,drop=TRUE]),
              is.numeric(SGCdata[,3L,drop=TRUE]), all(SGCdata[,3L,drop=TRUE]>=0),
              is.numeric(SGCdata[,4L,drop=TRUE]), 
              is.numeric(SGCdata[,5L,drop=TRUE]), all(SGCdata[,5L,drop=TRUE]>0),
              is.numeric(pars), 
              length(model)==1L, model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin))
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
        if (length(iSAR.Cycle)!=2L) {
            stop(paste("[NO=", NO[i], 
                       "]: should contain two 'SAR.Cycle'!", sep=""))
        } # end if.
        if (!("N" %in% iSAR.Cycle && "R" %in% iSAR.Cycle)) {
            stop(paste("[NO=", NO[i], 
                       "]: should contain 'SAR.Cycle' of 'N' and 'R'!", sep=""))
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
        scalingFactor <- fcn(regenerativeDose)/regenerativeSignal
        ###
        ScaledNaturalSignal[i,] <- scalingFactor*c(naturalSignal,naturalSignalError)       
    } # end for.
    ###
    rownames(ScaledNaturalSignal) <- paste("[NO=",NO,"]",sep="")
    colnames(ScaledNaturalSignal) <- c("Ltx", "seLtx")
    ###
    return(ScaledNaturalSignal)
} # end function scaleSGCN.default.
#####
