#####
lsNORM<- 
function(Data, model="gok", origin=FALSE, 
         maxiter=10, weight=TRUE, plot=TRUE) {
    UseMethod("lsNORM")
} ###	
### 2016.07.09.
lsNORM.default<- 
function(Data, model="gok", origin=FALSE, 
         maxiter=10, weight=TRUE, plot=TRUE) {
    ### Stop if not.
    stopifnot(ncol(Data)==5L, nrow(Data)>=5L,
              is.numeric(Data[,1L,drop=TRUE]), all(abs(Data[,1L]-round(Data[,1L]))<.Machine$double.eps^0.5),
              is.numeric(Data[,3L,drop=TRUE]), is.numeric(Data[,4L,drop=TRUE]), is.numeric(Data[,5L,drop=TRUE]),
              all(Data[,3L,drop=TRUE]>=0), all(Data[,5L,drop=TRUE]>0),
              length(model)==1L, is.character(model), model %in% c("line","exp","lexp","dexp","gok"),
              length(origin)==1L, is.logical(origin),
              length(maxiter)==1L, is.numeric(maxiter), maxiter>0, maxiter<=1e3, abs(maxiter-round(maxiter))<.Machine$double.eps^0.5,
              length(weight)==1L, is.logical(weight),
              length(plot)==1L, is.logical(plot))
    ###
    colnames(Data) <- c("Grain.NO","SAR.Cycle","Dose","Signal","Signal.Err")
    ###
    GrainNumber <- as.numeric(levels(factor(Data[,"Grain.NO",drop=TRUE])))
    n <- length(GrainNumber)
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
        if (sum(SarCyclei=="N")>1L) {
            stop(paste("Error: Grain.NO", GrainNumber[i],
                 " of Data should contain not more than one SAR.Cycle of N!", sep=""))
        } # end if.
    } # end for.
    ###
    ###
    ### Filter on Data (remove data with SAR.Cycle=N).
    Data <- Data[Data[,"SAR.Cycle",drop=TRUE]!="N",,drop=FALSE]
    GrainNumber <- as.numeric(levels(factor(Data[,"Grain.NO",drop=TRUE])))
    n <- length(GrainNumber)
    ###
    DataList <- vector(mode="list", length=n)
    for (i in seq(n)) {
        GrainIndex <- which(Data[,"Grain.NO",drop=TRUE]==GrainNumber[i])  
        DataList[[i]] <- Data[GrainIndex,,drop=FALSE]
    } # end for.
    ###
    ###
    originData3L <- cbind(c(unlist(sapply(DataList, function(x) x[,"Dose",drop=TRUE]))),
                          c(unlist(sapply(DataList, function(x) x[,"Signal",drop=TRUE]))),
                          c(unlist(sapply(DataList, function(x) x[,"Signal.Err",drop=TRUE]))))
    ###
    ###
    iter <- 0
    savedataList <- DataList 
    tol <- 1.0e-06
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
    repeat {
        Data3L <- cbind(c(unlist(sapply(DataList, function(x) x[,"Dose",drop=TRUE]))),
                        c(unlist(sapply(DataList, function(x) x[,"Signal",drop=TRUE]))),
                        c(unlist(sapply(DataList, function(x) x[,"Signal.Err",drop=TRUE]))))
        rsd_old <- sd(Data3L[,2L,drop=TRUE])/
                   mean(Data3L[,2L,drop=TRUE])
        ###
        ###
        res1 <- fitGrowth(Data3L, model=model, origin=origin,
                          weight=weight, plot=FALSE)
        pars <- res1$LMpars[,1L,drop=TRUE]
        ###
        ax <- 1.0e-05
        bx <- 1.0e+05
        ###
        for (i in seq(n)) {
            xx <- DataList[[i]][,"Dose",drop=TRUE]
            yy <- DataList[[i]][,"Signal",drop=TRUE]
            nd <- length(xx)
            SF <- 0
            fmin <- 0
            ###    
            res2 <- .Fortran("calcSF",as.double(ax),as.double(bx),as.double(xx),as.double(yy),
                             as.double(pars),as.integer(nd),as.integer(n2),as.integer(mdl),
                             SF=as.double(SF),fmin=as.double(fmin),PACKAGE="numOSL")
            ###
            DataList[[i]][,c("Signal","Signal.Err")] <- 
            DataList[[i]][,c("Signal","Signal.Err"),drop=FALSE]*res2$SF
        } # end for. 
        ###
        ###
        Data3L <- cbind(c(unlist(sapply(DataList, function(x) x[,"Dose",drop=TRUE]))),
                        c(unlist(sapply(DataList, function(x) x[,"Signal",drop=TRUE]))),
                        c(unlist(sapply(DataList, function(x) x[,"Signal.Err",drop=TRUE]))))
        rsd_new <- sd(Data3L[,2L,drop=TRUE])/
                   mean(Data3L[,2L,drop=TRUE])
        ###
        iter <- iter + 1L
        ###
        ###
        if (iter==1L) {
            LMpars1 <- res1$LMpars
            value1 <- res1$value
        } # end if
        ###
        ###
        if (iter==maxiter) break 
        if (abs(rsd_new-rsd_old)<=tol) break
    } # end repeat.
    ###
    SFs <- vector(length=n)
    for (i in seq(n)) {
        SFs[i] <- DataList[[i]][1L,"Signal"]/
                  savedataList[[i]][1L,"Signal"]
    } # end if.
    ###
    ###
    Data[,c("Signal","Signal.Err")] <- Data3L[,-1L,drop=FALSE]
    ###
    ### Fit the growth curve for the last time.
    res3 <- fitGrowth(Data3L, model=model, origin=origin,
                      weight=weight, plot=FALSE)
    ###
    output <- list("optData"=Data,
                   "sf"=SFs,
                   "iter"=iter,
                   "LMpars1"=LMpars1,
                   "value1"=value1,
                   "LMpars2"=res3$LMpars,
                   "value2"=res3$value)
    ###
    ###
    Plot3 <- function(Curvedata,pars,model,origin,
                      ylim,xlab,ylab,xaxt,legend) {
        ###
        plot(Curvedata[,-3L,drop=FALSE], type="p", pch=21, bg="black", cex=1.5,
             ylim=ylim, xlab=xlab, ylab=ylab, las=0, xaxt=xaxt, xaxs="r", 
             yaxs="i", cex.lab=1.5, cex.axis=1.5) 
        ###
        legend("topleft", legend=legend, yjust=2, ncol=1, cex=1.5, bty="o")
        ###
        x <- NULL
        cst <- ifelse(origin==TRUE, 0, pars[length(pars)])
        ###
        if (model=="line") {
            curve(pars[1L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="exp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="lexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*x+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="dexp") {
            curve(pars[1L]*(1.0-exp(-pars[2L]*x))+pars[3L]*(1.0-exp(-pars[4L]*x))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } else if (model=="gok") {
            curve(pars[1L]*(1.0-(1.0+pars[2L]*pars[3L]*x)^(-1.0/pars[3L]))+cst, 
                  type="l", add=TRUE, lwd=2, col="skyblue")
        } # end if.
        ###
        ###
        dose <- Curvedata[,1L,drop=TRUE]
        doseltx <- Curvedata[,2L,drop=TRUE]
        sdoseltx <- Curvedata[,3L,drop=TRUE]
        ###
        arrowIndex <- which(sdoseltx>0.05 & sdoseltx/doseltx>0.01)
        if (length(arrowIndex)>=1L) {
            arrows(x0=dose[arrowIndex], y0=doseltx[arrowIndex]-sdoseltx[arrowIndex]/2.0, 
                   x1=dose[arrowIndex], y1=doseltx[arrowIndex]+sdoseltx[arrowIndex]/2.0,
                   code=3, lwd=1, angle=90, length=0.05, col="black")
        } # end if.
        ###
        grid()
        box(lwd=2)
    } # end function Plot3.
    ###
    ###
    if (plot==TRUE) {
        layout(cbind(c(rep(1,9), 2, rep(3,9)),
                     c(rep(1,9), 2, rep(3,9))))
        ###
        ylim <- c(min(originData3L[,2L],0)*1.1, max(originData3L[,2L])*1.1)
        ###
        ### The first plot.
        par(mar=c(0,6.1,1.1,6.1))
        Curvedata <- originData3L 
        pars <- LMpars1[,1L,drop=TRUE]
        xlab <- ""
        ylab <- "Standardised OSL"
        xaxt <- "n"
        legend <- "Before LS-normalisation"
        Plot3(Curvedata, pars, model, origin, ylim, xlab, ylab, xaxt, legend)
        ###
        ### The second plot.
        par(mar=c(0,6.1,0,6.1))
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
        ###
        ### Then thrid plot.
        par(mar=c(4.1,6.1,0,6.1))
        Curvedata <- Data3L 
        pars <- res3$LMpars[,1L,drop=TRUE]
        xlab <- "Dose (Gy)"
        ylab <- "Normalised standardised OSL"
        xaxt <- "s"
        legend <- "After LS-normalisation"
        Plot3(Curvedata, pars, model, origin, ylim, xlab, ylab, xaxt, legend)
        ###
        ###
        par(mar=c(5,4,4,2)+0.1)
        layout(1L)
    } # end if. 
    ### 
    return(output)
} # end function lsNORM.default.
#####
