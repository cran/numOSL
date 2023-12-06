#####
dbED <-
function(EDdata, plot=TRUE, typ="pdf",
         from=NULL, to=NULL, step=NULL, nbin=15,
         pcolor="purple", psize=1.5, outfile=NULL) {
    UseMethod("dbED")
} #
### 2023.12.03.
dbED.default <-
function(EDdata, plot=TRUE, typ="pdf",
         from=NULL, to=NULL, step=NULL, nbin=15,
         pcolor="purple", psize=1.5, outfile=NULL) {
    ### Stop if not.
    stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=2L,
              all(EDdata[,2L,drop=TRUE]>=0),
              length(plot)==1L, is.logical(plot),
              length(typ)==1L, is.character(typ), typ %in% c("pdf","hist"),
              is.null(from) || is.numeric(from),
              is.null(to) || is.numeric(to),
              is.null(step) || is.numeric(step), 
              length(nbin)==1L, is.numeric(nbin),
              length(pcolor)==1L, is.character(pcolor),
              length(psize)==1L, is.numeric(psize),
              is.null(outfile) || is.character(outfile))
    ###
    if (!is.null(step)) {
        if (length(step)!=1L) stop("Error: step should be an one-element vector!")
        if (step<=0)  stop("Error: step should exceed zero!")
    } # end if.
    ### 
    if (!is.null(outfile)) {
        if (length(outfile)!=1L)  stop("Error: outfile should be an one-element vector!")
    } # end if.
    ###
    nED<-nrow(EDdata)
    ed1<-as.numeric(EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(EDdata[,2L,drop=TRUE])
    idx <- which(sed1<=0)
    if (length(idx)>=1L) {
        cat("NOTE: the ", idx, "-th De, zero error!\n") 
        sed1[idx] <- ed1[idx]*0.01
    } # end if.
    ###
    idx <-order(ed1)
    ed1 <- ed1[idx]
    sed1 <- sed1[idx]
    ###
    if (is.null(from)) {
        from<-ifelse(min(ed1)>=0,
                     min(ed1)*0.7,
                     min(ed1)*1.3)
    } # end if.
    if (is.null(to)) {
        to<-ifelse(max(ed1)>=0,
                   max(ed1)*1.3,
                   max(ed1)*0.7)
    } # end if.
    if (from>=to)  stop("Error: from must not exceed to!")
    ###
    ### Calculate weighted skewness and kurtosis.
    weight<-ed1/sed1
    meanED<-mean(ed1)
    sdED<-sd(ed1)
    medianED<-median(ed1)
    skewness<-sum(weight*((ed1-meanED)/sdED)^3L)/sum(weight)     
    Std.skewness<-sqrt(6L/nED)   
    kurtosis<-nED*(nED+1L)/(nED-1L)/(nED-2L)/(nED-3L)*
              sum(((ed1-meanED)/sdED)^4L)-
              3L*(nED-1L)^2L/(nED-2L)/(nED-3L)
    Std.kurtosis<-sqrt(24L/nED)
    ###
    ### Calculate weighted ED.
    weightED<-sum(ed1/sed1^2L)/sum(1.0/sed1^2L)
    sweightED <-1.0/sqrt(sum(1.0/sed1^2L))
    weightED<-c(weightED,sweightED)
    ###
    quantileED<-quantile(ed1, probs=c(0.05,0.1,0.15,
                         0.5,0.85,0.9,0.95))
    ###
    output<-list("weight.ED"=round(weightED, 3L),
                 "skewness"=round(c(skewness, Std.skewness),3L),
                 "kurtosis"=round(c(kurtosis, Std.kurtosis),3L),
                 "quantile.ED"=round(quantileED,3L))
    ###
    if (plot==TRUE) {
        opar <- par("new")
        on.exit(par("new"=opar))
        ###
        if (typ=="pdf") {
            if (is.null(step)) { 
                step<-max(diff(sort(ed1)))/10.0
                cat(paste("Default step size:", round(step,3L), "\n\n"))
            } # end if.
            ###
            spreadED<-seq(from=from, to=to, by=step)
            pdfMat<-matrix(nrow=length(spreadED), ncol=nED)
            ###
            for(i in seq(nED)) {
                pdfMat[,i]<-dnorm(x=spreadED, mean=ed1[i], 
                                  sd=sed1[i], log=FALSE)
            } # end if.
            pdfED<-rowSums(pdfMat)/sum(rowSums(pdfMat))
            ###
            if (!is.null(outfile))  {
                write.csv(cbind(spreadED, pdfED), 
                  file=paste(outfile,".csv",sep=""))
            } # end if.        
            ### 
            plot(spreadED, pdfED, mgp=c(2.5,1,0),
                 main="De Distribution", col="skyblue",
                 xlab="De (Gy)", ylab="Density", cex.lab=1.25,
                 xlim=c(from,to), type="l", lwd=3)
            xTicks<-axTicks(side=1L)
            maxYx<-spreadED[which.max(pdfED)]
            box(lwd=1)
            matpoints(spreadED, pdfMat/sum(rowSums(pdfMat)), 
                      type="l", lwd=0.9, lty=1, col="grey60")
        } else if (typ=="hist") {
            breaks<-pretty(ed1, n=nbin)
            HIST<-hist(ed1, breaks=breaks, 
                       main="De Distribution", border="grey50",
                       xlab="De (Gy)", col="skyblue",
                       xlim=c(from,to), mgp=c(2.5,1,0))
            xTicks<-axTicks(side=1L)
            maxYx<-HIST$mids[which.max(HIST$counts)]
            box(lwd=1)
        } # end if.
        ###
        ###
        par("new"=TRUE)
        plot(ed1, seq(nED), xlab="", ylab="", 
             xlim=c(from,to), type="p", 
             pch=23, cex=psize, bg=pcolor, 
             col=pcolor, xaxt="n", yaxt="n")
        ###
        ###
        arrowIndex <- which(sed1>0.05 & sed1/ed1>0.01)
        if (length(arrowIndex)>=1L) {
            arrows(x0=ed1[arrowIndex]-sed1[arrowIndex]/2.0, y0=(seq(nED))[arrowIndex],
                   x1=ed1[arrowIndex]+sed1[arrowIndex]/2.0, y1=(seq(nED))[arrowIndex],
                   code=3, lwd=psize, angle=90, length=0.05, col=pcolor)
        } # end if.
        ###
        if (typ=="pdf") {
            legend(ifelse(maxYx>median(xTicks),"left","right"), 
                   legend=c("Individual De", "Individual PDF", "Sum PDF", 
                            paste("N=",nED,sep=""),
                            paste("Mean=",round(meanED,2L)," (Gy)",sep=""),
                            paste("Median=",round(medianED,2L)," (Gy)",sep=""),
                            paste("SD=",round(sdED,2L)," (Gy)",sep="")), 
                            lty=c(NA,1,1,rep(NA,4)),col=c(pcolor,"grey60","skyblue",rep(NA,4)),
                            pch=c(23,rep(NA,6)), pt.bg=c(pcolor,rep(NA,6)), lwd=c(NA,2,2,rep(NA,4)),
                            cex=1, pt.cex=1.25, bty="n")
        } else {
            legend(ifelse(maxYx>median(xTicks),"left","right"), 
                   legend=c("Individual De", paste("N=",nED,sep=""),
                            paste("Mean=",round(meanED,2L)," (Gy)",sep=""),
                            paste("Median=",round(medianED,2L)," (Gy)",sep=""),
                            paste("SD=",round(sdED,2L)," (Gy)",sep="")), 
                            col=c(pcolor,rep(NA,4)),pch=c(23,rep(NA,4)),
                            cex=1, pt.bg=c(pcolor,rep(NA,4)), pt.cex=1.25, bty="n")

        } # end if.
        ###
    } # end if.
    ###
    return(output)
} # end function dbED.default.
#####           
