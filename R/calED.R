#####
calED<-
function(Curvedata, Ltx, model=c("line","exp","lexp","dexp"),
         origin=FALSE, nstart=100, upb=0.5, ErrorMethod=c("mc","sp"),
         nsim=1000, plot=TRUE, outfile=NULL) {
    UseMethod("calED")
} ###
### 2014.09.12.
calED.default<-
function(Curvedata, Ltx, model=c("line","exp","lexp","dexp"),
         origin=FALSE, nstart=100, upb=0.5, ErrorMethod=c("mc","sp"),
         nsim=1000, plot=TRUE, outfile=NULL) {
    ### Stop if not.
    stopifnot(ncol(Curvedata)==3L, 
              all(Curvedata[,1L,drop=TRUE]>=0),
              all(Curvedata[,3L,drop=TRUE]>0),
              (is.vector(Ltx) && length(Ltx)==2L)|| 
              (is.matrix(Ltx) && ncol(Ltx)==2L), all(Ltx>0), 
              all(model %in% c("line","exp","lexp","dexp")),
              length(origin)==1L, is.logical(origin),
              is.numeric(nstart), nstart>=10L, nstart<=5000L,
              length(upb)==1L, is.numeric(upb), upb>0, upb<=100,
              all(ErrorMethod %in% c("mc","sp")),
              is.numeric(nsim), nsim>=100L, nsim<=3000L,
              length(plot)==1L, is.logical(plot), 
              is.null(outfile) || is.character(outfile))
    ###
    dose<-as.numeric(Curvedata[,1L,drop=TRUE])
    doseltx<-as.numeric(Curvedata[,2L,drop=TRUE])
    sltx<-as.numeric(Curvedata[,3L,drop=TRUE])
    ndat<-nrow(Curvedata)
    ninltx<-ifelse(is.vector(Ltx), 1L, nrow(Ltx))
    n2<-if (model[1L]=="line") {
            1L+!origin
        } else if (model[1L]=="exp") {
            2L+!origin
        } else if (model[1L]=="lexp") {
            3L+!origin
        } else if (model[1L]=="dexp") {
            4L+!origin
        } # end if  
    if (model[1L]!="line" && max(Ltx)>max(doseltx)) {
        stop("Error: Ltx is too large!")
    } # end if
    ###npoints<-length(levels(factor(Curvedata[,1L,drop=TRUE])))
    if ( (model[1L]=="line" && ndat<n2) ||
         (model[1L]=="exp" && ndat<n2)  ||
         (model[1L]=="lexp" && ndat<n2) ||
         (model[1L]=="dexp" && ndat<n2) ) {
        stop("Error: data points is not enough!")
    } # end if
    ###
    if(is.vector(Ltx)) {
        inltx<-matrix(Ltx, nrow=1L) 
    } else {
        inltx<-Ltx
    } # end if
    outDose<-matrix(0, nrow=ninltx, ncol=2L)
    mcED<-matrix(0, nrow=nsim, ncol=ninltx)
    pars<-stdp<-vector(length=n2)
    model1<-if (model[1L]=="line") {
        0L } else if (model[1L]=="exp") {
        1L } else if (model[1L]=="lexp") {
        2L } else if (model[1L]=="dexp") {
        3L } # end if
    origin1<-ifelse(origin==TRUE, 0L, 1L)
    method<-ifelse(ErrorMethod[1L]=="sp", 0L, 1L)
    fvec1<-vector(length=ndat)
    fvalue<-0
    message<-0
    ###
    res<-.Fortran("calED",as.double(dose),as.double(doseltx),as.double(sltx),
                  as.integer(ndat),as.integer(ninltx),as.integer(n2),
                  as.double(inltx),outDose=as.double(outDose),mcED=as.double(mcED),
                  pars=as.double(pars),stdp=as.double(stdp),as.double(upb),
                  as.integer(model1),as.integer(origin1),as.integer(method),
                  as.integer(nstart),as.integer(nsim),fvec1=as.double(fvec1),
                  fvalue=as.double(fvalue),message=as.integer(message),PACKAGE="numOSL")
    if (res$message!=0) {
        stop("Error: fail in equivalent dose calculation!")
    } # end if
    ###
    LMpars<-cbind(res$pars,res$stdp)
    colnames(LMpars)<-c("Pars","Std.Pars")
    rownames(LMpars)<-(c("a","b","c","d","e"))[1L:n2]
    ###
    fit.value<-cbind(dose, doseltx, res$fvec1)
    colnames(fit.value)<-c("Redose", "Lx/Tx", "Fit.Lx/Tx")
    rownames(fit.value)<-paste("Redose", seq(ndat), sep="")
    ###
    ED<-matrix(res$outDose, ncol=2L)
    rownames(ED)<-paste("NO.", seq(ninltx), sep="")
    colnames(ED)<-c("ED", "Std.ED")
    ###
    if (ErrorMethod[1L]=="mc") {
        mcED<-matrix(res$mcED, ncol=ninltx)
        colnames(mcED)=paste("NO.", seq(ninltx), sep="")
    } else {
        mcED<-NULL
    } # end if
    ###
    output<-list("LMpars"=LMpars, 
                 "value"=res$fvalue,
                 "fit.value"=fit.value, 
                 "ED"=ED)
    ###
    if (ErrorMethod[1L]=="mc" && !is.null(outfile)) {
        write.csv(mcED, file=paste(outfile,".csv",sep=""))
    } # end if
    ###
    ###
    ###
    Plot1<-function(xvalue,yvalue,simED,Mainname) {
        ###
        lowerX<-min(min(dose,xvalue),0)*1.2
        upperX<-max(dose,xvalue)*1.2
        lowerY<-min(min(doseltx,yvalue),0)*1.2
        upperY<-max(doseltx,yvalue)*1.2
        ###
        par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
        plot(NA, NA, main=Mainname, xlab="Dose (Gy)", ylab="Standardised OSL",
             las=0, cex.lab=1, cex.main=1.25, xlim=c(lowerX,upperX),
             ylim=c(lowerY,upperY), xaxs="i", yaxs="i", lab=c(7,7,9))
        if (!is.null(simED) && all(yvalue>0)) {
            dmcED<-density(simED)
            dxy<-cbind(dmcED$x,dmcED$y)
            dxy[,2L]<-(dxy[,2L,drop=TRUE]-min(dxy[,2L,drop=TRUE]))/
                      (max(dxy[,2L,drop=TRUE])-min(dxy[,2L,drop=TRUE]))*
                      yvalue[1L]*0.8
            polygon(dxy, col="grey")
            rug(simED, quiet=TRUE)
        } # end if
        ###
        points(dose, doseltx, pch=21, cex=3, bg="white")
        ###
        x<-NULL
        if(origin==TRUE)  {
            if(model[1L]=="line") {
                curve(LMpars[1L,1L]*x, type="l", add=TRUE, 
                      lw=2, from=lowerX, to=upperX)
            } else if(model[1L]=="exp") {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x)),  
                     type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
            } else if(model[1L]=="lexp")  {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x, 
                     type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
            } else if(model[1L]=="dexp") {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+
                     LMpars[3L,1L]*(1.0-exp(-LMpars[4L,1L]*x)), 
                     type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
            } # end if
       } else {
           if(model[1L]=="line") {
               curve(LMpars[1L,1L]*x+LMpars[2L,1L], type="l", 
                     add=TRUE, lw=2, from=lowerX, to=upperX)
           } else if(model[1L]=="exp") {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L], 
                     type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
           } else if(model[1L]=="lexp")  {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x+LMpars[4L,1L], 
                     type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
           } else if(model[1L]=="dexp") {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+
                     LMpars[3L,1L]*(1.0-exp(-LMpars[4L,1L]*x))+LMpars[5L,1L], 
                     type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
           } # end if
       } # end if
       ###
       points(x=xvalue[1L], y=yvalue[1L], pch=23, cex=3, bg="grey")
       ###
       arrowsData<-Curvedata[Curvedata[,3L,drop=TRUE]>=1e-3,,drop=FALSE]
       options("warn"=-1)
       if (nrow(arrowsData)>=1L) {
           arrows(x0=arrowsData[,1L,drop=TRUE], 
                  y0=arrowsData[,2L,drop=TRUE]-arrowsData[,3L,drop=TRUE]/2L,
                  x1=arrowsData[,1L,drop=TRUE],
                  y1=arrowsData[,2L,drop=TRUE]+arrowsData[,3L,drop=TRUE]/2L,
                  code=3, lwd=2.5, angle=90, length=0.05, col="black")
       } # end if
       ###
       if (xvalue[2L]>=1e-3) {
           arrows(x0=xvalue[1L]-xvalue[2L]/2L, y0=yvalue[1L],
                  x1=xvalue[1L]+xvalue[2L]/2L, y1=yvalue[1L],
                  code=3, lwd=2.5, angle=90, length=0.05, col="black")
       } # ned if
       options("warn"=0)
       ###
       lines(x=c(0,xvalue[1L],xvalue[1L]),
             y=c(yvalue[1L],yvalue[1L],0),
             lty="dashed", lwd=2)    
       legend("topleft", legend=paste("ED=",round(xvalue[1L],2L), "+-",
              round(xvalue[2L],2L), " (Gy)", sep=""), yjust=2, ncol=1L,
              cex=par("cex"), bty="n")
        grid()
        box(lwd=2)
        par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
        ###
    } # end function Plot1
    ###
    ###
    ###
    if (plot==TRUE) {
        if (ninltx==1L) {
            Plot1(xvalue=ED[1L,,drop=TRUE],
                  yvalue=inltx[1L,,drop=TRUE],
                  simED=mcED[,1L,drop=TRUE],
                  Mainname="Growth Curve")
        } else {
           if (ninltx<=3L) {
               par(mfrow=c(1L,ninltx))
           } else if (ninltx==4L) {
               par(mfrow=c(2L,2L))
           } else {  
               par(mfrow=c(ifelse(ninltx%%4L==0L,ninltx/4L, 
                   ninltx%/%4L+1L), 4L))
           } # end if
           for (i in seq(ninltx)) {
               Plot1(xvalue=ED[i,,drop=TRUE],
                     yvalue=inltx[i,,drop=TRUE],
                     simED=mcED[,i,drop=TRUE],
                     Mainname=paste("NO.",i,sep=""))
           } # end for
           par(mfrow=c(1L,1L))          
       } # end if
   } # end if
   return(output)
} # end function calED.default
#####