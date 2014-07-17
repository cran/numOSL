#####################################################################################################################
###     The R package numOSL
###     R package for tackling basic numeric problems in optically stimulated luminescence dating.
###     Copyright (C) 2013 Peng Jun 
###
###     Written by Peng Jun.
###
###     This file is the R scripts (all R functions) for package numOSL.
###     URL: http://CRAN.R-project.org/package=numOSL
### 
###     The numOSL package is free software: you can redistribute it  
###     or modify it under the terms of the GNU General Public License as 
###     published by the Free Software Foundation, either version 3 of the 
###     License, or any later version.
###
###     Package numOSL is distributed in the hope that it will be 
###     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
###     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###     GNU General Public License for more details.
###
###     You should have received a copy of the GNU General Public License
###     along with numOSL.  If not, see <http://www.gnu.org/licenses/>.
###
###     If you wish to report bugs please contact the author:
###
###     Peng Jun
###     pengjun10@mails.ucas.ac.cn
###     University of Chinese Academy of Sciences.
############################################### FUNCTION RadialPlotter ###############################################
### ******************************************************************************************************************
#### Function RadialPlotter() is used to perform Galbraith's
#### statistical age models analysis. Models that can be 
###  analyzed include: CAM, FMM and MAM.
###
###     Author: Peng Jun, 2013.04.20, revised in 2013.07.26, 
###             revised in 2013.09.19, revised in 2013.10.09; revised in 2014.04.19.
###
### References: Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard 
###             Errors. Technometrics, 30 (3), pp. 271-281.
###
###             Galbraith, R.F., Roberts, R.G., 2012. Statistical aspects of equivalent dose and error calculation and  
###             display in OSL dating: An overview and some recommendations. Quaternary Geochronology, 11, pp. 1-27.
###
### Mandatory arguments---
###     EDdata: a data.frame, equivalent doses (two columns).
### ******************************************************************************************************************
RadialPlotter<-
function(EDdata, ncomp=0, addsigma=0, maxiter=500,
         maxcomp=6, algorithm=c("lbfgsb","port"),
         eps=.Machine$double.eps^0.5, plot=TRUE,
         pcolor="blue", psize=1.5, kratio=0.3,
         zscale=NULL, samplename=NULL)  {
  UseMethod("RadialPlotter")
} ###
### set default method for function RadialPlotter().
RadialPlotter.default<-
function(EDdata, ncomp=0, addsigma=0, maxiter=500,
         maxcomp=6, algorithm=c("lbfgsb","port"),
         eps=.Machine$double.eps^0.5, plot=TRUE,
         pcolor="blue", psize=1.5, kratio=0.3,
         zscale=NULL, samplename=NULL)  {
  ### stop if not. 
  ###
  ### For argument "EDdata".
  if(!is.data.frame(EDdata))    stop("Error: EDdata should be of type data.frame!")
  if(ncol(EDdata)!=2L)          stop("Error: EDdata should contain two columns!")
  if(nrow(EDdata)<5L)           stop("Error: number of equivalent dose values is too small!")
  if(!is.numeric(as.matrix(EDdata)))  
                                stop("Error: all elements in EDdata should be of type numeric!")
  if(any(!is.finite(unlist(unclass(EDdata)))))   
                                stop("Error: EDdata should not contain non-finite value!")
  if(any(EDdata[,1L,drop=TRUE]<=0))     
                                stop("Error: all equivalent dose values should larger than zero!")
  if(any(EDdata[,2L,drop=TRUE]<=0))     
                                stop("Error: std.error of equivalent dose should larger than zero!")
  ###
  ### For argument "maxcomp".
  if(!is.numeric(maxcomp))      stop("Error: maxcomp should be of type numeric!")
  if(length(maxcomp)!=1L)       stop("Error: maxcomp should be an one-element vector!")
  if(!is.finite(maxcomp))       stop("Error: maxcomp should not be a non-finite value!")
  if(abs(maxcomp-round(maxcomp))
     >=.Machine$double.eps^0.5) stop("Error: maxcomp should be an integer!")                                
  if(maxcomp<2L)                stop("Error: maxcomp should not smaller than 2!")
  if(maxcomp>=nrow(EDdata))     stop("Error: maxcomp should not exceed the number of equivalent doses!")
  if(maxcomp>9L)                stop("Error: maxcomp should not exceed 9!")
  ###
  ### For argument "ncomp".
  if(!is.numeric(ncomp))        stop("Error: ncomp should be of type numeric!")
  if(length(ncomp)!=1L)         stop("Error: ncomp should be an one-element vector!")
  if(!is.finite(ncomp))         stop("Error: ncomp should not be a non-finite value!")
  if(!ncomp %in% (-2L:maxcomp)) stop("Error: ncomp should be an integer ranging from -2 to maxcomp!")
  ###
  ### For argument "addsigma".
  if(!is.numeric(addsigma))     stop("Error: addsigma should be of type numeric!")
  if(length(addsigma)!=1L)      stop("Error: addsigma should be an one-element vector!")
  if(!is.finite(addsigma))      stop("Error: addsigma should not be a non-finite value!")
  if(addsigma<0)                stop("Error: addsigma should not smaller than 0!")
  ###
  ### For argument "maxiter".
  if(!is.numeric(maxiter))      stop("Error: maxiter should be of type numeric!")
  if(length(maxiter)!=1L)       stop("Error: maxiter should be an one-element vector!")
  if(!is.finite(maxiter))       stop("Error: maxiter should not be a non-finite value!")
  if(abs(maxiter-round(maxiter))
     >=.Machine$double.eps^0.5) stop("Error: maxiter should be an integer!")
  if(maxiter<10L)               stop("Error: maxiter should not less than 10!")
  if(maxiter>10000L)            stop("Error: maxiter should not larger than 10000!")
  ###
  ### For argument "algorithm".
  if(!is.character(algorithm))  stop("Error: algorithm should be of type character!")
  if(length(algorithm)>2L)      stop("Error: algorithm should contain at most two elements!")
  if(!all(algorithm %in% c("lbfgsb","port"))) 
                                stop("Error: only algorithm 'lbfgsb' or 'port' is available!")
  ###
  ### For argument "eps".
  if(!is.numeric(eps))          stop("Error: eps should be of type numeric!")
  if(length(eps)!=1L)           stop("Error: eps should be an one-element vector!")
  if(!is.finite(eps))           stop("Error: eps should not be a non-finite value!")
  if(eps<0)                     stop("Error: eps should be not smaller than 0!")
  ###
  ### For argument "plot".
  if(!is.logical(plot))         stop("Error: plot should be either TRUE or FALSE!")
  if(length(plot)!=1L)          stop("Error: plot should be an one-element vector!")
  ###
  ### For argument "pcolor".
  if(!is.character(pcolor))     stop("Error: pcolor should be of type character!")
  if(length(pcolor)!=1L)        stop("Error: pcolor should be an one-element vector!")
  ###
  ### For argument "psize".
  if(!is.numeric(psize))        stop("Error: psize should be of type numeric!")
  if(length(psize)!=1L)         stop("Error: psize should be an one-element vector!")
  if(!is.finite(psize))         stop("Error: psize should not be a non-finite value!")
  if(psize<0.1)                 stop("Error: psize is too small!")
  if(psize>5)                   stop("Error: psize is too large!")
  ###
  ### For argument "kratio".
  if(!is.numeric(kratio))       stop("Error: kratio should be of type numeric!")
  if(length(kratio)!=1L)        stop("Error: kratio should be an one-element vector!")
  if(!is.finite(kratio))        stop("Error: kratio should not be a non-finite value!")
  if(abs(kratio)>3)             stop("Error: magnitude of kratio is too large!")
  ###
  ### For argument "zscale".
  if(!is.null(zscale) &&
     !is.numeric(zscale))       stop("Error: zscale should be NULL or a numeric vector!")
  if(!is.null(zscale)) {
    if(any(!is.finite(zscale))) stop("Error: zscale should not contain non-finite value!")
  } # end if
  ###
  ### For argument "samplename".
  if(!is.null(samplename) && !is.character(samplename)) stop("Error: samplename should be NULL or be of type character!")
  if(!is.null(samplename) && length(samplename)!=1L)    stop("Error: samplename should be an one-element vector!")
  ###
  ###
  ### Set n, ED, sED.
  n<-nrow(EDdata) 
  ed<-drop(EDdata[,1L])
  error<-drop(EDdata[,2L])
  ###
  ###
  ### Function for radial plot drawing, this code is revised from Rex Galbraith. 
  ### *******************************************************************************
  RadialPlot<-function(Data,Pars,addsigma,zscale,samplename,kratio,pcolor,psize)  {
    z.i<-log(Data[,1L,drop=TRUE])
    se.i<-sqrt( (Data[,2L,drop=TRUE]/Data[,1L,drop=TRUE])^2L + addsigma^2L)
    ### If Pars=NULL, lines will not be plotted out.
    if(!is.null(Pars))  {
      Pars<-log(Pars)
    } # end if
    plot.area_ratio<-kratio
    z.0<-sum(z.i/se.i^2L)/sum(1.0/se.i^2L)
    x.i<-1.0/se.i
    y.i<-(z.i-z.0)/se.i
    h<-plot.area_ratio*(max(x.i)-min(x.i))/(max(y.i)-min(y.i))
    r.0<-max(sqrt(x.i^2L+y.i^2L*h^2L))
    if(is.null(zscale))  {
      mkest<-round(exp(pretty(z.i)), 0L)
      if(sd(z.i)/mean(z.i)<0.05)  {
        mkest<-round(exp(pretty(c(min(z.i)*0.98, max(z.i)*1.02))), 0L)
      } # end if
    } else  {
      mkest<-zscale
    } # end if
    mkr<-range(mkest)
    circ.x<-(0:200)/200
    circ.x<-mkr[1L]+circ.x*(mkr[2L]-mkr[1L])
    circ.z<-log(circ.x)
    zmkest<-log(mkest) 
    circ.x<-r.0/sqrt(1.0+h^2L*(circ.z-z.0)^2L)
    circ.y<-(circ.z-z.0)*circ.x
    xmk1<-r.0/sqrt(1.0+h^2L*(zmkest-z.0)^2L)
    xmk2<-xmk1*1.01
    xmk3<-1.01*1.04*r.0/sqrt(1.0+h^2L*(zmkest-z.0)^2L)
    ymk1<-(zmkest-z.0)*xmk1
    ymk2<-(zmkest-z.0)*xmk2
    ymk3<-(zmkest-z.0)*xmk3 
    may<-max(abs(y.i))
    yaxis.max<-if(may>12)   {
      may  } else if(may<3) {
      8    } else  {
      12   } # end if
    par(mfrow=c(1,1), oma=c(2,0,0,0), mar=c(5,4,4,4.5), xpd=TRUE, las=1, cex=1)
    plot(NA, NA, xlim=c(0,max(x.i)), ylim=c(-yaxis.max,yaxis.max), xlab="", ylab="",
         xaxt="n", yaxt="n", main=samplename, bty="n", typ="n", xaxs="i", yaxs="i")      
    lenpar<-length(Pars) 
    maxx<-which.max(circ.x)
    ### If Pars!=NULL, starting drawing lines out.
    if(!is.null(Pars))  {
      for(i in 1L:lenpar)  {
        if(Pars[i]-z.0>=0.0)  {
          Ci<-approx(circ.y[maxx:201L], circ.x[maxx:201L], ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i), rule=2, method="constant", f=0.5)$y
        }  else  {
          Ci<-approx(circ.y[1L:(maxx-1L)], circ.x[1L:(maxx-1L)], ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i), rule=2, method="constant",f=0.5)$y
        } # end if
       lines(c(0,Ci), c(0,(Pars[i]-z.0)*Ci), lty=1, col="black", lwd=1.5)
      } # end for
    } # end if
    lines(circ.x, circ.y, lwd=2)
    segments(xmk1, ymk1, xmk2, ymk2)
    text(xmk3, ymk3, round(mkest,2L), cex=1)
    text(max(x.i)*1.15,0,"De (Gy)", cex=1, srt=90)
    par(mfrow=c(1,1), oma=c(2,0,0,0.5), mar=c(5,4,4,4.5), xpd=TRUE, las=1, cex=1, new=TRUE)
    plot(NA, NA, xlim=c(0,max(x.i)), ylim=c(-yaxis.max, yaxis.max), yaxt="n", xaxt="n",
         xaxs="i", yaxs="i", ylab="Standardised Estimate", xlab="", type="n", bty="n")
    points(x.i, y.i, pch=21, col="black", cex=psize, bg=pcolor)
    axis(side=1, line=4, cex.axis=1, lwd=2)
    mtext(side=1, line=6, "Precision", cex=1)
    reticks.labels<-round(1/axTicks(side=1)*100,1L)
    reticks.values<-axTicks(side=1)
    axis(side=1, at=reticks.values[-1L], lwd=2, labels=reticks.labels[-1L], line=4, cex.axis=1, tck=0.02, padj=-4)
    mtext(side=1, line=1.5,"Relative Error (%)", cex=1)
    axis(side=2,at=c(-2,-1,0,1,2), lwd=2, labels=c("-2","","0","","2"), cex.axis=1)
    par(oma=c(0,0,0,0), xpd=FALSE,
        las=0, new=FALSE, mar=c(5,4,4,2)+0.1)
  } ### end function RadialPlot.
  ### ********************************************************************************
  ###
  ###
  tol<-.Machine$double.eps^0.5  
  ###
  ###   
  ### Function Rmam for port routines to do
  ### optimization for MAM3 or MAM4, programmed in 2013.07.25.
  ### ******************************************************
  Rmam<-function(EDdata,ncomp,addsigma,maxiter) {
    ### -1 for 3-parameter, -2 for 4-parameter.
    stopifnot(ncomp %in% c(-1L,-2L))
    ### add spread dispersion to relative standard error of Equivalent Dose.
    x<-sqrt((EDdata[,2L,drop=TRUE]/EDdata[,1L,drop=TRUE])^2L+addsigma^2L)
    ### change Equivalent Dose to log-scale.
    y<-log(drop(EDdata[,1L]))
    ###
    ### function minfunc to be minimized.
    ### **************************************************************************
    minfunc<-function(p)  {
      ### here ncomp is a gloable variable.
      if(ncomp==-1L) {
        ### for MAM3.
        u0<-(p[2L]/(p[3L])^2L+y/x^2L)/(1.0/(p[3L])^2L+1.0/x^2L)
        sigma0<-1.0/(sqrt(1.0/(p[3L])^2L+1.0/x^2L)) 
        prop1<-p[1L]/sqrt(2.0*pi)/x*exp(-(y-p[2L])^2L/2.0/x^2L)
        prop2<-(1.0-p[1L])/sqrt(2.0*pi*((p[3L])^2L+x^2L))* 
               (1.0-pnorm((p[2L]-u0)/sigma0))/(1.0-pnorm((p[2L]-p[2L])/p[3L]))*
                exp(-(y-p[2L])^2L/2.0/((p[3L])^2L+x^2L))
        return(-sum(log(prop1+prop2)))
      } else if(ncomp==-2L) {
        ### for MAM4.
        u0<-(p[3L]/(p[4L])^2L+y/x^2L)/(1.0/(p[4L])^2L+1.0/x^2L)
        sigma0<-1.0/(sqrt(1.0/(p[4L])^2L+1.0/x^2L)) 
        prop1<-p[1L]/sqrt(2.0*pi)/x*exp(-(y-p[2L])^2L/2.0/x^2L)
        prop2<-(1.0-p[1L])/sqrt(2.0*pi*((p[4L])^2L+x^2L))* 
               (1.0-pnorm((p[2L]-u0)/sigma0))/(1.0-pnorm((p[2L]-p[3L])/p[4L]))*
                exp(-(y-p[3L])^2L/2.0/((p[4L])^2L+x^2L))
        return(-sum(log(prop1+prop2)))
      } # end if
    } # end function minfunc
    ### ****************************************************************************
    ###
    ### function RmamErr to estimate parameters' Std.Err, programmed in 2013.07.25.
    ### ****************************************************************************
    RmamErr<-function(pars, ED, sED, addsigma) {
      npars<-length(pars)
      nED<-length(ED)
      spars<-vector(length=npars)
      ### tol<-.Machine$double.eps^0.3
      message<-0
      ### finite-difference approximation hessian matrix.
      res<-.Fortran("aperr", as.double(ED), as.double(sED), as.integer(nED),
                    as.double(pars), as.double(addsigma), spars=as.double(spars),
                    as.integer(npars), as.double(tol), message=as.integer(message), package="numOSL")
      if(res$message!=0) {
        return(NULL)
      } else {
        return(res$spars)
      } # end if    
    } # end function RmamErr
    ### ******************************************************************************
    ###
    ### set boundaries for parameters.
    if(ncomp==-1L) {
      lower<-c(1e-4, min(y), 1e-3)
      upper<-c(0.99, max(y), 5.0)
    } else if(ncomp==-2L) {
      lower<-c(1e-4, min(y), min(y), 1e-3)
      upper<-c(0.99, max(y), max(y), 5.0)
    } # end if
    ### find out three cluster using K-Means then sorting them by ascending order.
    kclus<-suppressWarnings(try(kmeans(x=y, centers=3L, iter.max=20L, nstart=100L), silent=TRUE))
    ### don't forget always doing an error checking.
    if(class(kclus)=="try-error") stop("Error: k-means fails in parameter initializing!")
    kclus<-sort(kclus$centers)
    ### set gamma, mu and sdy values.
    gama<-kclus[1L]
    mu<-c(kclus[2L], mean(kclus), mean(y))
    sdy<-sd(y)
    ### start the optimization by various initial values calling nlminb.
    cmaxlik<-maxlik<-1e30
    errorflag<-1
    bexist<-0
    if(ncomp==-1L) {
      ### for three parameters (MAM3).
      for(i in 1L:3L) {
        for(j in 1L:5L) {
          ini<-c(0.01*(5.0)^(i-1L), gama, sdy*0.4*j)
          res<-suppressWarnings(try(nlminb(start=ini,objective=minfunc,scale=1,
               control=list(iter.max=maxiter),lower=lower,upper=upper),silent=TRUE))
          if(class(res)!="try-error" && res$convergence==0) { 
            aperr<-RmamErr(res$par,drop(EDdata[,1L]),drop(EDdata[,2L]),addsigma=addsigma)
            if(!is.null(aperr) &&                 # std.error can be approximated
               res$objective<maxlik &&            # decreasing in minus maxlik
               all(abs(res$par-lower)>=1e-4) &&   # within low boundaries
               all(abs(res$par-upper)>=1e-4) ) {  # within up boundaries
              pars<-res$par
              error<-aperr
              maxlik<-res$objective
              errorflag<-0
            } # end if
            if(!is.null(aperr) &&                 # std.error can be approximated
               res$objective<cmaxlik &&           # decreasing in minus maxlik
               errorflag==1) {                    # might near to the boundaries
              cpars<-res$par
              cerror<-aperr
              cmaxlik<-res$objective
              bexist<-1
            } # end if
          } # end if
        } # end for
      } # end for
    } else if(ncomp==-2L) {
      ### for four parameters (MAM4).
      for(i in 1L:3L) {
        for(j in 1L:3L) {
          for(k in 1L:3L) {
            ini<-c(0.01*(5.0)^(i-1L), gama, mu[k], sdy*0.5*j)
            res<-suppressWarnings(try(nlminb(start=ini,objective=minfunc,scale=1,
                 control=list(iter.max=maxiter),lower=lower,upper=upper),silent=TRUE))
            if(class(res)!="try-error" && res$convergence==0) { 
               aperr<-RmamErr(res$par,drop(EDdata[,1L]),drop(EDdata[,2L]),addsigma=addsigma)
               if(!is.null(aperr) &&                 # std.error can be approximated
                  res$objective<maxlik &&            # decreasing in minus maxlik
                  all(abs(res$par-lower)>=1e-4) &&   # within low boundaries
                  all(abs(res$par-upper)>=1e-4) &&   # within up boundaries
                  res$par[2L]<=res$par[3L] ) {       # gamma must not larger than mu
                 pars<-res$par
                 error<-aperr
                 maxlik<-res$objective
                 errorflag<-0
               } # end if
              if(!is.null(aperr) &&                  # std.error can be approximated
                  res$objective<cmaxlik &&           # decreasing in minus maxlik
                  res$par[2L]<=res$par[3L] &&        # gamma must not larger than mu
                  errorflag==1) {                    # might near to the boundaries
                 cpars<-res$par
                 cerror<-aperr
                 cmaxlik<-res$objective
                 bexist<-1
              } # end if
            } # end if
          } # end for
        } # end for
      } # end for
    } # end if
    ### now output the results.
    if(errorflag==0) {
      ### reset parameters before output.
      if(ncomp==-1L) {
        pars[2L]<-exp(pars[2L])
        error[2L]<-pars[2L]*error[2L]
      } else if(ncomp==-2L) {
        pars[2L:3L]<-exp(pars[2L:3L])
        error[2L:3L]<-pars[2L:3L]*error[2L:3L]
      } # end if
      return(list("errorflag"=0,"maxlik"=-maxlik,
                  "pars"=cbind(pars,error),
                  "BIC"=ifelse(ncomp==-1L,
                   2.0*maxlik+3L*log(length(y)),
                   2.0*maxlik+4L*log(length(y)))))
    } else if(bexist==1){
      ### reset parameters before output.
      if(ncomp==-1L) {
        cpars[2L]<-exp(cpars[2L])
        cerror[2L]<-cpars[2L]*cerror[2L]
      } else if(ncomp==-2L) {
        cpars[2L:3L]<-exp(cpars[2L:3L])
        cerror[2L:3L]<-cpars[2L:3L]*cerror[2L:3L]
      } # end if
      return(list("errorflag"=1,"maxlik"=-cmaxlik,
                  "pars"=cbind(cpars,cerror),
                  "BIC"=ifelse(ncomp==-1L,
                   2.0*cmaxlik+3L*log(length(y)),
                   2.0*cmaxlik+4L*log(length(y)))))
    } else {
      return(NULL)
    } # end if
  } # end function Rmam
  ###
  ###
  ### ******************************************************
  ### Now calculate the commom age model based equivalent dose.
  z.i<-log(ed)
  se.i<-sqrt( (error/ed)^2L + addsigma^2L )
  commonED<-sum(z.i/se.i^2L)/sum(1.0/se.i^2L)
  scommonED <- 1.0/sqrt(sum(1.0/se.i^2L))
  ###
  ### Start do optimization.
  errorflag<-0
  maxlik<-BIC<-0
  loopcomp<-0
  BILI<-NULL
  ### For CAM and FMM analysis.
  if(ncomp %in% (0L:maxcomp))  {
    if (ncomp==0L)  {
      goodcomp<-0
      BILI<-matrix(0, nrow=maxcomp-1L, ncol=2L)
      ### Call Fortran subroutine FineComp() to pick out the appropriate number of component automatically.
      Results<-.Fortran("FineComp", as.double(ed), as.double(error), as.integer(n),
                        goodcomp=as.integer(goodcomp), as.double(addsigma), as.integer(maxiter),       
                        as.double(eps), as.integer(maxcomp), BILI=as.double(BILI), package="numOSL")
      BILI<-Results$BILI
      dim(BILI)<-c(maxcomp-1L, 2L)
      colnames(BILI)<-c("BIC","Maxlik")
      rownames(BILI)<-paste(rep("k=", maxcomp-1L), c(2L:maxcomp), sep ="") 
      ncomp<-Results$goodcomp
      loopcomp<-1L
    } # end if
    spars<-pars<-matrix(0,nrow=2L, ncol=ncomp)
    ### Call Fortran subroutine FineED1() to calculate parameters for central age model or finite mixture age model.
    Results<-.Fortran("FineED1", as.double(ed), as.double(error), as.integer(n), as.integer(ncomp),
                      as.double(addsigma), pars=as.double(pars), spars=as.double(spars), maxlik=as.double(maxlik),
                      BIC=as.double(BIC), as.integer(maxiter), as.double(eps), as.double(tol),
                      errorflag=as.integer(errorflag), package="numOSL")  
    ### Store errorflag, BIC, maxlik, ParsAndErrors.
    errorflag<-Results$errorflag
    BIC<-Results$BIC
    maxlik<-Results$maxlik
    ParsAndErrors<-cbind(matrix(Results$pars, byrow=TRUE, ncol=2L),
                         matrix(Results$spars, byrow=TRUE, ncol=2L))
    if(ncomp==1L)  { 
      ### Set matrix ParsAndErrors for CAM analysis.
      ParsAndErrors<-matrix(ParsAndErrors[,c(1L,3L,2L,4L)])
      rownames(ParsAndErrors)<-c("Sigma", "Std.Sigma", "ED", "Std.ED")
      colnames(ParsAndErrors)<-"CAM.Value"
      ### Plot or not (CAM).
      if(plot==TRUE)  {
        RadialPlot(Data=EDdata, Pars=ParsAndErrors[3L], addsigma=addsigma, zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } # end if
    }  else {  
      ### Set matrix ParsAndErrors for FMM analysis.
      ParsAndErrors<-ParsAndErrors[,c(1L,3L,2L,4L)][order(ParsAndErrors[,2L,drop=TRUE]), , drop=FALSE]
      colnames(ParsAndErrors)<-c("P","Std.P","ED","Std.ED") 
      rownames(ParsAndErrors) <- paste(rep("comp", ncomp), c(1L:ncomp), sep = "") 
      ### Plot or not (FMM).
      if(plot==TRUE) {
        ### Only plot out fmmED if FMM was fitting successfully.
        if(Results$errorflag==0)  {
          RadialPlot(Data=EDdata, Pars=ParsAndErrors[,3L,drop=TRUE], addsigma=addsigma, zscale=zscale, 
                     samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
        } else {
          ### Alternatively, just plot observtions out.
          RadialPlot(Data=EDdata, Pars=NULL, addsigma=addsigma, zscale=zscale, 
                     samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
        } # end if
      } # end if
    } # end if
  }  else  { 
    ### Set default algorithm.
    algorithm<-algorithm[1L]
    ### For MAM analysis. 
    if(algorithm=="lbfgsb") {
      ### For algorithm 'lbfgsb'.
      pars<-spars<-vector(length=2L-ncomp)
      ### Call Fortran subroutine MAM() to fit the minimum age models.
      Results<-.Fortran("MAM",as.double(ed), as.double(error), as.integer(n), pars=as.double(pars),
                        spars=as.double(spars), maxlik=as.double(maxlik), BIC=as.double(BIC), as.integer(2L-ncomp),
                        as.double(addsigma), as.integer(maxiter), as.double(tol), errorflag=as.integer(errorflag), package="numOSL")
      ### Store errorflag, BIC, maxlik, ParsAndErrors.
      errorflag<-Results$errorflag
      BIC<-Results$BIC
      maxlik<-Results$maxlik
      ParsAndErrors<-cbind(Results$pars, Results$spars)
    } else if(algorithm=="port") {
      ### For algorithm "port".
      Results<-Rmam(EDdata, ncomp=ncomp, addsigma=addsigma, maxiter=maxiter)
      ### Checking error of R function Rmam() and store errorflag, BIC, maxlik, ParsAndErrors.
      if(!is.null(Results)) {
        errorflag<-Results$errorflag
        BIC<-Results$BIC
        maxlik<-Results$maxlik
        ParsAndErrors<-Results$pars
      } else {
        ParsAndErrors<-matrix(0, nrow=2L-ncomp, ncol=2L) 
        errorflag<-1
        BIC<-maxlik<-0
      } # end if
    } # end if
    ### Set matrix ParsAndErrors for MAM.
    colnames(ParsAndErrors)<-c("Pars", "Std.Pars")
    if(ncomp==-1L)  {
      rownames(ParsAndErrors)<-c("P", "gamma", "sigma")
    } else if(ncomp==-2L) {
      rownames(ParsAndErrors)<-c("P", "gamma", "mu", "sigma")
    } # end if
    ### Plot or not (MAM).
    if(plot==TRUE)  {
      if(all(ParsAndErrors[,2L,drop=TRUE]<=0)==FALSE)  {
        ### Only plot mamED if MAM was fitted succesfully.
        RadialPlot(Data=EDdata, Pars=ParsAndErrors[2L,1L], addsigma=addsigma, zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } else {
        ### Alternatively, just plot observations out.
        RadialPlot(Data=EDdata, Pars=NULL, addsigma=addsigma, zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } # end if
    } # end if
  } # end if
  ###
  ### At this point the age model analysis have been terminated.
  ### organize the returned list.
  out<-list("errorflag"=errorflag,
            "commonED"=exp(commonED)*c(1L,scommonED),
            "ncomp"=ncomp, "maxcomp"=maxcomp, "loopcomp"=loopcomp,
            "pars"=ParsAndErrors, "BIC"=BIC, "maxlik"=maxlik, "BILI"=BILI)
  ### Set out to be of class 'RadialPlotter'.
  class(out)<-"RadialPlotter"
  ### Set out invisible
  invisible(out)
} # end function RadialPlotter.default.
###
### Set print method for object "RadialPlotter", programmed in 2013.07.25.
print.RadialPlotter<-function(x,...) {
  obj<-x
  ### Output the result to the terminal screen.
  cat("\n")
  cat("==================Results of RadialPlotter==================","\n\n")
  ### Error message.
  cat("Error message:",obj$errorflag,"\n\n")
  ### Common ED.
  cat("Common equivalent dose:",round(obj$commonED[1L],5L),"+-",round(obj$commonED[2L],5L),"\n\n")
  if(obj$ncomp==1L)  {
    ### Overdispersion in CAM.
    cat("Overdispersion and central dose are:","\n\n")
  }  else if (obj$ncomp %in% c(0L,2L:obj$maxcomp))  {
    ### FMM ED
    if(obj$errorflag==0) {
      if(obj$loopcomp==1L)  {
        ### Best number of components.
        cat("The best component number is estimated at:",obj$ncomp,"\n\n")
      } # end if
      ### Parameters for CAM, FMM.
      cat("Parameters and standard errors are:","\n")
    } else {
      cat("Error: singular model, some parameters' standard errors can not be estimated!","\n")
      cat("Note: in the Radial Plot, no characteristic dose was drawn!","\n\n")
    } # end if
  }  else  {
    ### MAM ED
    if(all(obj$pars[,2L,drop=TRUE]<=0)==FALSE)  {
      if(obj$errorflag==1) {
        ### Boundary checking.
        cat("Warning: at least one estimated parameter is near to the boundary!","\n\n")
      } # end if
      ### Parameters for MAM.
      cat("Parameters and standard errors are:","\n")
    } else  {
      ### Fails in MAM.
      cat("Error: singular model, some parameters' standard errors can not be estimated!","\n")
      cat("Note: in the Radial Plot, no characteristic dose was drawn!","\n\n")
    } # end if
  } # end if
  ### BIC, maxlik for CAM, FMM, MAM.
  if( (!obj$ncomp %in% (-2L:-1L) && obj$errorflag==0) || 
      (obj$ncomp %in% (-2L:-1L) && all(obj$pars[,2L,drop=TRUE]<=0)==FALSE) )  {
    cat("---------------------------------------","\n")
    print(round(as.data.frame(obj$pars),5L))
    cat("---------------------------------------","\n\n")
    cat("                      BIC value:",round(obj$BIC,5L),"\n\n")
    cat("Maximum logged likelihood value:",round(obj$maxlik,5L),"\n\n")
    ### Looped BIC and maxlik values for FMM.
    if(obj$loopcomp==1L)  {
      cat("BIC values and Maximum logged likelihood values are:","\n")
      cat("---------------------------","\n")
      print(round(as.data.frame(obj$BILI),5L))
      cat("---------------------------","\n\n")
    } # end if
  } # end if
  cat("===========================The End===========================","\n\n")
} # end function print.RadialPlotter.
################################################ END FUNCTION RadialPlotter ##########################################
###
################################################ FUNCTION calED ######################################################
### ******************************************************************************************************************
### Function calED() is used to calculate a number of equivalent dose values.
###
###     Author: Peng Jun, 2013.06.23; revised in 2013.07.26; revised in 2013.08.02; revised in 2013.09.18;
###             revised in 2013.11.12; revised in 2014.04.05.
###
### References: Roberts,H.M. and Duller,G.A.T., 2004. Standardised growth curves for optical dating 
###             of sediment using multiple-grain aliquots. Radiation Measurements 38, pp. 241-252.
###
###             Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from 
###             single aliquot regenerative dose measurements. Ancient TL 25, pp. 15-24.
###
### Mandatory arguments---
### Curvedata: a data.frame, three columns, data used for constructing a growth curve.
###       Ltx: a vector or a matrix, two columns that contains standlized signals and its standard errors.
### ******************************************************************************************************************
calED<-
function(Curvedata, Ltx, model=c("line","exp","line+exp"),
         origin=FALSE, nstart=100, upb=1, ErrorMethod=c("mc","sp"),
         nsim=1000, plot=TRUE, outfile=NULL)  {
   UseMethod("calED")
} ###
### Default method for function calED().
calED.default<-
function(Curvedata, Ltx, model=c("line","exp","line+exp"),
         origin=FALSE, nstart=100, upb=1, ErrorMethod=c("mc","sp"),
         nsim=1000, plot=TRUE, outfile=NULL)  {
  ### Stop if not.
  ###
  ### For argument "Curvedata".
  if(!is.data.frame(Curvedata))              stop("Error: Curvedata should be of type data.frame!")
  if(ncol(Curvedata)!=3L)                    stop("Error: Curvedata should contain three columns!")
  if(!is.numeric(as.matrix(Curvedata)))      stop("Error: all values in Curvedata should be of type numeric!")
  if(any(!is.finite(unlist(unclass(Curvedata)))))   
                                             stop("Error: Curvedata should not contain non-finite value!")
  if(any(Curvedata[,]<0))                    stop("Error: all values in Curvedata should not smaller than zero!")
  ###
  ### For argument "Ltx".
  if(!is.vector(Ltx) && !is.matrix(Ltx))     stop("Error: Ltx should be a vector or a matrix!")
  if(is.vector(Ltx)) {
    if(length(Ltx)!=2L)                      stop("Error: Ltx should be a two-element vector!")
    Ltx<-matrix(Ltx,nrow=1L)
  } # end if
  if(is.matrix(Ltx))  {
    if(nrow(Ltx)<1L)                         stop("Error: the row number of Ltx should not less than one!")
    if(nrow(Ltx)>20L)                        stop("Error: the row number of Ltx should not exceed 20!")
    if(ncol(Ltx)!=2L)                        stop("Error: Ltx should contain two columns!")
  } # end if
  if(!is.numeric(Ltx))                       stop("Error: all values in Ltx should be of type numeric!")
  if(any(!is.finite(Ltx)))                   stop("Error: Ltx should not contain non-finite value!")
  if(any(Ltx<=0))                            stop("Error: all values in Ltx should larger than zero!")
  if(model[1L]!="line") {
    if(max(Ltx)>max(Curvedata[,2L,drop=TRUE]))        
                                             stop("Error: some values in Ltx is too large!")
  } # end if
  ###
  ### For argument "model".
  if(!is.character(model))                   stop("Error: model should be of type character!")
  if(length(model)>3L)                       stop("Error: model should contain no more than three elements!")
  if(!all(model %in% 
      c("line","exp","line+exp")))           stop("Error: incorrect model, only 'line', 'exp' or 'line+exp' is available!")
  ###
  ### For argument "origin".
  if(!is.logical(origin))                    stop("Error: origin should be of type logical!")
  if(length(origin)!=1L)                     stop("Error: origin should be an one-element vector!")
  ###
  ### For argument "nstart".
  if(!is.numeric(nstart))                    stop("Error: nstart should be of type numeric!")
  if(length(nstart)!=1L)                     stop("Error: nstart should be an one-element vector!")
  if(!is.finite(nstart))                     stop("Error: nstart should not be a non-finite value!")
  if(abs(nstart-round(nstart))
     >=.Machine$double.eps^0.5)              stop("Error: nstart should be an integer!")
  if(nstart<5L)                              stop("Error: nstart should not below 5!")
  if(nstart>5000L)                           stop("Error: nstart should not exceed 5000!") 
  ###
  ### For argument "upb".         
  if(!is.numeric(upb))                       stop("Error: upb should be of type numeric!")
  if(length(upb)!=1L)                        stop("Error: upb should be an one-element vector!")
  if(!is.finite(upb))                        stop("Error: upb should not be a non-finite value!")
  if(upb<=0)                                 stop("Error: upb should larger than 0!")
  if(upb>100)                                stop("Error: upb should not exceed 100!")
  ###
  ### For argument "ErrorMethod".
  if(!is.character(ErrorMethod))             stop("Error: ErrorMethod should be of type character!")
  if(length(ErrorMethod)>2L)                 stop("Error: ErrorMethod should contain at most two elements!")
  if(!all(ErrorMethod %in% c("sp","mc")))    stop("Error: incorrect ErrorMethod, only 'sp' or 'mc' is available!")
  ###
  ### For argument "nsim".
  if(!is.numeric(nsim))                      stop("Error: nsim should be of type numeric!")
  if(length(nsim)!=1L)                       stop("Error: nsim should be an one-element vector!")
  if(!is.finite(nsim))                       stop("Error: nsim should not be a non-finite value!")
  if(abs(nsim-round(nsim))
     >=.Machine$double.eps^0.5)              stop("Error: nsim should be an integer!")
  if(nsim<100L)                              stop("Error: nsim should not smaller than 100!")
  if(nsim>3000L)                             stop("Error: nsim should not larger than 3000!")
  ###
  ### For argument "plot".
  if(!is.logical(plot))                      stop("Error: plot should be either TRUE or FALSE!")
  if(length(plot)!=1L)                       stop("Error: plot should be an one-element vector!")
  ###
  ###
  ### For argument "outfile".
  if(!is.null(outfile) && !is.character(outfile)) {
    stop("Error: outfile should be NULL or of type character!")
  } # end if
  if (!is.null(outfile)) {
    if(length(outfile)!=1L)  { 
        stop("Error: outfile should be an one-element vector!")
    } # end if
    if(ErrorMethod[1L]=="sp")  {
        cat('Note: ErrorMethod="sp", so nothing will be written to outfile!\n')
    } # end if
  } # end if
  ###
  ### Set default fitting model (linear).
  model<-model[1L]
  ### Set default ErrorMethod (monte carlo).
  ErrorMethod<-ErrorMethod[1L]
  ###
  ### Check if observations are enough for model fitting
  if(origin==TRUE) {
    if(model=="line" && length(levels(factor(Curvedata[,1L,drop=TRUE])))<1L) {
      stop("Error: fitting a linear model (origin) needs at least one paired independent observations!")
    } # end if
  } else {
    if(model=="line" && length(levels(factor(Curvedata[,1L,drop=TRUE])))<2L) {
      stop("Error: fitting a linear model (non-origin) needs at least two paired independent observations!")
    } # end if
  } # end if
  if(model=="exp" && length(levels(factor(Curvedata[,1L,drop=TRUE])))<3L) {
     stop("Error: fitting a exponential model needs at least three paired independent observations!")
  } # end if
  if(model=="line+exp" && length(levels(factor(Curvedata[,1L,drop=TRUE])))<4L) {
     stop("Error: fitting a linear+exponential model needs at least four paired independent observations!")
  } # end if
  ### 
  ### Parameters used in subroutine calED() or calED0().
  ### ===================================================
  Dose<-Curvedata[,1L,drop=TRUE]
  Doseltx<-as.matrix(Curvedata[,2L:3L])
  ndose<-nrow(Curvedata)
  ninltx<-nrow(Ltx)
  outDose<-matrix(0,nrow=ninltx, ncol=2L)
  mcED<-matrix(0,nrow=nsim,ncol=ninltx)
  npars<-if (origin==TRUE) {
             if(model=="line") {
                 1L } else if (model=="exp") {
                 2L } else if (model=="line+exp") {
                 3L } # end if
         } else {
             if(model=="line") {
                 2L } else if (model=="exp") {
                 3L } else if (model=="line+exp") {
                 4L } # end if
         } # end if
  pars<-parserrors<-vector(length=npars)
  predtval<-vector(length=ndose)
  value<- -99.0
  method<-ifelse(ErrorMethod=="sp",1L,2L)
  errorflag<-vector(length=2L)
  ### =====================================================
  ### Calculate equivalent doses.
  fFortran<-if(origin==FALSE) "calED" else "calED0"
  res<-.Fortran(fFortran, as.double(Dose), as.double(Doseltx), as.integer(ndose), as.integer(ninltx), as.double(Ltx), outDose=as.double(outDose),
                mcED=as.double(mcED), pars=as.double(pars), as.integer(npars), predtval=as.double(predtval), as.double(upb), as.integer(nstart),
                parserrors=as.double(parserrors), value=as.double(value), as.integer(method), as.integer(nsim), errorflag=as.integer(errorflag), package="numOSL")
  ### Error checking.
  if(res$errorflag[1L]!=123 || res$errorflag[2L]!=0) {
    stop("Error: fail in growth curve fitting, the model may be inappropriate, or upb needs to be modified!")
  } # end if
  ###  
  ###
  ### Set LMpars.
  LMpars<-cbind(res$pars, res$parserrors)
  colnames(LMpars)<-c("Pars", "Std.Pars")
  rownames(LMpars)<-(c("a", "b", "c", "d"))[1L:npars]    
  ###
  ### Set fit.value.
  fit.value<-cbind(Curvedata[,1L,drop=TRUE], Curvedata[,2L,drop=TRUE], res$predtval)
  colnames(fit.value)<-c("ReDose", "Lx/Tx", "Fit.Lx/Tx")
  rownames(fit.value)<-paste("ReDose.", 1L:ndose, sep="")
  ###
  ### Set ED values.
  ED<-matrix(res$outDose,ncol=2L)
  rownames(ED)<-paste("NO.", 1L:ninltx, sep="")
  colnames(ED)<-c("ED", "Std.ED")
  ###
  if (ErrorMethod=="mc") {
      mcED<-data.frame(matrix(res$mcED, ncol=ninltx))
      colnames(mcED)=paste("NO.", 1L:ninltx,sep="")
  } else {
      mcED<-NULL
  } # end if
  ###
  ### Results for output.
  output<-list("LMpars"=LMpars,"value"=res$value,
               "fit.value"=fit.value,"ED"=ED)
  ###
  if (ErrorMethod[1L]=="mc" && !is.null(outfile))  {
      write.csv(mcED,file=paste(outfile,".csv",sep=""))
  } # end if
  ###
  ###
  ### A function for plotting the results.
  PlotRegion<-function(xvalue,yvalue,simED,Mainname) {
      ###
      lowerX<-min(min(Curvedata[,1L,drop=TRUE],xvalue),0)
      upperX<-max(Curvedata[,1L,drop=TRUE],xvalue)*1.2
      lowerY<-min(min(Curvedata[,2L,drop=TRUE],yvalue),0)
      upperY<-max(Curvedata[,2L,drop=TRUE],yvalue)*1.2
      ### Set pars
      par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
      plot(NA, NA, main=Mainname, xlab="Dose (Gy)", ylab="Standardised OSL", las=0, cex.lab=1,
           cex.main=1.25, xlim=c(lowerX,upperX), ylim=c(lowerY,upperY), xaxs="i", yaxs="i", lab=c(7,7,9) )
      ###
      ### Add a density plot.
      if (!is.null(simED))  {
          dmcED<-density(simED)
          dxy<-cbind(dmcED$x,dmcED$y)
          dxy[,2L]<-(dxy[,2L,drop=TRUE]-min(dxy[,2L,drop=TRUE]))/(max(dxy[,2L,drop=TRUE])-min(dxy[,2L,drop=TRUE]))*yvalue[1L]*0.8
          polygon(dxy, col="grey")
          rug(simED, quiet=TRUE)
      } # end if
      ###
      ### Add ReDose as points.
      points(Curvedata[,1L:2L], pch=21, cex=3, bg="white")
      ###
      ### Add a fitted growth curve.
      x<-NULL
      if(origin==TRUE)  {
        if(model=="line") {
          curve(LMpars[1L,1L]*x, type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
        } else if(model=="exp") {
          curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x)), type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
        } else if(model=="line+exp")  {
          curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x, type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
        } # end if
      } else {
        if(model=="line") {
          curve(LMpars[1L,1L]*x+LMpars[2L,1L], type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
        } else if(model=="exp") {
          curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L], type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
        } else if(model=="line+exp")  {
          curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x+LMpars[4L,1L], type="l", add=TRUE, lw=2, from=lowerX, to=upperX)
        } # end if
      } # end if
      ###
      ### Add point (xvalue .VS. yvalue).
      points(xvalue[1L], yvalue[1L], pch=23, cex=3, bg="grey")
      ###
      ### Add error bars to Curvedata.
      arrowsData<-Curvedata[Curvedata[,3L,drop=TRUE]>=5e-3,,drop=FALSE]
      options("warn"=-1)
      if (nrow(arrowsData)>=1L)  {
          arrows(x0=arrowsData[,1L,drop=TRUE], y0=arrowsData[,2L,drop=TRUE]-arrowsData[,3L,drop=TRUE]/2L,
                 x1=arrowsData[,1L,drop=TRUE], y1=arrowsData[,2L,drop=TRUE]+arrowsData[,3L,drop=TRUE]/2L,
                 code=3, lwd=2.5, angle=90, length=0.05, col="black")
      } # end if
      ###
      ### Add an error bar to calculated xvalues.
      if (xvalue[2L]>=5e-3)  {
          arrows(x0=xvalue[1L]-xvalue[2L]/2L, y0=yvalue[1L], 
                 x1=xvalue[1L]+xvalue[2L]/2L, y1=yvalue[1L],
                 code=3, lwd=2.5, angle=90, length=0.05, col="black")
      } # end if
      options("warn"=0)
      ###
      ### Add a dash line. 
      lines(c(0, xvalue[1L], xvalue[1L]), c(yvalue[1L], yvalue[1L], 0), lty="dashed", lwd=2)
      ###
      ### Add a legend.
       legend("topleft", legend=paste("ED=", round(xvalue[1L],2L), " +- ", round(xvalue[2L],2L), " (Gy)", sep=""),
              yjust=2, ncol=1L, cex=par("cex"), bty="n")
      ###
      ### Add grid and box.
      grid()
      box(lwd=2)
      ### Reset plot(par) before leaving.
      par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
      ###
  } # end function PlotRegion.
  ###
  ###
  ### Plot the results or not.
  if (plot==TRUE)  {
      ### Set plot regions.
      if (ninltx==1L) {
          PlotRegion(xvalue=ED[1L,,drop=TRUE], yvalue=Ltx[1L,,drop=TRUE], simED=mcED[,1L], Mainname="Growth Curve")
      } else {
          if (ninltx<=3L)  {
              par(mfrow=c(1L,ninltx))
          } else if(ninltx==4L)  {
              par(mfrow=c(2L,2L))
          } else {
              par(mfrow=c(ifelse(ninltx%%4L==0L, ninltx/4L, ninltx%/%4L+1L), 4L))
          } # end if
          ###
          for (i in seq(ninltx)) {
              PlotRegion(xvalue=ED[i,,drop=TRUE], yvalue=Ltx[i,,drop=TRUE], simED=mcED[,i], Mainname=paste("NO.",i,sep="")) 
          } # end for
          ###
          par(mfrow=c(1L,1L))
      } # end if
  } # end if
  return(output)
} # end function calED.   
################################################# END FUNCTION calED #################################################
###
################################################## FUNCTION dbED #####################################################
### ******************************************************************************************************************
### Function dbED() is used to calculate statistical characters for equivalent dose values.
###
###     Author: Peng Jun, 2013.12.01.
###
### References: Galbraith, R., 2003. A simple homogeneity test for estimates of dose obtained using OSL.
###             Ancient TL, 21(2), pp.75-77.
###             
###             Bailey, R.M., Arnold, L.J., 2006. Statistical modelling of single grain quartz De 
###             distributions and an assessment of procedures for estimating burial dose.
###             Quaternary Science Reviews, 25 (19-20), pp.2475-2502.
###
###             Arnold, L.J., Bailey, R.M., Tucher, G.E., 2007. Statistical treatment of fluvial dose distributions
###             from southern Colorado arroyo deposits. Quaternary Geochronology, 2, pp.162-167.
###
###             Duller, G.A.T., 2008. Single-grain optical dating of Quaternary sediments: why aliquot size matters 
###             in luminescence dating. Boreas, 37 (4), pp. 589-612.
###
### Mandatory arguments---
###     EDdata: a data.frame, equivalent doses (two columns).
### *******************************************************************************************************************
dbED<-
function(EDdata, aliquot=c("sa","sg"), plot=TRUE, from=NULL,
         to=NULL, step=0.1, nbin=15, samplename=NULL)  {
  UseMethod("dbED")
} ###
### Set default method for function dbED.
dbED.default<-
function(EDdata, aliquot=c("sa","sg"), plot=TRUE, from=NULL,
         to=NULL, step=0.1, nbin=15, samplename=NULL)  {
  ### Stop if not
  ###
  ### For argument "EDdata".
  if(!is.data.frame(EDdata))    stop("Error: EDdata should be of type data.frame!")
  if(ncol(EDdata)!=2L)          stop("Error: EDdata should contain two columns!")
  if(!is.numeric(as.matrix(EDdata)))  
                                stop("Error: all elements in EDdata should be of type numeric!")
  if(any(!is.finite(unlist(unclass(EDdata)))))   
                                stop("Error: EDdata should not contain non-finite value!")
  if(any(EDdata[,1L,drop=TRUE]<=0))       
                                stop("Error: equivalent dose should larger than zero!")
  if(any(EDdata[,2L,drop=TRUE]<=0))       
                                stop("Error: std.error of equivalent dose should larger than zero!")
  if(nrow(EDdata)<5L)           stop("Error: at least five equivalent doses should be provided!")
  ###
  ### For argument "aliquot".
  if(!is.character(aliquot))    stop("Error: aliquot should be of type numeric!")
  if(length(aliquot)==1L) {
    if(!aliquot %in% c("sa","sg")) 
                                stop("Error: aliquot should be 'sa' or 'sg'!")
  } else if(length(aliquot)==2L) {
    if(!all(aliquot %in% c("sa","sg"))) 
                                stop("Error: only aliquot type of 'sa' and 'sg' are available!")
  } else {                      stop("Error: aliquot should not contain more than 2 element!")
  } # end if
  ### For argument "plot".
  if(!is.logical(plot))         stop("Error: plot should be of type logical!")
  if(length(plot)!=1L)          stop("Error: plot should be an one-element vector!")
  ###
  ### For argument "from".
  if(!is.numeric(from) &&
     !is.null(from))            stop("Error: from should be NULL or of type numeric!")
  if(!is.null(from)) {
    if(length(from)!=1L)        stop("Error: from should be an one-elemtn vector!")
    if(!is.finite(from))        stop("Error: from should not be a non-finite value!")
    if(from<0)                  stop("Error: from should not below zero!")
  } # end if
  ###
  ### For argument "to".
  if(!is.numeric(to) && 
     !is.null(to))              stop("Error: to should be NULL or of type numeric!")
  if(!is.null(to)) {
    if(length(to)!=1L)          stop("Error: to should be an one-elemtn vector!")
    if(!is.finite(to))          stop("Error: to should not be a non-finite value!")
    if(to<=0)                   stop("Error: to should larger than zero!")
  } # end if
  ###
  ### For argument "step".
  if(!is.numeric(step))         stop("Error: step should be of type numeric!")
  if(length(step)!=1L)          stop("Error: step should be an one-element vector!")
  if(!is.finite(step))          stop("Error: step should not be a non-finite value!")
  if(step<=1e-5)                stop("Error: step should larger than 1e-5!")
  ###
  ### For argument "nbin".
  if(!is.numeric(nbin))         stop("Error: nbin should be of type numeric!")
  if(length(nbin)!=1L)          stop("Error: nbin should be an one-element vector!")
  if(!is.finite(nbin))          stop("Error: nbin should not be a non-finite value!")
  if(abs(nbin-round(nbin))
     >=.Machine$double.eps^0.5) stop("Error: nbin should be an integer!")
  if(nbin<=2L)                  stop("Error: nbin should larger than 2!")
  ###
  ### For argument "samplename".
  if(!is.null(samplename) && !is.character(samplename)) 
                                stop("Error: samplename should be NULL or of type character!")
  if(!is.null(samplename)) {
    if(length(samplename)!=1L)  stop("Error: samplename should be an one-element vector!")
  } # end if
  ###
  ###
  nED<-nrow(EDdata)
  if(is.null(from)) {
    from<-min(EDdata[,1L,drop=TRUE])*0.7
  } # end if
  if(is.null(to)) {
    to<-max(EDdata[,1L,drop=TRUE])*1.1
  } # end if
  ###
  if(to<=from)                  stop("Error: to should larger than from!")
  ### 
  ### Plot or not
  if(plot==TRUE) {
    ###  
    layout(cbind(c(1L,1L,1L,2L,2L),c(1L,1L,1L,2L,2L)))
    par(mar=c(0,5.1,3.1,2.1))
    ###
    ### Sort equivalent dose values in increasing order.
    EDdata<-EDdata[order(EDdata[,1L,drop=TRUE],decreasing=FALSE),,drop=FALSE]
    ###
    spreadED<-seq(from=from, to=to, by=step)
    pdfMat<-matrix(nrow=length(spreadED), ncol=nED)
    for(i in 1L:nED) {
      pdfMat[,i]<-dnorm(x=spreadED, mean=EDdata[i,1L,drop=TRUE], sd=EDdata[i,2L,drop=TRUE], log=FALSE)
    } # end if
    pdfED<-rowSums(pdfMat)/sum(rowSums(pdfMat))
    ###
    ### Draw a probability density function plot.
    plot(spreadED, pdfED, main=ifelse(is.null(samplename),"De Distribution",samplename), ylab="Probability Density", 
         xlim=c(from,to), type="l", lwd=2, las=0, lab=c(7,7,9), xaxs="r", yaxs="r", xaxt="n", cex.main=1.5, cex.lab=1.5)
    xTicks<-axTicks(side=1L)
    rug(x=EDdata[,1L,drop=TRUE], side=1, lwd=1, quiet=TRUE)
    rug(x=EDdata[,1L,drop=TRUE], side=3, lwd=1, quiet=TRUE)
    ###
    ### Add equivalent doses as points to current plot.
    par("new"=TRUE)
    plot(EDdata[,1L,drop=TRUE], 1L:nED, xlab="",ylab="", xlim=c(from,to), type="p", las=0, 
         lab=c(7,7,9), xaxs="r", yaxs="r", pch=23, cex=3, bg="gray", xaxt="n", yaxt="n")
    axis(side=4, at=axTicks(side=2L), labels=axTicks(side=2L), las=0)
    par("new"=FALSE)
    ### 
    ### Add error bars to current plot.
    if(any(EDdata[,2L,drop=TRUE]>=1e-3)) {
      arrows(EDdata[EDdata[,2L,drop=TRUE]>=1e-3,1L,drop=TRUE]-EDdata[EDdata[,2L,drop=TRUE]>=1e-3,2L,drop=TRUE]/2L, (1L:nED)[EDdata[,2L,drop=TRUE]>=1e-3],
             EDdata[EDdata[,2L,drop=TRUE]>=1e-3,1L,drop=TRUE]+EDdata[EDdata[,2L,drop=TRUE]>=1e-3,2L,drop=TRUE]/2L, (1L:nED)[EDdata[,2L,drop=TRUE]>=1e-3],
             code=3, lwd=1.5, angle=90, length=0.05, col="black")
    } # end if
    box(lwd=1.5)
    par(mar=c(5.1,5.1,0,2.1))
    breaks<-pretty(EDdata[,1L,drop=TRUE], n=nbin)
    ###
    ### Draw a histogram.
    HIST<-hist(EDdata[,1L,drop=TRUE], breaks=breaks, main=NULL, xlab="De (Gy)", col="gray50",
               xlim=c(from,to), las=0, lab=c(7,7,9), xaxs="r", yaxs="i", xaxt="n", cex.main=1.5, cex.lab=1.5)
    histTicks<-axTicks(side=2L)
    axis(side=1, at=xTicks, labels=xTicks)
    axis(side=4, at=histTicks, labels=histTicks, las=0)
    box(lwd=1.5)
    legend(ifelse(HIST$mids[which.max(HIST$counts)]>median(xTicks),"topleft","topright"), legend=paste("N=",nED,sep=""), cex=1.5, bty="n")
    par(mar=c(5,4,4,2)+0.1)
    layout(1L)
  } # end if
  ###
  ###
  ### Calculate weighted skewness and kurtosis.
  weight<-EDdata[,1L,drop=TRUE]/EDdata[,2L,drop=TRUE]
  meanED<-mean(EDdata[,1L,drop=TRUE])
  sdED<-sd(EDdata[,1L,drop=TRUE])
  ### weighted skewness.
  skewness<-sum(weight*((EDdata[,1L,drop=TRUE]-meanED)/sdED)^3L)/sum(weight)
  STDskewness<-sqrt(6L/nED)
  ### kurtosis.
  kurtosis<-nED*(nED+1L)/(nED-1L)/(nED-2L)/(nED-3L)*sum(((EDdata[,1L,drop=TRUE]-meanED)/sdED)^4L)-
            3L*(nED-1L)^2L/(nED-2L)/(nED-3L)
  STDkurtosis<-sqrt(24L/nED)
  ###
  ### Significant tests of weighted skewness and kurtosis. 
  sigSkewness<-ifelse(abs(skewness)>2L*STDskewness, "Yes", "No")
  sigKurtosis<-ifelse(abs(kurtosis)>2L*STDkurtosis, "Yes", "No")
  ###
  ###
  ### Do a homogeneity test.
  resED<-EDdata[,2L,drop=TRUE]/EDdata[,1L,drop=TRUE]
  logED<-log(EDdata[,1L,drop=TRUE])
  weightED<-sum(logED/(resED)^2L)/sum(1.0/(resED)^2L)
  chisqValue<-sum((logED-weightED)^2L/(resED)^2L)
  pvalue<-pchisq(q=chisqValue, df=nED-1L, lower.tail = FALSE)
  homogeneity<-ifelse(pvalue>=0.01, "Yes", "No")
  ###
  ###
  ### Default aliquot (Single-aliquot).
  aliquot<-aliquot[1L]
  ### Calculate the overdispersion.
  Sigma<-RadialPlotter(EDdata,ncomp=1L,addsigma=0.0,maxcomp=3L,plot=FALSE)$pars[1L]
  ###
  ###
  ### Estimate a proposed age-model.
  proposal<-
  if(aliquot=="sg") {
    if(sigSkewness=="Yes") {
      if(Sigma>1.0) {
        "Lowest5%"
      } else {
        if(skewness>7.25) {
          "MAM3"
        } else {
          if(sigKurtosis=="Yes") {
            if(Sigma>0.4) {
              "MAM4"
            } else {
              if(kurtosis<0.0) {
                "Lowest5%"
              } else if(kurtosis>2.18) {
                if(Sigma>0.1) {
                  "MAM3"
                } else {
                  "MAM4"
                } # end if
              } else {
                if(Sigma>0.2) {
                  "MAM4"
                } else {
                  "MAM3"
                } # end if
              } # end if
            } # end if
          } else {
            if(kurtosis>STDkurtosis) {
              if(Sigma>0.1) {
                "MAM3"
              } else {
                "MAM4"
              } # end if
            } else {
              if(Sigma>0.1) {
                "Lowest5%"
              } else {
                "MAM4"
              } # end if
            } # end if
          } # end if
        } # end if
      } # end if
    } else {
      if(skewness/2L/STDskewness<0.5) {
        "CAM"
      } else {
        "Lowest5%"
      } # end if
    } # end if
  } else if(aliquot=="sa") {
    rstd<-sd(EDdata[,1L,drop=TRUE])/weighted.mean(x=EDdata[,1L,drop=TRUE],w=1.0/EDdata[,2L,drop=TRUE])
    if(skewness>STDskewness) {
      "MAM3"
    } else {
      if(kurtosis>STDkurtosis) {
        "MAM3"
      } else {
        if(Sigma>0.4) {
          "MAM3"
        } else {
          if(rstd>0.4) {
            "MAM3"
          } else {
            "CAM"
          } # end if
        } # end if 
      } # end if
    } # end if
  } # end if
  ###
  ###
  ### Summary equivalent dose values. 
  summaryED<-c(summary(EDdata[,1L,drop=TRUE]), "Std.dev"=sd(EDdata[,1L,drop=TRUE]), "OverDisp"=Sigma)
  ###
  ### Calculate mean of the lowest x% of equivalent dose.
  quantileED<-quantile(EDdata[,1L,drop=TRUE], probs=c(0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.7,0.9))
  percentED<-apply(as.matrix(quantileED), MARGIN=1L, function(x,y) x>=y, EDdata[,1L,drop=TRUE])
  Std.pED<-apply(percentED, MARGIN=2L, function(x,y) sd(y[x])/sqrt(sum(x)), EDdata[,1L,drop=TRUE])
  pED<-apply(percentED, MARGIN=2L, function(x,y) mean(y[x]), EDdata[,1L,drop=TRUE])
  ###
  list("summaryED"=round(summaryED,2L),
       "mlpED"    =round(rbind(pED,Std.pED),2L),
       "homogeneity"=data.frame("chisq.value"=chisqValue,"p.value"=pvalue,"homogeneity"=homogeneity,row.names=""),
       "skewness"=data.frame("Skewness"=round(skewness,2L),"Std.Skew"=round(STDskewness,2L),"significant"=sigSkewness,row.names=""),
       "kurtosis"=data.frame("Kurtosis"=round(kurtosis,2L),"Std.Kurt"=round(STDkurtosis,2L),"significant"=sigKurtosis,row.names=""),
       "proposal"=proposal)
  ###
} # end function dbED().
################################################ END FUNCTION dbED ################################################### 
####
################################################ FUNCTION decomp #####################################################
### ******************************************************************************************************************
### Function decomp() is used to decompose OSL signal
### curve, either type "CW" or "LM" can be analyzed.
###
###      Author: Peng Jun, 2013.06.05, revised in 2013.07.25, revised in 2013.09.18, 
###              revised in 2013.10.06, revised in 2013.11.17.
###
###  References: Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
###              applications. Geochronometria 13, 135–141.
###
###              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
###              OSL decay curves. Radiation Measurements 41, 886-891.
###
###              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
###              luminescence components in different quartz samples: implications for dose measurement. 
###              Radiation Measurements, 37 (4-5), pp. 441-449.
###
### Mandatory arguments---
###     Sigdata: a dataframe, OSL decay curve (two columns).
### ******************************************************************************************************************
decomp<-
function(Sigdata, ncomp=3, typ=c("cw","lm"),
         control.args=list(), transf=FALSE,
         LEDpower=60, LEDwavelength=470, plot=TRUE,
         xlog=TRUE, lwd=3, samplename=NULL, outfile=NULL) {
  UseMethod("decomp")
} ###
### default method for function decomp().
decomp.default<-
function(Sigdata, ncomp=3, typ=c("cw","lm"),
         control.args=list(), transf=FALSE,
         LEDpower=60, LEDwavelength=470, plot=TRUE,
         xlog=TRUE, lwd=3, samplename=NULL, outfile=NULL) {
  ### stop if not.
  ###
  ### For argument "Sigdata".
  if(!is.data.frame(Sigdata))             stop("Error: Sigdata should be of type data.frame!")
  if(ncol(Sigdata)<2L)                    stop("Error: Sigdata should contain at least two columns!")
  if(ncol(Sigdata)>2L && typ[1L]!="cw")   stop("Error: only signal table of type CW-OSL can be analyzed currently!")
  if(!is.numeric(as.matrix(Sigdata)))           
                                          stop("Error: all elements in Sigdata should be of type numeric!")
  if(any(!is.finite(unlist(unclass(Sigdata)))))   
                                          stop("Error: Sigdata should not contain non-finite value!")
  if(any(Sigdata[,1L,drop=TRUE]<0))       stop("Error: the first column of Sigdata [time] should larger than zero!")
  if(any(Sigdata[,-1L]<0))                stop("Error: all signal values should larger than zero!")
  ###
  ### For argument "ncomp".
  if(!is.numeric(ncomp))                  stop("Error: ncomp should be of type numeric!")
  if(length(ncomp)!=1L)                   stop("Error: ncomp should be an one-element vector!")
  if(!is.finite(ncomp))                   stop("Error: ncomp should not be a non-finite value!")
  if(!ncomp %in% (1L:7L) )                stop("Error: ncomp should be a integer in the space [1,7]!")
  ###
  ### For argument "typ".
  if(!is.character(typ))                  stop("Error: typ should be of type character!")
  if(length(typ)==1L) {
    if(!typ %in% c("cw","lm"))            stop("Error: typ should be 'cw' or 'lm'!")
    if(ncol(Sigdata)==2L) {
      if(typ=="cw" && which.max(Sigdata[,2L,drop=TRUE])>5L)      
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
      if(typ=="lm" && which.max(Sigdata[,2L,drop=TRUE])<=5L)     
                                          stop("Error: incorrect typ, signal may not be of type LM-OSL!")
    } else if(ncol(Sigdata)>2L) {
      if(which.max(apply(Sigdata[,-1L],MARGIN=1L,mean))>5L)        
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
    } # end if
  } else if(length(typ)==2L) {
    if(!all(typ %in% c("cw","lm")))       stop("Error: incorrect typ, please choose one between 'cw' and 'lm'!")
    if(ncol(Sigdata)==2L) {
      if(typ[1L]=="cw" && which.max(Sigdata[,2L,drop=TRUE])>5L)      
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
      if(typ[1L]=="lm" && which.max(Sigdata[,2L,drop=TRUE])<=5L)     
                                          stop("Error: incorrect typ, signal may not be of type LM-OSL!")
    } else if(ncol(Sigdata)>2L) {
      if(which.max(apply(Sigdata[,-1L],MARGIN=1L,mean))>5L)        
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
    } # end if
  } else {                                stop("Error: typ should contain no more than two elements!")
  } # end if
  ###
  ### For argument "control.args".
  if(class(control.args)!="list")         stop("Error: control.args should be a list!")
  if(length(control.args)>5L)             stop("Error: number of parameters in control.args should not exceed 5!")
  if(!all(names(control.args) %in% list("factor","f","cr","maxiter","tol")))              
                                          stop("Error: incorrect parameter in control.args!")
  ###
  ### For argument "transf".
  if(!is.logical(transf))                 stop("Error: transf should be of type logical!")
  if(length(transf)!=1L)                  stop("Error: transf should be an one-element vector!")
  ###
  ### For argument "LEDpower".
  if(!is.numeric(LEDpower))               stop("Error: LEDpower should be of type numeric!")
  if(length(LEDpower)!=1L)                stop("Error: LEDpower should be an one-element vector!")
  if(!is.finite(LEDpower))                stop("Error: LEDpower should not be a non-finite value!")
  if(LEDpower<=0)                         stop("Error: LEDpower should larger than zero!")
  ###
  ### For argument "LEDwavelength".
  if(!is.numeric(LEDwavelength))          stop("Error: LEDwavelength should be of type numeric!")
  if(length(LEDwavelength)!=1L)           stop("Error: LEDwavelength should be an one-element vector!")
  if(!is.finite(LEDwavelength))           stop("Error: LEDwavelength should not be a non-finite value!")
  if(LEDwavelength<=0)                    stop("Error: LEDwavelength should larger than zero!")
  ###
  ### For argument "plot".
  if(!is.logical(plot))                   stop("Error: plot should be either TRUE or FALSE!")
  if(length(plot)!=1L)                    stop("Error: plot should be an one-element vector!")
  ###
  ### For argument "xlog".
  if(!is.logical(xlog))                   stop("Error: xlog should be TRUE or FALSE!")
  if(length(xlog)!=1L)                    stop("Error: xlog should be an one-element vector!")
  ###
  ### For argument "lwd".
  if(!is.numeric(lwd))                    stop("Error: lwd should be of type numeric!")
  if(length(lwd)!=1L)                     stop("Error: lwd should be an one-element vector!")
  if(!is.finite(lwd))                     stop("Error: lwd should not be a non-finite value!")
  if(lwd<=0)                              stop("Error: lwd should larger than zero!")
  ###
  ### For argument "samplename".
  if(!is.null(samplename) && !is.character(samplename))           
                                          stop("Error: samplename should be NULL or of type character!")
  if(!is.null(samplename)) {
    if(length(samplename)!=1L)            stop("Error: sample name should be an one-element vector!")
  } # end if
  ###
  ### For argument "outfile".
  if(!is.null(outfile) && !is.character(outfile))             
                                          stop("Error: outfile should be NULL or of type character!")
  if(!is.null(outfile)) {
    if(ncol(Sigdata)>2L)                  stop("Error: signal table is not allowed to be written out!")
    if(length(outfile)!=1L)               stop("Error: outfile should be an one-element vector!")
  } # end if
  ###
  ###
  ### Public parameters for Fortran   
  ### subroutine decomp() and fitlm().
  if(ncol(Sigdata)==2L) {
    ### Decompose signal curve.
    ntim<-nrow(Sigdata)
    tim<-drop(Sigdata[,1L])
    sig<-drop(Sigdata[,2L])
  } else if(ncol(Sigdata)>2L) {
    ### Decompose signal table (CW).
    ntim<-nrow(Sigdata)*(ncol(Sigdata)-1L)
    Sigdata<-Sigdata[,c(1L,order(Sigdata[1L,-1L],decreasing=TRUE)+1L)]
    tim<-rep(Sigdata[,1L,drop=TRUE],ncol(Sigdata)-1L)
    sig<-c(Sigdata[,-1L],recursive=TRUE)
  } # end if 
  pars<-Stdpars<-vector(length=2L*ncomp)
  value<--99.0
  predtval<-vector(length=ntim)
  ### Set default typ (CW).
  typ<-typ[1L]
  ###
  ### Do a type dependent fitting.
  if(typ=="cw") {
    ### For OSL signal of type "CW".
    ### default arguments for differential evolution.
    ### ****************************************************
    args<-list(factor=20L,f=0.5,cr=0.9,maxiter=1000L,tol=0.1)   
    args[names(control.args)]<-control.args
    factor<-args$factor
    f<-args$f
    cr<-args$cr
    maxiter<-args$maxiter
    tol<-args$tol
    ### For argument "factor".
    if(!is.numeric(factor))       stop("Error: control.args: factor should be of type numeric!")
    if(length(factor)!=1L)        stop("Error: control.args: factor should be an one-element vector!")
    if(!is.finite(factor))        stop("Error: control.args: factor should not be a non-finite value!")
    if(abs(factor-round(factor))
       >=.Machine$double.eps^0.5) stop("Error: control.args: factor should be an integer!")
    if(factor<3L)                 stop("Error: control.args: factor is too small!") 
    if(factor>100L)               stop("Error: control.args: factor is too large!")  
    ### For argument "f".           
    if(!is.numeric(f))            stop("Error: control.args: f should be of type numeric!")
    if(length(f)!=1L)             stop("Error: control.args: f should be an one-element vector!")
    if(!is.finite(f))             stop("Error: control.args: f should not be a non-finite value!")
    if(f<=0.0 || f>1.2)           stop("Error: control.args: incorrect f!")
    ### For argument "cr".
    if(!is.numeric(cr))           stop("Error: control.args: cr should be of type numeric!")
    if(length(cr)!=1L)            stop("Error: control.args: cr should be an one-element vector!")
    if(!is.finite(cr))            stop("Error: control.args: cr should not be a non-finite value!")
    if(cr<=0.0 || cr>1.0)         stop("Error: control.args: incorrect cr!")
    ### For argument "maxiter".
    if(!is.numeric(maxiter))      stop("Error: control.args: maxiter should be of type numeric!")
    if(length(maxiter)!=1L)       stop("Error: control.args: maxiter should be an one-element vector!")
    if(!is.finite(maxiter))       stop("Error: control.args: maxiter should not be a non-finite value!")
    if(abs(maxiter-round(maxiter))
       >=.Machine$double.eps^0.5) stop("Error: control.args: maxiter should be an integer!")
    if(maxiter<10L)               stop("Error: control.args: maxiter is too small!")
    if(maxiter>10000L)            stop("Error: control.args: maxiter is too large!") 
    ### For argument "tol".
    if(!is.numeric(tol))          stop("Error: control.args: tol should be of type numeric!")
    if(length(tol)!=1L)           stop("Error: control.args: tol should be an one-element vector!")
    if(!is.finite(tol))           stop("Error: control.args: tol should not be a non-finite value!")
    if(tol<0.0)                   stop("Error: control.args: incorrect tol!")
    ###
    ### ****************************************************
    errorflag<-vector(length=3L)
    ### Call Fortran subroutine decomp().
    res<-.Fortran("decomp", as.integer(ncomp), as.double(tim), as.double(sig),
                  as.integer(ntim), pars=as.double(pars), Stdpars=as.double(Stdpars),
                  value=as.double(value), transf=as.integer(transf), predtval=as.double(predtval),
                  as.integer(factor), as.double(f), as.double(cr), as.integer(maxiter), as.double(tol),
                  errorflag=as.integer(errorflag), package="numOSL")
    ### Error checking.
    if(all(res$pars<0.0))  {
       stop(paste("Error: signal can not be decomposed to",ncomp,"components, or control.args need modification!"))
    } # end if
  } else if(typ=="lm") {
    ### For OSL signal of type "LM".
    errorflag<-0
    ### Call Fortran subroutine fitlm().
    res<-.Fortran("fitlm", as.integer(ncomp), as.double(tim), as.double(sig), as.integer(ntim),
                  pars=as.double(pars), Stdpars=as.double(Stdpars), value=as.double(value), predtval=as.double(predtval),
                  transf=as.integer(transf), errorflag=as.integer(errorflag), package="numOSL")
    ### Error checking.
    if(all(res$pars<0.0)) {
      stop(paste("Error: signal can not be decomposed to",ncomp,"components!"))
    } # end if
  } # end if
  ### print(res$errorflag)
  ###
  ### Decide Photoionisation cross-section (cm2) (may not be ture!).
  h<-6.62606957e-34
  ny<-299792458/(LEDwavelength/10^9L)
  E<-h*ny
  LEDpower<-LEDpower/1000.0
  ###
  ### Reshpae parameters for output.
  pars<-cbind(res$pars[1L:ncomp], res$Stdpars[1L:ncomp], 
              res$pars[-(1L:ncomp)], res$Stdpars[-(1L:ncomp)],
              res$pars[-(1L:ncomp)]/LEDpower*E)
  pars<-pars[order(pars[,3L], decreasing=TRUE), , drop=FALSE]
  colnames(pars)<-c("Ithn", "Std.Ithn", "Lamda", "Std.Lamda", "Pcs")
  rownames(pars)<-paste("Comp.", 1L:ncomp, sep="")
  ###
  ### Reset parameters and errorflag.
  if(typ=="cw") {
    if(all(pars[,c(2L,4L)]<0.0)) {
      pars[,c(2L,4L)]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } else if(typ=="lm") {
    if(all(pars[,c(2L,4L)]<0.0)) {
      pars[,c(2L,4L)]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } # end if
  ###
  ### Calculate signal values for each component.
  if(transf==FALSE) {
    CompSig<-apply(cbind(pars[,1L], pars[,3L]), MARGIN=1L,
             function(x) if(typ=="cw") x[1L]*exp(-x[2L]*tim) else 
             x[1L]*(tim/max(tim))*exp(-x[2L]*tim^2L/2.0/max(tim)))
    SigProp<-(pars[,1L]/pars[,3L])/sum(pars[,1L]/pars[,3L])
  } else if(transf==TRUE) {
    CompSig<-apply(cbind(pars[,1L], pars[,3L]), MARGIN=1L,
             function(x) if(typ=="cw") x[1L]*x[2L]*exp(-x[2L]*tim) else 
             x[1L]*x[2L]*(tim/max(tim))*exp(-x[2L]*tim^2L/2.0/max(tim)))
    SigProp<-pars[,1L]/sum(pars[,1L])
  } # end if
  CompSig<-cbind(res$predtval, CompSig)
  colnames(CompSig)<-c("Fit.Signal", paste("Comp.", 1L:ncomp, sep=""))
  ###
  ### Plot out result out or not ?
  if(plot==TRUE) {
    ### Set layout.
    layout(cbind(c(1L,1L,1L,2L),c(1L,1L,1L,2L)))
    ### Set plot pars.
    par(bg="grey95",
        mar=c(0,5.1,3.1,1.1))
    ### Add a scatter plot (Time .VS. Signal).
    plot(tim, sig, main=samplename, log=ifelse(xlog==TRUE,"x",""), las=0, cex.main=1.5,
         lab=c(7,7,9), ylim=c(-max(sig)*0.01, max(sig)*1.01), cex.lab=1.5, xaxt="n", ylab="Photon Counts",
         xaxs="r", yaxs="i", type="p", pch=21, cex=ifelse(typ=="cw",1.5,1.25), bg="white", col="black") 
    XaxisCentral<-median(axTicks(side=1L))
    ### Set colors
    colors<-c("blue", "red", "green", "deepskyblue", "purple", "orange", "brown")
    ###
    ### Lines Time .VS. Fitted values.
    if(ncol(Sigdata)==2L) {
      ### For signal curve
      x<-seq(min(tim),max(tim),by=(max(tim)-min(tim))/ntim/100L)
    } else {
      ### For signal table (CW).
      x<-seq(min(tim),max(tim),by=(max(tim)-min(tim))/ntim/2L)
    } # end if
    ### If the model has been transfomed?
    if(transf==FALSE) {
      if(typ=="cw") {
        lines(x, eval(parse(text=paste("pars[",1:ncomp,",1]*exp(-pars[",1:ncomp,",3]*x)",collapse="+",sep=""))), 
              lwd=lwd, col="black", lty="solid")
      } else if(typ=="lm") {
        lines(x,eval(parse(text=paste("pars[",1:ncomp,",1]*(x/max(tim))*exp(-pars[",1:ncomp,",3]*x^2/2/max(tim))",collapse="+",sep=""))), 
              lwd=lwd, col="black", lty="solid")
      } # end if
    } else if(transf==TRUE) {
      if(typ=="cw") {
        lines(x, eval(parse(text=paste("pars[",1:ncomp,",1]*pars[",1:ncomp,",3]*exp(-pars[",1:ncomp,",3]*x)",collapse="+",sep=""))), 
              lwd=lwd, col="black", lty="solid")
      } else if(typ=="lm") {
        lines(x,eval(parse(text=paste("pars[",1:ncomp,",1]*pars[",1:ncomp,",3]*(x/max(tim))*exp(-pars[",1:ncomp,",3]*x^2/2/max(tim))",collapse="+",sep=""))), 
              lwd=lwd, col="black", lty="solid")
      } # end if
    } # end if
    ### 
    ### Lines Time .VS. Component signal (1 to ncomp).
    for(i in 1L:ncomp) {
      if(transf==FALSE) {
        if(typ=="cw") {
          curve(pars[i,1L]*exp(-pars[i,3L]*x), lwd=lwd, col=colors[i], lty= "solid", add=TRUE)
        } else if(typ=="lm") {
          curve(pars[i,1L]*(x/max(tim))*exp(-pars[i,3L]*x^2L/2.0/max(tim)), lwd=lwd, col=colors[i], lty= "solid", add=TRUE)
        } # end if
      } else if(transf==TRUE) {
        if(typ=="cw") {
          curve(pars[i,1L]*pars[i,3L]*exp(-pars[i,3L]*x), lwd=lwd, col=colors[i], lty= "solid", add=TRUE)
        } else if(typ=="lm") {
          curve(pars[i,1L]*pars[i,3L]*(x/max(tim))*exp(-pars[i,3L]*x^2L/2.0/max(tim)), lwd=lwd, col=colors[i], lty= "solid", add=TRUE)
        } # end if
      } # end if
    } # end for
    ###
    ### Add a legend to the plot.
    legend(ifelse(typ=="cw", "topright", ifelse(tim[which.max(sig)]>XaxisCentral, "topleft","topright")), 
           legend=c("Fitted.Curve", paste("Comp.", 1L:ncomp," (",round(SigProp*100,2L),"%)", sep="")),
           col=c("black", colors[1L:ncomp]), pch=c(21,rep(NA,ncomp)), lty="solid",
           yjust=2, ncol=1, cex=1.5, bty="o", lwd=lwd, pt.bg="white")
  ### Add grid and box to the plot.
  grid(equilogs=FALSE)
  box(lwd=2L)
  ###
  ### Add residuals as a new plot.
  par(mar=c(5.1,5.1,0,1.1))
  plot(tim, sig-CompSig[,1L], log=ifelse(xlog==TRUE,"x",""), las=0, lab=c(7,7,9), 
       ylim=c(min(sig-CompSig[,1L])*1.01, max(sig-CompSig[,1L])*1.01),xlab="Stimulated Time (s)", cex.lab=1.5,
       ylab="Residuals", xaxs="r", yaxs="i", type="o", pch=21, cex=ifelse(typ=="cw",0.75,0.5), bg="black", col="gray") 
  abline(h=0)
  box(lwd=2L)
  ### Reset pars brefore leaving.
  par(bg="transparent",
      mar=c(5,4,4,2)+0.1)
  layout(1L)
  } # end if
  ###
  ### Prepare results for output.
  out<-list("Comp.Signal"=CompSig, "pars"=pars,
            "value"=res$value, "errorflag"=errorflag)
  ###
  ### If wirte fitted signal values to a file or not ?
  if(!is.null(outfile)) {
    write.csv(CompSig, file=paste(outfile, ".csv",sep=""))
  } # end if
  ###
  ### Set output to be invisible!
  invisible(out)
} # end function decomp
############################################### END FUNCTION decomp ###################################################
###
################################################ FUNCTION decompc #####################################################
### *******************************************************************************************************************
### Function decompc() is used to decompose OSL signal curve (with a constant subtracted),
### either type "CW" or "LM" can be analyzed.
###
###      Author: Peng Jun, 2013.12.16; revised in 2014.01.04.
###
###  References: Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
###              applications. Geochronometria 13, 135–141.
###
###              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
###              OSL decay curves. Radiation Measurements 41, 886-891.
###
###              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
###              luminescence components in different quartz samples: implications for dose measurement. 
###              Radiation Measurements, 37 (4-5), pp. 441-449.
###
### Mandatory arguments---
###     Sigdata: a dataframe, OSL decay curve (two columns).
### *******************************************************************************************************************
decompc<-
function(Sigdata, ncomp=2, typ=c("cw","lm"),
         control.args=list(), LEDpower=60, LEDwavelength=470,
         plot=TRUE, xlog=TRUE, lwd=3, samplename=NULL, outfile=NULL) {
  UseMethod("decompc")
} ###
### default method for function decompc().
decompc.default<-
function(Sigdata, ncomp=2, typ=c("cw","lm"),
         control.args=list(), LEDpower=60, LEDwavelength=470,
         plot=TRUE, xlog=TRUE, lwd=3, samplename=NULL, outfile=NULL) {
  ### stop if not
  ###
  ### For argument "Sigdata".
  if(!is.data.frame(Sigdata))             stop("Error: Sigdata should be of type data.frame!")
  if(ncol(Sigdata)<2L)                    stop("Error: Sigdata should contain at least two columns!")
  if(ncol(Sigdata)>2L && typ[1L]!="cw")   stop("Error: only signal table of type CW-OSL can be analyzed currently!")
  if(!is.numeric(as.matrix(Sigdata)))           
                                          stop("Error: all elements in Sigdata should be of type numeric!")
  if(any(!is.finite(unlist(unclass(Sigdata)))))   
                                          stop("Error: Sigdata should not contain non-finite value!")
  if(any(Sigdata[,1L,drop=TRUE]<0))       stop("Error: the first column of Sigdata [time] should larger than zero!")
  if(any(Sigdata[,-1L]<0))                stop("Error: all signal values should larger than zero!")
  ###
  ### For argument "ncomp".
  if(!is.numeric(ncomp))                  stop("Error: ncomp should be of type numeric!")
  if(length(ncomp)!=1L)                   stop("Error: ncomp should be an one-element vector!")
  if(!is.finite(ncomp))                   stop("Error: ncomp should not be a non-finite value!")
  if(!ncomp %in% (1L:7L) )                stop("Error: ncomp should be a integer in the space [1,7]!")
  ###
  ### For argument "typ".
  if(!is.character(typ))                  stop("Error: typ should be of type character!")
  if(length(typ)==1L) {
    if(!typ %in% c("cw","lm"))            stop("Error: typ should be either 'cw' or 'lm'!")
    if(ncol(Sigdata)==2L) {
      if(typ=="cw" && which.max(Sigdata[,2L,drop=TRUE])>5L)      
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
      if(typ=="lm" && which.max(Sigdata[,2L,drop=TRUE])<=5L)     
                                          stop("Error: incorrect typ, signal may not be of type LM-OSL!")
    } else if(ncol(Sigdata)>2L) {
      if(which.max(apply(Sigdata[,-1L],MARGIN=1L,mean))>5L)        
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
    } # end if
  } else if(length(typ)==2L) {
    if(!all(typ %in% c("cw","lm")))       stop("Error: incorrect typ, please choose one between 'cw' and 'lm'!")
    if(ncol(Sigdata)==2L) {
      if(typ[1L]=="cw" && which.max(Sigdata[,2L,drop=TRUE])>5L)      
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
      if(typ[1L]=="lm" && which.max(Sigdata[,2L,drop=TRUE])<=5L)     
                                          stop("Error: incorrect typ, signal may not be of type LM-OSL!")
    } else if(ncol(Sigdata)>2L) {
      if(which.max(apply(Sigdata[,-1L],MARGIN=1L,mean))>5L)        
                                          stop("Error: incorrect typ, signal may not be of type CW-OSL!")
    } # end if
  } else {                                stop("Error: typ should contain no more than two elements!")
  } # end if
  ###
  ### For argument "control.args".
  if(class(control.args)!="list")         stop("Error: control.args should be a list!")
  if(length(control.args)>5L)             stop("Error: number of parameters in control.args should not exceed 5!")
  if(!all(names(control.args) %in% list("factor","f","cr","maxiter","tol")))              
                                          stop("Error: incorrect parameter in control.args!")
  ###
  ### For argument "LEDpower".
  if(!is.numeric(LEDpower))               stop("Error: LEDpower should be of type numeric!")
  if(length(LEDpower)!=1L)                stop("Error: LEDpower should be an one-element vector!")
  if(!is.finite(LEDpower))                stop("Error: LEDpower should not be a non-finite value!")
  if(LEDpower<=0)                         stop("Error: LEDpower should larger than zero!")
  ###
  ### For argument "LEDwavelength".
  if(!is.numeric(LEDwavelength))          stop("Error: LEDwavelength should be of type numeric!")
  if(length(LEDwavelength)!=1L)           stop("Error: LEDwavelength should be an one-element vector!")
  if(!is.finite(LEDwavelength))           stop("Error: LEDwavelength should not be a non-finite value!")
  if(LEDwavelength<=0)                    stop("Error: LEDwavelength should larger than zero!")
  ###
  ### For argument "plot".
  if(!is.logical(plot))                   stop("Error: plot should be TRUE or FALSE!")
  if(length(plot)!=1L)                    stop("Error: plot should be an one-element vector!")
  ###
  ### For argument "xlog".
  if(!is.logical(xlog))                   stop("Error: xlog should be TRUE or FALSE!")
  if(length(xlog)!=1L)                    stop("Error: xlog should be an one-element vector!")
  ###
  ### For argument "lwd".
  if(!is.numeric(lwd))                    stop("Error: lwd should be of type numeric!")
  if(length(lwd)!=1L)                     stop("Error: lwd should be an one-element vector!")
  if(!is.finite(lwd))                     stop("Error: lwd should not be a non-finite value!")
  if(lwd<=0)                              stop("Error: lwd should larger than zero!")
  ###
  ### For argument "samplename".
  if(!is.null(samplename) && !is.character(samplename))           
                                          stop("Error: samplename should be NULL or of type character!")
  if(!is.null(samplename)) {
    if(length(samplename)!=1L)            stop("Error: sample name should be an one-element vector!")
  } # end if
  ###
  ### For argument "outfile".
  if(!is.null(outfile) && !is.character(outfile))             
                                          stop("Error: outfile should be NULL or of type character!")
  if(!is.null(outfile)) {
    if(ncol(Sigdata)>2L)                  stop("Error: signal table is not allowed to be written out!")
    if(length(outfile)!=1L)               stop("Error: outfile should be an one-element vector!")
  } # end if
  ###
  ###
  ### Public parameters for Fortran  
  ### subroutine decomp() and fitlm().
  if(ncol(Sigdata)==2L) {
    ### Decompose signal curve.
    ntim<-nrow(Sigdata)
    tim<-drop(Sigdata[,1L])
    sig<-drop(Sigdata[,2L])
  } else if(ncol(Sigdata)>2L) {
    ### Decompose signal table (CW).
    ntim<-nrow(Sigdata)*(ncol(Sigdata)-1L)
    Sigdata<-Sigdata[,c(1L,order(Sigdata[1L,-1L],decreasing=TRUE)+1L)]
    tim<-rep(Sigdata[,1L,drop=TRUE],ncol(Sigdata)-1L)
    sig<-c(Sigdata[,-1L],recursive=TRUE)
  } # end if 
  pars<-Stdpars<-vector(length=2L*ncomp+1L)
  value<--99.0
  predtval<-vector(length=ntim)
  ### Set default typ (CW)
  typ<-typ[1L]
  ###
  ### Do a type dependent fitting.
  if(typ=="cw") {
    ### For OSL signal of type "CW".
    ### default arguments for differential evolution.
    ### ****************************************************
    args<-list(factor=20L,f=0.5,cr=0.9,maxiter=1000L,tol=0.1)   
    args[names(control.args)]<-control.args
    factor<-args$factor
    f<-args$f
    cr<-args$cr
    maxiter<-args$maxiter
    tol<-args$tol
    ### For argument "factor".
    if(!is.numeric(factor))       stop("Error: control.args: factor should be of type numeric!")
    if(length(factor)!=1L)        stop("Error: control.args: factor should be an one-element vector!")
    if(!is.finite(factor))        stop("Error: control.args: factor should not be a non-finite value!")
    if(abs(factor-round(factor))
       >=.Machine$double.eps^0.5) stop("Error: control.args: factor should be an integer!")
    if(factor<3L)                 stop("Error: control.args: factor is too small!") 
    if(factor>100L)               stop("Error: control.args: factor is too large!")  
    ### For argument "f".          
    if(!is.numeric(f))            stop("Error: control.args: f should be of type numeric!")
    if(length(f)!=1L)             stop("Error: control.args: f should be an one-element vector!")
    if(!is.finite(f))             stop("Error: control.args: f should not be a non-finite value!")
    if(f<=0.0 || f>1.2)           stop("Error: control.args: incorrect f!")
    ### For argument "cr".
    if(!is.numeric(cr))           stop("Error: control.args: cr should be of type numeric!")
    if(length(cr)!=1L)            stop("Error: control.args: cr should be an one-element vector!")
    if(!is.finite(cr))            stop("Error: control.args: cr should not be a non-finite value!")
    if(cr<=0.0 || cr>1.0)         stop("Error: control.args: incorrect cr!")
    ### For argument "maxiter".
    if(!is.numeric(maxiter))      stop("Error: control.args: maxiter should be of type numeric!")
    if(length(maxiter)!=1L)       stop("Error: control.args: maxiter should be an one-element vector!")
    if(!is.finite(maxiter))       stop("Error: control.args: maxiter should not be a non-finite value!")
    if(abs(maxiter-round(maxiter))
       >=.Machine$double.eps^0.5) stop("Error: control.args: maxiter should be an integer!")
    if(maxiter<10L)               stop("Error: control.args: maxiter is too small!")
    if(maxiter>10000L)            stop("Error: control.args: maxiter is too large!") 
    ### For argument "tol".
    if(!is.numeric(tol))          stop("Error: control.args: tol should be of type numeric!")
    if(length(tol)!=1L)           stop("Error: control.args: tol should be an one-element vector!")
    if(!is.finite(tol))           stop("Error: control.args: tol should not be a non-finite value!")
    if(tol<0.0)                   stop("Error: control.args: incorrect tol!")
    ###
    ### ****************************************************
    errorflag<-vector(length=3L)
    ### Call Fortran subroutine decomp_C().
    res<-.Fortran("decomp_C", as.integer(ncomp), as.double(tim), as.double(sig),
                  as.integer(ntim), pars=as.double(pars), Stdpars=as.double(Stdpars),
                  value=as.double(value), predtval=as.double(predtval), as.integer(factor),
                  as.double(f), as.double(cr), as.integer(maxiter), as.double(tol),
                  errorflag=as.integer(errorflag), package="numOSL")
    ### Error checking.
    if(all(res$pars<0.0))  {
       stop(paste("Error: signal can not be decomposed to",ncomp,"components, or control.args need modification!"))
    } # end if
  } else if(typ=="lm") {
    ### For OSL signal of type "LM".
    errorflag<-0
    ### Call Fortran subroutine fitlm().
    res<-.Fortran("fitlm_C", as.integer(ncomp), as.double(tim), as.double(sig), as.integer(ntim),
                  pars=as.double(pars), Stdpars=as.double(Stdpars), value=as.double(value), 
                  predtval=as.double(predtval), errorflag=as.integer(errorflag), package="numOSL")
    ### Error checking
    if(all(res$pars<0.0)) {
      stop(paste("Error: signal can not be decomposed to",ncomp,"components!"))
    } # end if
  } # end if
  ### print(res$errorflag)
  ###
  ### Decide Photoionisation cross-section (cm2) (may not be ture!).
  h<-6.62606957e-34
  ny<-299792458/(LEDwavelength/10^9L)
  E<-h*ny
  LEDpower<-LEDpower/1000.0
  ###
  ### Reshpae parameters for output.
  pars<-cbind(res$pars[1L:ncomp], res$Stdpars[1L:ncomp], 
              res$pars[(ncomp+1L):(2L*ncomp)], res$Stdpars[(ncomp+1L):(2L*ncomp)],
              res$pars[(ncomp+1L):(2L*ncomp)]/LEDpower*E)
  pars<-pars[order(pars[,3L], decreasing=TRUE), , drop=FALSE]
  colnames(pars)<-c("Ithn", "Std.Ithn", "Lamda", "Std.Lamda", "Pcs")
  rownames(pars)<-paste("Comp.", 1L:ncomp, sep="")
  constant<-c("Constant"=res$pars[2L*ncomp+1L],"Std.Constant"=res$Stdpars[2L*ncomp+1L])
  ###
  ### Reset parameters and errorflag.
  if(typ=="cw") {
    if(all(pars[,c(2L,4L)]<0.0)) {
      pars[,c(2L,4L)]<-NA
      constant[2L]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } else if(typ=="lm") {
    if(all(pars[,c(2L,4L)]<0.0)) {
      pars[,c(2L,4L)]<-NA
      constant[2L]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } # end if
  ###
  ### Calculate signal values for each component
  CompSig<-apply(cbind(pars[,1L], pars[,3L]), MARGIN=1L,
                 function(x) if(typ=="cw") x[1L]*x[2L]*exp(-x[2L]*tim) else 
                 x[1L]*x[2L]*(tim/max(tim))*exp(-x[2L]*tim^2L/2.0/max(tim)))
  SigProp<-pars[,1L]/sum(pars[,1L])
  ###
  if(typ=="cw") {
    CompSig<-cbind(res$predtval, CompSig, constant[1L])
  } else if (typ=="lm") {
    CompSig<-cbind(res$predtval, CompSig, constant[1L]*tim/max(tim))
  } # end if
  colnames(CompSig)<-c("Fit.Signal", paste("Comp.", 1L:ncomp, sep=""), "Constant")
  ###
  ### Plot out result out or not ?
  if(plot==TRUE) {
    ### Set layout
    layout(cbind(c(1L,1L,1L,2L),c(1L,1L,1L,2L)))
    ### Set plot pars.
    par(bg="grey95",
        mar=c(0,5.1,3.1,1.1))
    ### Add a scatter plot (Time .VS. Signal).
    plot(tim, sig, main=samplename, log=ifelse(xlog==TRUE,"x",""), las=0, cex.main=1.5,
         lab=c(7,7,9), ylim=c(-max(sig)*0.01, max(sig)*1.01), cex.lab=1.5, xaxt="n", ylab="Photon Counts",
         xaxs="r", yaxs="i", type="p", pch=21, cex=ifelse(typ=="cw",1.5,1.25), bg="white", col="black") 
    XaxisCentral<-median(axTicks(side=1L))
    ### Set colors.
    colors<-c("blue", "red", "green", "deepskyblue", "purple", "orange", "brown")
    ###
    ### Lines Time .VS. Fitted values.
    if(ncol(Sigdata)==2L) {
      ### For signal curve
      x<-seq(min(tim),max(tim),by=(max(tim)-min(tim))/ntim/100L)
    } else {
      ### For signal table (CW).
      x<-seq(min(tim),max(tim),by=(max(tim)-min(tim))/ntim/2L)
    } # end if
    ### 
    ###
    if(typ=="cw") {
      lines(x, eval(parse(text=paste("pars[",1:ncomp,",1]*pars[",1:ncomp,",3]*exp(-pars[",1:ncomp,",3]*x)",collapse="+",sep="")))+constant[1L], 
            lwd=lwd, col="black", lty="solid")
    } else if(typ=="lm") {
      lines(x,eval(parse(text=paste("pars[",1:ncomp,",1]*pars[",1:ncomp,",3]*(x/max(tim))*exp(-pars[",1:ncomp,",3]*x^2/2/max(tim))",collapse="+",sep="")))+constant[1L]*x/max(tim), 
            lwd=lwd, col="black", lty="solid")
    } # end if
    ### 
    ### Add constant as back counts.
    if(typ=="cw") {
      points(tim, rep(constant[1L],ntim), type="l", lty="dashed", lwd=lwd)
    } else if(typ=="lm") {
      points(tim, constant[1L]*tim/max(tim), type="l", lty="dashed", lwd=lwd)
    } # end if
    ###
    ### Lines Time .VS. Component signal (1 to ncomp).
    for(i in 1L:ncomp) {
      if(typ=="cw") {
        curve(pars[i,1L]*pars[i,3L]*exp(-pars[i,3L]*x), lwd=lwd, col=colors[i], lty="solid", add=TRUE)
      } else if(typ=="lm") {
        curve(pars[i,1L]*pars[i,3L]*(x/max(tim))*exp(-pars[i,3L]*x^2L/2.0/max(tim)), lwd=lwd, col=colors[i], lty="solid", add=TRUE)
      } # end if
    } # end if
    ###
    ### Add a legend to the plot.
    legend(ifelse(typ=="cw", "topright", ifelse(tim[which.max(sig)]>XaxisCentral, "topleft","topright")), 
           legend=c("Fitted.Curve", paste("Comp.", 1L:ncomp," (",round(SigProp*100,2L),"%)", sep=""),"Constant"),
           col=c("black", colors[1L:ncomp], "black"), pch=c(21,rep(NA,ncomp+1L)), lty=c(rep("solid",ncomp+1L),"dashed"),
           yjust=2, ncol=1, cex=1.5, bty="o", lwd=lwd, pt.bg="white")
  ### Add grid and box to the plot.
  grid(equilogs=FALSE)
  box(lwd=2L)
  ###
  ### Add residuals as a new plot.
  par(mar=c(5.1,5.1,0,1.1))
  plot(tim, sig-CompSig[,1L], log=ifelse(xlog==TRUE,"x",""), las=0, lab=c(7,7,9), 
       ylim=c(min(sig-CompSig[,1L])*1.01, max(sig-CompSig[,1L])*1.01),xlab="Stimulated Time (s)", cex.lab=1.5,
       ylab="Residuals", xaxs="r", yaxs="i", type="o", pch=21, cex=ifelse(typ=="cw",0.75,0.5), bg="black", col="gray") 
  abline(h=0)
  box(lwd=2L)
  ### Reset pars brefore leaving.
  par(bg="transparent",
      mar=c(5,4,4,2)+0.1)
  layout(1L)
  } # end if
  ###
  ### Prepare results for output.
  out<-list("Comp.Signal"=CompSig, "pars"=pars, "constant"=constant,
            "value"=res$value, "errorflag"=errorflag)
  ###
  ### If wirte fitted signal values to a file or not ?
  if(!is.null(outfile)) {
    write.csv(CompSig, file=paste(outfile, ".csv",sep=""))
  } # end if
  ###
  ### Set output to be invisible!
  invisible(out)
} # end function decompc.
############################################## END FUNCTION decompc ###################################################
###
################################################## FUNCTION fastED ####################################################
### *******************************************************************************************************************
### Function fastED() is used to calculate a fast-component equivalent dose value.
### 
###     Author: Peng Jun, 2013.11.21; revised in 2013.12.17; revised in 2014.01.04; revised in 2014.04.06.
###
### References: Murray, A.S., Wintle, A.G., 2000. Luminescence dating of quartz using improved single-aliquot
###             regenerative-dose protocol. Radiation Measurements, 32, pp.57-73.
###
###             Li, S.H., Li, B., 2006. Dose measurement using the fast component of LM-OSL signals from quartz.
###             Radiation Measurements, 41, pp.534-541.
###
### Mandatory arguments---
###    Sigdata: a data.frame, OSL decay curve (more than five columns).  
###     Redose: a vector, regenerative doses used for dose response curve construction.    
### *******************************************************************************************************************
fastED<-
function(Sigdata, Redose, ncomp=2, constant=TRUE,
         control.args=list(), typ="cw", nstart=100, upb=1,
         ErrorMethod=c("mc","sp"), model=NULL, origin=NULL) {
  UseMethod("fastED")
} #
### Set default method for function fastED().
fastED.default<-
function(Sigdata, Redose, ncomp=2, constant=TRUE,
         control.args=list(), typ="cw", nstart=100, upb=1,
         ErrorMethod=c("mc","sp"), model=NULL, origin=NULL) {
  ### Stop if not.
  ###
  ### For argument "Sigdata".
  if(!is.data.frame(Sigdata))                 stop("Error: Sigdata should be of type data.frame!")
  if(ncol(Sigdata)<5L)                        stop("Error: Sigdata should contain at least 5 columns!")
  if(ncol(Sigdata) %% 2L != 1L)               stop("Error: Sigdata should contain odd number of columns!")
  if(!is.numeric(as.matrix(Sigdata)))         stop("Error: all value in Sigdata should be of type numeric!")
  if(any(!is.finite(unlist(unclass(Sigdata)))))   
                                              stop("Error: Sigdata should not contain non-finite value!")
  if(any(Sigdata[,1L,drop=TRUE]<0))           stop("Error: all time values should larger than zero!")
  if(any(abs(Sigdata[,-1L]-round(Sigdata[,-1L]))>=.Machine$double.eps^0.5))
                                              stop("Error: all photon counts (signal) should be of type integer!")
  if(any(Sigdata[,-1L]<0))                    stop("Error: all photon counts (signal) should not below zero!")
  ###
  ### For argument "Redose".
  if(!is.numeric(Redose))                     stop("Error: Redose should be a numeric vector!")
  if(any(!is.finite(Redose)))                 stop("Error: Redose should not contain non-finite value!")
  if(any(Redose<0))                           stop("Error: all Redose values should be not smaller than 0!")
  if(!length(which(abs(Redose)<=.Machine$double.eps^0.5)) %in% c(0L,1L))         
                                              stop("Error: only one zero-Redose is allowed!")
  if(length(Redose)!=(ncol(Sigdata)-3L)/2L)   stop("Error: mismatch between number of Redose values and number of decay curves!")
  if(!(length(Redose)-length(as.numeric(levels(factor(Redose))))) %in% c(0L,1L))
                                              stop("Error: only one repeated-Redose is allowed!")
  ###
  ### For argument "ncomp".
  if(!is.numeric(ncomp))                      stop("Error: ncomp should be a numeric vector!")
  if(length(ncomp)!=1L)                       stop("Error: ncomp should be an one-element vector!")
  if(!is.finite(ncomp))                       stop("Error: ncomp should not be a non-finite value!")
  if(!ncomp %in% (2L:4L))                     stop("Error: ncomp should be an integer in the space [2,4]!")
  ###
  ### For argument "constant".
  if(!is.logical(constant))                   stop("Error: constant should be of type logical!")
  if(length(constant)!=1L)                    stop("Error: constant should be an one-element vector!")
  ###
  ### For argument "control.args".
  if(class(control.args)!="list")             stop("Error: control.args should be a list!")
  if(length(control.args)>5L)                 stop("Error: number of parameters in control.args should not exceed 5!")
  if(!all(names(control.args) %in% 
     list("factor","f","cr",
          "maxiter","tol")))                  stop("Error: incorrect names of parameters in control.args!")
  ###
  ### For argument "typ".
  if(!is.character(typ))                      stop("Error: typ should be of type character!")
  if(length(typ)!=1L)                         stop("Error: typ should be an one-element vector!")
  if(typ!="cw")                               stop("Error: only CW-OSL can be analyzed currently!")
  ###
  ### For argument "nstart".
  if(!is.numeric(nstart))                     stop("Error: nstart should be of type numeric!")
  if(length(nstart)!=1L)                      stop("Error: nstart should be an one-element vector!")
  if(!is.finite(nstart))                      stop("Error: nstart should not be a non-finite value!")
  if(abs(nstart-round(nstart))
     >=.Machine$double.eps^0.5)               stop("Error: nstart should be an integer!")
  if(nstart<5L)                               stop("Error: nstart should not be less than 5!")
  if(nstart>5000L)                            stop("Error: nstart should not exceed 5000!") 
  ###
  ### For argument "upb".      
  if(!is.numeric(upb))                        stop("Error: upb should be of type numeric!")
  if(length(upb)!=1L)                         stop("Error: upb should be an one-element vector!")
  if(!is.finite(upb))                         stop("Error: upb should not be an non-finite value!")
  if(upb<=0)                                  stop("Error: upb should larger than 0!")
  if(upb>100)                                 stop("Error: upb should not exceed 100!")
  ###
  ### For argument "ErrorMethod".
  if(!is.character(ErrorMethod))              stop("Error: ErrorMethod should be of type character!")
  if(length(ErrorMethod)>2L)                  stop("Error: ErrorMethod should contain at most two elements!")
  if(!all(ErrorMethod %in% c("sp","mc")))     stop("Error: incorrect ErrorMethod, only 'sp' and 'mc' are available!")
  ###
  ### For argument "model".
  if(!is.null(model) && !is.character(model)) stop("Error: model should be NULL or of type character!")
  if(!is.null(model)) {
    if(length(model)!=1L)                     stop("Error: model should be an one-element vector!")
    if(!model%in%c("line","exp","line+exp"))  stop("Error: incorrect name of model!")
  } # end if
  ### For argument "origin".
  if(!is.logical(origin) && !is.null(origin)) stop("Error: origin should be either NULL or a logical variable!")
  if(!is.null(origin)) {
    if(length(origin)!=1L)                    stop("Error: origin should be an one-element vector!")
  } # end if
  ###
  ###
  ##################################
  ### Set the function for plotting.
  PlotRegion<-function(tim,sig,pars,cvalue,samplename) {
    ###
    par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
    ###
    ### Add a scatter plot (Time .VS. Signal).
    plot(tim, sig, main=samplename, log="x", las=0, cex.main=1.25, lab=c(7,7,9),
         ylim=c(-max(sig)*0.01, max(sig)*1.01), cex.lab=1, xlab="Stimulated time (s)", 
         ylab="Photon Counts", xaxs="r", yaxs="i", type="p", pch=21, cex=1.5, bg="white", col="black") 
    ###
    ### Set colors
    colors<-c("blue", "red", "green", "deepskyblue")
    ###
    ### Draw lines with Time .VS. Fitted values.
    x<-seq(min(tim), max(tim), by=(max(tim)-min(tim))/length(tim)/100L)
    lines(x, eval(parse(text=paste("pars[",1L:nrow(pars),",1L]*pars[",1L:nrow(pars),",3L]*exp(-pars[",1L:nrow(pars),",3L]*x)",collapse="+",sep="")))+cvalue, 
          lwd=3, col="black", lty="solid")
    ### 
    ### Draw lines with Time .VS. Component signal (from 1 to ncomp).
    for(i in 1L:nrow(pars)) {
      curve(pars[i,1L]*pars[i,3L]*exp(-pars[i,3L]*x), lwd=3, col=colors[i], lty="solid", add=TRUE)
    } # end for
    ###
    ### Add a legend.
    if(cvalue>0.0) {
      abline(h=cvalue, lty="dashed", lwd=3)
      legend("topright",legend=c("Fitted.Curve", paste("Comp.", 1L:nrow(pars), sep=""),"Constant"), col=c("black", colors[1L:nrow(pars)], "black"), 
             pch=c(21,rep(NA,nrow(pars)+1L)), lty=c(rep("solid",nrow(pars)+1L),"dashed"), yjust=2, ncol=1L, cex=par("cex"), bty="o", lwd=3, pt.bg="white")
    } else {
      legend("topright",legend=c("Fitted.Curve", paste("Comp.", 1L:nrow(pars), sep="")), col=c("black", colors[1L:nrow(pars)]), 
             pch=c(21,rep(NA,nrow(pars))), lty="solid", yjust=2, ncol=1L, cex=par("cex"), bty="o", lwd=3, pt.bg="white")
    } # end if
    ###
    ### Add grid and box to the plot.
    grid(equilogs=FALSE)
    box(lwd=2L)
    ### Reset plot pars before leaving.
    par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
    ###
  } # end function PlotRegion.
  ###########################
  ###
  ###
  ### Set dimensions of ararys.
  x<-ncol(Sigdata)-1L
  fLtx<-matrix(nrow=x,ncol=2L)
  pars<-vector("list",length=x)
  ###
  ### Fit natural decay signal.
  if(constant==FALSE) {
      res<-decomp(Sigdata[,c(1L,2L)], ncomp=ncomp, typ=typ, control.args=control.args, transf=TRUE, plot=FALSE)
  } else {
      res<-decompc(Sigdata[,c(1L,2L)], ncomp=ncomp, typ=typ, control.args=control.args, plot=FALSE)
  } # end if
  ### Stop if the natural decay curve can not be decomposed.
  if(res$errorflag!=0) {
      stop("Error: standard errors of parameters of natural decay curve cannot be estimated!")
  } # end if
  ###
  ### Set the plot regions.
  par(mfrow=c(ifelse(x%%4L==0L,x/4L+1L,x%/%4L+1L),4L))
  ### Store fitted pars of natural decay curve.
  pars[[1L]]<-res$pars[,-5L,drop=FALSE]
  ### Store natural Lx (the fastest component) and its std.error.
  fLtx[1L,]<-res$pars[which.max(res$pars[,3L,drop=TRUE]),c(1L,2L),drop=TRUE]
  ### Plot natural decay signal.
  PlotRegion(tim=Sigdata[,1L,drop=TRUE], sig=Sigdata[,2L,drop=TRUE], pars=res$pars, 
             cvalue=ifelse(constant==FALSE,0.0,res$constant[1L]), samplename="Natural")
  ###
  ###
  ### Analyze regenerative decay curves.
  for(i in 3L:ncol(Sigdata)) {
      if(constant==FALSE) {
          res<-try(decomp(Sigdata[,c(1L,i)], ncomp=ncomp, typ=typ, control.args=control.args, transf=TRUE, plot=FALSE), silent=TRUE)
      } else {
          res<-try(decompc(Sigdata[,c(1L,i)], ncomp=ncomp, typ=typ, control.args=control.args, plot=FALSE), silent=TRUE)
      } # end if
    ###
    ### Set name for each decay curve.
    samplename<-ifelse(i==3L,"Test[Natural]",ifelse(i%%2L==0L,paste("Redose.",i/2L-1L,sep=""), paste("Test[Redose.",(i-1L)/2L-1L,"]",sep="")))
    ### Error checking, and plot the results.
    if(class(res)!="try-error" && res$errorflag==0) {
      ### Store pars of Lx or Tx and their std.error.
      pars[[i-1L]]<-res$pars[,-5L,drop=FALSE]
      ### Store Lx or Tx (the fastest component) and their std.error.
      fLtx[i-1L,]<-res$pars[which.max(res$pars[,3L,drop=TRUE]),c(1L,2L),drop=TRUE]
      ### Plot the ith decay curve.
      PlotRegion(tim=Sigdata[,1L,drop=TRUE], sig=Sigdata[,i,drop=TRUE], pars=res$pars, 
                 cvalue=ifelse(constant==FALSE,0.0,res$constant[1L]), samplename=samplename)
    } else {
      cat(paste("Note: fail in decomposing the ",i-1L,"th decay curve!\n", sep=""))
      par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
      plot(x=Sigdata[,1L,drop=TRUE], y=Sigdata[,i,drop=TRUE], main=samplename, log="x", las=0, cex.main=1.25, lab=c(7,7,9),
           ylim=c(-max(Sigdata[,i,drop=TRUE])*0.01, max(Sigdata[,i,drop=TRUE])*1.2), cex.lab=1, xlab="Stimulated time (s)",
           ylab="Photon Counts", xaxs="r", yaxs="i", type="p", pch=21, cex=1.5, bg="white", col="black")
      grid(equilogs=FALSE)
      box(lwd=2L)
      par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
    } # end if
  } # end for
  ###
  ###
  ### Set decayRateMat (a matrix that consists of a series of decay rates).
  ### (remove decay rates of 0-regenerative dose decay curve from decayRateMat if possible).
  if(length(which(abs(Redose)<=.Machine$double.eps^0.5))==1L)  {
      zeroDoseindex<-2L*which(abs(Redose)<=.Machine$double.eps^0.5)+1L
      if(zeroDoseindex %in% which(sapply(pars,function(x) is.null(x))==FALSE) ) {
          decayRateMat<-matrix(unlist(sapply(pars[-zeroDoseindex],function(x) x[,3L,drop=TRUE])), ncol=ncomp, byrow=TRUE)
      } else {
          decayRateMat<-matrix(unlist(sapply(pars,function(x) x[,3L,drop=TRUE])), ncol=ncomp, byrow=TRUE)
      } # end if
  } else {
      decayRateMat<-matrix(unlist(sapply(pars,function(x) x[,3L,drop=TRUE])), ncol=ncomp, byrow=TRUE)
  } # end if
  ### 
  ###
  ### Set the fast-component decay rates, and calculate an uniformity index.
  fastRate<-c("mean"=mean(decayRateMat[,1L,drop=TRUE]),"Std.mean"=sd(decayRateMat[,1L,drop=TRUE])/sqrt(nrow(decayRateMat)))
  WithinVar<-sum(apply(decayRateMat, MARGIN=1L, function(x) sum((x-mean(x))^2L)) )/nrow(decayRateMat)/(ncol(decayRateMat)-1L)
  BetweenVar<-sum(apply(decayRateMat, MARGIN=1L, function(x,y) (mean(x)-mean(y))^2L, decayRateMat) )/(nrow(decayRateMat)-1L)
  uniformity<-WithinVar/BetweenVar
  ###
  ###
  ### Reset Curvedata and set standardised OSLs (for constructing the growth curve).
  fLtx[,2L]<-fLtx[,2L,drop=TRUE]/fLtx[,1L,drop=TRUE]
  fLxTx<-fLtx[(1L:nrow(fLtx))%%2L==1L,1L,drop=TRUE]/fLtx[(1L:nrow(fLtx))%%2L==0L,1L,drop=TRUE]
  sfLxTx<-fLxTx*sqrt( (fLtx[(1L:nrow(fLtx))%%2L==1L,2L,drop=TRUE])^2L + (fLtx[(1L:nrow(fLtx))%%2L==0L,2L,drop=TRUE])^2L )
  Curvedata<-data.frame("Redose"=Redose,"OSL"=fLxTx[-1L],"std.OSL"=sfLxTx[-1L])
  ### Remove NAs from Curvedata (for constructing the growth curve).
  Curvedata<-Curvedata[complete.cases(Curvedata),,drop=FALSE]
  ###
  ###
  ### Set the natrual standardlized OSL.
  NatureLxTx<-c(fLxTx[1L],sfLxTx[1L]) 
  if(any(!is.finite(NatureLxTx))) { 
    stop("Error: fail in calculating standardised natural OSL!")
  } # end if
  ###
  ###
  ### Calculate the recycling ratio with the repeated regenerative dose (remove it from Curvedata if possible).
  if(length(Curvedata[,1L,drop=TRUE])-length(as.numeric(levels(factor(Curvedata[,1L,drop=TRUE]))))==1L) {
    RepeatIndex<-apply(as.matrix(as.numeric(levels(factor(Curvedata[,1L,drop=TRUE])))), 1L, 
                       function(x,y) which(abs(x-y)<=.Machine$double.eps^0.5), Curvedata[,1L,drop=TRUE])
    RepeatIndex<-unlist(RepeatIndex[sapply(RepeatIndex,length)==2L])
    RecycleRatio<-Curvedata[,2L,drop=TRUE][RepeatIndex[2L]]/Curvedata[,2L,drop=TRUE][RepeatIndex[1L]]
    if(RecycleRatio>1.2 || RecycleRatio<0.8) {
      Curvedata<-Curvedata[-RepeatIndex[2L],,drop=FALSE]
      cat("Note: recycling ratio is too large (small) and the repeated dose was removed from the growth curve!\n")
    } # end if
  } else {
    cat("Note: recycling ratio is not available!\n")
    RecycleRatio<-NULL
  } # end if
  ###
  ###
  ### Calculate recuperation (remove it from Curvedata if possible).
  if(length(which(abs(Curvedata[,1L,drop=TRUE])<=.Machine$double.eps^0.5))==1L) {
    Recuperation<-Curvedata[which(abs(Curvedata[,1L,drop=TRUE])<=.Machine$double.eps^0.5),2L,drop=TRUE]/NatureLxTx[1L]
    if(Recuperation>=0.2) {
      Curvedata<-Curvedata[-which(abs(Curvedata[,1L,drop=TRUE])<=.Machine$double.eps^0.5),,drop=FALSE]
      cat("Note: recuperation is too large so the 0-regenerative dose was removed from the growth curve!\n")
    } # end if
  } else {
    cat("Note: recuperation is not available!\n")
    Recuperation<-NULL
  } # end if
  ###
  ###
  ### Fit the growth curve.
  Models<-if(is.null(model)) c("line","exp","line+exp") else model
  Origins<-if(is.null(origin)) c(TRUE,FALSE) else origin
  minResidual<-1e30
  ### Find the most approporate fitting model automatically (minimum sum of squared residuals).
  for(i in Models) {
    for(j in Origins) {
      ### Try to call function calED().
      res<-try(calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=i, origin=j, nstart=nstart, upb=upb,
                     ErrorMethod=ErrorMethod, nsim=100L, plot=FALSE), silent=TRUE)
      ### Error checking.
      if(class(res)!="try-error") {
        if(res$value<minResidual)  {
          model<-i
          origin<-j
          minResidual<-res$value
          OK<-1L
        } # end if
      } # end if
    } # end for
  } # end for
  ###
  if (exists("OK")) {
      res<-calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=model, origin=origin, 
                 nstart=nstart, upb=upb, ErrorMethod=ErrorMethod, nsim=1000L, plot=TRUE)
  } else {
      par(bg="grey95", mgp=c(2,1,0), mar=c(3,3,2,1)+0.1)
      plot(Curvedata[,c(1L,2L)], main="Growth Curve", pch=21, cex=3, bg="white", xlab="Dose (Gy)", ylab="Standardised OSL",
           las=0, cex.main=1.25, cex.lab=1, xlim=c(0,max(Curvedata[,1L,drop=TRUE])*1.2), ylim=c(0,max(Curvedata[,2L,drop=TRUE])*1.2), 
           xaxs="i", yaxs="i", lab=c(7,7,9))
      grid(equilogs=FALSE)
      box(lwd=2L)
      par(bg="transparent", mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)
      stop("Error: fail in fitting the fast-component growth curve!")
  } # end if
  ### Reset par("mfrow") before leaving.
  par(mfrow=c(1L,1L))
  ### Output results.
  list("pars"=pars, "decayRate"=list("decayRateMat"=decayRateMat,"uniformity"=uniformity, "fastRate"=fastRate),
       "Curvedata"=Curvedata, "Ltx"=NatureLxTx, "model"=paste(model,"(origin=",origin,")",sep=""), "LMpars"=res$LMpars,
       "value"=res$value, "ED"=res$ED, "RecyclingRatio"=RecycleRatio, "Recuperation"=Recuperation)
  ###
} # end function fastED().
################################################ END FUNCTION fastED #################################################
###
### Default function for reporting S3 class object "mcAgeModels".
##################################################### FUNCTION reportSAM #############################################
### Author: Peng Jun, 2014.03.14; revised in 2014.04.12.
reportSAM<-
function(obj, burn=NULL, thin=NULL, plot=TRUE, outfile=NULL,...)  {
    UseMethod("reportSAM")
} #
reportSAM.default<-
function(obj, burn=NULL, thin=NULL, plot=TRUE, outfile=NULL,...)  {
  ### For argument "obj".
  if (class(obj)!="mcAgeModels")  {
      stop("Error: obj should be an object of class 'mcAgeModels'!")
  } # end if
  if (length(obj)!=6L) {
      stop("Error: incorrect object obj!")
  } # end if
  if (!all(names(obj)==c("EDdata","addsigma","model","iflog","nsim","chains"))) {
      stop("Error: incorrect members in object obj!")
  } # end if
  if (obj$nsim<100L) {
      stop("Error: the chains might be too short for a meaningful report!")
  } # end if
  ###
  ### For arugment "burn".
  if(!is.null(burn) && !is.numeric(burn))  {
      stop("Error: burn should be either NULL or be of type numeric!")
  } # end if 
  if(!is.null(burn))  {
      if(length(burn)!=1L)  {
          stop("Error: burn should be an one-element vector!")
      } # end if
      if(!is.finite(burn))  {
          stop("Error: burn should not be a non-finite value!")
      } # end if
      if(abs(burn-round(burn))>=.Machine$double.eps^0.5)  {
          stop("Error: burn should be an integer!")
      } # end if
      if(burn<0L)  {
          stop("Error: burn should not below zero!")
      } # end if 
      if(burn>=obj$nsim)  {
          stop("Error: burn should be less than the row number of the chains!")
      } # end if
  } else {
      ### Default number of burn.
      burn<-floor(obj$nsim/5L)
  } # end if
  ###
  ### For argument "thin".
  if(!is.null(thin) && !is.numeric(thin))  { 
      stop("Error: thin should either be NULL or of type numberic!")
  } # end if
  if(!is.null(thin))  {
      if(length(thin)!=1L)  {
          stop("Error: thin should be an one-element vector!")
      } # end if
      if(!is.finite(thin))  {
          stop("Error: thin should not be a non-finite value!")
      } # end if
      if(abs(thin-round(thin))>=.Machine$double.eps^0.5)  {
          stop("Error: thin should be an integer!")
      } # end if
      if(thin<1L) {
          stop("Error: thin should not less than one!")
      } # end if
      if( (obj$nsim-burn)/thin<16)  {
          stop("Error: resultant chains are too short after the burning and thining!")
      } # end if
  } else {
      ### Fefault number of thin.
      thin<-5L
  } # end if
  ###
  ### For argument "plot".
  if(!is.logical(plot))  {
      stop("Error: plot should either be TRUE or FALSE!")
  } # end if
  if(length(plot)!=1L)  {
      stop("Error: plot should be an one-element vector!")
  } # end if
  ###
  ### For argument "outfile"
  if(!is.null(outfile) && !is.character(outfile))  {
      stop("Error: outfile should either be NULL or of type character!")
  } # end if
  if (!is.null(outfile))  {
      if(length(outfile)!=1L) {
          stop("Error: outfile should be an one-element vector!")
      } # end if
  } # end if
  ###
  ###
  ### A handy function for calculating logged likelihood values.
  calLoglik<-function(y,x,pars,model)  {
      if(model=="CAM")  {
          Mu<-ifelse(obj$iflog==TRUE, log(pars[1L]), pars[1L])
          Sigmma<-ifelse(obj$iflog==TRUE, pars[2L], pars[2L]*pars[1L])
          ###
          Loglik<-1.0/sqrt(2.0*pi)/sqrt(x^2L+Sigmma^2L)*
                  exp(-(y-Mu)^2L/2.0/(x^2L+Sigmma^2L))
          ###
          return( sum(log(Loglik)) )
          ###
          ###
      } else if (model=="MAM3") {
          P<-pars[1L]
          Gamma<-ifelse(obj$iflog==TRUE, log(pars[2L]), pars[2L])
          Sigmma<-pars[3L]
          Mu0<-(Gamma/Sigmma^2L+y/x^2L)/(1.0/Sigmma^2L+1.0/x^2L)
          Sigmma0<-1.0/sqrt(1.0/Sigmma^2L+1.0/x^2L)
          ###
          part1<-P/sqrt(2.0*pi)/x*exp(-(y-Gamma)^2L/2.0/x^2L)
          part2<-(1.0-P)/sqrt(2.0*pi)/sqrt(Sigmma^2L+x^2L)*
                 (1.0-pnorm((Gamma-Mu0)/Sigmma0))*2.0*
                 exp(-(y-Gamma)^2L/2.0/(Sigmma^2L+x^2L))
          Loglik<-part1+part2
          ###
          return( sum(log(Loglik)) )
          ###
          ###
      } else if (model=="MAM4") {
          P<-pars[1L]
          Gamma<-ifelse(obj$iflog==TRUE, log(pars[2L]), pars[2L])
          Mu<-ifelse(obj$iflog==TRUE, log(pars[3L]), pars[3L])
          Sigmma<-pars[4L]
          Mu0<-(Mu/Sigmma^2L+y/x^2L)/(1.0/Sigmma^2L+1.0/x^2L)
          Sigmma0<-1.0/sqrt(1.0/Sigmma^2L+1.0/x^2L)
          ###
          part1<-P/sqrt(2.0*pi)/x*exp(-(y-Gamma)^2L/2.0/x^2L)
          part2<-(1.0-P)/sqrt(2.0*pi)/sqrt(Sigmma^2L+x^2L)*
          (1.0-pnorm((Gamma-Mu0)/Sigmma0))/(1.0-pnorm((Gamma-Mu)/Sigmma))*
          exp(-(y-Mu)^2L/2.0/(Sigmma^2L+x^2L))
          Loglik<-part1+part2
          ###
          return( sum(log(Loglik)) )
          ###
          ###
      } else {
          Ps<-pars[1L:(length(pars)/2L)]
          Mus<-if (obj$iflog==TRUE) { 
                  log( pars[(length(pars)/2L+1L):length(pars)] )
               } else {
                   pars[(length(pars)/2L+1L):length(pars)]
               } # end if
          ###
          Loglik<-Ps[1L]/sqrt(2.0*pi)/x*exp(-(y-Mus[1L])^2L/2.0/x^2L)
          for(i in 2L:(length(pars)/2L)) {
              Loglik<-Loglik+Ps[i]/sqrt(2.0*pi)/x*exp(-(y-Mus[i])^2L/2.0/x^2L)
          } # end if
          ###
          return( sum(log(Loglik)) )
      } # end if
  } # end function calLoglik
  ###
  ###
  ### Reshape the chains with arguments "burn" and "thin".
  if (burn>0L)  {
      chains<-obj$chains[-(1L:burn),,drop=FALSE]
  } else {
      chains<-obj$chains
  } # end if
  ###
  chains<-chains[seq(from=1L,to=obj$nsim-burn,by=thin),,drop=FALSE]
  Pars<-apply(chains,MARGIN=2L,mean)
  Std.Pars<-apply(chains,MARGIN=2L,sd)
  ###
  ### Calculate the maximum logged likelihood value.
  if (obj$iflog==TRUE)  {
      yyy<-obj$EDdata[,1L,drop=TRUE]
      xxx<-obj$EDdata[,2L,drop=TRUE]
      xxx<-sqrt( (xxx/yyy)^2L + (obj$addsigma)^2L)
      yyy<-log(yyy)
  } else {
      yyy<-obj$EDdata[,1L,drop=TRUE]
      xxx<-obj$EDdata[,2L,drop=TRUE]
      xxx<-sqrt( xxx^2L + (obj$addsigma)^2L)
  } # end if
  maxlik<-try(calLoglik(yyy,xxx,Pars,obj$model),silent=TRUE)
  if (class(maxlik)=="try-error")  {
      cat("Warning: logged likelihood value (maxlik) is not available!\n")
      maxlik<-NULL
  } # end if 
  ###
  ###
  ### Set quantiles for each variable.
  Probs<-t( apply(chains,MARGIN=2L,quantile, probs=c(0.025,0.25,0.5,0.75,0.975)) )
  ### Write out the chains or not.
  if(!is.null(outfile)) {
    write.csv(chains,file=paste(outfile,".csv",sep=""))
  } # end if
  ###
  ### Plot or not ? 
  if (plot==TRUE) {
      ###
      par(mfrow=c(ncol(chains),3L))
      par(mgp=c(2,1,0),
          mar=c(3,3,2,1)+0.1)
      namesPars<-colnames(chains)
      for (i in seq(ncol(chains)))  {
          DS<-density(chains[,i,drop=TRUE])
          plot(DS, main=paste("Density of ",namesPars[i]), ylab="")
          polygon(DS, col="grey")
          rug(chains[,i,drop=TRUE], quiet=TRUE)
          plot(chains[,i,drop=TRUE], type="l", main=paste("Trace of ",namesPars[i],sep=""), xlab="Iterations", ylab="")
          Autc<-acf(chains[,i,drop=TRUE], lag.max=30L, plot=FALSE)
          plot(Autc$lag, Autc$acf, main=paste("Autocorrelation of ",namesPars[i],sep=""), xlab="Lag", ylab="", type="h")
          abline(h=0)
      } # end for
      par(mfrow=c(1,1))
      par(mgp=c(3,1,0),
          mar=c(5,4,4,2)+0.1)
  } # end if
  ###
  list("pars"=round(data.frame("Pars"=Pars,"Std.Pars"=Std.Pars),5L),
       "quantile"=round(Probs,5L), "maxlik"=maxlik)
} # end function reportSAM.
##################################################### END FUNCTION reportSAM ###########################################
###
################################################### FUNCTION mcFMM #####################################################
### ********************************************************************************************************************
#### Function mcFMM() is used to perform Galbraith's statistical age models analysis with the Slice Sampling method.
#### Models that can be analyzed include: CAM and FMM(1-4).
###
###     Author: Peng Jun, 2014.03.14; revised in 2014.04.12.
###
### References: Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard 
###             Errors. Technometrics, 30 (3), pp. 271-281.
###
###             Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
###             vol. 31, no. 3, pp. 705-767.
###
### Mandatory arguments---
###     EDdata: a data.frame, equivalent doses (two columns).
### ******************************************************************************************************************
mcFMM<-
function(EDdata, ncomp=1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
    UseMethod("mcFMM")
} # end function mcFMM
### Default method for function mcFMM().
mcFMM.default<-
function(EDdata, ncomp=1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
  ###
  ### For argument "EDdata".
  if(!is.data.frame(EDdata))  {
      stop("Error: EDdata should be of type data.frame!")
  } # end if
  if(ncol(EDdata)!=2L)  {
      stop("Error: EDdata should contain two columns!")
  } # end if
  if(!is.numeric(as.matrix(EDdata))) {
      stop("Error: elements in EDdata should be all of type numeric!")
  } # end if
  if(any(!is.finite(unlist(unclass(EDdata)))))  {
      stop("Error: EDdata should not contain non-finite value!")
  } # end if
  if(any(EDdata[,2L,drop=TRUE]<=0))  {
      stop("Error: all std.errors in EDdata should be larger than zero!")
  } # end if
  ###
  ### For argument "ncomp".
  if(!is.numeric(ncomp))  {
      stop("Error: ncomp should be of type numeric!")
  } # end if
  if(length(ncomp)!=1L)  {
      stop("Error: ncomp should be an one-element vector!")
  } # end if
  if(!is.finite(ncomp))  {
      stop("Error: ncomp should not be a non-finite value!")
  } # end if
  if(!ncomp %in% (1L:4L))  {
      stop("Error: ncomp should be an integer in the space [1,4]!")
  } # end if
  ###
  ### For argument "addsigma".
  if(!is.numeric(addsigma))  {
      stop("Error: addsigma should be of type numeric!")
  } # end if
  if(length(addsigma)!=1L)  {
      stop("Error: addsigma should be an one-element vector!")
  } # end if
  if(!is.finite(addsigma))  {
      stop("Error: addsigma should not be a non-finite value!")
  } # end if
  if(addsigma<0)  {
      stop("Error: addsigma should not be smaller than zero!")
  } # end if
  ###
  ### For argument "iflog".
  if(!is.logical(iflog))  {
      stop("Error: iflog should be either TRUE or FALSE!")
  } # end if
  if(length(iflog)!=1L)  {
      stop("Error: iflog should be an one-element vector!")
  } # end if
  if(any(EDdata[,1L,drop=TRUE]<=0) && iflog==TRUE)  {
      stop("Error: minus (zero) ED values could not be logged, please set iflog to be FALSE!")
  } # end if
  ### 
  ### For argument "nsim".
  if(!is.numeric(nsim))  {
      stop("Error: nsim should be of type numeric!")
  } # end if
  if(length(nsim)!=1L)  {
      stop("Error: nsim should be an one-element vector!")
  } # end if
  if(!is.finite(nsim))  { 
      stop("Error: nsim should not be a non-finite value!")
  } # end if
  if(abs(nsim-round(nsim))>=.Machine$double.eps^0.5)  {
      stop("Error: nsim should be an integer!")
  } # end if
  if(nsim<100L || nsim>2e5)  {
       stop("Error: nsim should be in the space [1e2,2e5]!")
  } # end if
  ###
  ###
  ###
  ### For argument "inis".
  if(class(inis)!="list")  {
      stop("Error: inis should be a list!")
  } # end if
  if(length(inis)>2L*ncomp)  {
      stop("Error: incorrect number of parameters in inis!")
  } # end if
  if(ncomp==1L && !all(names(inis) %in% c("mu","sigma")) )  {
      stop("Error: incorrect names of parameters (CAM) in inis!")
  } # end if
  if(ncomp>1L && !all(names(inis) %in% c(paste("p",1L:ncomp,sep=""),paste("mu",1L:ncomp,sep=""))))  {
      stop("Error: incorrect names of parameters (FMM) in inis!")
  } # end if
  ###
  ### Set ranges of Mus.
  rangeED<-range(EDdata[,1L,drop=TRUE])
  if(all(EDdata[,1L,drop=TRUE]>0))  {
      lowerMus<-rangeED[1L]*0.999
      upperMus<-rangeED[2L]*1.001
  } else if (all(EDdata[,1L,drop=TRUE]<=0))  {
      lowerMus<-rangeED[1L]*1.001
      upperMus<-rangeED[2L]*0.999
  } else {
      lowerMus<-rangeED[1L]*1.001
      upperMus<-rangeED[2L]*1.001
  } # end if
  ### 
  if (ncomp==1L)  {
      ### Default inis for simulation of CAM.
      args.inis<-list("mu"=mean(EDdata[,1L,drop=TRUE]),"sigma"=0.5)
      upperSigma<-ifelse(iflog==TRUE,5.0,var(EDdata[,1L,drop=TRUE]))
  } else {
      if (ncomp==2L)  {
          args.inis<-list("p1"=1.0,"mu1"=mean(EDdata[,1L,drop=TRUE]),
                          "p2"=1.0,"mu2"=mean(EDdata[,1L,drop=TRUE]))
      } # end if
      ###
      if (ncomp==3L)  {
          ### Default inis for simulation of FMM3.
          args.inis<-list("p1"=1.0,"mu1"=mean(EDdata[,1L,drop=TRUE]),
                          "p2"=1.0,"mu2"=mean(EDdata[,1L,drop=TRUE]),
                          "p3"=1.0,"mu3"=mean(EDdata[,1L,drop=TRUE]))
      } # end if
      ###
      if (ncomp==4L)  {
          ### Default inis for simulation of FMM3.
          args.inis<-list("p1"=1.0,"mu1"=mean(EDdata[,1L,drop=TRUE]),
                          "p2"=1.0,"mu2"=mean(EDdata[,1L,drop=TRUE]),
                          "p3"=1.0,"mu3"=mean(EDdata[,1L,drop=TRUE]),
                          "p4"=1.0,"mu4"=mean(EDdata[,1L,drop=TRUE]))
      } # end if
  } # end if
  ###
  ###
  ### Pass specified parameters to list "args.inis".
  args.inis[names(inis)]<-inis
  ### For all arguments in "inis".
  if(any(!sapply(args.inis,is.numeric))) {
      stop("Error: all members in inis should be of type numeric!")
  } # end if
  if(any(sapply(args.inis,length)!=1L)) {
      stop("Error: all members in inis should be one-element vector!")
  } # end if
  if(any(!sapply(args.inis,is.finite))) {
      stop("Error: all members in inis should not be non-finite value!")
  } # end if
  if(ncomp==1L) {
      if(args.inis["mu"]<=lowerMus || args.inis["mu"]>=upperMus)  {
          stop(paste("Error: mu value should be initialized in the space (",lowerMus," ,",upperMus,")!",sep=""))
      } # end if
      if(args.inis["sigma"]<=0 || args.inis["sigma"]>=upperSigma)  {
          stop(paste("Error: sigma should be initialized in the space (",0," ,",upperSigma,")!",sep=""))
      } # end if
  } else {
      if(any(args.inis[paste("p",1L:ncomp,sep="")]<=0)) {
          stop("Error: all p values in inis should be larger than zero!")
      } # end if
      if(any(args.inis[paste("mu",1L:ncomp,sep="")]<=lowerMus |
             args.inis[paste("mu",1L:ncomp,sep="")]>=upperMus))  {
          stop(paste("Error: all mu values should be initialized in the space (",lowerMus," ,",upperMus,")!",sep=""))
      } # end if
  } # end if
  ###
  ###
  ### For argument "control.args".
  if(class(control.args)!="list")  {
      stop("Error: control.args should be a list!")
  } # end if
  if(length(control.args)>3L)  {
      stop("Error: control.args should contain three parameters (w, m, nstart)!")
  } # end if
  if(!all(names(control.args) %in% c("w","m","nstart")))  {
      stop("Error: incorrect names of parameter in control.args!")
  } # end if
  ### Default parameters in list "args.control".
  args.control<-list(w=1.0,m=-100.0,nstart=1L)
  ### Pass specified parameters to list "args.control".
  args.control[names(control.args)]<-control.args
  ### Set w, m, nstart
  w<-args.control$w
  m<-args.control$m
  nstart<-args.control$nstart
  ### For arguments in "control.args".
  if(!is.numeric(w) || !is.numeric(m) || !is.numeric(nstart))  {
      stop("Error: all members in control.args should be of type numeric!")
  } # end if
  if(length(w)!=1L || length(m)!=1L || length(nstart)!=1L)  {
      stop("Error: all members in conrol.args should be one-element vector!")
  } # end if
  if(!is.finite(w) || !is.finite(m) || !is.finite(nstart))  {
      stop("Error: all members in control.args should not be non-finite value!")
  } # end if
  if(w<1e-2)  {
      stop("Error: w is too small, the sampling will be very inefficient!")
  } # end if
  if(w>1e3)  {
      stop("Error: w is too large (not exceed 1e3)!")
  } # end if
  if(m>1e9)  {
     stop("Error: m should not exceed 1e9!")
  } # end if
  if (abs(nstart-round(nstart))>=.Machine$double.eps^0.5)  {
      stop("Error: nstart should be an integer!")
  } # end if
  if(nstart<=0L || nstart>1e6)  {
      stop("Error: nstart should be in the space [1,1e6]!")
  } # end if
  ###
  ###
  nED<-nrow(EDdata)
  ED<-EDdata[,1L,drop=TRUE]
  Error<-EDdata[,2L,drop=TRUE]
  iflag<-0 
  chains<-matrix(nrow=nsim,ncol=2L*ncomp)
  if(ncomp==1L) {
      inis<-unlist(args.inis,use.names=FALSE)
  } else {
      inis<-matrix(unlist(args.inis,use.names=FALSE),ncol=ncomp)
  } # end if
  subFortran<- ifelse(ncomp==1L, "mcCAM", paste("mcFMM",ncomp,sep=""))
  ###
  res<-.Fortran(subFortran, as.integer(nED), as.integer(nsim), as.double(ED), as.double(Error),
                as.double(addsigma), as.double(inis), as.integer(iflog), as.integer(nstart),
                as.double(w), as.double(m), chains=as.double(chains), iflag=as.integer(iflag),
                NAOK=TRUE, package="numOSL") 
  ### Error checking.
  if(res$iflag!=0)  {
      Count<-sum(res$chains[1L:nsim]>0)
      if (Count<=1L) {
          stop("Error: unreasonable initials, chains break at the beginning of the simulation!")
      } # end if
      cat(paste("Warning: the chains crashed down at the ",Count," th simulation!\n",sep=""))
  } # end if
  ###
  ###
  chains<-matrix(res$chains,ncol=2L*ncomp)
  ### Reshape the chains.
  if(res$iflag!=0)  {
      chains<-chains[1L:(Count-1L), ,drop=FALSE]
  } # end if
  ###
  if (ncomp==1L) {
      colnames(chains)<-c("mu","sigma")
      if (iflog==TRUE) {
          chains[,1L]<-exp(chains[,1L,drop=TRUE])
      } else {
          chains[,2L]<-chains[,2L,drop=TRUE]/chains[,1L,drop=TRUE]
      } # end if     
  } else {
      ### Sort (increasing) the chains before output for the FMM models.
      chains<-t( apply(array(t(as.matrix(chains)),dim=c(ncomp,2L,nrow(chains))), 
                       MARGIN=3L, function(x)  x[order(x[,2L,drop=TRUE]),] ) )
      colnames(chains)<-c(paste("p",1L:ncomp,sep=""),paste("mu",1L:ncomp,sep=""))
      if (iflog==TRUE) {
          chains[,(ncomp+1L):(2L*ncomp)]<-exp( chains[,(ncomp+1L):(2L*ncomp)] )
      } # end if
  } # end if
  ###
  ###
  ### Set the output.
  out<-list("EDdata"=EDdata, "addsigma"=addsigma, "model"=ifelse(ncomp==1L,"CAM",paste("FMM",ncomp,sep="")), 
            "iflog"=iflog, "nsim"=ifelse(res$iflag==0,nsim,Count-1L), "chains"=data.frame(chains))
  class(out)<-"mcAgeModels"
  invisible(out)
} # end function mcFMM.      
################################################### END FUNCTION mcFMM ###################################################
###
################################################### FUNCTION mcMAM #######################################################
### **********************************************************************************************************************
#### Function mcMAM() is used to perform Galbraith's statistical age models analysis with the Slice Sampling method.
#### Models that can be analyzed include: MAM3 and MAM4.
###
###     Author: Peng Jun, 2014.03.14; revised in 2014.04.12.
###
### References: Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard 
###             Errors. Technometrics, 30 (3), pp. 271-281.
###
###             Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
###             vol. 31, no. 3, pp. 705-767.
###
### Mandatory arguments---
###     EDdata: a data.frame, equivalent doses (two columns).
### **********************************************************************************************************************
mcMAM<-
function(EDdata, ncomp=-1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
    UseMethod("mcMAM")
} # end function mcMAM
### Set default method for function mcMAM().
mcMAM.default<-
function(EDdata, ncomp=-1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
  ###
  ### For argument "EDdata".
  if(!is.data.frame(EDdata))  {
      stop("Error: EDdata should be of type data.frame!")
  } # end if
  if(ncol(EDdata)!=2L)  {
      stop("Error: EDdata should contain two columns!")
  } # end if
  if(!is.numeric(as.matrix(EDdata)))  {
      stop("Error: elements in EDdata should be all of type numeric!")
  } # end if
  if(any(!is.finite(unlist(unclass(EDdata)))))  {
      stop("Error: EDdata should not contain non-finite value!")
  } # end if
  if(any(EDdata[,2L,drop=TRUE]<=0))  {
      stop("Error: all std.errors in EDdata should be larger than zero!")
  } # end if
  ###
  ### For argument "ncomp".
  if(!is.numeric(ncomp))  {
      stop("Error: ncomp should be of type numeric!")
  } # end if
  if(length(ncomp)!=1L)  {
      stop("Error: ncomp should be an one-element vector!")
  } # end if
  if(!is.finite(ncomp))  {
      stop("Error: ncomp should not be a non-finite value!")
  } # end if
  if(!ncomp %in% c(-1L,-2L) )  {
      stop("Error: ncomp should be either -1 (MAM3) or -2 (MAM4)!")
  } # end if
  ###
  ### For argument "addsigma".
  if(!is.numeric(addsigma))  {
      stop("Error: addsigma should be of type numeric!")
  } # end if
  if(length(addsigma)!=1L)  {
      stop("Error: addsigma should be an one-element vector!")
  } # end if
  if(!is.finite(addsigma))  {
      stop("Error: addsigma should not be a non-finite value!")
  } # end if
  if(addsigma<0)  {
      stop("Error: addsigma should not be smaller than zero!")
  } # end if
  ###
  ### For argument "iflog".
  if(!is.logical(iflog))  {
      stop("Error: iflog should be either TRUE or FALSE!")
  } # end if
  if(length(iflog)!=1L)  {
      stop("Error: iflog should be an one-element vector!")
  } # end if
  if(any(EDdata[,1L,drop=TRUE]<=0) && iflog==TRUE)  {
      stop("Error: minus (zero) ED values could not be logged, please set iflog to be FALSE!")
  } # end if
  ###
  ### For argument "nsim".
  if(!is.numeric(nsim))  {
      stop("Error: nsim should be of type numeric!")
  } # end if
  if(length(nsim)!=1L)  {
      stop("Error: nsim should be an one-element vector!")
  } # end if
  if(!is.finite(nsim))  {
      stop("Error: nsim should not be a non-finite value!")
  } # end if
  if (abs(nsim-round(nsim))>=.Machine$double.eps^0.5)  {
      stop("Error: nsim should be an integer!")
  } # end if
  if(nsim<100L || nsim>2e5)  {
      stop("Error: nsim should be in the space [1e2,2e5]!")
  } # end if
  ###
  ###
  ### For argument "inis".
  if(class(inis)!="list")  {
      stop("Error: inis should be a list!")
  } # end if
  if(length(inis)>2L-ncomp)  {
      stop("Error: incorrect number of parameters in inis!")
  } # end if
  if(ncomp==-1L && !all(names(inis) %in% c("p","gamma","sigma")) )  {
      stop("Error: incorrect names of parameters (MAM3) in inis!")
  } # end if
  if(ncomp==-2L && !all(names(inis) %in% c("p","gamma","mu","sigma")) )  {
      stop("Error: incorrect names of parameters (MAM4) in inis!")
  } # end if
  ###
  ### Set the boundary of gamma (mu).
  rangeED<-range(EDdata[,1L,drop=TRUE])
  if(all(EDdata[,1L,drop=TRUE]>0))  {
      lowerGamma<-rangeED[1L]*0.999
      upperGamma<-rangeED[2L]*1.001
  } else if (all(EDdata[,1L,drop=TRUE]<=0))  {
      lowerGamma<-rangeED[1L]*1.001
      upperGamma<-rangeED[2L]*0.999
  } else {
      lowerGamma<-rangeED[1L]*1.001
      upperGamma<-rangeED[2L]*1.001
  } # end if
  ###
  if (ncomp==-1L)  {
      ### Default inis for simulation of MAM3.
      args.inis<-list("p"=0.5,"gamma"=mean(EDdata[,1L,drop=TRUE]),"sigma"=0.5)
  } else if (ncomp==-2L) {
      ### Default inis for simulation of MAM4.
      args.inis<-list("p"=0.5,"gamma"=min(EDdata[,1L,drop=TRUE]),"mu"=mean(EDdata[,1L,drop=TRUE]),"sigma"=0.5)
  } # end if
  ### Pass specified parameters to list "args.inis".
  args.inis[names(inis)]<-inis
  if(any(!sapply(args.inis,is.numeric)) )  {
      stop("Error: all members in inis should be of type numeric!")
  } # end if
  if(any(sapply(args.inis,length)!=1L))  {
      stop("Error: all members in inis should be one-element vector!")
  } # end if
  if(any(!sapply(args.inis,is.finite))) {
      stop("Error: all members in inis should not be non-finite value!")
  } # end if
  if(args.inis["p"]<=0 || args.inis["p"]>=1)  {
      stop("Error: p should be initialized in the space (0, 1)!")
  } # end if
  if(args.inis["gamma"]<=lowerGamma || args.inis["gamma"]>=upperGamma)  {
      stop(paste("Error: gamma should be initialized in the space (",lowerGamma," ,",upperGamma,")!",sep=""))
  } # end if
  if(ncomp==-2L && (args.inis["mu"]<=lowerGamma || args.inis["mu"]>=upperGamma) )  {
      stop(paste("Error: mu should be initialized in the space (",lowerGamma," ,",upperGamma,")!",sep=""))
  } # end if
  ###
  upperSigma<-ifelse(iflog==TRUE,5.0,var(EDdata[,1L,drop=TRUE]))
  if(args.inis["sigma"]<=0 || args.inis["sigma"]>=upperSigma)  {
      stop(paste("Error: sigma should be initialized in the space (0,",upperSigma,")!",sep=""))
  } # end if
  ###
  ### For argument "control.args".
  if(class(control.args)!="list")  {
      stop("Error: control.args should be a list!")
  } # end if
  if(length(control.args)>3L)  {
      stop("Error: control.args should contain at most three parameters (w, m, nstart)!")
  } # end if
  if(!all(names(control.args) %in% c("w","m","nstart")))  {
      stop("Error: incorrect names of parameters in control.args!")
  } # end if
  ### Default parameters in list "args.control".
  args.control<-list(w=1.0,m=-100.0,nstart=1L)
  ### Pass specified parameters to list "args.control".
  args.control[names(control.args)]<-control.args
  ### Set w, m, nstart
  w<-args.control$w
  m<-args.control$m
  nstart<-args.control$nstart
  ### For argument in "control.args".
  if(!is.numeric(w) || !is.numeric(m) || !is.numeric(nstart))  {
      stop("Error: all members in control.args should be of type numeric!")
  } # end if
  if(length(w)!=1L || length(m)!=1L || length(nstart)!=1L)  {
      stop("Error: all members in control.args should be one-element vector!")
  } # end if
  if(!is.finite(w) || !is.finite(m) || !is.finite(nstart))  {
      stop("Error: all members in control.args should not be non-finite value!")
  } # end if
  if(w<1e-2)  {
      stop("Error: w is too small, the sampling will be very inefficient!")
  } # end if
  if(w>1e3)  {
      stop("Error: w is too large (not exceed 1e3)!")
  } # end if
  if(m>1e9)  {
     stop("Error: m should not exceed 1e9!")
  } # end if
  if (abs(nstart-round(nstart))>=.Machine$double.eps^0.5)  {
      stop("Error: nstart should be an integer!")
  } # end if
  if(nstart<=0L || nstart>1e6)  {
      stop("Error: nstart should be in the space [1,1e6]!")
  } # end if
  ###
  ###
  ### Arguments for Fortran subroutine mcMAM3() and mcMAM4().
  nED<-nrow(EDdata)
  ED<-EDdata[,1L,drop=TRUE]
  Error<-EDdata[,2L,drop=TRUE]
  iflag<-0
  inis<-unlist(args.inis,use.names=FALSE)
  chains<-matrix(nrow=nsim, ncol=2L-ncomp)
  subFortran<-ifelse(ncomp==-1L,"mcMAM3","mcMAM4")
  ###
  res<-.Fortran(subFortran, as.integer(nED), as.integer(nsim), as.double(ED), as.double(Error),
                as.double(addsigma), as.double(inis), as.integer(iflog), as.integer(nstart),
                as.double(w), as.double(m), chains=as.double(chains), iflag=as.integer(iflag),
                NAOK=TRUE, package="numOSL")
  ### Error checking.
  if(res$iflag!=0)  {
      Count<-sum(res$chains[1L:nsim]>0)
      if (Count<=1L) {
          stop("Error: unreasonable initials, chains break at the beginning of the simulation!")
      } # end if
      cat(paste("Warning: the chains crashed down at the ",Count," th simulation!\n",sep=""))
  } # end if
  ###
  ###
  ### Reshape the chains.
  chains<-data.frame(matrix(res$chains,ncol=2L-ncomp,dimnames=list(NULL,names(args.inis))))
  if (res$iflag!=0)  {
      chains<-chains[1L:(Count-1L),,drop=FALSE]
  } # end if
  ###
  if (iflog==TRUE) {
      if (ncomp==-1L)  {
          chains[,2L]<-exp(chains[,2L,drop=TRUE])   
      } else if (ncomp==-2L)  {
          chains[,2L:3L]<-exp(chains[,2L:3L])
      } # end if
  } # end if
  ### 
  ### Set the output.
  out<-list("EDdata"=EDdata, "addsigma"=addsigma, "model"=ifelse(ncomp==-1L,"MAM3","MAM4"),
            "iflog"=iflog, "nsim"=ifelse(res$iflag==0,nsim,Count-1L), "chains"=chains)
  class(out)<-"mcAgeModels"
  invisible(out)  
} # end function mcMAM.
#################################################### END FUNCTION mcMAM ############################################
