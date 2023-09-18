#####
reportMC<-
function(obj, burn=1e4, thin=5,
         plot=TRUE, outfile=NULL,...) {
    UseMethod("reportMC")
} #
### 2023.09.10.
reportMC.default<-
function(obj, burn=1e4, thin=5,
         plot=TRUE, outfile=NULL,...) {
    ### Stop if not.
    stopifnot(class(obj)=="mcAgeModels",
              is.list(obj), length(obj)==6L,
              all(names(obj)==c("EDdata","addsigma",
              "model","iflog","nsim","chains")),
              length(burn)==1L, is.numeric(burn), 
              burn>=0L, burn<obj$nsim,
              length(thin)==1L, is.numeric(thin), 
              thin>=1L, thin<(obj$nsim-burn)/100L,
              length(plot)==1L, is.logical(plot),
              is.null(outfile) || is.character(outfile))
    ###
    ### R function calculates the logged maximum likelihhod.
    calLoglik<-function(y,x,pars,model)  {
        if(model=="CAM")  {
            Mu<-ifelse(obj$iflog==TRUE,log(pars[1L]),pars[1L])
            Sigmma<-pars[2L]
            Loglik<-1.0/sqrt(2.0*pi)/sqrt(x^2L+Sigmma^2L)*
                    exp(-(y-Mu)^2L/2.0/(x^2L+Sigmma^2L))
            maxlik <- sum(log(Loglik))
            bic <- -2.0*maxlik+2.0*log(length(x))
        } else if (model %in% c("MAM3","MXAM3")) {
            P<-pars[1L]
            if (model=="MAM3") {
                Gamma<-ifelse(obj$iflog==TRUE,log(pars[2L]),pars[2L])
            } else if (model=="MXAM3") {
                Gamma<-ifelse(obj$iflog==TRUE,-log(pars[2L]),-pars[2L])
            } # end if.    
            Sigmma<-pars[3L]
            Mu0<-(Gamma/Sigmma^2L+y/x^2L)/(1.0/Sigmma^2L+1.0/x^2L)
            Sigmma0<-1.0/sqrt(1.0/Sigmma^2L+1.0/x^2L)
            part1<-P/sqrt(2.0*pi)/x*exp(-(y-Gamma)^2L/2.0/x^2L)
            part2<-(1.0-P)/sqrt(2.0*pi)/sqrt(Sigmma^2L+x^2L)*
                   (1.0-pnorm((Gamma-Mu0)/Sigmma0))*2.0*
                   exp(-(y-Gamma)^2L/2.0/(Sigmma^2L+x^2L))
            Loglik<-part1+part2
            maxlik <- sum(log(Loglik)) 
            bic <- -2.0*maxlik+3.0*log(length(x))
        } else if (model %in% c("MAM4","MXAM4")) {
            P<-pars[1L]
            if (model=="MAM4") {
                Gamma<-ifelse(obj$iflog==TRUE,log(pars[2L]),pars[2L])
                Mu<-ifelse(obj$iflog==TRUE,log(pars[3L]),pars[3L])
            } else if (model=="MXAM4") {
                Gamma<-ifelse(obj$iflog==TRUE,-log(pars[2L]),-pars[2L])
                Mu<-ifelse(obj$iflog==TRUE,-log(pars[3L]),-pars[3L])
            } # end if.
            Sigmma<-pars[4L]
            Mu0<-(Mu/Sigmma^2L+y/x^2L)/(1.0/Sigmma^2L+1.0/x^2L)
            Sigmma0<-1.0/sqrt(1.0/Sigmma^2L+1.0/x^2L)
            part1<-P/sqrt(2.0*pi)/x*exp(-(y-Gamma)^2L/2.0/x^2L)
            part2<-(1.0-P)/sqrt(2.0*pi)/sqrt(Sigmma^2L+x^2L)*
            (1.0-pnorm((Gamma-Mu0)/Sigmma0))/
            (1.0-pnorm((Gamma-Mu)/Sigmma))*
            exp(-(y-Mu)^2L/2.0/(Sigmma^2L+x^2L))
            Loglik<-part1+part2
            maxlik <- sum(log(Loglik)) 
            bic <- -2.0*maxlik+4.0*log(length(x))
        } else {
            Ps<-pars[1L:(length(pars)/2L)]
            Mus<-if (obj$iflog==TRUE) { 
                log( pars[(length(pars)/2L+1L):length(pars)] )
            } else {
                pars[(length(pars)/2L+1L):length(pars)]
            } # end if
            Loglik<-Ps[1L]/sqrt(2.0*pi)/x*
                    exp(-(y-Mus[1L])^2L/2.0/x^2L)
            for(i in 2L:(length(pars)/2L)) {
                Loglik<-Loglik+Ps[i]/sqrt(2.0*pi)/x*
                        exp(-(y-Mus[i])^2L/2.0/x^2L)
            } # end if
            maxlik <- sum(log(Loglik)) 
            bic <- -2.0*maxlik+(length(pars)-1.0)*log(length(x))
        } # end if
        return(c(maxlik, bic))
    } # end function calLoglik
    ###
    ### Burn-in. 
    if (burn>0L) {
        chains<-obj$chains[-seq(burn),,drop=FALSE]
    } else {
        chains<-obj$chains
    } # end if
    ###
    ### Thinning.
    chains<-chains[seq(from=1L,to=obj$nsim-burn,by=thin),,drop=FALSE]
    Pars<-apply(chains,MARGIN=2L,mean)
    Std.Pars<-apply(chains,MARGIN=2L,sd)
    ###
    funPars2 <- function(x) {
        ddd <- density(x)
        ddd$x[which.max(ddd$y)]
    } # end function funPars2.
    Pars2 <- apply(chains, MARGIN=2L, funPars2)
    Pars2[1:(length(Pars2)/2)] <- Pars2[1:(length(Pars2)/2)]/sum(Pars2[1:(length(Pars2)/2)])
    ###
    ed1<-as.numeric(obj$EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(obj$EDdata[,2L,drop=TRUE])
    ###
    if (obj$iflog==TRUE)  {
        yyy<-ed1
        xxx<-sed1
        xxx<-sqrt((xxx/yyy)^2L+(obj$addsigma)^2L)
        if (!obj$model %in% c("MXAM3","MXAM4")) {
            yyy<-log(yyy)
        } else {
            yyy<- -log(yyy)
        } # end if.
    } else {
        yyy<-ed1
        xxx<-sed1
        xxx<-sqrt(xxx^2L+(obj$addsigma)^2L)
        if (obj$model %in% c("MXAM3","MXAM4")) {
            yyy<- -yyy
        } # end if.
    } # end if
    ### Calculate the maxliklihood value.
    maxlikbic1<-try(calLoglik(yyy,xxx,Pars,obj$model),silent=TRUE)
    if (inherits(maxlikbic1,what="try-error")==TRUE)  {
        cat("Warning: 'maxlik' and 'bic' cannot be calculated using the mean parameters!\n")
        maxlik1<-bic1<-NULL
    } else {
        maxlik1 <- maxlikbic1[1L]
        bic1 <- maxlikbic1[2L]
    } # end if.
    maxlikbic2<-try(calLoglik(yyy,xxx,Pars2,obj$model),silent=TRUE)
    if (inherits(maxlikbic2,what="try-error")==TRUE)  {
        cat("Warning: 'maxlik' and 'bic' cannot be calculated using the mode parameters!\n")
        maxlik2<-bic2<-NULL
    } else {
        maxlik2 <- maxlikbic2[1L]
        bic2 <- maxlikbic2[2L]
    } # end if.
    ###
    ###
    Probs<-t(apply(chains,MARGIN=2L,quantile, 
             probs=c(0.025,0.25,0.5,0.75,0.975)))
    ###
    if(!is.null(outfile)) {
        write.csv(chains,file=paste(outfile,".csv",sep=""))
    } # end if
    ###
    ###
    output<-list("pars"=cbind("Mean.Pars"=Pars,"Sd.Pars"=Std.Pars,"Mode.pars"=Pars2),
                 "quantile"=Probs, "maxlik"=c(maxlik1,maxlik2), "bic"=c(bic1,bic2))
    ###
    if (plot==TRUE) {
        ###
        opar <- par("mfrow", "mgp", "mar")
        on.exit(par(opar))
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
            plot(chains[,i,drop=TRUE], type="l", main=paste("Trace of ",
                 namesPars[i],sep=""), xlab="Iterations", ylab="")
            Autc<-acf(chains[,i,drop=TRUE], lag.max=30L, plot=FALSE)
            plot(Autc$lag, Autc$acf, main=paste("Auto.Cor of ",
                 namesPars[i],sep=""), xlab="Lag", ylab="", type="h")
            abline(h=0)
        } # end for
        ###
    } # end if
    ###
    return(output)
} # end function reportMC.default.
#####
