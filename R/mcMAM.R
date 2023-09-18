#####
mcMAM<-
function(EDdata, ncomp=-1, addsigma=0, iflog=TRUE, 
         nsim=5e4, inis=list(), control.args=list()) {
    UseMethod("mcMAM")
} #
### 2023.09.10.
mcMAM.default<-
function(EDdata, ncomp=-1, addsigma=0, iflog=TRUE, 
         nsim=5e4, inis=list(), control.args=list()) {
    ### Stop if not.
    stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=5L,
              all(EDdata[,2L,drop=TRUE]>0),
              length(ncomp)==1L, ncomp %in% c(-1L,-2L,-3L,-4L),
              length(addsigma)==1L, is.numeric(addsigma),
              length(iflog)==1L, is.logical(iflog),
              length(nsim)==1L, is.numeric(nsim), nsim>=100L, nsim<=2e5,
              is.list(inis), is.list(control.args),
              all(names(inis) %in% c("p","gamma","mu","sigma")), 
              all(names(control.args) %in% c("w","m","nstart")))
    ###
    ed1<-as.numeric(EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(EDdata[,2L,drop=TRUE])
    nED<-nrow(EDdata)
    if (ncomp %in% c(-1L,-2L))  {
        inverse <- FALSE 
    } else {
        inverse <- TRUE
    } # end if.
    ###
    if (iflog==TRUE && any(ed1<=0)) {
        stop("Error: minus(zero) ED canot be logged!")
    } # end if
    ### Bounds of dose population.
    rangeED<-range(ed1)
    if(all(ed1>0))  {
        lowerGamma<-rangeED[1L]*0.999
        upperGamma<-rangeED[2L]*1.001
    } else if (all(ed1<=0))  {
        lowerGamma<-rangeED[1L]*1.001
        upperGamma<-rangeED[2L]*0.999
    } else {
        lowerGamma<-rangeED[1L]*1.001
        upperGamma<-rangeED[2L]*1.001
    } # end if
    ###
    ### Default inis.
    args.inis<-list("p"=0.5,"gamma"=quantile(ed1,probs=0.25,names=FALSE),
                    "mu"=quantile(ed1,probs=0.5,names=FALSE),
                    "sigma"=ifelse(iflog==TRUE,0.3,3.0))
    args.inis[names(inis)]<-inis
    stopifnot(args.inis[["p"]]>0.0, args.inis[["p"]]<1.0,
              args.inis[["gamma"]]>lowerGamma, 
              args.inis[["gamma"]]<upperGamma,
              args.inis[["mu"]]>lowerGamma, 
              args.inis[["mu"]]<upperGamma,
              args.inis[["sigma"]]>0.0,
              args.inis[["sigma"]]<ifelse(iflog==TRUE,5.0,50.0))
    ###
    ### Default arguments of slice sampling.
    args.control<-list(w=1, m=-100, nstart=1L)
    args.control[names(control.args)]<-control.args
    stopifnot(args.control[["w"]]>=1e-2,
              args.control[["w"]]<=1e3,
              args.control[["m"]]<=1e9,
              args.control[["nstart"]]>=1L,
              args.control[["nstart"]]<=1000L)
    ###
    w<-args.control$w
    m<-args.control$m
    nstart<-args.control$nstart
    inis<-if (ncomp==-1L) {
        c(args.inis[["p"]],
          args.inis[["gamma"]],
          args.inis[["sigma"]])
    } else {
        c(args.inis[["p"]],
          args.inis[["gamma"]],
          args.inis[["mu"]],
          args.inis[["sigma"]])
    } # end if
    iflag<-0
    chains<-matrix(0,nrow=nsim,ncol=ifelse(ncomp %in% c(-1L,-3L),3L,4L))
    ###
    if (ncomp %in% c(-1L,-3L)) {
        res <- .Fortran("mcMAM3",as.integer(nED),as.integer(nsim),as.double(ed1), 
                        as.double(sed1),as.double(addsigma),as.double(inis), 
                        as.integer(iflog),as.integer(nstart),as.double(w), 
                        as.double(m),chains=as.double(chains),as.integer(inverse),
                        iflag=as.integer(iflag),PACKAGE="numOSL")
    } else if (ncomp %in% c(-2L,-4L)) {
        res <- .Fortran("mcMAM4",as.integer(nED),as.integer(nsim),as.double(ed1), 
                        as.double(sed1),as.double(addsigma),as.double(inis), 
                        as.integer(iflog),as.integer(nstart),as.double(w), 
                        as.double(m),chains=as.double(chains),as.integer(inverse), 
                        iflag=as.integer(iflag),PACKAGE="numOSL")
    } # end if.
    ###
    if (res$iflag!=0) {
        niter<-sum(res$chains[seq(nsim)]>0)
        stop(paste("Error: the simulation failed at the ", 
                   niter+1L,"th iteration!",sep=""))
    } # end if
    ###
    if (ncomp==-1L)  {
        myNAMES <-  c("Prop","MAM3.De","Sigma")
    } else if (ncomp==-2L)  {
        myNAMES <- c("Prop","MAM4.De","Mu","Sigma")
    } else if (ncomp==-3L) {
        myNAMES <-  c("Prop","MXAM3.De","Sigma")
    } else if (ncomp==-4L) {
        myNAMES <- c("Prop","MXAM4.De","Mu","Sigma")
    } # end if.
    chains<-matrix(res$chains,ncol=ifelse(ncomp %in% c(-1L,-3L),3L,4L),
                   dimnames=list(NULL,myNAMES))
    ###
    if (iflog==TRUE) {
        if (ncomp==-1L) {
            chains[,2L]<-exp(chains[,2L])
            model <- "MAM3"
        } else if (ncomp==-3L) {
            chains[,2L]<-exp(-chains[,2L])
            model <- "MXAM3"
        } else if (ncomp==-2L) {
            chains[,2L:3L]<-exp(chains[,2L:3L])
            model <- "MAM4"
        } else if (ncomp==-4L) {
            chains[,2L:3L]<-exp(-chains[,2L:3L])
            model <- "MXAM4"
        } # end if
    } else {
        if (ncomp==-1L) {
            model <- "MAM3"
        } else if (ncomp==-2L) {
            model <- "MAM4"
        } else if (ncomp==-3L) {
            chains[,2L]<- -chains[,2L]
            model <- "MXAM3"
        } else if (ncomp==-4L) {
            chains[,2L:3L]<- -chains[,2L:3L]
            model <- "MXAM4"
        } # end if
        
    } # end if.
    ###
    output<-list("EDdata"=EDdata, 
                 "addsigma"=addsigma, 
                 "model"=model,
                 "iflog"=iflog, 
                 "nsim"=nsim, 
                 "chains"=chains)
    class(output)<-"mcAgeModels"
    invisible(output)
    ###
} # end function mcMAM.
#####
