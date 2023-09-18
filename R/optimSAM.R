#######
    optimSAM <- function(EDdata, model="cam", addsigma=0, iflog=TRUE, maxcomp=6) {

        UseMethod("optimSAM")

    } # end function optimSAM.
    ### 2023.09.10.
    optimSAM.default <- function(EDdata, model="cam", addsigma=0, iflog=TRUE, maxcomp=6) { 

        stopifnot(ncol(EDdata)==2, nrow(EDdata)>=5,
                  all(EDdata[,2,drop=TRUE]>0),
                  length(model)==1, is.character(model),
                  length(addsigma)==1, addsigma>=0,
                  length(iflog)==1, is.logical(iflog),
                  length(maxcomp)==1, maxcomp%in%(2:13))

        ###
        ed <- as.numeric(EDdata[,1,drop=TRUE])
        sed <- as.numeric(EDdata[,2,drop=TRUE])

        if (iflog==TRUE && any(ed<=0)) {

            iflog <- FALSE

            cat("Note: zero or negative De values, set iflog=FALSE!\n")

        } # end if.

        ###
        fmm_name <- c("fmm0","fmm1","fmm2","fmm3","fmm4",
                      "fmm5","fmm6","fmm7","fmm8","fmm9")
        mam_name <- c("mam3","mam4","mxam3","mxam4")
        
        ###
        if (!model %in% c("com", "cam", fmm_name, mam_name)) 
            stop("Error: incorrect model!")

        ### Optimize the common age model.
        ###-----------------------------------------------------------
        comED <- function(ed, sed, addsigma, iflog) {

            if (iflog==TRUE) {

                sigma2 <- (sed/ed)^2+addsigma^2
                d <- log(ed)

            } else {

                sigma2 <- sed^2+addsigma^2
                d <- ed

            } # end if.

            ###
            w <- 1.0/sigma2
            v1 <- sum(w*d)/sum(w)
            v2 <- sqrt(1.0/sum(w))

            ###
            maxlik <- sum(log(1.0/sqrt(2.0*pi)/sqrt(sigma2)*exp(-(v1-d)^2/2.0/sigma2)))
            bic <- -2.0*maxlik + log(length(ed))

            ###
            if (iflog==TRUE) {

                pars <- matrix(c(exp(v1),exp(v1)*v2), nrow=1)

            } else {

                pars <- matrix(c(v1,v2), nrow=1)

            } # end if.

            ###
            rownames(pars) <- "COM.De"
            colnames(pars) <- c("Pars", "Se.Pars")
            
            ###
            output <- list("pars"=pars, "bic"=bic, "maxlik"=maxlik)
            return(output)

        } # end funciton comED.
        ###

        ### Optimize the central age model. 
        ###-----------------------------------------------------------
        camED <- function(ed, sed, addsigma=0, iflog=TRUE) {

            if (iflog==TRUE) {

                sz <- sqrt((sed/ed)^2+addsigma^2)
                z <- log(ed)

            } else {

                sz <- sqrt(sed^2+addsigma^2)
                z <- ed

            } # end if.

            ###
            sigma <- 0.1
            wz <- 1.0/(sigma^2+sz^2)
            mu <- sum(z*wz)/sum(wz)

            ###
            repeat {

                mu1 <- sum(z*wz)/sum(wz)
                sigma1 <- sigma*sqrt(sum(wz^2*(z-mu1)^2/sum(wz)))
                wz <- 1.0/(sigma1^2+sz^2)

                ###
                if (abs(mu-mu1)+abs(sigma-sigma1)<1.0e-8) break
                
                ###
                sigma <- sigma1
                mu <- mu1

            } # end repeat.

            ###
            maxlik <- sum(log(1.0/sqrt(2.0*pi)*sqrt(wz)*exp(-(z-mu)^2*wz/2.0)))
            bic <- -2.0*maxlik + 2.0*log(length(ed))

            ###
            se_mu <- 1/sqrt(sum(wz))
            se_sigma <- 1.0/sqrt(2.0*mu*sum(wz^2))
            
            ###
            if (iflog==TRUE) {

                mu <- exp(mu)
                se_mu <- mu*se_mu

            } # end if.

            ###
            pars <- cbind(c(mu,sigma),c(se_mu,se_sigma)) 
            rownames(pars) <- c("CAM.De","CAM.OD")
            colnames(pars) <- c("Pars", "Se.Pars")

            ###
            output <- list("pars"=pars, "maxlik"=maxlik, "bic"=bic)
            return(output)

        } # end function camED.  
        ###camED(EDdata$gl1[,1], EDdata$gl1[,2], addsigma=0, iflog=TRUE)


        ### Calculate the standard errors of parameters  
        ### of the finite mixture age model.
        ###---------------------------------------------------------------------
        apfmmstd <- function(p, mu, z, s) {

            ###
            w <- 1.0/s^2

            ### 
            ndat <- length(z)
            ncomp <- length(p)

            ###
            pdf <- matrix(nrow=ndat, ncol=ncomp)
            for (j in 1:ncomp) {

                pdf[,j] <- p[j]*sqrt(w)*exp(-0.5*w*(z-mu[j])^2)

            } # end for.

            ###
            prob <- pdf/rowSums(pdf)

            ###
            aaMat <- bbMat <- matrix(0, nrow=ndat, ncol=ncomp)
            for (j in 1:ncomp) {

                aaMat[,j] <- w*(z-mu[j])
                bbMat[,j] <- -w + (w*(z-mu[j]))^2

            } # end for.
            
            ###
            aMat <- matrix(0, nrow=ncomp-1, ncol=ncomp-1)  
            for (i in 1:(ncomp-1)) {

                for (j in 1:(ncomp-1)) {

                    aMat[i,j] <- sum((prob[,i]/p[i]-prob[,ncomp]/p[ncomp])*
                                     (prob[,j]/p[j]-prob[,ncomp]/p[ncomp]))

                } # end for.

            } # end for.

            ###
            bMat <- matrix(0,nrow=ncomp-1,ncol=ncomp)
            iMat <- diag(ncomp)
            for (i in 1:(ncomp-1)) {

                for (j in 1:ncomp) {

                    bMat[i,j] <- sum(prob[,j]*aaMat[,j]*(prob[,i]/p[i]-
                                     prob[,ncomp]/p[ncomp]-iMat[i,j]/p[i]+
                                     iMat[ncomp,j]/p[ncomp]))
                
                } # end for.

            } # end for.

            ###
            cMat <- matrix(0,nrow=ncomp,ncol=ncomp)
            for (i in 1:ncomp) {

                for (j in 1:ncomp) {

                    cMat[i,j] <- sum(prob[,i]*prob[,j]*aaMat[,i]*aaMat[,j]-
                                     iMat[i,j]*bbMat[,i]*prob[,i])

                } # end for.

            } # end for.

            ###
            covar <- matrix(nrow=2*ncomp-1,ncol=2*ncomp-1)
            covar[1:(ncomp-1),1:(ncomp-1)] <- aMat
            covar[1:(ncomp-1),ncomp:(2*ncomp-1)] <- bMat
            covar[ncomp:(2*ncomp-1),1:(ncomp-1)] <-  t(bMat)
            covar[ncomp:(2*ncomp-1),ncomp:(2*ncomp-1)] <- cMat

            ###
            invcovar <- try(solve(covar,tol=.Machine$double.xmin),silent=TRUE)

            ###
            if (inherits(invcovar,what="try-error")==FALSE) {

                if (all(diag(invcovar)>=0.0) && sum(invcovar[1:(ncomp-1),1:(ncomp-1)])>=0.0) {

                    sep <- c(sqrt(diag(invcovar[1:(ncomp-1),1:(ncomp-1),drop=FALSE])), 
                             sqrt(sum(invcovar[1:(ncomp-1),1:(ncomp-1)])))
                    semu <- sqrt(diag(invcovar[ncomp:(2*ncomp-1),ncomp:(2*ncomp-1),drop=FALSE]))

                } else {

                    stop("Error: negative value in the diagonal of the inverse Hessian matrix!")

                } # end if.

            } else {

                stop("Error: the Hessian matrix is a singular matrix!")

            } # end if.

            ###
            list("sep"=sep, "semu"=semu)

        } # end function apfmmstd.
            
        ###
        ### Optimize the finite mixture age model.
        ###-----------------------------------------------------------------
        fmmED <- function(ed, sed, ncomp, addsigma=0, iflog=TRUE) {

            ###
            if (ncomp==1)    {

                output <- comED(ed=ed, sed=sed, addsigma=addsigma, iflog=iflog)
                pars <- cbind(matrix(c(1.0,0.0),nrow=1), output$pars)
                rownames(pars) <- "Comp1"
                colnames(pars) <- c("Prop","Se.Prop","De","Se.De") 
                output$pars <- pars

                ###
                return(output)

            } # end if.

            ###
            if (iflog==TRUE) {

                sed <- sqrt((sed/ed)^2+addsigma^2)
                ed <- log(ed)

            } else {

                sed <- sqrt(sed^2+addsigma^2)

            } # end if.

            ###
            w <- 1.0/sed^2

            ### 
            ndat <- length(ed)
            
            ###
            pf <- matrix(nrow=ndat, ncol=ncomp)
            ###muvec <- min(ed) +(max(ed)-min(ed))/(ncomp+3)*((1:(ncomp+4))-1)
            muvec <- min(ed) +(max(ed)-min(ed))/(ncomp+4)*((1:(ncomp+5))-1)
            maxval <- -1.0e20

            ###
            info <- 1
            for (kk in 1:6) {

                inip <- rep(1.0/ncomp, ncomp)
                inimu <- muvec[kk:(kk+ncomp-1)]
                ###
                repeat {

                    for (j in 1:ncomp) {

                        pf[,j] <- inip[j]*sqrt(w)*exp(-0.5*w*(ed-inimu[j])^2)

                    } # end for.

                    ###
                    if (any(!is.finite(pf))) 
                        stop("Error: Improper initial parameters!")

                    ###
                    pp <- pf/rowSums(pf)
                    wp <- w*pp
                    pv <- pp*w*ed

                    ###
                    p1 <- colSums(pp)/ndat
                    mu1 <- colSums(pv)/colSums(wp)

                    ###
                    if (sum(abs(inip-p1))+sum(abs(inimu-mu1))<=1.0e-8) break

                    ###
                    inip <- p1
                    inimu <- mu1

                } # end repeat.

                ###
                maxlik <- sum(log(1.0/sqrt(2.0*pi)*rowSums(pf)))
                bic <- -2.0*maxlik + (2.0*ncomp-1.0)*log(ndat)
                SSEE <- try(apfmmstd(p=inip, mu=inimu, z=ed, s=sed), silent=TRUE)

                ###
                if (maxlik>maxval && inherits(SSEE,what="try-error")==FALSE) {

                    cp <- inip
                    csep <- SSEE$sep
                    cmu <- inimu
                    csemu <- SSEE$semu
                    cmaxlik <- maxlik
                    cbic <- bic
                    info <- 0
                    maxval <- maxlik
                
                } # end if.

            } # end for.

            ###
            if (info==1) stop("Error: optimization of FMM failed!")
            
            ###
            p <- cp
            sep <- csep
            mu <- cmu
            semu <- csemu
            maxlik <- cmaxlik
            bic <- cbic
            
            ###
            if (iflog==TRUE) {

                mu <- exp(mu)
                semu <- mu*semu
              
            } # end if.

            ###
            idx <- order(mu)
            pars <- cbind(p[idx],sep[idx],mu[idx],semu[idx])
            colnames(pars) <- c("Prop","Se.Prop","FMM.De","Se.FMM.De") 
            rownames(pars)<-paste(rep("Comp", ncomp), seq(ncomp), sep="") 

            ###
            output <- list("pars"=pars, "maxlik"=maxlik, "bic"=bic)
            return(output)
        } # end function fmmED. 
        ###fmmED(EDdata$al3[,1], EDdata$al3[,2], ncomp=3, iflog=TRUE,addsigma=0.05)

        ###
        ### Optimize the minimum age model.
        ###----------------------------------------------------------------------------
        mamED <- function(ed, sed, ncomp, addsigma=0, iflog=TRUE) {

            ###
            if (iflog==TRUE) {

                x <- sqrt((sed/ed)^2+addsigma^2)
                y <- log(ed)

            } else {
             
                x <- sqrt(sed^2+addsigma^2)
                y <- ed

            } # end if.

            ###
            ndat <- length(y)
       
            ### 
            ### Function to be minimized.
            minfunc <- function(p)  {

                if (ncomp %in% c(-1L,-3L)) {

                    ### For MAM3 and MXAM3.
                    u0 <- (p[2L]/(p[3L])^2L+y/x^2L)/(1.0/(p[3L])^2L+1.0/x^2L)
                    sigma0 <- 1.0/(sqrt(1.0/(p[3L])^2L+1.0/x^2L)) 
                    prop1 <- p[1L]/sqrt(2.0*pi)/x*exp(-(y-p[2L])^2L/2.0/x^2L)

                    ###
                    if (ncomp==-1L) {

                        prop2 <- (1.0-p[1L])/sqrt(2.0*pi*((p[3L])^2L+x^2L))* 
                                  2.0*(1.0-pnorm((p[2L]-u0)/sigma0))*
                                  exp(-(y-p[2L])^2L/2.0/((p[3L])^2L+x^2L))

                    } else {

                        prop2 <- (1.0-p[1L])/sqrt(2.0*pi*((p[3L])^2L+x^2L))* 
                                  2.0*(pnorm((p[2L]-u0)/sigma0))*
                                  exp(-(y-p[2L])^2L/2.0/((p[3L])^2L+x^2L))

                    } # end if.

                    ###
                    sumlog <- -sum(log(prop1+prop2))
                    if (is.finite(sumlog)) return(sumlog) else return(1.0e20)

                } else if (ncomp %in% c(-2L,-4L)) {

                    ### For MAM4 and MXAM4.
                    u0 <- (p[3L]/(p[4L])^2L+y/x^2L)/(1.0/(p[4L])^2L+1.0/x^2L)          
                    sigma0 <- 1.0/(sqrt(1.0/(p[4L])^2L+1.0/x^2L)) 
                    prop1 <- p[1L]/sqrt(2.0*pi)/x*exp(-(y-p[2L])^2L/2.0/x^2L)
 
                    ###
                    if (ncomp==-2L) {

                        prop2 <- (1.0-p[1L])/sqrt(2.0*pi*((p[4L])^2L+x^2L))* 
                                 (1.0-pnorm((p[2L]-u0)/sigma0))/(1.0-pnorm((p[2L]-p[3L])/p[4L]))*
                                 exp(-(y-p[3L])^2L/2.0/((p[4L])^2L+x^2L))

                    } else {

                        prop2 <- (1.0-p[1L])/sqrt(2.0*pi*((p[4L])^2L+x^2L))* 
                                 (pnorm((p[2L]-u0)/sigma0))/(pnorm((p[2L]-p[3L])/p[4L]))*
                                 exp(-(y-p[3L])^2L/2.0/((p[4L])^2L+x^2L))

                    } # end if.

                    ###
                    sumlog <- -sum(log(prop1+prop2))
                    if (is.finite(sumlog)) return(sumlog) else return(1.0e20)

                } # end if.
           
            } # end function minfunc.

            ###
            ### Set boundaries.
            if (ncomp %in% c(-1L,-3L)) {

                lower <- c(1e-4, min(y), 1e-3)
                upper <- c(1-1e-4, max(y), ifelse(iflog==TRUE,5.0,50.0))

            } else if (ncomp %in% c(-2L,-4L)) {

                lower <- c(1e-4, min(y), min(y), 1e-3)
                upper <- c(1-1e-4, max(y), max(y), ifelse(iflog==TRUE,5.0,50.0))

            } # end if
  
            ### 
            kclus <- kmeans(x=y, centers=3L, iter.max=50L, nstart=100L)
            kclus <- sort(kclus$centers)
        
            ###
            ### Set gamma, mu and sdy.
            gama <- kclus[1L]
            mu <- c(kclus[2L], mean(kclus), mean(y))
            sdy <- sd(y)
        
            ### 
            cmaxlik <- maxlik <- 1.0e20
            errorflag <- 1
            bexist <- 0

            ###
            ### Do optimization with various initial values.
            if (ncomp %in% c(-1L,-3L)) {

                ### For MAM3 and MXAM3.
                for(i in seq(3L)) {

                    for(j in seq(5L)) {

                        for (k in seq(3L)) {

                            ###
                            ini <- c(0.01*(5.0)^(i-1L), mu[k], sdy*0.4*j)

                            ###
                            res <- suppressWarnings(try(optim(par=ini, fn=minfunc, gr=NULL, 
                                   method="L-BFGS-B", control=list(maxit=1000), lower=lower,
                                   upper=upper, hessian=TRUE), silent=TRUE))
                       
                            ###
                            if (inherits(res, what="try-error")==FALSE && res$convergence==0) { 
                            
                                ###
                                min.obj <- res$value
                            
                                ###
                                stdp <- suppressWarnings(try(sqrt(diag(solve(res$hessian,
                                        tol=.Machine$double.xmin))),silent=TRUE))
                              
                                ###
                                if (inherits(stdp,what="try-error")==FALSE && all(is.finite(stdp)) && 
                                    min.obj<maxlik && all(abs(res$par-lower)>=1e-5) && 
                                    all(abs(res$par-upper)>=1e-5)) {

                                    pars <- res$par
                                    error <- stdp
                                    maxlik <- min.obj
                                    errorflag <- 0

                                } # end if
                            
                                ###
                                if (inherits(stdp,what="try-error")==FALSE && 
                                    all(is.finite(stdp)) && min.obj<cmaxlik && errorflag==1) {  
                  
                                    cpars <- res$par
                                    cerror <- stdp
                                    cmaxlik <- min.obj
                                    bexist <- 1

                                } # end if
                            
                            } # end if

                        } # end if
                    
                    } # end for
               
                } # end for
            
            } else if (ncomp %in% c(-2L,-4L)) {
            
                ### For MAM4 and MXAM4.
                for (i in seq(3L)) {

                    for (j in seq(3L)) {

                        for (k in seq(3L)) {
                        
                            ###
                            if (ncomp==-2L) {

                                ### MAM4.
                                ini <- c(0.01*(5.0)^(i-1L), gama, mu[k], sdy*0.5*j)

                            } else if (ncomp==-4L) {

                                ### MXAM4.
                                ini <- c(0.01*(5.0)^(i-1L), mu[k], gama, sdy*0.5*j)

                            } # end if.

                            ###
                            res <- suppressWarnings(try(optim(par=ini, fn=minfunc, gr=NULL, 
                                   method="L-BFGS-B", control=list(maxit=1000), lower=lower,
                                   upper=upper, hessian=TRUE), silent=TRUE))
                            
                            ###
                            if (inherits(res,what="try-error")==FALSE && res$convergence==0) { 
                            
                                ###
                                min.obj <- res$value
                                
                                ###
                                stdp <- suppressWarnings(try(sqrt(diag(solve(res$hessian))),silent=TRUE))
                                
                                ###
                                if (inherits(stdp,what="try-error")==FALSE && all(is.finite(stdp))  &&   
                                    min.obj<maxlik && all(abs(res$par-lower)>=1e-5) && 
                                    all(abs(res$par-upper)>=1e-5)) { 
     
                                    pars <- res$par
                                    error <- stdp
                                    maxlik <- min.obj
                                    errorflag <- 0

                                } # end if
                            
                                ###
                                if (inherits(stdp,what="try-error")==FALSE && all(is.finite(stdp)) &&                
                                    min.obj<cmaxlik && errorflag==1) {  

                                    cpars <- res$par
                                    cerror <- stdp
                                    cmaxlik <- min.obj
                                    bexist <- 1

                                } # end if
                            
                            } # end if
                        
                        } # end for
                   
                    } # end for
               
                } # end for
            
            } # end if
        
            ###
            ### Output the results.
            if (errorflag==0) {
            
                ###
                maxlik <- -maxlik

                ### Reset parameters before output.
                if (ncomp %in% c(-1L,-3L)) {

                    if (iflog==TRUE) {

                        pars[2L] <- exp(pars[2L])
                        error[2L] <- pars[2L]*error[2L]

                    } # end if.
                
                    ###
                    bic <- -2.0*maxlik+3L*log(ndat)
                
                } else if (ncomp %in% c(-2L,-4L)) {

                    if (iflog==TRUE) {

                        pars[2L:3L] <- exp(pars[2L:3L])
                        error[2L:3L] <- pars[2L:3L]*error[2L:3L]

                    } # end if.

                    ###
                    bic <- -2.0*maxlik+4L*log(ndat)

                } # end if
            
                ###
                pars <- cbind(pars, error)
            
            } else if (bexist==1) {

                ###
                maxlik <- -cmaxlik

                ### Reset parameters before output.
                if (ncomp %in% c(-1L,-3L)) {

                    if (iflog==TRUE) {

                        cpars[2L] <- exp(cpars[2L])
                        cerror[2L] <- cpars[2L]*cerror[2L]

                    } # end if.

                    ###
                    bic <- -2.0*maxlik+3L*log(ndat)

                } else if (ncomp %in% c(-2L,-4L)) {

                    if (iflog==TRUE) {

                        cpars[2L:3L] <- exp(cpars[2L:3L])
                        cerror[2L:3L] <- cpars[2L:3L]*cerror[2L:3L]

                    } # end if.

                    ###
                    bic <- -2.0*maxlik+4L*log(ndat)

                } # end if.
         
                ###
                pars <- cbind(cpars, cerror)
            
            } else {

               stop("Error: optimization of MAM/MXAM failed!")

            } # end if.
        
            ###
            colnames(pars) <- c("Pars", "Se.Pars")
            if (ncomp==-1L) {

                rownames(pars) <- c("Prop", "MAM3.De", "Sigma")

            } else if (ncomp==-3L) {

                rownames(pars) <- c("Prop", "MXAM3.De", "Sigma")

            } else if (ncomp==-2L) {

                rownames(pars) <- c("Prop", "MAM4.De", "Mu", "Sigma")

            } else if (ncomp==-4L) {

                rownames(pars) <- c("Prop", "MXAM4.De", "Mu", "Sigma")

            } # end if.
            
            ###
            output <- list("pars"=pars, "maxlik"=maxlik, "bic"=bic)
            return(output)

        } # end function mamED.
        #mamED(EDdata$al3[,1], EDdata$al3[,2], ncomp=-1, addsigma=0, iflog=TRUE)
 
        ###
        ### START THE DRIVE FUNCTION.
        ###----------------------------------------------------------------------------------------------
        if (model=="com") {

            SAM <- comED(ed=ed, sed=sed, addsigma=addsigma, iflog=iflog)
            return(SAM)

        } else if (model=="cam")  { 

            SAM <- camED(ed=ed, sed=sed, addsigma=addsigma, iflog=iflog)
            return(SAM)

        } else if (model %in% fmm_name) {

            if (model=="fmm0") {
                
                minval <- 1.0e20
                lbmat <- matrix(NA,nrow=maxcomp,ncol=3)
 
                ###
                info <- 1
                for (i in 1:maxcomp) {

                    ###
                    lbmat[i,1] <- i

                    ###
                    SAM <- try(fmmED(ed=ed, sed=sed, ncomp=i,   
                           addsigma=addsigma, iflog=iflog), silent=TRUE)
                    
                    ###
                    if (inherits(SAM,what="try-error")==FALSE) {

                        ###
                        lbmat[i,2:3] <- c(SAM$maxlik, SAM$bic) 

                        ###
                        if (SAM$bic<minval) {

                            SAM1 <- SAM       
                            minval <- SAM$bic
                            info <- 0

                        } # end if.

                    } # end if.

                } # end for.
                       
                ###
                if (info==1) stop("Error: automatically find the optimal FMM failed!")

                ###
                colnames(lbmat) <- c("ncomp","Maxlik","BIC")
                SAM1$lbmat <- lbmat

                ###
                return(SAM1)

            } else {

                ncomp <- as.numeric(substr(model,start=4,stop=4))              
                SAM <- fmmED(ed=ed, sed=sed, ncomp=ncomp, addsigma=addsigma, iflog=iflog)
                return(SAM)

            } # end if.

        } else if (model %in% mam_name) {

            ncomp <- -seq(mam_name)[mam_name==model]
            SAM <- mamED(ed=ed, sed=sed, ncomp=ncomp, addsigma=addsigma, iflog=iflog)
            return(SAM)

        } # end if.

    } # end function optimSAM.default.
    ###optimSAM(EDdata$al3, model="fmm0", addsigma=0, iflog=FALSE, maxcomp=6)
