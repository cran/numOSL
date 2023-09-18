#######
    sensSAM <- function(EDdata, model, sigmaVEC=NULL, iflog=TRUE, 
                        maxcomp=8, plot=TRUE, outfile=NULL) {

        UseMethod("sensSAM")

    } # end function sensSAM.
    ### 2023.09.10.
    sensSAM.default <- function(EDdata, model, sigmaVEC=NULL, iflog=TRUE, 
                                maxcomp=8, plot=TRUE, outfile=NULL) {
        ###
        stopifnot(ncol(EDdata)==2, nrow(EDdata)>=5,
                  all(EDdata[,2,drop=TRUE]>0),
                  length(model)==1, is.character(model),
                  is.null(sigmaVEC) || (is.numeric(sigmaVEC) && all(sigmaVEC>=0)),
                  length(iflog)==1, is.logical(iflog),
                  length(maxcomp)==1, maxcomp%in%(2:13),
                  length(plot)==1, is.logical(plot),
                  is.null(outfile) || (length(outfile)==1L && is.character(outfile)))

        ###
        fmm_name <- c("fmm0","fmm1","fmm2","fmm3","fmm4",
                      "fmm5","fmm6","fmm7","fmm8","fmm9")
        mam_name <- c("mam3","mam4","mxam3","mxam4")
        
        ###
        if (!model %in% c("com", "cam", fmm_name, mam_name)) 
            stop("Error: incorrect model!")

        ###
        if (iflog==TRUE && any(EDdata[,1,drop=TRUE]<=0)) iflog <- FALSE

        ###
        if (is.null(sigmaVEC)) {

            if (iflog==TRUE) {

                sigmaVEC <- seq(from=0,to=0.3,by=0.01)

            } else {

                xer <- mean(as.numeric(EDdata[,2,drop=TRUE]))*0.3
                sigmaVEC <- round(seq(from=0,to=xer,by=xer/29.0),3)

            } # end if.

        } # end if.

        ###
        nsigma <- length(sigmaVEC) 
        PARS <- vector(length=nsigma, mode="list")
        names(PARS) <- paste("addsigma=",sigmaVEC,sep="")

        ###
        if (model=="com") {

            nx <- 5
            cln <- c("addsigma","Com.De","Se.Com.De","Maxlik","BIC")

        } else if (model=="cam") {

            nx <- 7
            cln <- c("addsigma","CAM.De","Se.CAM.De",
                     "CAM.OD","Se.CAM.OD","Maxlik","BIC")

        } else if (model=="mam3") {

            nx <- 9
            cln <- c("addsigma","Prop","Se.Prop","MAM3.De","Se.MAM3.De",
                     "Sigma","Se.Sigma","Maxlik","BIC")

        } else if (model=="mxam3") {

            nx <- 9
            cln <- c("addsigma","Prop","Se.Prop","MXAM3.De","Se.MXAM3.De",
                     "Sigma","Se.Sigma","Maxlik","BIC")

        } else if (model=="mam4") {

            nx <- 11
            cln <- c("addsigma","Prop","Se.Prop","MAM4.De","Se.MAM4.De",
                     "Mu","Se.Mu","Sigma","Se.Sigma","Maxlik","BIC")

        } else if (model=="mxam4") {

            nx <- 11
            cln <- c("addsigma","Prop","Se.Prop","MXAM4.De","Se.MXAM4.De",
                     "Mu","Se.Mu","Sigma","Se.Sigma","Maxlik","BIC")
              
        } else {

            ###
            ncomp <- as.numeric(substr(model,start=4,stop=4))

            ###
            if (ncomp>0) {

                nx <- 3 + ncomp*4
                cln <- c("addsigma",paste(c("Prop","Se.Prop"),rep(1:ncomp,each=2),sep=""),
                         paste(c("FMM.De","Se.FMM.De"),rep(1:ncomp,each=2),sep=""),"Maxlik","BIC") 

            } else {
             
                nx <- 4
                cln <- c("addsigma", "ncomp", "Maxlik", "BIC")

            } # end if.  
            
        } # end if.

        ###
        mat <- matrix(nrow=nsigma, ncol=nx)
        colnames(mat) <- cln
        
        ###
        if (!model %in% c("com","cam","mam3","mxam3","mam4","mxam4")) {

            PropDeD <- matrix(NA,nrow=nsigma, ncol=4)

        } # end if
        
        ###
        pb <- txtProgressBar(min=1, max=nsigma, initial=1, style=3)

        ###
        for (i in 1:nsigma) {

            SAMX <- try(optimSAM(EDdata, model=model, addsigma=sigmaVEC[i], 
                                 iflog=iflog, maxcomp=maxcomp), silent=TRUE)

            ###
            if (inherits(SAMX,what="try-error")==FALSE) {

                ###
                PARS[[i]] <- SAMX$pars

                if (model=="com") {

                    mat[i,] <- c(sigmaVEC[i], SAMX$pars[1,,drop=TRUE], 
                                 SAMX$maxlik, SAMX$bic)

                } else if (model %in% c("cam","mam3","mxam3","mam4","mxam4")) {

                    mat[i,] <- c(sigmaVEC[i], c(t(SAMX$pars)), SAMX$maxlik, SAMX$bic)

                } else {

                    if (model=="fmm0") {

                        mat[i,] <- c(sigmaVEC[i], nrow(SAMX$pars), SAMX$maxlik, SAMX$bic)

                    } else {

                        mat[i,] <- c(sigmaVEC[i], c(t(SAMX$pars[,1:2])), 
                                     c(t(SAMX$pars[,3:4])), SAMX$maxlik, SAMX$bic)

                    } # end if.

                    ###
                    PropDeD[i,] <- SAMX$pars[which.max(SAMX$pars[,1]),,drop=TRUE]

                } # end if.

            } else {

                PARS[[i]] <- NA
                mat[i,] <- c(sigmaVEC[i], rep(NA,nx-1))

            } # end if.

            ###
            setTxtProgressBar(pb,i)

        } # end for.

        ###
        close(pb)

        ###
        if (plot==TRUE) {

            opar <- par("mfrow", "mgp", "mar")
            on.exit(par(opar))

            ###
            if (iflog==TRUE) {

                mylab <- expression(paste("Sigmab (",sigma[b],")"))

            } else {

                mylab <- expression(paste("Sigmab (",sigma[b],") (Gy)"))

            } # end if.

            ###
            if (model %in% c("com","fmm1")) {

                layout(rbind(c(1,2), c(3,4)),respect=TRUE)
                par("mgp"=c(2.5,1,0))
                par("mar"=c(4.1, 4.1, 3.1, 2.1))

                ###  
                offset <- ifelse(model=="com",0,2)

                ###    
                plot(mat[,1], round(mat[,2+offset],2), type="o", pch=21, bg="skyblue", 
                     xlab=mylab, ylab=ifelse(model=="com","COM [De] (Gy)","FMM [De1] (Gy)"),
                     col="skyblue", cex.axis=1.2, cex.lab=1.2, 
                     ylim=c(min(mat[,2+offset]-mat[,3+offset]/2.0,na.rm=TRUE),
                     max(mat[,2+offset]+mat[,3+offset]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,2+offset]-mat[,3+offset]/2.0, 
                                 x1=mat[,1], y1=mat[,2+offset]+mat[,3+offset]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))
          
                ###
                par("mar"=c(5.1, 4.1, 2.1, 2.1))
                plot(mat[,1], mat[,4+offset], type="o", pch=21, bg="purple", col="purple",
                     xlab=mylab, ylab="Likelihood value",cex.axis=1.2, cex.lab=1.2)

                ###
                plot(mat[,1], mat[,5+offset], type="o", pch=21, bg="red", col="red",
                     xlab=mylab, ylab="BIC value", cex.axis=1.2, cex.lab=1.2)

            } else if (model=="cam") {

                layout(rbind(c(1,2), c(3,4)),respect=TRUE)
                par("mgp"=c(2.5,1,0))
                par("mar"=c(4.1, 4.1, 3.1, 2.1))
                
                ###
                plot(mat[,1], round(mat[,2],2), type="o", pch=21, bg="skyblue", 
                     xlab=mylab, ylab="CAM [De] (Gy)", col="skyblue", cex.axis=1.2, cex.lab=1.2,
                     ylim=c(min(mat[,2]-mat[,3]/2.0,na.rm=TRUE),max(mat[,2]+mat[,3]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,2]-mat[,3]/2.0, 
                                 x1=mat[,1], y1=mat[,2]+mat[,3]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ###
                plot(mat[,1], mat[,4], type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=ifelse(iflog==TRUE, "CAM [OD]", "CAM [OD] (Gy)"),
                     ylim=c(min(mat[,4]-mat[,5]/2.0,na.rm=TRUE),max(mat[,4]+mat[,5]/2.0,na.rm=TRUE)),
                     cex.axis=1.2, cex.lab=1.2)

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,4]-mat[,5]/2.0, 
                                 x1=mat[,1], y1=mat[,4]+mat[,5]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))
          
                ###
                plot(mat[,1], mat[,6], type="o", pch=21, bg="purple", col="purple",
                     xlab=mylab, ylab="Likelihood value",cex.axis=1.2, cex.lab=1.2)

                ###
                plot(mat[,1], mat[,7], type="o", pch=21, bg="red", col="red",
                     xlab=mylab, ylab="BIC value",cex.axis=1.2, cex.lab=1.2)   

            } else if (model %in% c("mam3","mxam3")) {

                layout(rbind(c(1,2),c(3,4),c(5,6)),respect=TRUE)
                par("mgp"=c(2.5,1,0))
                par("mar"=c(4.1, 4.1, 2.1, 1.1))
                
                ### Figure 1.
                ###------------------------------------------------------------------------
                plot(mat[,1], round(mat[,2],2), type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=ifelse(model=="mam3","MAM3 [Prop]", "MXAM3 [Prop]"), 
                     cex.axis=1.2, cex.lab=1.2, ylim=c(min(mat[,2]-mat[,3]/2.0,na.rm=TRUE),
                     max(mat[,2]+mat[,3]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,2]-mat[,3]/2.0, 
                                 x1=mat[,1], y1=mat[,2]+mat[,3]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 2.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,4], type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=ifelse(model=="mam3","MAM3 [De] (Gy)","MXAM3 [De] (Gy)"), cex.axis=1.2,
                     cex.lab=1.2, ylim=c(min(mat[,4]-mat[,5]/2.0,na.rm=TRUE),max(mat[,4]+mat[,5]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,4]-mat[,5]/2.0, 
                                 x1=mat[,1], y1=mat[,4]+mat[,5]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 3.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,6], type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=paste(ifelse(model=="mam3","MAM3 ","MXAM3 "),
                     ifelse(iflog==TRUE, "[Sigma]", "[Sigma] (Gy)"),sep=""),
                     cex.axis=1.2, cex.lab=1.2, ylim=c(min(mat[,6]-mat[,7]/2.0,na.rm=TRUE),
                     max(mat[,6]+mat[,7]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,6]-mat[,7]/2.0, 
                                 x1=mat[,1], y1=mat[,6]+mat[,7]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))
          
                ### Figure 4.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,8], type="o", pch=21, bg="purple", col="purple",
                     xlab=mylab, ylab="Likelihood value",cex.axis=1.2, cex.lab=1.2)

                ### Figure 5.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,9], type="o", pch=21, bg="red", col="red",
                     xlab=mylab, ylab="BIC value",cex.axis=1.2, cex.lab=1.2)   

            } else if (model %in% c("mam4","mxam4")) {
                
                layout(rbind(c(1,2),c(3,4),c(5,6)),respect=TRUE)
                par("mgp"=c(2.5,1,0))
                par("mar"=c(4.1, 4.1, 2.1, 1.1))
                
                ### Figure 1.
                ###------------------------------------------------------------------------
                plot(mat[,1], round(mat[,2],2), type="o", pch=21, bg="skyblue", 
                     xlab=mylab, ylab=ifelse(model=="mam4","MAM4 [Prop]","MXAM4 [Prop]"), 
                     col="skyblue",cex.axis=1.2, cex.lab=1.2, ylim=c(min(mat[,2]-mat[,3]/2.0,na.rm=TRUE),
                     max(mat[,2]+mat[,3]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,2]-mat[,3]/2.0, 
                                 x1=mat[,1], y1=mat[,2]+mat[,3]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 2.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,4], type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=ifelse(model=="mam4","MAM4 [De] (Gy)","MXAM4 [De] (Gy)"),cex.axis=1.2, 
                     cex.lab=1.2, ylim=c(min(mat[,4]-mat[,5]/2.0,na.rm=TRUE),max(mat[,4]+mat[,5]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,4]-mat[,5]/2.0, 
                                 x1=mat[,1], y1=mat[,4]+mat[,5]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 3.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,6], type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=ifelse(model=="mam4","MAM4 [Mu] (Gy)","MXAM4 [Mu] (Gy)"), 
                     cex.axis=1.2, cex.lab=1.2, ylim=c(min(mat[,6]-mat[,7]/2.0,na.rm=TRUE),
                     max(mat[,6]+mat[,7]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,6]-mat[,7]/2.0, 
                                 x1=mat[,1], y1=mat[,6]+mat[,7]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 4.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,8], type="o", pch=21, bg="skyblue", col="skyblue",
                     xlab=mylab, ylab=paste(ifelse(model=="mam4","MAM4 ","MAM4 "),
                     ifelse(iflog==TRUE, "[Sigma]", "[Sigma] (Gy)"),sep=""), cex.axis=1.2, cex.lab=1.2,
                     ylim=c(min(mat[,8]-mat[,9]/2.0,na.rm=TRUE),max(mat[,8]+mat[,9]/2.0,na.rm=TRUE)))

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=mat[,8]-mat[,9]/2.0, 
                                 x1=mat[,1], y1=mat[,8]+mat[,9]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))
          
                ### Figure 5.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,10], type="o", pch=21, bg="purple", col="purple",
                     xlab=mylab, ylab="Likelihood value", cex.axis=1.2, cex.lab=1.2)

                ### Figure 6.
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,11], type="o", pch=21, bg="red", col="red",
                     xlab=mylab, ylab="BIC value", cex.axis=1.2, cex.lab=1.2)   

            } else if (model=="fmm0") {

                layout(rbind(c(1,2),c(3,4),c(5,6)),respect=TRUE)
                par("mgp"=c(2.5,1,0))
                par("mar"=c(4.1, 4.1, 2.1, 1.1))

                ###
                idx <- which.min(mat[,4])

                ### Figure 1.
                ###------------------------------------------------------------------------
                plot(mat[,1], PropDeD[,1], type="o", pch=21, bg="skyblue", xlab=mylab,  
                     ylab="FMM [Prop_D]", col="skyblue", cex.axis=1.2, cex.lab=1.2,
                     ylim=c(min(PropDeD[,1]-PropDeD[,2]/2.0,na.rm=TRUE),
                            max(PropDeD[,1]+PropDeD[,2]/2.0,na.rm=TRUE)))
                abline(v=mat[idx,1], lwd=2, col="green")

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=PropDeD[,1]-PropDeD[,2]/2.0, 
                                 x1=mat[,1], y1=PropDeD[,1]+PropDeD[,2]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 2.
                ###------------------------------------------------------------------------
                plot(mat[,1], PropDeD[,3], type="o", pch=21, bg="skyblue", xlab=mylab,  
                     ylab="FMM [De_D] (Gy)", col="skyblue", cex.axis=1.2, cex.lab=1.2,
                     ylim=c(min(PropDeD[,3]-PropDeD[,4]/2.0,na.rm=TRUE),
                            max(PropDeD[,3]+PropDeD[,4]/2.0,na.rm=TRUE)))
                abline(v=mat[idx,1], lwd=2, col="green")

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=PropDeD[,3]-PropDeD[,4]/2.0, 
                                 x1=mat[,1], y1=PropDeD[,3]+PropDeD[,4]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))
                
                ### Figure 3.  
                ###------------------------------------------------------------------------  
                plot(mat[,1], mat[,2], type="o", pch=21, bg="brown", xlab=mylab,  
                     ylab="FMM [k]", col="brown", cex.axis=1.2, cex.lab=1.2)
                abline(v=mat[idx,1], lwd=2, col="green")
                
                
                ### Figure 4.  
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,3], type="o", pch=21, bg="purple", col="purple",
                     xlab=mylab, ylab="Likelihood value",cex.axis=1.2, cex.lab=1.2)
                abline(v=mat[idx,1], lwd=2, col="green")

                ### Figure 5.  
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,4], type="o", pch=21, bg="red", col="red",
                     xlab=mylab, ylab="BIC value", cex.axis=1.2, cex.lab=1.2)
                abline(v=mat[idx,1], lwd=2, col="green")

                ### Figure 5.  
                ###------------------------------------------------------------------------
                par("mar"=c(1,1,1,1))
                plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
                legend("center", legend=c(paste("Sigmab range: [",round(min(sigmaVEC),2)," ,",round(max(sigmaVEC),2),"]",sep=""), 
                       paste("Optimal Sigmab: ",round(mat[idx,1],2),ifelse(iflog==TRUE,""," (Gy)"),sep=""),
                       paste("Optimal k: ",mat[idx,2],sep=""),paste("Minimized BIC: ", round(mat[idx,4],2),sep=""), 
                       paste("Dominant Prop: ",round(PropDeD[idx,1],2),"+-",round(PropDeD[idx,2],2),sep=""),
                       paste("Dominant De: ",round(PropDeD[idx,3],2),"+-",round(PropDeD[idx,4],2)," (Gy)",sep="")),
                       yjust=2, ncol=1,bty="n")

            } else {

                layout(rbind(c(1,2), c(3,4)),respect=TRUE)
                par("mgp"=c(2.5,1,0))
                par("mar"=c(4.1, 4.1, 3.1, 2.1))

                ###
                idx <- which.min(mat[,nx])

                ### Figure 1.
                ###------------------------------------------------------------------------
                plot(mat[,1], PropDeD[,1], type="o", pch=21, bg="skyblue", xlab=mylab,  
                     ylab="Dominant Prop", col="skyblue", cex.axis=1.2, cex.lab=1.2,
                     ylim=c(min(PropDeD[,1]-PropDeD[,2]/2.0,na.rm=TRUE),
                            max(PropDeD[,1]+PropDeD[,2]/2.0,na.rm=TRUE)))
                ###abline(v=mat[idx,1], lwd=2, col="green")

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=PropDeD[,1]-PropDeD[,2]/2.0, 
                                 x1=mat[,1], y1=PropDeD[,1]+PropDeD[,2]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))

                ### Figure 2.
                ###------------------------------------------------------------------------
                plot(mat[,1], PropDeD[,3], type="o", pch=21, bg="skyblue", xlab=mylab,  
                     ylab="Dominant De (Gy)", col="skyblue", cex.axis=1.2, cex.lab=1.2,
                     ylim=c(min(PropDeD[,3]-PropDeD[,4]/2.0,na.rm=TRUE),
                            max(PropDeD[,3]+PropDeD[,4]/2.0,na.rm=TRUE)))
                ###abline(v=mat[idx,1], lwd=2, col="green")

                ###
                suppressWarnings(arrows(x0=mat[,1], y0=PropDeD[,3]-PropDeD[,4]/2.0, 
                                 x1=mat[,1], y1=PropDeD[,3]+PropDeD[,4]/2.0, code=3, 
                                 lwd=1, angle=90, length=0.05, col="skyblue"))
                         
                ### Figure 3.  
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,nx-1], type="o", pch=21, bg="purple", col="purple",
                     xlab=mylab, ylab="Likelihood value",cex.axis=1.2, cex.lab=1.2)
                ###abline(v=mat[idx,1], lwd=2, col="green")
                

                ### Figure 4.  
                ###------------------------------------------------------------------------
                plot(mat[,1], mat[,nx], type="o", pch=21, bg="red", col="red",
                     xlab=mylab, ylab="BIC value", cex.axis=1.2, cex.lab=1.2)
                ###abline(v=mat[idx,1], lwd=2, col="green")  

            } # end if.

        } # end if.
        
        ###
        if (!is.null(outfile)) write.csv(mat, file=paste(outfile, ".csv",sep=""))

        ###
        output <- list("pars"=PARS,"mat"=mat)   
        return(invisible(output))

    } # end function sensSAM.default.
    ###sensSAM(EDdata$al3, model="mam3", iflog=TRUE)
    ###
