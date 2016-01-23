subroutine calED(dose,ltx,sltx,ndat,ninltx,n2,inltx,&
                 outDose,mcED,pars,stdp,upb,model,uw,&
                 method,nstart,nsim,fvec1,fmin,message)
!--------------------------------------------------------------------
! Subroutine calED() is used for calculating  
! ED values and assess their standard errors.
!--------------------------------------------------------------------
!        dose(ndat):: input, real values, dose values.
!         ltx(ndat):: input, real values, Lx/Tx values
!        sltx(ndat):: input, real values, errors of Lx/Txs.
!              ndat:: input, integer, number of data points.
!            ninltx:: input, integer, number of ED values.
!                n2:: input, integer, number of pars (>=1).
!   inltx(ninltx,2):: input, real values, natural Lx/Txs and errors.
! outDose(ninltx,2):: output, real values, estimated EDs and errors.
! mcED(nsim,ninltx):: output, real values, simulated ED values.
!          pars(n2):: output, real values, pars of growth curve.
!          stdp(n2):: output, real values, std errors of pars.
!               upb:: input, real value, upper bound of b value.
!             model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp.
!                uw:: input, integer, 0=un-weighted, 1=weighted.
!            method:: input, integer, 0=sp, 1=mc.
!            nstart:: input, integer, number of trails.
!              nsim:: input, integer, number of simulations.
!       fvec1(ndat):: output, real values, fitted Lx/Tx values.
!              fmin:: output, real value, minimized objective.
!           message:: output, integer, 0=success, 1=fail.
!--------------------------------------------------------------------
! Author:: Peng Jun, 2014.10.01; revised in 2016.01.21.
!--------------------------------------------------------------------
! Dependence:: subroutine linefit; subroutine inipars;---------------
!              subroutine lmfit; subroutine interpolate;-------------
!              subroutine r8vec_normal.------------------------------
!--------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, ninltx, n2, model,&
                                  uw, method, nstart, nsim
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat), sltx(ndat),&
                                  upb, inltx(ninltx,2)
    real   (kind=8), intent(out):: outDose(ninltx,2), mcED(nsim,ninltx),&
                                   pars(n2), stdp(n2), fvec1(ndat), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: i, j, info, seed, mcCount, niter
    real   (kind=8):: wght1(ndat), minCond, ran(2), outp(3),&
                      locp(5), cpars(n2), cstdp(n2),&
                      cfvec1(ndat), cfmin, cond, grad,&
                      maxDose, ifmin, avg, low(ninltx), up(ninltx),&
                      low1(ninltx), up1(ninltx), sumdose, sumdose2,& 
                      simltx(ndat), mcOSL(1), mcDose,&
                      Xm, Ym, aaa, bbb, ccc, upDose,&
                      DDD, DDD1
    !
    outDose = -99.0
    mcED = -99.0
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    message = 1
    !
    if (uw==0) then
        ! Un-weighted.
        wght1 = 1.0
    else if (uw==1) then
        ! Weighted.
        wght1 = sltx
    end if
    !
    if (model==0) then
        ! Linear model.
        call linefit(dose,ltx,wght1,ndat,pars,&
                     stdp,n2,fvec1,fmin,info)
        if (info==0) message = 0
        if (message/=0) return
    else if (model==1 .or. model==2 .or. model==3) then
        ! Non-linear model.
        minCond = 1.0D+20
        maxDose = maxval(dose)
        ! Initialization.
        loopA: do i=1, nstart
            call random_number(ran)
            ran = ran*upb
            call inipars(ran(1),ran(2),model,n2,&
                         dose,ltx,wght1,ndat,outp,info)
            if (info/=0) cycle loopA
            !
            locp = 0.0
            if (model==1 .or. model==2) then
                locp(1) = outp(1)
                locp(2) = ran(1)
                locp(3) = outp(2)
                locp(4) = outp(3)
            else if (model==3) then
                locp(1) = outp(1)
                locp(2) = ran(1)
                locp(3) = outp(2)
                locp(4) = ran(2)
                locp(5) = outp(3)
            end if 
            cpars = locp(1:n2)
            !
            call lmfit(dose,ltx,wght1,ndat,cpars,cstdp,n2,&
                       model,cfvec1,cfmin,cond,info)
            if (info/=0) cycle loopA
            !
            ! Check the saturating level.
            locp = 0.0
            locp(1:n2) = cpars
            if (model==1 .or. model==3) then
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                       locp(3)*locp(4)*exp(-locp(4)*maxDose)
            end if
            if (model==2) then
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                       locp(3)
            end if
            if (grad<1.0D-07) cycle loopA
            !
            if (cond<minCond) then
                pars = cpars
                stdp = cstdp
                fvec1 = cfvec1
                fmin = cfmin
                minCond = cond
                message = 0
            end if
        end do loopA
        !
        if (message/=0) return
    end if
    !
    !
    ! Calculate quivalent dose values. 
    if (model==0) then
        locp = 0.0
        locp(1:n2) = pars
        outDose(:,1) = (inltx(:,1)-locp(2))/locp(1)
    else if (model==1 .or. model==2 .or. model==3) then 
        !
        if (model==1 .or. model==3) then
            locp = 0.0
            locp(1:n2) = pars
            !
            if (model==1) then
                aaa = locp(1)
                bbb = locp(2)
                ccc = locp(3)
            else if (model==3) then
                if (locp(2)>locp(4)) then
                    aaa = locp(3)
                    bbb = locp(4)
                    ccc = locp(5)
                else 
                    aaa = locp(1)
                    bbb = locp(2)
                    ccc = locp(5)
                end if
            end if
            !
            DDD = 1.0D+00/bbb
            !
            Xm = -log(1.0D-03/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc
            !
            if (maxval(inltx(:,1))>Ym) then
                message = 1
                return
            end if
            !
        end if
        ! 
        if (model==1 .or. model==3) then
            upDose = Xm
        else 
            upDose = 5.0*maxDose
        end if
        !
        do i=1, ninltx
            call interpolate(-10.0D+00,upDose,inltx(i,1),&
                             pars,n2,model,outDose(i,1),ifmin)
            ! Check the quality of interpolation.
            if (ifmin>1.0D-06) then
                message = 1
                return
            end if
        end do
        !
    end if
    !
    !
    ! Assess the standard errors  
    ! of equivalent dose values.
    if (method==0) then
        ! Simple transformation.
        avg = sum((ltx-fvec1)**2)/real(ndat)
        low = inltx(:,1) - sqrt((inltx(:,2))**2+avg)
        up = inltx(:,1) + sqrt((inltx(:,2))**2+avg)
        !
        if (model==0) then 
            locp = 0.0
            locp(1:n2) = pars
            low1 = (low-locp(2))/locp(1)
            up1 = (up-locp(2))/locp(1)
        else if (model==1 .or. model==2 .or. model==3) then
            do i=1, ninltx
                !
                call interpolate(-10.0D+00,upDose,low(i),&
                                 pars,n2,model,low1(i),ifmin)
                if (ifmin>1.0D-06) then
                    message = 1
                    return
                end if
                !
                call interpolate(-10.0D+00,upDose,up(i),&
                                 pars,n2,model,up1(i),ifmin)
                if (ifmin>1.0D-06) then
                    message = 1
                    return
                end if
                !
            end do
        end if
        !
        outDose(:,2) = (up1-low1)/2.0
    else if (method==1) then
        ! Monte Carlo simulation.
        seed = 123456789
        do i=1, ninltx
            mcCount = 0
            sumdose = 0.0
            sumdose2 = 0.0
            niter = 0
            Innerloop: do
                ! Record the number of iterations.
                niter = niter + 1
                if (niter>200*nsim) then
                    message = 1
                    return
                end if
                ! Simulate standardised OSLs.
                do j=1, ndat
                    call r8vec_normal(1,ltx(j),sltx(j),seed,simltx(j))
                end do
                !
                if (model==0) then
                    call linefit(dose,simltx,wght1,ndat,cpars,&
                                 cstdp,n2,cfvec1,cfmin,info)
                    if (info/=0) cycle Innerloop
                else if (model==1 .or. model==2 .or. model==3) then
                    call random_number(ran)
                    ran = ran*upb
                    call inipars(ran(1),ran(2),model,n2,dose,&
                                 simltx,wght1,ndat,outp,info)
                    if (info/=0) cycle Innerloop
                    !
                    locp = 0.0
                    if (model==1 .or. model==2) then
                        locp(1) = outp(1)
                        locp(2) = ran(1)
                        locp(3) = outp(2)
                        locp(4) = outp(3)
                    else if (model==3) then
                        locp(1) = outp(1)
                        locp(2) = ran(1)
                        locp(3) = outp(2)
                        locp(4) = ran(2)
                        locp(5) = outp(3)
                    end if 
                    cpars = locp(1:n2)
                    !
                    call lmfit(dose,simltx,wght1,ndat,cpars,cstdp,&
                               n2,model,cfvec1,cfmin,cond,info)
                    if (info/=0) cycle Innerloop
                    !
                    ! Check the saturating level.
                    locp = 0.0
                    locp(1:n2) = cpars
                    if (model==1 .or. model==3) then
                        grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                               locp(3)*locp(4)*exp(-locp(4)*maxDose)
                        if (model==1) then
                            DDD1 = 1.0D+00/locp(2)
                        else if (model==3) then
                            DDD1 = 1.0D+00/min(locp(2),locp(4))
                        end if
                        if (DDD1<0.5*DDD) cycle Innerloop
                    end if
                    if (model==2) then
                        grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                               locp(3)
                    end if
                    if (grad<1.0D-07) cycle Innerloop
                end if
                !
                ! Simulate natural standardised OSL.
                call r8vec_normal(1,inltx(i,1),inltx(i,2),seed,mcOSL(1))
                !
                if (model==0) then
                    locp = 0.0
                    locp(1:n2) = cpars
                    mcDose = (mcOSL(1)-locp(2))/locp(1)
                else if (model==1 .or. model==2 .or. model==3) then 
                    call interpolate(-10.0D+00,upDose,mcOSL(1),&
                                     cpars,n2,model,mcDose,ifmin)
                    ! Check the quality of interpolation.
                    if (ifmin>1.0D-06) cycle Innerloop
                end if
                ! Record the values.
                mcCount = mcCount + 1
                sumdose = sumdose + mcDose
                sumdose2 = sumdose2 + mcDose**2
                mcED(mcCount,i) = mcDose
                if (mcCount==nsim) exit Innerloop
            end do Innerloop
            !
            outDose(i,2) = sqrt((real(nsim)*sumdose2-sumdose**2)/&
                                 real(nsim)/real(nsim-1))
        end do
    end if
    !
    return
end subroutine calED
