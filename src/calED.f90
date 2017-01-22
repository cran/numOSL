subroutine calED(dose,ltx,sltx,ndat,n2,inltx,&
                 outDose,mcED,pars,stdp,model,uw,&
                 nsim,fvec1,fmin,saturateDose,& 
                 acceptRate,message)
!--------------------------------------------------------------------
! Subroutine calED is used for calculating  
! ED values and assess their standard errors.
!--------------------------------------------------------------------
!        dose(ndat):: input, real values, dose values.
!         ltx(ndat):: input, real values, Lx/Tx values
!        sltx(ndat):: input, real values, errors of Lx/Txs.
!              ndat:: input, integer, number of data points.
!                n2:: input, integer, number of pars (>=1).
!        inltx(1,2):: input, real values, natural Lx/Tx and error.
!      outDose(1,2):: output, real values, estimated ED and error.
!        mcED(nsim):: output, real values, simulated ED values.
!          pars(n2):: output, real values, pars of growth curve.
!          stdp(n2):: output, real values, std errors of pars.
!             model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp.
!                uw:: input, integer, 0=un-weighted, 1=weighted.
!              nsim:: input, integer, number of simulations.
!       fvec1(ndat):: output, real values, fitted Lx/Tx values.
!              fmin:: output, real value, minimized objective.
!      saturateDose:: output, real value, saturate dose.
!        acceptRate:: output, real vlaue, accept rate in MC.
!           message:: output, integer, error indicator,
!                     0=success,
!                     1=fail in DRC fitting,
!                     2=natural OSL saturated,
!                     3=fail in ED calculating,
!                     4=fail in ED error estimation using MC.
!--------------------------------------------------------------------
! Author:: Peng Jun, 2017.01.13.
!--------------------------------------------------------------------
! Dependence:: subroutine linefit; 
!              subroutine inipars
!              subroutine lmfit; 
!              subroutine interpolate;
!              subroutine r8vec_normal;
!              subroutine calDEXPxm. 
!--------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2, model, uw, nsim
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat),&
                                  sltx(ndat), inltx(1,2)
    real   (kind=8), intent(out):: outDose(1,2), mcED(nsim),&
                                   pars(n2), stdp(n2),&
                                   fvec1(ndat), fmin,&
                                   saturateDose, acceptRate
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: i, j, info, seed, mcCount, niter, ookk
    real   (kind=8):: wght1(ndat), minValue, ran(2), outp(3),&
                      locp(5), cpars(n2), cstdp(n2),cfvec1(ndat),&
                      cfmin, grad, maxDose, ifmin, sumdose, sumdose2,& 
                      simltx(ndat), mcOSL(1), mcDose, Xm, Ym, aaa, bbb,& 
                      ccc, ddd, eee, inib(24), maxSig
    !
    outDose = -99.0
    mcED = -99.0
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    saturateDose = -99.0
    acceptRate = 0.0 
    message = 0
    !
    if (uw==0) then
        ! Un-weighted.
        wght1 = 1.0
    else if (uw==1) then
        ! Weighted.
        wght1 = sltx
    end if
    !
    !
    ! Fit growth curve.
    if (model==0) then
        !
        ! Linear model.
        call linefit(dose,ltx,wght1,ndat,pars,&
                     stdp,n2,fvec1,fmin,info)
        ! Check fit.
        if (info/=0) then
            message = 1
            return
        end if
        !
    else 
        !
        ! Non-linear model.
        ookk = 1
        minValue = 1.0D+20
        maxDose = maxval(dose)
        !
        do i=1, 12
            inib(2*i-1:2*i) = (10.0)**(i-11)*(/1.0, 5.0/)
        end do
        !
        ! Initialization.
        if (model==1 .or. model==2) then
            !
            ! EXP or LEXP.
            loopA: do i=1, 24
                !
                ran(1) = inib(i)
                ran(2) = 0.0
                !
                call inipars(ran(1),ran(2),model,n2,&
                             dose,ltx,wght1,ndat,outp,info)
                if (info/=0) cycle loopA
                !
                locp = 0.0
                locp(1) = outp(1)
                locp(2) = ran(1)
                locp(3) = outp(2)
                locp(4) = outp(3)
                !
                cpars = locp(1:n2)
                !
                call lmfit(dose,ltx,wght1,ndat,cpars,cstdp,&
                           n2,model,cfvec1,cfmin,info)
                ! Check fit.
                if (info/=0) cycle loopA
                !
                ! Check saturating level.
                locp = 0.0
                locp(1:n2) = cpars
                if (model==1) then
                    grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)
                else if (model==2) then
                    grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+locp(3)
                end if
                !
                if (grad<1.0D-07) cycle loopA
                !
                if (cfmin<minValue) then
                    pars = cpars
                    stdp = cstdp
                    fvec1 = cfvec1
                    fmin = cfmin
                    minValue = cfmin
                    ookk = 0
                end if
            end do loopA
            !
        else if (model==3) then
            !
            ! DEXP.
            loopB: do i=1, 24
                loopC: do j=i, 24
                    !
                    ran(1) = inib(i)
                    ran(2) = inib(j)
                    !
                    call inipars(ran(1),ran(2),model,n2,&
                                 dose,ltx,wght1,ndat,outp,info)
                    if (info/=0) cycle loopC
                    !
                    locp = 0.0
                    locp(1) = outp(1)
                    locp(2) = ran(1)
                    locp(3) = outp(2)
                    locp(4) = ran(2)
                    locp(5) = outp(3)
                    !
                    cpars = locp(1:n2)
                    !
                    call lmfit(dose,ltx,wght1,ndat,cpars,cstdp,&
                               n2,model,cfvec1,cfmin,info)
                    if (info/=0) cycle loopC
                    !
                    ! Check saturating level.
                    locp = 0.0
                    locp(1:n2) = cpars
                    grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                           locp(3)*locp(4)*exp(-locp(4)*maxDose)
                    !
                    if (grad<1.0D-07) cycle loopC
                    !
                    if (cfmin<minValue) then
                        pars = cpars
                        stdp = cstdp
                        fvec1 = cfvec1
                        fmin = cfmin
                        minValue = cfmin
                        ookk = 0
                    end if
                end do loopC
            end do loopB
            !
        end if
        !
        if (ookk/=0) then
            message = 1
            return
        end if
        !
    end if
    !
    !
    !
    ! Calculate quivalent dose. 
    if (model==0) then
        !
        ! Linear model.
        locp = 0.0
        locp(1:n2) = pars
        outDose(1,1) = (inltx(1,1)-locp(2))/locp(1)
        !
    else
        !
        ! Non-linear model.
        maxSig = inltx(1,1)
        !
        locp = 0.0
        locp(1:n2) = pars
        !
        if (model==1) then
            !
            ! EXP.
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            Xm = -log(1.0D-05/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc
            !
        else if (model==2) then
            !
            ! LEXP.
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            if (ccc<1.0D-05) Xm = -log((1.0D-05-ccc)/aaa/bbb)/bbb
            if (ccc>=1.0D-05) Xm = 1.0D+05
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc*Xm + ddd
            !
        else if (model==3) then
            !
            ! DEXP. 
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            eee = locp(5)
            !
            call calDEXPxm(-50.0D+00,1.0D+05,1.0D-05,locp(1:4),Xm)
            !
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc*(1.0D+00-exp(-ddd*Xm)) + eee
            !
        end if
        !  
        saturateDose = Xm
        !
        if (maxSig>=Ym) then
            message = 2
            return
        end if
        !
        call interpolate(-50.0D+00,saturateDose,inltx(1,1),pars,&
                         n2,model,outDose(1,1),ifmin)
        ! Check quality of interpolation.
        if (ifmin>1.0D-02) then
            message = 3
            return
        end if
        !
    end if
    !
    !
    ! Assess standard errors of equivalent dose values.
    ! Monte Carlo simulation.
    seed = 123456789
    !
    mcCount = 0
    sumdose = 0.0
    sumdose2 = 0.0
    niter = 0
    !
    Innerloop: do
        !
        ! Record the number of iterations.
        niter = niter + 1
        if (niter>100*nsim) then
            message = 4
            return
        end if
        !
        ! Simulate standardised regenerative OSLs.
        do j=1, ndat
            call r8vec_normal(1,ltx(j),sltx(j),seed,simltx(j))
        end do
        !
        !
        if (model==0) then
            !
            ! Linear model.
            !
            ! Fit random growth curve.
            call linefit(dose,simltx,wght1,ndat,cpars,&
                         cstdp,n2,cfvec1,cfmin,info)
            ! Check fit.
            if (info/=0) cycle Innerloop
            !
        else
            !
            ! Non-linear model.
            !
            ! Fit random growth curve.
            cpars = pars
            call lmfit(dose,simltx,wght1,ndat,cpars,&
                       cstdp,n2,model,cfvec1,cfmin,info)
            ! Check fit.
            if (info==1) cycle Innerloop
            !
            ! Check saturating level.
            locp = 0.0
            locp(1:n2) = cpars
            if (model==1) then
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)
            else if (model==2) then
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+locp(3)
            else if (model==3) then
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                       locp(3)*locp(4)*exp(-locp(4)*maxDose)
            end if
            if (grad<1.0D-07) cycle Innerloop
            !
        end if
        !
        ! Simulate standardised natural  OSL.
        call r8vec_normal(1,inltx(1,1),inltx(1,2),seed,mcOSL(1))
        !
        ! Calculate simulated equivalent dose.
        if (model==0) then
            !
            ! Linear model.
            locp = 0.0
            locp(1:n2) = cpars
            mcDose = (mcOSL(1)-locp(2))/locp(1)
            !
        else 
            !
            ! Non-linear model.
            call interpolate(-50.0D+00,saturateDose,mcOSL(1),cpars,&
                             n2,model,mcDose,ifmin)
            ! Check quality of interpolation.
            if (ifmin>1.0D-02) cycle Innerloop
            !
            if (abs(mcDose)>5.0*abs(outDose(1,1))) cycle Innerloop
            !
        end if
        !
        ! Record values.
        mcCount = mcCount + 1
        sumdose = sumdose + mcDose
        sumdose2 = sumdose2 + mcDose**2
        mcED(mcCount) = mcDose
        !
        if (mcCount==nsim) exit Innerloop
        !
    end do Innerloop
    !
    !
    outDose(1,2) = sqrt((real(nsim)*sumdose2-sumdose**2)/&
                         real(nsim)/real(nsim-1))
    !
    if (outDose(1,2) .ne. outDose(1,2)) then
        message = 4
        return
    end if 
    !
    acceptRate = real(nsim)/real(niter)*100.0
    !
    return
end subroutine calED
