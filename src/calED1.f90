subroutine calED1(dose,ltx,sltx,ndat,n2,inltx,outDose,&
                  errMethod,mcED,pars,stdp,uw,nsim,fvec1,&
                  fmin,saturateDose,acceptRate,message)
!--------------------------------------------------------------------
! Subroutine calED1 is used for calculating  
! ED values and assess their standard errors.
!--------------------------------------------------------------------
!        dose(ndat):: input, real values, dose values.
!         ltx(ndat):: input, real values, Lx/Tx values
!        sltx(ndat):: input, real values, errors of Lx/Txs.
!              ndat:: input, integer, number of data points.
!                n2:: input, integer, number of pars (>=1).
!        inltx(1,2):: input, real values, natural Lx/Txs and errors.
!      outDose(1,2):: output, real values, estimated EDs and errors.
!         errMethod:: input, integer, 0=SP, 1=MC.
!        mcED(nsim):: output, real values, simulated ED values.
!          pars(n2):: output, real values, pars of growth curve.
!          stdp(n2):: output, real values, std errors of pars.
!                uw:: input, integer, 0=un-weighted, 1=weighted.
!              nsim:: input, integer, number of simulations.
!       fvec1(ndat):: output, real values, fitted Lx/Tx values.
!              fmin:: output, real value, minimized objective.
!      saturateDose:: output, real value, saturate dose.
!        acceptRate:: output, real value, accept rate of MC.
!           message:: output, integer, error indicator,
!                     0=success, 
!                     1=fail in DRC fitting,
!                     2=natural OSL saturated,
!                     3=fail in ED calculating,
!                     4=fail in ED error estimation using SP (or MC).
!--------------------------------------------------------------------
! Author:: Peng Jun, 2017.04.01.
!--------------------------------------------------------------------
! Dependence:: subroutine inipars;
!              subroutine lmfit1; 
!              subroutine interpolate;
!              subroutine r8vec_normal.
!--------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2, uw, nsim, errMethod
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat),&
                                  sltx(ndat), inltx(1,2)
    real   (kind=8), intent(out):: outDose(1,2), mcED(nsim),&
                                   pars(n2), stdp(n2),&  
                                   saturateDose, acceptRate,&
                                   fvec1(ndat), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: i, j, info, seed, mcCount, niter, ookk
    integer(kind=4), parameter:: model=7
    real   (kind=8):: wght1(ndat), minValue, ran(2), outp(3),&
                      locp(4), cpars(n2), cstdp(n2), rcv(11),&
                      cfvec1(ndat), cfmin, maxDose, ifmin,& 
                      sumdose, sumdose2, simltx(ndat), mcOSL(1), mcDose,& 
                      Xm, Ym, aaa, bbb, ccc, ddd, inib(24), limSig, derivative,&
                      avgDev, spErr, spltx1, spltx2, spDose1, spDose2
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
    ookk = 1
    minValue = 1.0D+20
    maxDose = maxval(dose)
    !
    do i=1, 12
        inib(2*i-1:2*i) = (10.0)**(i-11)*(/1.0, 5.0/)
    end do
    !
    rcv = (/((10.0)**(i-6), i=1, 11)/)
    !
    ! Initialization.
    loopA: do i=1, 24
        loopB: do j=1, 11
            !
            ran(1) = inib(i)
            ran(2) = rcv(j)
            !
            call inipars(ran(1),ran(2),model,n2,&
                         dose,ltx,wght1,ndat,outp,info)
            if (info/=0) cycle loopB
            !
            locp = 0.0
            locp(1) = outp(1)
            locp(2) = ran(1)
            locp(3) = ran(2)
            locp(4) = outp(2)
            !
            cpars = locp(1:n2)
            !
            call lmfit1(dose,ltx,wght1,ndat,cpars,cstdp,&
                        n2,cfvec1,cfmin,info)
            ! Check fit.
            if (info/=0) cycle loopB
            !
            if (cfmin<minValue) then
                pars = cpars
                stdp = cstdp
                fvec1 = cfvec1
                fmin = cfmin
                minValue = cfmin
                ookk = 0
            end if
        end do loopB
    end do loopA
    !
    if (ookk/=0) then
        message = 1
        return
    end if
    !
    !
    ! Calculate equivalent dose.
    locp = 0.0
    locp(1:n2) = pars
    !   
    aaa = locp(1)
    bbb = locp(2)
    ccc = locp(3)
    ddd = locp(4)
    !
    derivative = 1.0D-06
    Xm = (derivative)**(-ccc/(1.0+ccc))*(1.0-(derivative)**(ccc/(1.0+ccc))*&
         (1.0/aaa/bbb)**(ccc/(1.0+ccc)))*(1.0/aaa/bbb)**(-ccc/(1.0+ccc))/bbb/ccc
    Ym = aaa*(1.0-(1.0+bbb*ccc*Xm)**(-1.0/ccc)) + ddd
    !
    saturateDose = Xm
    !
    if (inltx(1,1)>=Ym) then
        message = 2
        return
    end if
    !
    call interpolate(-50.0D+00,saturateDose,inltx(1,1),pars,&
                     n2,model,outDose(1,1),ifmin)
    ! Check quality of interpolation.
    if (ifmin>1.0D-03) then
        message = 3
        return
    end if
    !
    !
    ! Assess standard error of equivalent dose.
    if (errMethod==0) then
        !
        ! Simple transformation.
        avgDev = sqrt(sum((ltx-fvec1)**2))/real(ndat)
        spErr = sqrt(avgDev**2+(inltx(1,2))**2)
        !
        spltx1 = inltx(1,1) - spErr
        spltx2 = inltx(1,1) + spErr
        !
        if (spltx2>=Ym) then
            message = 4
            return
        end if
        !
        call interpolate(-50.0D+00,saturateDose,spltx1,pars,&
                         n2,model,spDose1,ifmin)
        ! Check quality of interpolation.
        if (ifmin>1.0D-03) then
            message = 4
            return
        end if
        !
        call interpolate(-50.0D+00,saturateDose,spltx2,pars,&
                         n2,model,spDose2,ifmin)
        ! Check quality of interpolation.
        if (ifmin>1.0D-03) then
            message = 4
            return
        end if
        !
        outDose(1,2) = (spDose2-spDose1)/2.0
        !
    else if (errMethod==1) then
        !
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
            ! Fit random growth curve.
            cpars = pars
            call lmfit1(dose,simltx,wght1,ndat,cpars,&
                        cstdp,n2,cfvec1,cfmin,info)
            ! Check fit. 
            if (info==1) cycle Innerloop
            !
            locp = 0.0
            locp(1:n2) = cpars
            !
            ! Simulate natural standardised OSL.
            call r8vec_normal(1,inltx(1,1),inltx(1,2),seed,mcOSL(1))
            !
            limSig = locp(1) + locp(4)
            if (mcOSL(1)>0.999*limSig) cycle Innerloop
            !
            ! Calculate simulated equivalent dose.
            call interpolate(-50.0D+00,saturateDose,mcOSL(1),cpars,&
                             n2,model,mcDose,ifmin)
            ! Check quality of interpolation.
            if (ifmin>1.0D-03) cycle Innerloop
            if (mcDose>0.999*saturateDose) cycle Innerloop
            if (abs(mcDose)>5.0*abs(outDose(1,1))) cycle Innerloop
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
    end if
    !
    return
end subroutine calED1
