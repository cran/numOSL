subroutine calSGCED_fort(n2,inltx,pars,outDose,errMethod,avgDev,mcED,&
                         saturateDose,model,nsim,acceptRate,message)
!------------------------------------------------------------------------------
! Subroutine calSGCED_fort is used for calculating ED values
! using the SGC method using known growth curve parameters.
!------------------------------------------------------------------------------
!                n2:: input, integer, number of pars (>=1).
!        inltx(1,2):: input, real values, natural OSL and errors.
!          pars(n2):: input, real values, pars of growth curve.
!      outDose(1,2):: output, real values, estimated EDs and errors.
!         errMethod:: input, integer, 0=SP, 1=MC.
!            avgDev:: input, real value, average fit error.
!        mcED(nsim):: output, real values, simulated ED values.
!      saturateDose:: output, real value, saturate dose. 
!             model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp, 7=gok.
!              nsim:: input, integer, number of simulations.
!        acceptRate:: output, real value, accept rate of Monte Carlo simulation.
!           message:: output, integer, error indicator,
!                     0=success,
!                     1=natural OSL saturated,
!                     2=fail in ED caculating,
!                     3=fail in ED error estimation using SP (or MC).
!------------------------------------------------------------------------------
! Author:: Peng Jun, 2017.04.04. 
!------------------------------------------------------------------------------
! Dependence:: subroutine interpolate; subroutine r8vec_normal;
!              subroutine calDEXPxm.
!------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n2, model, nsim, errMethod
    real   (kind=8), intent(in):: inltx(1,2), pars(n2), avgDev
    real   (kind=8), intent(out):: outDose(1,2), mcED(nsim),& 
                                   saturateDose, acceptRate
    integer(kind=4), intent(out):: message
    ! Local variables.
    real   (kind=8):: locp(5), aaa, bbb, ccc, ddd, eee,& 
                      Xm, Ym, ifmin, sumdose, sumdose2,& 
                      mcOSL(1), mcDose, limSig, derivative,&
                      spErr, spltx1, spltx2, spDose1, spDose2
    integer(kind=4):: seed, mcCount, niter
    !
    outDose = -99.0
    mcED = -99.0
    saturateDose = -99.0
    acceptRate = 0.0 
    message = 0
    !
    locp = 0.0
    locp(1:n2) = pars
    !
    derivative = 1.0D-06
    !
    ! Calculate SGC equivalent doses.
    if (model==0) then
        !
        ! Linear model.
        saturateDose = 1.0D+05
        !
        outDose(1,1) = (inltx(1,1)-locp(2))/locp(1)
        !
    else 
        !
        ! Non-linear model.
        if (model==1) then
            ! EXP.
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            Xm = -log(derivative/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc
        else if (model==2) then
            ! LEXP.
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            if (ccc < derivative) Xm = -log((derivative-ccc)/aaa/bbb)/bbb
            if (ccc >= derivative) Xm = 1.0D+05
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) +ccc*Xm + ddd
        else if (model==3) then
            ! DEXP.
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            eee = locp(5)
            !
            call calDEXPxm(-50.0D+00,1.0D+05,derivative,locp(1:4),Xm)
            !
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc*(1.0D+00-exp(-ddd*Xm)) + eee
        else if (model==7) then
            ! GOK.
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            Xm = (derivative)**(-ccc/(1.0+ccc))*(1.0-(derivative)**(ccc/(1.0+ccc))*&
                 (1.0/aaa/bbb)**(ccc/(1.0+ccc)))*(1.0/aaa/bbb)**(-ccc/(1.0+ccc))/bbb/ccc 
            Ym = aaa*(1.0-(1.0+bbb*ccc*Xm)**(-1.0/ccc)) + ddd
        end if
        !
        saturateDose = Xm
        !
        if (inltx(1,1)>=Ym) then
            message = 1
            return
        end if
        !
        call interpolate(-50.0D+00,saturateDose,inltx(1,1),pars,&
                         n2,model,outDose(1,1),ifmin)  
        ! Check quality of interpolation.
        if (ifmin>1.0D-03) then
            message = 2
            return
        end if
        !
    end if
    !
    !
    ! Assess standard errors of SGC equivalent doses.
    if (errMethod==0) then
        !
        ! Simple transformation.
        spErr = sqrt(avgDev**2+(inltx(1,2))**2)
        !
        spltx1 = inltx(1,1) - spErr
        spltx2 = inltx(1,1) + spErr
        !
        if (model==0) then
            spDose1 = (spltx1-locp(2))/locp(1)
            spDose2 = (spltx2-locp(2))/locp(1)
        else 
            if (spltx2>=Ym) then
                message = 3
                return
            end if
            !
            call interpolate(-50.0D+00,saturateDose,spltx1,pars,&
                             n2,model,spDose1,ifmin)  
            ! Check quality of interpolation.
            if (ifmin>1.0D-03) then
                message = 3
                return
            end if
            !
            call interpolate(-50.0D+00,saturateDose,spltx2,pars,&
                             n2,model,spDose2,ifmin)  
            ! Check quality of interpolation.
            if (ifmin>1.0D-03) then
                message = 3
                return
            end if
            !
        end if
        !
        outDose(1,2) = (spDose2-spDose1)/2.0
        !
    else if (errMethod==1) then
        seed = 123456789
        mcCount = 0
        sumdose = 0.0
        sumdose2 = 0.0
        niter = 0
        !
        Innerloop: do 
            niter = niter + 1
            if (niter>100*nsim) then
                message = 3
                return
            end if
            !
            ! Simulate standardised natural OSLs.
            call r8vec_normal(1,inltx(1,1),inltx(1,2),seed,mcOSL(1))
            !
            if (model==1) then
                limSig = locp(1) + locp(3)
                if (mcOSL(1)>0.999*limSig) cycle Innerloop
            else if (model==3) then
                limSig = locp(1) + locp(3) + locp(5)
                if (mcOSL(1)>0.999*limSig) cycle Innerloop
            else if (model==7) then
                limSig = locp(1) + locp(4)
                if (mcOSL(1)>0.999*limSig) cycle Innerloop
            end if
            !
            ! Calculate simulated equivalent dose.
            if (model==0) then
                ! Linear model.
                mcDose = (mcOSL(1)-locp(2))/locp(1)
            else 
                ! Non-linear model.
                call interpolate(-50.0D+00,saturateDose,mcOSL(1),pars,&
                                 n2,model,mcDose,ifmin)
                ! Check quality of interpolation.
                if (ifmin>1.0D-03) cycle Innerloop
                if (mcDose>0.999*saturateDose) cycle Innerloop
                if (abs(mcDose)>5.0*abs(outDose(1,1))) cycle Innerloop
            end if
            !
            ! Record values. 
            mcCount = mcCount + 1
            sumdose = sumdose + mcdose
            sumdose2 = sumdose2 + mcdose**2
            mcED(mcCount) = mcDose
            !
            if (mcCount==nsim) exit Innerloop
            !
        end do Innerloop
        !
        outDose(1,2) = sqrt((real(nsim)*sumdose2-sumdose**2)/&
                             real(nsim)/real(nsim-1))
        !
        if (outDose(1,2) .ne. outDose(1,2)) then
            message = 3
            return
        end if 
        !
        acceptRate = real(nsim)/real(niter)*100.0
        !
    end if
    !
    return
end subroutine calSGCED_fort
