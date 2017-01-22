subroutine calSGCED(n2,inltx,pars,outDose,mcED,saturateDose,&
                    model,nsim,acceptRate,message)
!------------------------------------------------------------------------------
! Subroutine calSGCED is used for calculating ED values
! using the SGC method using known growth curve parameters.
!------------------------------------------------------------------------------
!                n2:: input, integer, number of pars (>=1).
!        inltx(1,2):: input, real values, natural OSL and errors.
!          pars(n2):: input, real values, pars of growth curve.
!      outDose(1,2):: output, real values, estimated EDs and errors.
!        mcED(nsim):: output, real values, simulated ED values.
!      saturateDose:: output, real value, saturate dose. 
!             model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp, 7=gok.
!              nsim:: input, integer, number of simulations.
!        acceptRate:: output, real value, accept rate of Monte Carlo simulation.
!           message:: output, integer, error indicator,
!                     0=success,
!                     1=natural OSL saturated,
!                     2=fail in ED caculating,
!                     3=fail in ED error estimation.
!------------------------------------------------------------------------------
! Author:: Peng Jun, 2017.01.13.
!------------------------------------------------------------------------------
! Dependence:: subroutine interpolate; subroutine r8vec_normal.
!------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n2, model, nsim
    real   (kind=8), intent(in):: inltx(1,2), pars(n2)
    real   (kind=8), intent(out):: outDose(1,2), mcED(nsim),& 
                                   saturateDose, acceptRate
    integer(kind=4), intent(out):: message
    ! Local variables.
    real   (kind=8):: locp(5), aaa, bbb, ccc, ddd, eee,& 
                      Xm, Ym, ifmin, maxSig, sumdose,& 
                      sumdose2, mcOSL(1), mcDose
    integer(kind=4):: seed, mcCount, niter
    !
    outDose = -99.0
    mcED = -99.0
    saturateDose = -99.0
    acceptRate = 0.0 
    message = 0
    !
    ! Calculate SGC equivalent doses.
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
        locp = 0.0
        locp(1:n2) = pars
        !
        if (model==1) then
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            Xm = -log(1.0D-05/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc
        else if (model==2) then
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            if (ccc<1.0D-05) Xm = -log((1.0D-05-ccc)/aaa/bbb)/bbb
            if (ccc>=1.0D-05) Xm = 1.0D+05
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) +ccc*Xm + ddd
        else if (model==3) then
            if (locp(2)>locp(4)) then
                aaa = locp(3) 
                bbb = locp(4)
                ccc = locp(1)
                ddd = locp(2)
                eee = locp(5)
            else 
                aaa = locp(1)
                bbb = locp(2)
                ccc = locp(3)
                ddd = locp(4)
                eee = locp(5)
            end if
            Xm = -log(1.0D-05/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc*(1.0D+00-exp(-ddd*Xm)) + eee
        else if (model==7) then
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            Xm = (1.0D-05)**(-ccc/(1.0+ccc))*(1.0-(1.0D-05)**(ccc/(1.0+ccc))*&
                 (1.0/aaa/bbb)**(ccc/(1.0+ccc)))*(1.0/aaa/bbb)**(-ccc/(1.0+ccc))/bbb/ccc 
            Ym = aaa*(1.0-(1.0+bbb*ccc*Xm)**(-1.0/ccc)) + ddd
        end if
        !
        saturateDose = Xm
        maxSig = inltx(1,1)
        !
        if (maxSig>=Ym) then
            message = 1
            return
        end if
        !
        call interpolate(-50.0D+00,saturateDose,inltx(1,1),pars,&
                         n2,model,outDose(1,1),ifmin)  
        ! Check quality of interpolation.
        if (ifmin>1.0D-02) then
            message = 2
            return
        end if
        !
    end if
    !
    !
    ! Assess standard errors of SGC equivalent doses.
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
        ! Calculate simulated equivalent dose.
        if (model==0) then
            ! Linear model.
            locp = 0.0
            locp(1:n2) = pars
            mcDose = (mcOSL(1)-locp(2))/locp(1)
        else 
            ! Non-linear model.
            call interpolate(-50.0D+00,saturateDose,mcOSL(1),pars,&
                             n2,model,mcDose,ifmin)
            ! Check quality of interpolation.
            if (ifmin>1.0D-02) cycle Innerloop
            !
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
    return
end subroutine calSGCED
