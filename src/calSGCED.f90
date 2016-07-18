subroutine calSGCED(ninltx,n2,inltx,pars,&
                    outDose,model,message,whichErr)
!---------------------------------------------------------------------------------
! Subroutine calSGCED is used for calculating ED values
! using the SGC method using known growth curve parameters.
!---------------------------------------------------------------------------------
!            ninltx:: input, integer, number of ED values.
!                n2:: input, integer, number of pars (>=1).
!   inltx(ninltx,2):: input, real values, natural OSL and errors.
!          pars(n2):: input, real values, pars of growth curve.
! outDose(ninltx,2):: output, real values, estimated EDs and errors.
!             model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp, 7=gok.
!           message:: output, integer, error indicator,
!                     0=success,
!                     1=natural OSL saturated,
!                     2= fail in ED caculating,
!                     3= fail in ED error estimation using sp.
!          whichErr:: output, integer, the index of inltx in which calculation failed.
!---------------------------------------------------------------------------------
! Author:: Peng Jun, 2016.07.16.
!---------------------------------------------------------------------------------
! Dependence:: subroutine interpolate.
!---------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: ninltx, n2, model
    real   (kind=8), intent(in):: inltx(ninltx,2), pars(n2)
    real   (kind=8), intent(out):: outDose(ninltx,2)
    integer(kind=4), intent(out):: message, whichErr
    ! Local variables.
    real   (kind=8):: locp(5), aaa, bbb, ccc, ddd, eee, Xm, Ym,& 
                      upDose, ifmin, low(ninltx), up(ninltx),&
                      low1(ninltx), up1(ninltx), maxSig
    integer(kind=4):: i
    !
    outDose = -99.0
    message = 0
    whichErr = 0
    !
    ! Calculate SGC equivalent doses.
    if (model==0) then
        !
        ! Linear model.
        locp = 0.0
        locp(1:n2) = pars
        outDose(:,1) = (inltx(:,1)-locp(2))/locp(1)
        !
    else 
        !
        ! Non-linear model.
        !
        locp = 0.0
        locp(1:n2) = pars
        !
        if (model==1) then
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            Xm = -log(1.0D-04/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc
        else if (model==2) then
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            if (ccc<1.0D-04) Xm = -log((1.0D-04-ccc)/aaa/bbb)/bbb
            if (ccc>=1.0D-04) Xm = 1.0D+05
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
            Xm = -log(1.0D-04/aaa/bbb)/bbb
            Ym = aaa*(1.0D+00-exp(-bbb*Xm)) + ccc*(1.0D+00-exp(-ddd*Xm)) + eee
        else if (model==7) then
            aaa = locp(1)
            bbb = locp(2)
            ccc = locp(3)
            ddd = locp(4)
            Xm = (0.0001)**(-ccc/(1.0+ccc))*(1.0-(0.0001)**(ccc/(1.0+ccc))*&
                 (1.0/aaa/bbb)**(ccc/(1.0+ccc)))*(1.0/aaa/bbb)**(-ccc/(1.0+ccc))/bbb/ccc 
            Ym = aaa*(1.0-(1.0+bbb*ccc*Xm)**(-1.0/ccc)) + ddd
        end if
        !
        do i=1, ninltx
            maxSig = inltx(i,1) + inltx(i,2)
            !
            if (maxSig>=Ym) then
                message = 1
                whichErr = i
                return
            end if
        end do
        !
        upDose = Xm
        !
        do i=1, ninltx
            call interpolate(-50.0D+00,upDose,inltx(i,1),pars,&
                             n2,model,outDose(i,1),ifmin)  
            ! Check quality of interpolation.
            if (ifmin>1.0D-02) then
                message = 2
                whichErr = i
                return
            end if
        end do
        !
    end if
    !
    !
    ! Assess standard errors of SGC equivalent doses.
    low = inltx(:,1) - inltx(:,2)
    up = inltx(:,1) + inltx(:,2)
    !
    if (model==0) then
        !
        ! Linear model.
        locp = 0.0
        locp(1:n2) = pars
        low1 = (low-locp(2))/locp(1)
        up1 = (up-locp(2))/locp(1)
        !
    else 
        !
        ! Non-linear model. 
        do i=1, ninltx
            call interpolate(-50.0D+00,upDose,low(i),pars,&
                             n2,model,low1(i),ifmin)
            ! Check quality of interpolation.
            if (ifmin>1.0D-02) then
                whichErr = i
                message = 3
                return
            end if
            !
            call interpolate(-50.0D+00,upDose,up(i),pars,&
                             n2,model,up1(i),ifmin)
            ! Check quality of interpolation.
            if (ifmin>1.0D-02) then
                whichErr = i
                message = 3
                return
            end if
        end do
        !
    end if
    !
    outDose(:,2) = (up1-low1)/2.0
    !
    return
end subroutine calSGCED
