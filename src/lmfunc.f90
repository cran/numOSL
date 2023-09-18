subroutine lmfunc(nd,n2,pars,fvec,&
                  iflag,xd,yd,syd,model)
!-----------------------------------------------------
! Subroutine lmfunc is used for calculating
! the residual vector of a given model.
! ----------------------------------------------------
!         nd:: input, integer, number of points.
!         n2:: input, integer, number ([2,15]) of pars.
!   pars(n2):: input, real vlaues, pars.
!   fvec(nd):: output, real values, residuals.
!      iflag:: integer.
!     xd(nd):: input, real values, observations X.
!     yd(nd):: input, real values, observations Y.
!    syd(nd):: input, real values, weight of Y.
!      model:: input, integer: 1=exp;
!                              2=lexp;
!                              3=dexp;
!                              4=CW-OSL;
!                              5=LM-OSL.
!-----------------------------------------------------
! Author:: Peng Jun, 2023.08.30.
!-----------------------------------------------------
! Dependence:: NO.
!-----------------------------------------------------
    ! Arguments.
    integer:: nd, n2, iflag, model
    real(kind(1.0d0)):: pars(n2), fvec(nd),& 
                        xd(nd), yd(nd), syd(nd)
    ! Local variables.
    real(kind(1.0d0)):: xx(15), xd1(nd), xd2(nd)
    integer:: i
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    if (model==1) then
        ! Exponential model.
        fvec = xx(1)*(1.0-exp(-xx(2)*xd))+&
               xx(3)
        fvec = (fvec-yd)/syd
        return
    end if 
    !
    if (model==2) then
         ! Exponential plus linear model.
         fvec = xx(1)*(1.0-exp(-xx(2)*xd))+&
                       xx(3)*xd+xx(4)
         fvec = (fvec-yd)/syd
         return
    end if
    !
    if (model==3) then
        ! Double saturating model.
        fvec = xx(1)*(1.0-exp(-xx(2)*xd))+&
               xx(3)*(1.0-exp(-xx(4)*xd))+&
               xx(5)
        fvec = (fvec-yd)/syd
        return
    end if
    !
    if (model==4) then
        ! CW-OSL decay curve.
        if (mod(n2,2)==0) then
            fvec = 0.0
            do i=1, n2/2
                fvec = fvec + xx(i)*xx(i+n2/2)*&
                       exp(-xx(i+n2/2)*xd)
            end do
        else if (mod(n2,2)==1) then
            fvec = xx(n2)
            do i=1, (n2-1)/2
                fvec = fvec + xx(i)*xx(i+(n2-1)/2)*&
                       exp(-xx(i+(n2-1)/2)*xd)
            end do
        end if
        fvec = (fvec-yd)/syd
        return
    end if
    !
    if (model==5) then
        ! LM-OSL decay curve.
        xd1 = xd/xd(nd)
        xd2 = xd**2/xd(nd)/2.0
        if (mod(n2,2)==0) then
            fvec = 0.0
            do i=1, n2/2
                fvec = fvec + xx(i)*xd1*xx(i+n2/2)*&
                       exp(-xx(i+n2/2)*xd2)
            end do 
        else if (mod(n2,2)==1) then 
            fvec = xx(n2)*xd1
            do i=1, (n2-1)/2
                fvec = fvec + xx(i)*xd1*xx(i+(n2-1)/2)*&
                       exp(-xx(i+(n2-1)/2)*xd2) 
            end do  
        end if
        fvec = (fvec-yd)/syd
        return
    end if
    !
end subroutine lmfunc          
