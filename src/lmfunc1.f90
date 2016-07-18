subroutine lmfunc1(nd,n2,pars,fvec,iflag,&
                   xd,yd,syd,lower,upper)
!-----------------------------------------------------
! Subroutine lmfunc1 is used for calculating
! the residual vector of a GOK model.
! ----------------------------------------------------
!         nd:: input, integer, number of points.
!         n2:: input, integer, number (3-4) of pars.
!   pars(n2):: input, real vlaues, pars.
!   fvec(nd):: output, real values, residuals.
!      iflag:: integer.
!     xd(nd):: input, real values, observations X.
!     yd(nd):: input, real values, observations Y.
!    syd(nd):: input, real values, weight of Y.
!  lower(n2):: input, real values, lower bounds.
!  upper(n2):: input, rea values, upper bounds.
!-----------------------------------------------------
! Author:: Peng Jun, 2016.06.09.
!-----------------------------------------------------
! Dependence:: NO.
!-----------------------------------------------------
    ! Arguments.
    integer(kind=4):: nd, n2, iflag
    real   (kind=8):: pars(n2), fvec(nd),& 
                      xd(nd), yd(nd), syd(nd),&
                      lower(n2), upper(n2)
    ! Local variables.
    real   (kind=8):: xx(4)
    integer(kind=4):: i
    !
    ! Bound constraints.
    do i=1, n2
        if (pars(i)<lower(i))  then
            pars(i) = lower(i)
        else if (pars(i)>upper(i)) then
            pars(i) = upper(i)
        end if
    end do
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    fvec = xx(1)*(1.0-(1.0+xx(2)*xx(3)*xd)**&
           (-1.0/xx(3)))+xx(4)
    fvec = (fvec-yd)/syd
    !
    return
    !
end subroutine lmfunc1
