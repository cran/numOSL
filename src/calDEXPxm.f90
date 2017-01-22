subroutine calDEXPxm(ax,bx,Ym,pars,Xm)
!-----------------------------------------------------
! Subroutine calDexpXmis used for calculating
! the suturation dose of a DEXP model.
!-----------------------------------------------------
!       ax:: input, real value, lower limit.
!       bx:: input, real value, upper limit.
!       Ym:: input, real value, Y value.
!  pars(4):: input, real vlaues, parameters a, b, c, d.
!       xm:: output, real value, X value.
!-----------------------------------------------------
! Author:: Peng Jun, 2016.12.31.
!-----------------------------------------------------
! Dependence:: None.
!-----------------------------------------------------
    implicit none
    ! Arguments.
    real   (kind=8), intent(in):: ax, bx, Ym, pars(4)
    real   (kind=8), intent(out):: Xm
    !
    integer(kind=4), parameter:: ntry=1000000
    real   (kind=8):: vstep, x, y
    integer(kind=4):: i, niter
    !
    vstep = (bx-ax)/real(ntry-1)
    !
    niter = 0
    !
    do i=1, ntry
        x = ax +  real(i-1)*vstep
        y = pars(1)*pars(2)*exp(-pars(2)*x)+pars(3)*pars(4)*exp(-pars(4)*x)
        if (y<=Ym) exit
        niter = niter + 1
    end do

    if (niter<ntry)  then
        Xm=x
    else 
        Xm=bx
    end if
    !
    return
end subroutine calDEXPxm
