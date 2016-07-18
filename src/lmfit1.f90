subroutine lmfit1(xd,yd,syd,nd,pars,stdp,&
                  n2,fvec1,fmin,message)
!----------------------------------------------------
! Subroutine lmfit1() is used for fitting a GOK
! model using the Levenberg-Marquardt algorithm.
!----------------------------------------------------
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
!   syd(nd):: input, real values, weigth of Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, parameters.
!  stdp(n2):: input, real vlaues, error of pars.
!        n2:: input, integer, number of pars (>=2).
! fvec1(nd):: output, real vlaues, fitted Y.
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, error indicator,
!             0: success,
!             1: fail,
!             2: stdp<=0.
!------------------------------------------------------
! Author:: Peng Jun, 2016.07.07.
!------------------------------------------------------
! Dependence:: subroutine lmdif1_bd; 
!              subroutine lmfunc1;
!              subroutine inverse_sym.
!------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: nd, n2
    real   (kind=8), intent(in):: xd(nd), yd(nd), syd(nd)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: stdp(n2), fvec1(nd), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, ifault, i
    integer(kind=4), parameter:: model=7
    real   (kind=8), parameter:: tol=1.490116D-08
    real   (kind=8):: fvec(nd), hess(n2,n2),&
                      diag(n2), avgdv,& 
                      lower(n2), upper(n2)
    external:: lmfunc1
    !
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    ! 
    ! Set lower and upper bounds.
    lower = 1.0e-10
    upper = 1.0e10
    if (n2==4) lower(n2) = -1.0e10
    !
    call lmdif1_bd(lmfunc1,nd,n2,pars,fvec,tol,& 
                   info,xd,yd,syd,lower,upper,hess)
    !
    if (info==1 .or. info==2 .or. info==3) then
        message = 0
    else 
        message =1
        return
    end if
    !
    fvec1 = yd + fvec*syd
    fmin = sum(fvec**2)
    !
    if (nd==n2) then
        avgdv = fmin
    else if (nd>n2) then
        avgdv = fmin / real(nd-n2)
    end if
    !
    call inverse_sym(hess,n2,ifault)
    if (ifault/=0) then
        message = 1
        return
    end if
    !
    do i=1, n2
        diag(i) = hess(i,i) * avgdv
    end do
    if (any(diag<=0.0)) then
        message = 2
        return
    end if
    !
    stdp = sqrt(diag)
    !
    return
end subroutine lmfit1
