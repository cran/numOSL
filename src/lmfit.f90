subroutine lmfit(xd,yd,syd,nd,pars,stdp,n2,&
                 model,fvec1,fmin,message)
!----------------------------------------------------
! Subroutine lmfit() is used for fitting a non-linear
! model using the Levenberg-Marquardt algorithm.
!----------------------------------------------------
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
!   syd(nd):: input, real values, weigth of Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, parameters.
!  stdp(n2):: input, real vlaues, error of pars.
!        n2:: input, integer, number of pars (>=2).
!     model:: input, integer: 1=exp;
!                             2=lexp;
!                             3=dexp;
!                             4=CW-OSL;
!                             5=LM-OSL.
! fvec1(nd):: output, real vlaues, fitted Y.
!      fmin:: output, real value, minimized objective.
!   message:: output, error indicator, 
!             0: success, 
!             1: fail, 
!             2: stdp<=0.
!------------------------------------------------------
! Author:: Peng Jun, 2023.08.30.
!------------------------------------------------------
! Dependence:: subroutine lmdif1; 
!              subroutine lmfunc;
!              subroutine inverse_sym.
!------------------------------------------------------
    implicit none
    ! Arguments.
    integer, intent(in):: nd, n2, model
    real(kind(1.0d0)), intent(in):: xd(nd), yd(nd), syd(nd)
    real(kind(1.0d0)), intent(inout):: pars(n2)
    real(kind(1.0d0)), intent(out):: stdp(n2), fvec1(nd),fmin
    integer, intent(out):: message
    ! Local variables.
    integer:: info, ifault, i
    real(kind(1.0d0)), parameter:: tol=1.490116D-08
    real(kind(1.0d0)):: fvec(nd), hess(n2,n2),&
                        diag(n2), avgdv
    external:: lmfunc
    !
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    ! 
    call lmdif1(lmfunc,nd,n2,pars,fvec,& 
                tol,info,xd,yd,syd,model,hess)
    !
    if (info==1 .or. info==2 .or. info==3) then
        message = 0
    else 
        message =1
        return
    end if
    !
    if (model==1) then
        if (n2==2 .and. any(pars<=0.0)) then
            message = 1
            return
        end if
        if (n2==3 .and. any(pars(1:n2-1)<=0.0)) then
            message = 1
            return
        end if
    end if
    !
    if (model==2) then
        if (n2==3 .and. any(pars<=0.0)) then
            message = 1
            return
        end if
        if (n2==4 .and. any(pars(1:n2-1)<=0.0)) then
            message = 1
            return
        end if
    end if
    !
    if (model==3) then
        if (n2==4 .and. any(pars<=0.0)) then
            message = 1
            return
        end if
        if (n2==5 .and. any(pars(1:n2-1)<=0.0)) then
            message = 1
            return
        end if
    end if
    !
    if (model==4 .or. model==5) then
        if (any(pars<=0.0)) then
            message = 1
            return
        end if
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
end subroutine lmfit
