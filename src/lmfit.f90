subroutine lmfit(xd,yd,syd,nd,pars,stdp,n2,&
                 model,fvec1,fmin,cond,message)
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
!      cond:: output, real vlaue, logged condition(Inf).
!   message:: output, integer, 0=success, 1=fail.
!------------------------------------------------------
! Author:: Peng Jun, 2014.10.02.
!------------------------------------------------------
! Dependence:: subroutine lmdif1; subroutine lmfunc;
!              subroutine numHess, subroutine inverse.
!------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: nd, n2, model
    real   (kind=8), intent(in):: xd(nd), yd(nd), syd(nd)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: stdp(n2), fvec1(nd),& 
                                   fmin, cond
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, errorflag, singular, i
    real   (kind=8), parameter:: tol=1.490116D-08
    real   (kind=8):: fvec(nd), hess(n2,n2), diag(n2),&
                      cond1, cond2, avgdv
    external:: lmfunc
    !
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    cond = -99.0
    ! 
    call lmdif1(lmfunc,nd,n2,pars,fvec,& 
                tol,info,xd,yd,syd,model)
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
    call numHess(xd,yd,syd,nd,model,&
                 pars,n2,hess,errorflag)
    if (errorflag/=0) then
        message = 1
        return
    end if
    !
    cond1 = maxval(sum(abs(hess),dim=2))
    !
    call inverse(hess,n2,singular)
    if (singular/=0) then
        message = 1
        return
    end if
    !
    cond2 = maxval(sum(abs(hess),dim=2))
    !
    cond = log(cond1*cond2) 
    !
    do i=1, n2
        diag(i) = hess(i,i) * avgdv
    end do
    if (any(diag<0.0)) then
        message = 1
        return
    end if
    !
    stdp = sqrt(diag)
    !
    return
end subroutine lmfit
