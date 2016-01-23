subroutine linefit(xd,yd,syd,nd,pars,stdp,&
                   n2,fvec1,fmin,message)
!-----------------------------------------------------
! Subroutine linefit() is used for fitting
! a line of the formula y=a*x or y=a*x+b.
!-----------------------------------------------------
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observation Y.
!   syd(nd):: input, real vlaues, weight of Y.
!        nd:: input, integer, number of points.
!  pars(n2):: output, real values, parameters.
!  stdp(n2):: output, real values, errors of pars.
!        n2:: input, integer, number ([1,2]) of pars.
! fvec1(nd):: output, real vlaues, fitted Y.
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, 0=success, 1=fail.
!-----------------------------------------------------
! Author:: Peng Jun, 2014.10.02.
!-----------------------------------------------------
! Dependence:: subroutine numHess, subroutine inverse.
!-----------------------------------------------------
    ! Arguments.
    integer(kind=4), intent(in):: nd, n2
    real   (kind=8), intent(in):: xd(nd), yd(nd),&
                                  syd(nd)
    real   (kind=8), intent(out):: pars(n2), stdp(n2),&
                                   fvec1(nd), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    real   (kind=8):: wght(nd), xx(2),&
                      hess(n2,n2), diag(n2), avgdv
    integer(kind=4):: i, errorflag, singular
    integer(kind=4), parameter:: model=0
    !
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    message = 0 
    !
    wght = 1.0/syd**2
    xx = 0.0
    !
    if (n2==1) then
        if (sum(wght*xd**2)==0.0) then
            message = 1
            return
        end if
        xx(1) = sum(wght*xd*yd) / sum(wght*xd**2)
    else if (n2==2) then
        if (sum(wght)*sum(wght*xd**2)==&
            sum(wght*xd)*sum(wght*xd)) then
            message = 1
            return
        end if
        xx(1) = (sum(wght)*sum(wght*xd*yd)-&
                 sum(wght*xd)*sum(wght*yd))/&
                (sum(wght)*sum(wght*xd**2)-&
                 sum(wght*xd)*sum(wght*xd))
        xx(2) = (sum(wght*yd)-xx(1)*sum(wght*xd))/& 
                 sum(wght)
    end if
    !
    fvec1 = xx(1)*xd + xx(2)
    fmin = sum(wght*(yd-fvec1)**2)
    !
    if (nd==n2) then
        avgdv = fmin
    else if (nd>n2) then
        avgdv = fmin / real(nd-n2)
    end if
    !
    pars(1:n2) = xx(1:n2)
    !
    if (pars(1)<=0.0) then
        message = 1
        return
    end if
    !
    call numHess(xd,yd,syd,nd,model,&
                 pars,n2,hess,errorflag)
    if (errorflag/=0) then
        message = 1
        return
    end if
    !
    call inverse(hess,n2,singular)
    if (singular/=0) then
        message = 1
        return
    end if
    !
    do i=1, n2
        diag(i) = hess(i,i) * avgdv
    end do
    if (any(diag<0.0)) then
        message = 1
        return
    end if
    stdp = sqrt(diag)
    !
    return
end subroutine linefit