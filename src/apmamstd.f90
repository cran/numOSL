subroutine apmamstd(ed,sed,ndat,pars,&
                    stdp,np,iflag)
!--------------------------------------------------
! Subroutine apmamstd() is used for estimating
! the standard errors of a minimum age model 
! using finite difference approximation.
!--------------------------------------------------
! ed(ndat):: input, real values, logged EDs.
! sed(ndat:: input, real values, relative stdEDs.
!     ndat:: input, integer, number of data points.
! pars(np):: input, real vlaues, parameters.
! stdp(np):: output, estimated standard errors.
!       np:: input, integer, number ([3,4]) of pars.
!    iflag:: output, integer, 0=success, 1=fail.
!--------------------------------------------------
! Author:: Peng Jun ,2014.09.30.
!--------------------------------------------------
! Dependence:: subroutine numHess;-----------------
!              subroutine inverse.-----------------
!--------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, np
    real   (kind=8), intent(in):: ed(ndat), sed(ndat),& 
                                  pars(np)
    real   (kind=8), intent(out):: stdp(np)
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    integer(kind=4):: i, errorflag, singular
    real   (kind=8):: hess(np,np), diag(np), syd(ndat)
    integer(kind=4), parameter:: model=6
    !
    iflag = 0
    stdp = -99.0
    !
    syd = 1.0
    call numHess(sed,ed,syd,ndat,model,&
                 pars,np,hess,errorflag)
    if (errorflag/=0) then
        iflag = 1
        return
    end if
    !
    call inverse(hess,np,singular)
    if (singular/=0) then
        iflag = 1
        return
    end if
    !
    do i=1, np
        diag(i) = hess(i,i)
    end do
    if (any(diag<0.0)) then
        iflag = 1
        return
    end if
    !
    stdp = sqrt(diag)
    !
    return
end subroutine apmamstd      
