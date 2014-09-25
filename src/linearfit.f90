subroutine linearFit(dose,ltx,sltx,ndat,pars,&
                     stdp,n2,fvec1,fvalue,info)
!-------------------------------------------------------------------
! Subroutine LinearFit() is used for fitting a
! linear growth curve by weighted linear fitting.
!-------------------------------------------------------------------
!  dose(ndat):: input, real values, the dose values.
!   ltx(ndat):: input, real values, the Lx/Tx values.
!  sltx(ndat):: input, real vlaues, the errors of Lx/Tx values.
!        ndat:: input, integer, number of data points.
!    pars(n2):: output, real values, estimated parameters.
!    stdp(n2):: output, real values, standard errors of parameters.
!          n2:: input, integer, fitting type, 1=origin, 2=non-origin.
! fvec1(ndat):: output, real values, fitted Lx/Tx values.
!      fvalue:: output, real value, the minimized chi-square value.
!        info:: output, integer, 0=success, 1=fail.
!-------------------------------------------------------------------
! Author:: Peng Jun, 2014.09.04.
!-------------------------------------------------------------------
! Dependence:: subroutine numHess; subroutine inverse;--------------
!              inner function fcn.----------------------------------
!-------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2
    real    (kind=8), intent(in):: dose(ndat), ltx(ndat), sltx(ndat)
    real    (kind=8), intent(out):: pars(n2), stdp(n2),&
                                    fvec1(ndat), fvalue
    integer(kind=4), intent(out):: info
    ! Local variables.
    real   (kind=8):: wght(ndat), locp(2), hess(n2,n2), diag(n2)
    integer(kind=4):: i, errorflag, singular
    ! 
    ! Weight by inverse variance.
    wght = 1.0/sltx**2
    !
    if (n2==1) then
        locp(1) = sum(wght*dose*ltx) / sum(wght*dose**2)
        fvec1 = locp(1)*dose
        fvalue = sum(wght*(ltx-fvec1)**2)
    else if (n2==2) then
        locp(1) = (sum(wght)*sum(wght*dose*ltx)-&
                   sum(wght*dose)*sum(wght*ltx)) /&
                  (sum(wght)*sum(wght*dose**2)-&
                   sum(wght*dose)*sum(wght*dose))
        locp(2) = (sum(wght*ltx)-locp(1)*sum(wght*dose)) / sum(wght)
        fvec1 = locp(1)*dose + locp(2)
        fvalue = sum(wght*(ltx-fvec1)**2)
    end if
    !
    pars(1:n2) = locp(1:n2)
    stdp = -99.0
    info = 0
    !
    if (pars(1)<=0.0) then
        info = 1
        return
    end if
    !
    call numHess(fcn,pars,n2,hess,errorflag)
    if (errorflag/=0) then
        info = 1
        return
    end if
    !
    call inverse(hess,n2,singular)
    if (singular/=0) then
        info = 1
        return
    end if
    !
    do i=1, n2
        diag(i) = hess(i,i)
    end do
    if (any(diag<0.0)) then
         info = 1
         return
    end if
    stdp = sqrt(diag)
    !
    return
    !
    contains
        ! Calculate the sqrt of the sum of 
        ! the squared weighted residuals.
        function fcn(x)
            implicit none
            real   (kind=8):: x(n2), fcn
            real   (kind=8):: locx(2)
            !
            locx(1:n2) = x(1:n2)
            !
            if (n2==1) then
                fcn = sqrt(sum(wght*(ltx-locx(1)*dose)**2))
            else if (n2==2) then
                fcn = sqrt(sum(wght*(ltx-locx(1)*dose-locx(2))**2))
            end if
            !
            return
        end function fcn
end subroutine linearFit  
