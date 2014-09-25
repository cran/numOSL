subroutine camED(ed,sed,ndat,addsigma,&
                 pars,stdp,maxlik,bic)
!------------------------------------------------------
! Subroutine camED() is used for estimating
! parameters of a central age model.
!------------------------------------------------------
!  ed(ndat):: input, real values, logged EDs.
! sed(ndat):: input, real values, errors of logged EDs.
!      ndat:: input, integer, number of data points.
!  addsigma:: input, real value, additional error.
! pars(2,1):: output, real vlaues, estimated pars.
! stdp(2,1):: output, real values, estimated stdpars.
!    maxlik:: output, real vlaue, maximum likelihood.
!       bic:: output, real value, BIC value.
!-------------------------------------------------------
! Author:: Peng Jun, 2014.09.16.
!-------------------------------------------------------
! Dependence:: NO.
!-------------------------------------------------------
! Reference:: Galbraith, RF, Green PF, 1990. Estimating 
!             the component ages in a finite mixture. 
!             Nuclear Tracks and Radiation Measurements, 
!             17, page 197-206.
!-------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat
    real   (kind=8), intent(in):: ed(ndat), sed(ndat), addsigma
    real   (kind=8), intent(out):: pars(2,1), stdp(2,1),&
                                   maxlik, bic
    ! Local variables.
    real   (kind=8):: z(ndat), sz(ndat),&
                      wz(ndat), newp(2,1)
    real   (kind=8), parameter:: eps=1.0D-08,&
    PI=3.141592653589793238462643383279502884197D+00
    integer(kind=4), parameter:: maxiter=10000
    integer(kind=4):: i
    !
    z = ed
    sz = sqrt(sed**2+addsigma**2)
    !
    pars(1,1) = 0.1
    wz = 1.0/((pars(1,1))**2+sz**2)
    pars(2,1) = sum(z*wz)/sum(wz)
    !
    loopA: do i=1, maxiter
        newp(2,1) = sum(z*wz)/sum(wz)
        newp(1,1) = pars(1,1)*&
                    sqrt(sum((wz**2)*(z-newp(2,1))**2/sum(wz)))
        wz = 1.0/((newp(1,1))**2+sz**2)
        !
        if (abs(pars(1,1)-newp(1,1))+&
            abs(pars(2,1)-newp(2,1))<eps) then
            exit loopA
        else 
            pars = newp
        end if
        !
    end do loopA
    !
    maxlik = sum(log(1.0/sqrt(2.0*PI)*sqrt(wz)*&
                 exp(-(z-pars(2,1))**2*wz/2.0)))
    !
    bic = -2.0*maxlik - 2.0*log(real(ndat))
    !
    stdp(1,1) = 1.0/sqrt(2.0*pars(1,1)*sum(wz**2))
    pars(2,1) = exp(pars(2,1))
    stdp(2,1) = pars(2,1)/sqrt(sum(wz))
    !
    return
end subroutine camED         
