subroutine CAM(ED,Error,n,spreadsigma,maxiter,&
               eps,pars,spars,maxlik,BIC)
!-----------------------------------------------------------------------------------------------------------------------
! CAM() is used to do Central Age Model based equivalent dose calculation, using logged equivalent dose for analyzing.
! ======================================================================================================================
! ED(n),input                   :: real values, the unlogged equivalent dose values. 
!
! Error(n), input               :: real values, the assocaited absolute error of dquivalent dose.
!
! n, input                      :: integer, the size of dquivalent doses.
!
! spreadsigma, input,           :: real value, the spread dispersion to the relative error of dquivalent dose.
!                                
! maxiter, input                :: integer, the allowed maximum iterative number.
!
! eps, input                    :: real value, the maximum tolerance for stopping the iteration.
!
! pars(2,1), output             :: real values, the estimated parameters of central age model.
!
! spars(2,1), output            :: real values, the standard errors of estimated parameters.
!
! maxlik, output                :: real value, the maximum logged-likelihood value.
!
! BIC, output                   :: real value, the BIC value.
! ======================================================================================================================
! Author:: Peng Jun, 2013.03.04, revised in 2013.03.05.
!
! Dependence:: No
!
! References:: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!              Standard Errors.Techno-metrics, 30, page 271-281.
!
!              Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!              finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!-----------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),               intent(in)::n
  integer(kind=4),               intent(in)::maxiter
  real   (kind=8),               intent(in)::eps
  real   (kind=8),               intent(in)::spreadsigma
  real   (kind=8),dimension(n),  intent(in)::ED
  real   (kind=8),dimension(n),  intent(in)::Error
  real   (kind=8),dimension(2,1),intent(out)::pars
  real   (kind=8),dimension(2,1),intent(out)::spars
  real   (kind=8),               intent(out)::maxlik
  real   (kind=8),               intent(out)::BIC
  !
  ! Local variables
  real   (kind=8),dimension(n)::z,sz,wz
  real   (kind=8),dimension(2,1)::newpars
  integer(kind=4)::i
  !
  ! Change data to log-scale
  z=dlog(ED)
  sz=dsqrt((Error/ED)**2+spreadsigma**2)
  ! Guess initial sigma
  pars(1,1)=0.1D+00
  !
  wz=1.0D+00/((pars(1,1))**2+sz**2)
  ! Guess initial Central Dose value
  pars(2,1)=sum(z*wz)/sum(wz)
  !
  ! Estimate sigma and Central Dose iteratively
  do i=1,maxiter
    !
    newpars(2,1)=sum(z*wz)/sum(wz)
    newpars(1,1)=pars(1,1)*dsqrt(sum((wz**2)*(z-newpars(2,1))**2/sum(wz)))
    wz=1.0D+00/((newpars(1,1))**2+sz**2)
    !
    if (abs(pars(1,1)-newpars(1,1))+&
        abs(pars(2,1)-newpars(2,1))<eps)  then
       goto 100
     else
       pars=newpars
     end if
     !
   end do
  !
  ! Calculate maximum logged likelihood value
  100 maxlik=0.5D+00*sum(dlog(wz))-0.5D+00*sum(wz*(z-pars(2,1))**2)
  ! Calculate the BIC value
  BIC=-2.0D+00*maxlik-2.0D+00*dlog(real(n,kind=8))
  ! Transform Central Dose
  pars(2,1)=dexp(pars(2,1))
  ! Transform Std.Error of Central Dose
  spars(2,1)=pars(2,1)/dsqrt(sum(wz))
  ! Transform Std.Error of Sigma
  spars(1,1)=1.0D+00/dsqrt(2.0D+00*pars(1,1)*sum(wz**2))
  !
  return
end subroutine CAM  
