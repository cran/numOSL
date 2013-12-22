subroutine interpolate2(Dose,ltx,pars,npars,lowb,upb,value)
!-------------------------------------------------------------------------------------------------------------------
! interpolate2() is a subroutine used for interpolating a Dose value from a Dose-Response Curve (origin).
! ==================================================================================================================
! 
! npars,        input:: integer, length of parameters (1, or 2, or 3).
!
! pars(npars),  input:: real values, characteristic values for a Dose-Response curve.
!
! Dose,        output:: real value, calculated equivalent dose value.
!
! value,       output:: real value, a minimized object value.
!
! ltx,          input:: real value, the standardlized signal value from which a Dose is to be estimated.
!
! lowb,         input:: rea lvalue, low boundary of a interval from which the interpolation to take place.
!
! upb,          input:: real value, up boundary of a interval from which the interpolation to take place.
! ==================================================================================================================
! Author:: Peng Jun, 2013.09.20.
!
! Dependence:: function fmin2.
!--------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                 intent(in)::npars
  real   (kind=8),                 intent(in)::ltx
  real   (kind=8),                 intent(in)::lowb,upb
  real   (kind=8),dimension(npars),intent(in)::pars
  real   (kind=8),                 intent(out)::Dose
  real   (kind=8),                 intent(out)::value
  !
  ! local variables
  real   (kind=8),parameter::gtol=1.490116e-08  ! .Machine$double.eps^0.5 in R
  real   (kind=8),dimension(3)::cpars
  real   (kind=8)::fmin2
  !
  cpars=0.0D+00
  cpars(1:npars)=pars
  !
  ! initialize Dose and value
  Dose=0.0D+00 
  value=0.0D+00
  ! 
  ! calculate Dose
  Dose=fmin2(lowb,upb,gtol,ltx,cpars,npars)
  !
  ! calculate value
  if(npars==1) then
    value=(cpars(1)*Dose-ltx)**2
  else if(npars==2) then
    value=(cpars(1)*(1.0D+00-dexp(-cpars(2)*Dose))-ltx)**2
  else if(npars==3) then
    value=(cpars(1)*(1.0D+00-dexp(-cpars(2)*Dose))+cpars(3)*Dose-ltx)**2
  end if
  return
end subroutine interpolate2
