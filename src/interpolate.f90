subroutine interpolate(Dose,ltx,pars,npars,lowb,upb,value)
!-------------------------------------------------------------------------------------------------------------------
! Interpolate() is a subroutine used for interpolating a Dose value from a Dose-Response Curve. (Non-origin).
! ==================================================================================================================
!
! npars,        input:: integer, length of parameters.
!
! pars(npars),  input:: real values, characteristic values for a Dose-Response curve.
!
! Dose,        output:: real value, calculated Dose value.
!
! value,       output:: real value, a minimized object.
!
! ltx,          input:: real value, the standardlized signal value from which a Dose is to be estimated.
!
! lowb,         input:: rea lvalue, low boundary of a interval from which the interpolation to take place.
!
! upb,          input:: real value, up boundary of a interval from which the interpolation to take place.
! ===================================================================================================================
! Author:: Peng Jun, 2013.06.22.
!
! Dependence:: function fmin.
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
  real   (kind=8),dimension(4)::cpars
  real   (kind=8)::fmin
  !
  cpars=0.0D+00
  cpars(1:npars)=pars
  !
  ! initialize Dose and value
  Dose=0.0D+00 
  value=0.0D+00
  ! 
  ! calculate Dose
  Dose=fmin(lowb,upb,gtol,ltx,cpars,npars)
  !
  ! calculate value
  if(npars==2) then
    value=(cpars(1)*Dose+cpars(2)-ltx)**2
  else if(npars==3) then
    value=(cpars(1)*(1.0D+00-dexp(-cpars(2)*Dose))+cpars(3)-ltx)**2
  else if(npars==4) then
    value=(cpars(1)*(1.0D+00-dexp(-cpars(2)*Dose))+cpars(3)*Dose+cpars(4)-ltx)**2
  end if
  return
end subroutine interpolate
