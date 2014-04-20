subroutine inipars(bvalue,model,ndat,xdat,&
                   ydat,outpars,tol,info)
!--------------------------------------------------------------------------------------------
! Subroutine inipars() is used to initialize parameters for a growth curve with the 
! formla y=k*(1-exp(-b*x))+c or y=k*(1-exp(-b*x))+c*x+d, with a privde b values and some
! paired observations, k, c (or k, c, d) can be estimated using a Linear Algebra method.
! ===========================================================================================
! 
! model,                 input:: integer, model=3 for y=k*(1-exp(-b*x))+c;
!                                         model=4 for y=k*(1-exp(-b*x))+c*x+d.
!
! ndat,                  input:: integer, length of the observations.
!
! bvalue,                input:: integer, initial b value.
!
! tol,                   input:: real value, maximum tolerance for identify linear independent.
!
! xdat(ndat),            input:: real values, paired observations(x).
!
! ydat(ndat),            input:: real values, paired observations(y).
!
! outpars(model-1),     output:: real values, estimated parameters (k and c, or k, c and d).
!
! info,                 output:: integer value, error message: 
!                                1) if vectors are linear indenpendent, info=0;
!                                2) if vectors are linear denpendent,   info=1.
! ===========================================================================================
! Author:: Peng Jun, 2013.08.01; revised in 2014.04.02.
!
! Dependence:: subroutine GJordan.
!---------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                   intent(in)::model
  integer(kind=4),                   intent(in)::ndat
  real   (kind=8),                   intent(in)::bvalue
  real   (kind=8),                   intent(in)::tol
  real   (kind=8),dimension(ndat),   intent(in)::xdat, ydat
  real   (kind=8),dimension(model-1),intent(out)::outpars
  integer(kind=4),                   intent(out)::info
  ! Local variables.
  real   (kind=8),dimension(ndat,model-1)   ::coeff
  real   (kind=8),dimension(ndat,3)         ::coeff3
  real   (kind=8),dimension(model-1,model-1)::coeffmatrix
  real   (kind=8),dimension(model-1,1)      ::coeffy
  real   (kind=8),dimension(ndat,1)         ::cydat
  !
  info=0
  coeff3=-99.0D+00
  ! Return -99.0 if error appears.
  outpars=-99.0D+00
  !
  if(model==3) then
    coeff3(:,1)=1.0D+00-dexp(-bvalue*xdat)
    coeff3(:,2)=1.0D+00
  else if(model==4) then
    coeff3(:,1)=1.0D+00-dexp(-bvalue*xdat)
    coeff3(:,2)=xdat
    coeff3(:,3)=1.0D+00
  end if
  ! 
  cydat(:,1)=ydat
  coeff(:,1:model-1)=coeff3(:,1:model-1)
  ! Calculate coefficients.
  coeffmatrix=matmul(transpose(coeff),coeff)
  coeffy=matmul(transpose(coeff),cydat)
  call GJordan(coeffmatrix,coeffy,model-1,1,info,tol)
  ! Error checking.
  if(info/=0) return
  !
  outpars=coeffy(:,1)
  return
end subroutine inipars
