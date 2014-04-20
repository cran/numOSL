subroutine inipars0(bvalue,model,ndat,xdat,&
                    ydat,outpars,tol,info)
!---------------------------------------------------------------------------------------
! Subroutine inipars0() is used to initialize parameters of a growth curve with the 
! formla y=k*(1-exp(-b*x)) or y=k*(1-exp(-b*x))+c*x, with a privde b values and some
! paired observations, k, (or k, c) can be estimated using a Linear Algebra method.
! ======================================================================================
! 
! model,                 input:: integer, model=2 for y=k*(1-exp(-b*x));
!                                         model=3 for y=k*(1-exp(-b*x))+c*x.
!
! ndat,                  input:: integer, length of the observations.
!
! bvalue,                input:: integer, initial b value.
!
! tol,                   input:: real value, maximum tolerance for identify linear independent.
!
! xdat(ndat),ydat(ndat), input:: real values, paired observations.
!
! outpars(model-1),     output:: real values, estimated parameters (k, or k and c).
!
! info,                 output:: integer value, error message:
!                                1) if work vectors are linear indenpendent, info=0; 
!                                2) if work vectors are linear denpendent,   info=1.
! ==============================================================================================
! Author:: Peng Jun, 2013.09.20; revised in 2014.04.02.
!
! Dependence:: subroutine GJordan.
!-----------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                    intent(in)::model
  integer(kind=4),                    intent(in)::ndat
  real   (kind=8),                    intent(in)::bvalue
  real   (kind=8),                    intent(in)::tol
  real   (kind=8),dimension(ndat),    intent(in)::xdat, ydat
  real   (kind=8),dimension(model-1),intent(out)::outpars
  integer(kind=4),                   intent(out)::info
  ! Local variables.
  real   (kind=8),dimension(ndat,2)::coeff
  real   (kind=8),dimension(2,2)   ::coeffmatrix
  real   (kind=8),dimension(2,1)   ::coeffy
  real   (kind=8),dimension(ndat,1)::cydat
  !
  info=0
  ! Return -99.0 if error appears.
  outpars=-99.0D+00
  !
  if(model==2) then
    outpars=sum((1.0D+00-dexp(-bvalue*xdat))*ydat)/&
            sum((1.0D+00-dexp(-bvalue*xdat))**2)
  else if(model==3) then
    coeff(:,1)=1.0D+00-dexp(-bvalue*xdat)
    coeff(:,2)=xdat
    !
    cydat(:,1)=ydat
    ! Calculate coefficients.
    coeffmatrix=matmul(transpose(coeff),coeff)
    coeffy=matmul(transpose(coeff),cydat)
    call GJordan(coeffmatrix,coeffy,2,1,info,tol)
    ! Error checking.
    if(info/=0) return
    !
    outpars=coeffy(:,1)
  end if
  ! 
  return
end subroutine inipars0
