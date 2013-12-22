subroutine inipars(bvalue,model,ndat,xdat,&
                   ydat,outpars,tol,info)
!--------------------------------------------------------------------------------------------
! Subroutine inipars() is used to initialize parameters for dose-response model with the 
! formla y=k*(1-exp(-b*x))+c or y=k*(1-exp(-b*x))+c*x+d, user need to privde b values and 
! paired observations, then k, c (or k, c, d) can be estimated using Linear Algebra method.
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
! info(2),              output:: integer values, error messages: 
!                                1.1) if no error appears in array allocation, info(1)=0;
!                                1.2) if error appears in array allocation, info(1)=1;
!                                2.1) if vectors are linear indenpendent, info(2)=0;
!                                2.2) if vectors are linear denpendent, info(2)=1.
! ===========================================================================================
! Author:: Peng Jun, 2013.08.01.
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
  integer(kind=4),                   intent(out)::info(2)
  ! Local variables
  real   (kind=8),dimension(:,:),allocatable::coeff
  real   (kind=8),dimension(:,:),allocatable::coeffmatrix
  real   (kind=8),dimension(:,:),allocatable::coeffy
  real   (kind=8),dimension(ndat,1)::cydat
  !
  info=0
  ! Return -99.0 if error appears
  outpars=-99.0D+00
  if(model==3) then
    allocate(coeff(ndat,2),coeffmatrix(2,2),coeffy(2,1),stat=info(1))
    if(info(1)/=0) return
    coeff(:,1)=1.0D+00-dexp(-bvalue*xdat)
    coeff(:,2)=1.0D+00
  else if(model==4) then
    allocate(coeff(ndat,3),coeffmatrix(3,3),coeffy(3,1),stat=info(1))
    if(info(1)/=0) return
    coeff(:,1)=1.0D+00-dexp(-bvalue*xdat)
    coeff(:,2)=xdat
    coeff(:,3)=1.0D+00
  end if
  ! Calculate coefficients
  cydat(:,1)=ydat
  coeffmatrix=matmul(transpose(coeff),coeff)
  coeffy=matmul(transpose(coeff),cydat)
  call GJordan(coeffmatrix,coeffy,model-1,1,info(2),tol)
  if(info(2)/=0) return
  outpars=coeffy(:,1)
  deallocate(coeff,coeffmatrix,coeffy)
  ! Now return
  return
end subroutine inipars
