subroutine aperr(ED,Error,nED,pars,sigma,spars,npars,tol,message)
!--------------------------------------------------------------------------------------------------------------------
! Subroutine aperr() is used to approximate parameters' standard errors in a minimum age model.
! ===================================================================================================================
!
! nED,           input:: integer, length of equivalent dose values.
!
! npars,         input:: integer, number of parameters, either 3 or 4.
!
! ED(nED),       input:: real values, equivalent dose values (unlogged).
!
! Error(nED),    input:: real values, absolute errors of equivalent dose values.
!
! pars(npars),   input:: real values, parameters used to estimate corresponded std.errors of parameters.
!
! spars(npars), output:: real values, estimated std.errors.
!
! tol,           input:: real values, the allowed maximum tolerance for hessian maxtrix to be diagnosed as non-sigular.
!
! message,      output:: integer, error message generated during the calculation:
!                        1.1) if standard errors can be estimated, message=0;
!                        1.2) if standard errors can not be estimated, message=1.
! ====================================================================================================================
!     Author:: Peng Jun, 2013.07.24.
!
! Dependence:: subroutine fdhessian; subroutine inverse; subroutine diag.
!---------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                 intent(in)::nED
  integer(kind=4),                 intent(in)::npars
  real   (kind=8),dimension(nED),  intent(in)::ED
  real   (kind=8),dimension(nED),  intent(in)::Error
  real   (kind=8),dimension(npars),intent(in)::pars
  real   (kind=8),                 intent(in)::sigma
  real   (kind=8),                 intent(in)::tol
  real   (kind=8),dimension(npars),intent(out)::spars
  integer(kind=4),                 intent(out)::message
  !
  ! Arguments for subroutine fdhessian
  real   (kind=8), parameter:: minAbsPar=0.0D+00
  real   (kind=8),dimension(npars,npars)::hessian
  real   (kind=8),dimension(npars)::grad
  integer(kind=4),dimension(4)::errorflag
  real   (kind=8)::value
  !
  ! Arguments for subroutine inverse
  integer(kind=4)::ifault
  !
  ! Local variables
  real   (kind=8),dimension(nED)::sED,sError
  !
  ! Return all output to be -99.0 if error appears
  spars=-99.0D+00
  message=0
  ! 
  ! Add spread sigma value to relative errors of EDs
  sError=dsqrt((Error/ED)**2+sigma**2)
  ! Transform ED data to log-scale
  sED=dlog(ED)
  !
  ! Call fdhessian to approximate hessian matrix 
  ! of minimum age models with provided parameters
  call fdhessian(pars,sED,sError,npars,nED,tol,&
                 minAbsPar,hessian,grad,value,errorflag)
  ! Error checking of fdhessian
  if(any(errorflag/=0))   then
    message=1
    return
  end if
  !
  ! Inverse the hessian matrix to estimate the
  ! standard errors of provided parameters 
  call inverse(hessian,npars,ifault,tol)
  !
  ! Check if error presents during inversing hessian
  if(ifault/=0)  then
    message=1
    return
  end if
  !
  ! Extract diagnal elements of the 
  ! inversed hessian matrix
  call diag(hessian,npars,spars)
  !
  ! Check if any diagnal elment of 
  ! inversed hessian is below zero
  if(any(spars<0.0D+00)) then
    message=1
    return
  end if
  ! Further transform spars
  spars=dsqrt(spars)
  ! 
  return
end subroutine aperr
