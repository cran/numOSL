subroutine lmfit1(xdat,ydat,ndat,pars,npars,calErr,&
                  stderror,predtval,value,tol,info)
!---------------------------------------------------------------------------------------------------------------------------
! lmfit1() is used to fit a dose-response curve, linear, Exponential or linear+Exponential models are available (non-origin)
! ==========================================================================================================================
!
! xdat(ndat),           input:: real values, the OSL dose date (De[Gy]).
!
! ydat(ndat),           input:: real values, the standardlized OSL signal data (Lx/Tx).
!
! ndat,                 input:: integer, the length of dose or signal data.
!
! pars(npars),   input/output:: real values, the initial guess pars, will be overwritten for output.
!
! npars,                input:: integer, the length of parameters.
!
! calErr,               input:: logical value, calculate pars' std.error or not.
!
! stderror(npars),     output:: real values, estimated std.error of parameters.
!
! predtval(ndat),      output:: real values, fitted values correspond to ydat.
!
! value,               output:: real value, sum of square of residuals.
!
! tol,                  input:: real value, real value,  Termination occurs when the algorithm estimates either
!                               that the relative error in the sum of squares is at most TOL or that
!                               the relative error between X and the solution is at most TOL.
!
! info(2),             output:: integer, error message generated during the calling of lmfit1:
!                                1.1) if parameters can be estimated, info(1)=123;
!                                1.2) if parameters can not be estimated, info(1) will be one in (0,4,5,6,7);
!                                2.1) if parameters' std.errors can be estimated, info(2)=0;
!                                2.2) if fail in calling lmhess() to approximate the hessian matrix, info(2)=1;
!                                2.3) if the hessian matrix can not be inversed, info(2)=2;
!                                2.4) if any diagnal elements of inversed hessian matrix is below zero, info(2)=3.
! ============================================================================================================================
! Note:: if npars=2, a linear model of the form: y=ax+b will be fitted, where ndat>=2;
!        if npars=3, a Exponential model of the form: y=a(1-exp(-bx))+c  will be fitted, where ndat>=3;
!        if npars=4, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx+d will be fitted, where ndat>=4.
!
! Author:: Peng Jun, 2013.05.26, revised in 2013.08.04.
!
! Dependence:: subroutine growfunc; subroutine lmder1; subroutine lmhess; subroutine inverse; subroutine diag.
!-----------------------------------------------------------------------------------------------------------------------------
  
  implicit none
  integer(kind=4),                 intent(in)::ndat
  integer(kind=4),                 intent(in)::npars
  real   (kind=8),                 intent(in)::tol
  logical,                         intent(in)::calErr
  real   (kind=8),dimension(ndat), intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat), intent(out)::predtval
  real   (kind=8),                 intent(out)::value
  integer(kind=4),dimension(2),    intent(out)::info
  !
  ! Variables for subroutine lmhess()
  real   (kind=8),parameter::lmtol=2.220446D-16      ! maximum tolerance for singular matrix diagnosation (.Machine$double.eps in R)
  real   (kind=8),parameter::minAbsPar=0.0D+00       ! used in lmhess 
  real   (kind=8),dimension(npars,npars)::hessian    ! hessian matrix obtained with finite-difference approximation
  real   (kind=8),dimension(npars)::gradient         ! gradient obtained with finite-difference approximation
  integer(kind=4),dimension(5)::hesserror            ! error message generated in subroutine lmhess
  real   (kind=8)::hessvalue                         ! value in subroutine lmhess
  !
  ! Variables for subroutine inverse()
  integer(kind=4)::inverror                          ! error message generated in subroutine inverse
  ! Variables for subroutine growfunc() and lmder1()
  integer(kind=4)::ldfjac,iflag
  real   (kind=8),dimension(ndat)::fvec              ! fitted residual in vector form
  real   (kind=8),dimension(ndat,npars)::fjac        ! jacobian matrix
  integer(kind=4)::lmerr
  ! 
  ldfjac=ndat
  ! Default stderrors if encoutering any error during the calculation
  stderror=-99.0D+00
  ! Specify iflag to be 1 to calculate fvec
  iflag=1
  !
  ! Use initial parameters to caculate fvec, 
  ! which will be used in subroutine lmder1()
  call growfunc(xdat,ydat,ndat,npars,pars,&
                fvec,fjac,ldfjac,iflag)
  !
  ! Optimize initial pars using Levenberg-Marquadt
  ! method and return parameters(pars)
  call lmder1(growfunc,ndat,npars,pars,fvec,&
              fjac,ldfjac,tol,lmerr,xdat,ydat)
  ! Calculate sum of squre of residual (value) 
  value=sum((fvec)**2)
  ! Calculate fitted values correspond to ydat (predtval)
  predtval=fvec+ydat
  !
  ! check if any error appears during calling lmder1(), 
  if(lmerr==1 .or. lmerr==2 .or. lmerr==3 ) then
    ! Successed
    info(1)=123
  else 
    ! Failed
    info(1)=lmerr
  end if  
  ! If don't ask to calculate Std.pars, return
  if(calErr .eqv. .FALSE.) return
  ! Estimate parameters' standard errors
  call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
              hessian,gradient,hessvalue,hesserror,npars)
  ! Check if any error appears when calling lmhess()
  if(any(hesserror/=0))  then
    info(2)=1
    return
  end if
  !
  ! Set hessian to be its inverse 
  call inverse(hessian,npars,inverror,lmtol)
  ! Check if any error appears during calling inverse()
  if(inverror==1)  then
    info(2)=2
    return
  end if
  ! Extract diagnal elements from the inversed hessian 
  ! matrix and storing them in arrary stderror
  call diag(hessian,npars,stderror)
  ! Check if any elements in stderror is below zero
  if(any(stderror<0)) then
    info(2)=3
    return
  end if
  ! Reset stderror
  stderror=dsqrt(stderror)
  ! Terminate and return
  return
end subroutine lmfit1
