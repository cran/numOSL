subroutine lmfit2(xdat,ydat,ndat,pars,npars,calErr,&
                  stderror,predtval,value,tol,info)
!---------------------------------------------------------------------------------------------------------------------
! lmfit1() is used to fit a dose-response cuerve in OSL dating, linear,  
! Exponential or linear+Exponential models are available. (origin).
! ====================================================================================================================
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
! calErr,               input:: logical, calculate parameters' std.error or not.
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
!                                1) if lmfit1 can be called successfully, info(1)=123;
!                                2) if subroutine lmder1 is called with error, info(1) will be one in (0,4,5,6,7);
!                                3) if pars' std.errors can be estimated, info(2)=0;
!                                4) if error when calling lmhess to approximate model's hessian matrix, info(2)=1;
!                                5) if hessian can not be inversed, info(2)=2;
!                                6) if any diagnal elements of inversed hessian matrix is below zero, info(2)=3.
!
! Note:: if npars=1, a linear model of the form: y=ax will be fitted, where ndat>=1;
!        if npars=2, a Exponential model of the form: y=a(1-exp(-bx))  will be fitted, where ndat>=2;
!        if npars=3, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx will be fitted, where ndat>=3.
! ======================================================================================================================
! Author:: Peng Jun, 2013.09.20.
!
! Dependence:: subroutine growfunc2; subroutine lmder1; subroutine lmhess; subroutine inverse; subroutine diag.
!-----------------------------------------------------------------------------------------------------------------------
  
  implicit none
  integer(kind=4),                 intent(in)::ndat
  integer(kind=4),                 intent(in)::npars
  real   (kind=8),                 intent(in)::tol
  real   (kind=8),dimension(ndat), intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat), intent(out)::predtval
  real   (kind=8),                 intent(out)::value
  integer(kind=4),                 intent(out)::info(2)
  logical, intent(in):: calErr
  ! Variables for subroutine lmhess
  real   (kind=8),parameter::lmtol=2.220446D-16   ! maximum tolerance for singular matrix diagnosation (.Machine$double.eps in R)
  real   (kind=8),parameter::minAbsPar=0.0D+00    ! used in lmhess 
  real   (kind=8),dimension(npars,npars)::hessian ! hessian matrix obtained with finite-difference approximation
  real   (kind=8),dimension(npars)::gradient      ! gradient obtained with finite-difference approximation
  integer(kind=4),dimension(5)::hesserror         ! error message generated in subroutine lmhess
  real   (kind=8)::hessvalue                      ! value in subroutine lmhess
  ! variables for subroutine inverse
  integer(kind=4)::inverror                       ! error message generated in subroutine inverse
  ! Variables for subroutine growfunc2 and lmder1
  integer(kind=4)::ldfjac,iflag
  real   (kind=8),dimension(ndat)::fvec           ! fitted residual in vector form
  real   (kind=8),dimension(ndat,npars)::fjac     ! jacobian matrix
  integer(kind=4)::lmerr
  ! 
  info=0
  ldfjac=ndat
  ! Default stderrors if error appears
  stderror=-99.0D+00
  ! specify iflag to be 1 to calculate fvec
  iflag=1
  !
  ! Using initial pars to caculate fvec, 
  ! which will be used in subroutine lmder1
  call growfunc2(xdat,ydat,ndat,npars,pars,&
                 fvec,fjac,ldfjac,iflag)
  !
  ! Optimizing initial pars using Levenberg-Marquadt
  ! method and return (pars)
  call lmder1(growfunc2,ndat,npars,pars,fvec,&
              fjac,ldfjac,tol,lmerr,xdat,ydat)
  ! Calculate sum of squre of residual (value) 
  value=sum((fvec)**2)
  ! calculate fitted values correspond to ydat (predtval)
  predtval=fvec+ydat
  !
  ! Check if any error appears when calling lmder1, 
  ! set info(1)=lmerr if so, else reset info(1) to 123 
  if(lmerr==1 .or. lmerr==2 .or. lmerr==3 ) then
    info(1)=123
  else 
    info(1)=lmerr
  end if  
  ! If need not calculate Std.pars, return now
  if(calErr .eqv. .FALSE.) return
  ! Estimate pars' standard errors
  call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
              hessian,gradient,hessvalue,hesserror,npars+5)
  ! Check if any error appears when calling lmhess, 
  ! if so, reset info(2) to be 1 and return
  if(any(hesserror/=0))  then
    info(2)=1
    return
  end if
  !
  ! Set hessian to be its inverse 
  call inverse(hessian,npars,inverror,lmtol)
  ! Check if any error appears when calling inverse, 
  ! if so, reset info(2) to be 2 and return
  if(inverror==1)  then
    info(2)=2
    return
  end if
  ! Extract diagnal elements from inversed hessian 
  ! matrix and storing it in arrary stderror
  call diag(hessian,npars,stderror)
  ! Check if any elements in stderror is 
  ! below zero, if so, set info(2) to be 3
  if(any(stderror<0)) info(2)=3
  ! Rest stderror
  stderror=dsqrt(stderror)
  ! Now return
  return
end subroutine lmfit2
