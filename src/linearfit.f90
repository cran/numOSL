subroutine linearfit(xdat,ydat,ndat,pars,npars,&
                     stderror,predtval,value,info)
!------------------------------------------------------------------------------------------
! linearfit() is used to fit a growth curve of type linear (for both origin and non-origin). 
! =========================================================================================
!
! xdat(ndat),           input:: real values, the dose values (De[Gy]).
!
! ydat(ndat),           input:: real values, the signal values (Lx/Tx).
!
! ndat,                 input:: integer, the length of dose values.
!
! pars(npars),         output:: real values, the estimated parameters.
!
! npars,                input:: integer, the length of parameters (either 1 or 2).
!
! stderror(npars),     output:: real values, estimated std.errors.
!
! predtval(ndat),      output:: real values, fitted values correspond to ydat.
!
! value,               output:: real value, sum of squared residuals.
!
! info(2),             output:: integer, error message generated during the calling:
!                        1.1) if successed in estimating parameters, info(1)=123;
!                        1.2) if error generated in calculating parameters, info(1)=0;
!                        2.1) if parameters' std.errors can be estimated successfully, info(2)=0;
!                        2.2) if fail in calling lmhess to approximate model's hessian matrix, info(2)=1;
!                        2.3) if fail in inverse the hessian maxtrix, info(2)=2;
!                        2.4) if any diagnal elements of inversed hessian matrix is below zero, info(2)=3.
! ===========================================================================================================
! Note:: if npars=1, a linear model of the form: y=ax will be fitted, where ndat>=1, (origin).
!        if npars=2, a linear model of the form: y=ax+b will be fitted, where ndat>=2, (non-origin).
!
! Author:: Peng Jun, 2013.09.21; revised in 2014.04.03.
!
! Dependence:: subroutine lmhess; subroutine inverse; subroutine diag.
!------------------------------------------------------------------------------------------------------------
  
  implicit none
  integer(kind=4),                 intent(in)::ndat
  integer(kind=4),                 intent(in)::npars
  real   (kind=8),dimension(ndat), intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(out)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat), intent(out)::predtval
  real   (kind=8),                 intent(out)::value
  integer(kind=4),                 intent(out)::info(2)
  !
  ! Variables for subroutine lmhess
  real   (kind=8),parameter::lmtol=2.220446D-16   !.Machine$double.eps in R        
  real   (kind=8),parameter::minAbsPar=0.0D+00    ! used in lmhess 
  real   (kind=8),dimension(npars,npars)::hessian ! hessian matrix obtained with finite-difference approximation
  real   (kind=8),dimension(npars)::gradient      ! gradient obtained with finite-difference approximation
  integer(kind=4),dimension(5)::hesserror         ! error message generated in subroutine lmhess
  real   (kind=8)::hessvalue                      ! value in subroutine lmhess
  ! Variables for subroutine inverse
  integer(kind=4)::inverror                       ! error message generated in subroutine inverse
  ! Local variables
  real   (kind=8),dimension(ndat)::xxdat
  integer(kind=4)::lmhessmodel
  !
  pars=-99.0D+00
  stderror=-99.0D+00
  predtval=-99.0D+00
  value=-99.0D+00
  info(1)=123
  info(2)=0
  ! Calculate pars
  if(npars==1) then
    if(sum(xdat**2)<=lmtol) then
      info(1)=0
      return
    end if
    pars(1)=sum(xdat*ydat)/sum(xdat**2)
    predtval=pars(1)*xdat
    value=sum((ydat-predtval)**2)
  else if(npars==2) then
    if(sum(xdat**2)<=lmtol) then
      info(1)=0
      return
    end if
    xxdat=xdat-sum(xdat)/real(ndat,kind=8)
    pars(1)=sum(xxdat*ydat)/sum(xxdat**2)
    pars(2)=(sum(ydat)-sum(xdat)*pars(1))/real(ndat,kind=8)
    predtval=pars(1)*xdat+pars(2)
    value=sum((ydat-predtval)**2)
  end if
  !
  ! Set lmhessmodel
  if(npars==1) then
    lmhessmodel=6
  else if(npars==2) then
    lmhessmodel=2
  end if
  ! Calculate stderror
  call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
              hessian,gradient,hessvalue,hesserror,lmhessmodel)
  ! Check if any error appears when calling lmhess, 
  ! if so, reset info to be 1 and return
  if(any(hesserror/=0))  then
    info(2)=1
    return
  end if
  !
  ! Set hessian to be its inverse 
  call inverse(hessian,npars,inverror,lmtol)
  ! Check if any error appears when calling inverse, 
  ! if so, reset info to be 2 and return
  if(inverror/=0)  then
    info(2)=2
    return
  end if
  ! Extract diagnal elements from inversed hessian 
  ! matrix and storing it in arrary stderror
  call diag(hessian,npars,stderror)
  ! Check if any elements in stderror is 
  ! below zero, if so, set info to be 3 and return
  if(any(stderror<0)) then 
    info(2)=3
  end if
  ! Rset stderror
  stderror=dsqrt(stderror)
  ! Terminate and return
  return
end subroutine linearfit
