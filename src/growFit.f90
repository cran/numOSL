subroutine growFit(xdat,ydat,ndat,pars,npars,origin,&
                   stderror,predtval,value,tol,info)
!----------------------------------------------------------------------------------------
! growFit() is used to fit a dose-response curve, linear,Exponential,linear+Exponential
! models are available (that pass the origin or not can be set).
! =======================================================================================
!
! xdat(ndat),        input:: real values, the dose values (De[Gy]).
!
! ydat(ndat),        input:: real values, the standardlized signal values (Lx/Tx).
!
! ndat,              input:: integer, the length of dose values (2,3 for origin and 3,4 for non-origin).
!
! origin,            input:: logical, pass the origin or not.
!
! pars(npars),input/output:: real values, the initials, will be overwritten for output.
!
! npars,             input:: integer, the length of parameters.
!
! stderror(npars),  output:: real values, estimated std.error of parameters.
!
! predtval(ndat),   output:: real values, fitted values correspond to ydat.
!
! value,            output:: real value, sum of squared residuals.
!
! tol,               input:: real value, termination occurs when the algorithm estimates either
!                            that the relative error in the sum of squares is at most TOL or that
!                            the relative error between X and the solution is at most TOL.
!
! info(2),          output:: integer, error message generated during the calling of lmfit1:
!                   1) if parameters can be estimated, info(1)=123, or else info(1) in (0,4,5,6,7);
!                   2) if parameters' std.errors can be estimated, info(2)=0, or else info(2) in (1,2,3).
! =============================================================================================================
! Note:: for origin=FALSE:
!        if npars=3, a Exponential model of the form: y=a(1-exp(-bx))+c  will be fitted, where ndat>=3;
!        if npars=4, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx+d will be fitted, where ndat>=4.
!        and for origin=TRUE:
!        if npars=2, a Exponential model of the form: y=a(1-exp(-bx))  will be fitted, where ndat>=2;
!        if npars=3, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx will be fitted, where ndat>=3.
!
! Author:: Peng Jun, 2013.05.26, revised in 2013.08.04; revised in 2014.04.02.
!
! Dependence:: subroutine growfunc; subroutine growfunc0(); subroutine lmder1; 
!              subroutine lmhess; subroutine inverse; subroutine diag.
!--------------------------------------------------------------------------------------------------------------
  
  implicit none
  integer(kind=4),                    intent(in)::ndat
  integer(kind=4),                    intent(in)::npars
  real   (kind=8),                    intent(in)::tol
  logical,                            intent(in)::origin
  real   (kind=8),dimension(ndat),    intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),  intent(out)::stderror
  real   (kind=8),dimension(ndat),   intent(out)::predtval
  real   (kind=8),                   intent(out)::value
  integer(kind=4),dimension(2),      intent(out)::info
  !
  ! Local variables.
  real   (kind=8),parameter             ::lmtol=1.490116e-08  !.Machine$double.eps^0.5 in R  
  real   (kind=8),parameter             ::minAbsPar=0.0D+00     
  real   (kind=8),dimension(npars,npars)::hessian  
  real   (kind=8),dimension(npars)      ::gradient       
  integer(kind=4),dimension(5)          ::hesserror          
  real   (kind=8)                       ::hessvalue                       
  integer(kind=4)                       ::inverror                        
  integer(kind=4)                       ::iflag
  real   (kind=8),dimension(ndat)       ::fvec            
  real   (kind=8),dimension(ndat,npars) ::fjac      
  integer(kind=4)                       ::lmerr
  ! 
  info=0
  ! Default stderrors.
  stderror=-99.0D+00
  ! Specify iflag to be 1 to calculate fvec.
  iflag=1
  !
  ! Use initial parameters to caculate fvec, 
  ! which will be used in subroutine lmder1()
  if (origin .eqv. .false.) then
    call growfunc(xdat,ydat,ndat,npars,pars,&
                  fvec,fjac,ndat,iflag)
  else 
    call growfunc0(xdat,ydat,ndat,npars,pars,&
                   fvec,fjac,ndat,iflag)
  end if
  !
  ! Optimize initial pars using Levenberg-Marquadt
  ! method and return parameters(pars)
  if (origin .eqv. .false.) then
    call lmder1(growfunc,ndat,npars,pars,fvec,&
                fjac,ndat,tol,lmerr,xdat,ydat)
  else 
    call lmder1(growfunc0,ndat,npars,pars,fvec,&
                fjac,ndat,tol,lmerr,xdat,ydat)
  end if
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
  ! 
  ! Estimate parameters' standard errors
  if (origin .eqv. .false.) then
    call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
                hessian,gradient,hessvalue,hesserror,npars)
  else 
    call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
                hessian,gradient,hessvalue,hesserror,npars+5)
  end if
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
end subroutine growFit
