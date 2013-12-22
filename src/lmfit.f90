subroutine lmfit(xdat,ydat,ndat,pars,npars,typ,transf,&
                 stderror,predtval,value,tol,info)
!--------------------------------------------------------------------------------------------------------------------------------
! Specifying a fitting model, lmfit() will fit the model using Levenberg-Marquadt method.
! typ=1 has the formula I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t), where k=1:7;
! typ=2 has the formula I(t)=a1*(t/max(t))*exp(-b1*t^2)/2/max(t))+...+
!                            ak*(t/max(t))*exp(-bk*t^2)/2/max(t)), where k=1:7.
! ===============================================================================================================================
!
! ndat,                input:: integer, length of xdat or y dat.
!
! npars,               input:: integer, dimension of parameters (or length of pars).
!
! typ,                 input:: integer, 1 for fitting 'cw', 2 for fitting 'lm'.
!
! transf,              output:: logical, whether transform se(a[i]) to se(a[i]/b[i]) or not.
!
! xdat(ndat),          input:: real values, xdat (or independent variables x).
!
! ydat(ndat),          input:: real values, ydat (or dependent variables y).
!
! pars(npars),  input/output:: real values, initial guess values for parameters to be optimized,
!                              overwritten to be final results for output.
!
! stderror(npars),    output:: real values, estimated standard errors for parameters.
!
! predtval(ndat),     output:: real values, fited values correspond to ydat.
!
! value,              output:: real value, final total residual error.
!
! tol,                 input:: real value,  termination occurs when the algorithm estimates either that the relative error in the  
!                              sum of squares is at most TOL or that the relative error between X and the solution is at most TOL.
!
! info,               output:: error message generated during the calculation:
!                              1.1) successful work, info=1;
!                              1.2) error when calling lmder1, info=2;
!                              1.3) at least one estimated parameter is below zero, info=3;
!                              1.4) error when calling lmhess to approximate model's hessian matrix, info=4;
!                              1.5) error when attempt to inverse the approximated hessian matrix, info=5;
!                              1.6) if diagnal elements of inversed hessian matrix has at least one minus value, info=6.
! ================================================================================================================================
! Author:: Peng Jun, 2013.05.21, revised in 2013.05.22, revised in 2013.07.24, revised in 2013.10.05.
!
! Dependence:: subroutine lmfunc; subroutine lmfunc1; subroutine lmder1; subroutine lmhess; subroutine inverse; subroutine diag.
!---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ndat
  integer(kind=4),intent(in)::npars
  integer(kind=4),intent(in)::typ
  real   (kind=8),dimension(ndat),intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat),intent(out)::predtval
  real   (kind=8),intent(in)::tol
  logical,        intent(in)::transf
  real   (kind=8),intent(out)::value
  integer(kind=4),intent(out)::info
  !
  ! Variables for subroutine lmhess
  real   (kind=8),parameter::lmtol=2.220446D-16   ! used for singular matrix diagnose in subroutine (.Machine$double.eps in R)
                                                  ! lmhess and subroutine inverse
  real   (kind=8),parameter::minAbsPar=0.0D+00    ! for lmhess use
  real   (kind=8),dimension(npars,npars)::hessian ! hessian matrix by finite-difference approximation
  real   (kind=8),dimension(npars)::gradient      ! gradient by finite-difference approximation
  integer(kind=4),dimension(5)::hesserror         ! error message generated in lmhess
  ! Variables for subroutine inverse
  integer(kind=4)::inverror                       ! for singular matrix diagnose in subroutine inverse
  ! Local variables
  integer(kind=4)::ldfjac,iflag
  real   (kind=8),dimension(ndat)::fvec           ! fitted residual in vector form
  real   (kind=8),dimension(ndat,npars)::fjac     ! jacobian matrix
  integer(kind=4)::i
  !
  ldfjac=ndat
  !
  ! Default returned stderrors if error appears
  stderror=-99.0D+00
  ! Specify iflag to be 1 to calculate fvec
  iflag=1
  !
  ! Using initial pars to caculate fvec, which will be used in subroutine lmder1
  if(typ==1) then
    ! For 'cw'
    call lmfunc(xdat,ydat,ndat,npars,pars,fvec,fjac,ldfjac,iflag)
  else if(typ==2) then
    ! For 'lm'
    call lmfunc1(xdat,ydat,ndat,npars,pars,fvec,fjac,ldfjac,iflag)
  end if
  !
  ! Optimizing initial pars using Levenberg-Marquadt method
  ! and return pars and info for output
  if(typ==1) then
    ! For 'cw'
    call lmder1(lmfunc,ndat,npars,pars,fvec,&
                fjac,ldfjac,tol,info,xdat,ydat)
  else if(typ==2) then
    ! For 'lm'
    call lmder1(lmfunc1,ndat,npars,pars,fvec,&
                fjac,ldfjac,tol,info,xdat,ydat)
  end if
  !
  ! Calculate sum of squre of residual 
  value=sum((fvec)**2)
  ! Calculate fitted values correspond to ydat
  predtval=fvec+ydat
  !
  ! Check if any error appears when calling lmder1, and reset info 
  ! to 2 and return if so, else reset info to 1 and continue
  if(info==1 .or. info==2 .or. info==3) then
    info=1
  else 
    info=2
    return
  end if  
  !
  ! Check if any estimated parameter is below zero,
  ! if so, reset info to be 3 and return
  if(any(pars<0.0D+00)) then
    info=3
    return
  end if
  ! Estimate pars' standard errors with subroutine lmhess
  if(typ==1) then
    call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
                hessian,gradient,value,hesserror,1)
  else if(typ==2) then
    call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
                hessian,gradient,value,hesserror,5)
  end if
  ! Reset value after calling subroutine lmhess
  value=value**2
  ! Check if any error appears when calling lmhess, 
  ! if so, reset info to be 4 and return
  if(any(hesserror/=0))  then
    info=4
    return
  end if
  !
  ! Set hessian to be its inverse form
  call inverse(hessian,npars,inverror,lmtol)
  ! Check if any error appears when calling inverse, 
  ! if so, reset info to be 5 and return
  if(inverror==1)  then
    info=5
    return
  end if
  ! Extract diagnal elements from inversed hessian 
  ! matrix and save it in arrary stderror
  call diag(hessian,npars,stderror)
  ! check if any elements in stderror is 
  ! below zero, if so, set info to be 6 and return
  if(any(stderror<0)) then
    info=6
    return
  end if
  ! Rest stderror
  stderror=dsqrt(stderror)
  !
  ! transform se(a[i]) to se(a[i]/b[i]) or not?
  if(transf .eqv. .true.) then
    ! If transf=true, then std.err of ithn (a1,a2,...) need to be resetted
    do i=1,npars/2
      stderror(i)=pars(i)/pars(i+npars/2)*dsqrt(hessian(i,i)/pars(i)**2+&
                  hessian(i+npars/2,i+npars/2)/pars(i+npars/2)**2-&
                   2.0D+00*hessian(i,i+npars/2)/pars(i)/pars(i+npars/2))
    end do
  end if
  !
  return
end subroutine lmfit
