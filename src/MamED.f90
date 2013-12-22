subroutine MamED(ED,Error,nED,pars,spars,value,npars,&
                 sigma,maxiter,tol,bound,message)
!--------------------------------------------------------------------------------------------------------------------------
! MamED() attempts to estimate the values of the parameters in Minimum Age Models
! by using the L-BFGS-B method, insteading for a globle minimum, the L-BFGS-B routine
! results in a local minimum, and sometimes even end with nonsense, suggesting what
! a ill-condition problem is the mimimum age models! 
! The results rely strongly on the initial guess values, but in this subroutine, by
! adjusting the sigma, converge can be achieved in most situations.
! Note that this subroutine do logged equivalent dose analysis, which means that all the
! equivalent dose data must be larger than 0.
! =========================================================================================================================
! ED(nED), input            :: real values, the equivalent dose, must be of un-logged.
!
! Error(nED), input         :: real values, the associated absolute errors for the equivalent dose.
!
! nED, input                :: integer, the size of the equivalent dose data.
!
! npars, input              :: integer,the dimension of the problems, 3 for MAM3, 4 for MAM4.
!
! pars(npars), input/output :: real values, the estimated parameters.
!
! spars(npars), output      :: real values, the estimated standard errors for parameters.
!
! value , output            :: real value, the estimated minus logged maximum logged likelihood value.
!
! sigma, input              :: real value, the added spread to the equivalent data.
!
! maxiter, input            :: integer, the allowed maximum iterative numbers.
!
! tol, input                :: real value, the allowed maximum tolerance for the hessian to be non-sigular.
!
! bound(2,npars), input     :: real values, the low and up boundary for L-BFGS-B method.
!
! message(5),output         :: integer values, the error message generated in the analysis:
!                              1.1) if gradient estimation in l-bfgs-b using subroutine gradient() success, message(1)=0;
!                              1.2) if gradient estimation in l-bfgs-b using subroutine gradient() fail, message(1)=1;
!                              2.1) if hessian matrix can be approximated and inversed, message(2)=0;
!                              2.2) if hessian matrix can not be approximated and inversed, message(2)=1;
!                              3.1) if parameters' standard errors can be approximated, message(3)=0;
!                              3.2) if parameters' standard errors can not be approximated, message(3)=1;
!                              4.1) if no estimated parameter is near the boundary, message(4)=0;
!                              4.2) if any estimated parameter is near the boundary, message(4)=1;
!                              5.1) if for MAM4, estimated value gama<=mu, message(5)=0;
!                              5.2) if for MAM4, estimated value gama>mu, message(5)=1.
! =============================================================================================================================                           
! Author :: Peng Jun, 2013.03.15, revised in 2013.04.01.
! 
! Reference :: Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of
!              single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design
!              and statistical models. Archaeometry, 41, pp. 339-364.
!
! Dependence :: subroutine setulb; subroutine fdhessian; subroutine gradient;
!               subroutine inverse; subroutine diag.
!------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                   intent(in)::nED
  integer(kind=4),                   intent(in)::npars
  integer(kind=4),                   intent(in)::maxiter
  real   (kind=8),                   intent(in)::tol
  real   (kind=8),                   intent(in)::sigma
  real   (kind=8),dimension(2,npars),intent(in)::bound
  real   (kind=8),dimension(nED),    intent(in)::ED
  real   (kind=8),dimension(nED),    intent(in)::Error
  real   (kind=8),                   intent(out)::value
  integer(kind=4),dimension(5),      intent(out)::message
  real   (kind=8),dimension(npars),  intent(inout)::pars
  real   (kind=8),dimension(npars),  intent(out)::spars
  !
  ! Variables for subroutine setulb
  integer(kind=4),parameter::m=5, iprint=-1
  real   (kind=8),parameter::factr=1.0D+07, pgtol=0.0D+00
  logical,dimension(4)::lsave
  integer(kind=4):: nbd(npars),&
                    iwa(3*npars), isave(44)
  real   (kind=8):: f, l(npars), &
                    u(npars), g(npars), dsave(29), &
                    wa(2*m*npars+5*npars+11*m*m+8*m) 
  character(len=60)::task,csave 
  !
  ! Variables for subroutine fdhessian
  real   (kind=8),parameter::minAbsPar=0.0D+00     
  real   (kind=8)::grad(npars)
  real   (kind=8)::hessian(npars,npars)
  integer(kind=4)::errorflag(4)
  !
  ! Local variables for subroutine gradient
  integer(kind=4),parameter::iter=6
  integer(kind=4)::failure
  !
  ! Local variables for subroutine inverse
  integer(kind=4)::ifault
  !
  ! Local variables 
  real   (kind=8)::sED(nED),sError(nED) 
  !
  ! Transform ED data to log-scale and add spread sigma value
  sError=dsqrt((Error/ED)**2+sigma**2)
  sED=dlog(ED)
  !
  ! Set bounds on the pars
  nbd(1:npars)=2
  !
  l=bound(1,:)
  u=bound(2,:)
  !
  ! Start do optimizing using L-BFGS-B method
  task='START'
  spars=0.0D+00
  message=0
  !
  100 continue
  call setulb(npars,m,pars,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint, &
              csave,lsave,isave,dsave)
  ! Check if the maximum iteration has been reached
  if(isave(30)==maxiter)  goto 200
  !
  if (task(1:2) .eq. 'FG') then
    ! Now provide the gradient and value for the
    ! problem to be optimized
    call gradient(pars,sED,sError,npars,iter,&
                  nED,grad,value,failure)
    ! check if Inf or NaN presents in 
    ! gradient or value
    if( failure==1 )  then
      message(1)=1  
      return
    end if
    !
    g=grad
    f=value
    !
    goto 100
    !
  end if
  !
  if (task(1:5) .eq. 'NEW_X') then 
    !
    if (isave(34) .ge. 99)  then
      task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
    end if
    !
    if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f)))  then
      task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
    end if
    !
    goto 100
    !
  end if
  ! At this point the optimization is over
  200 continue
  !
  ! Estimate the standard errors of parameters
  call fdhessian(pars,sED,sError,npars,nED,tol, &
                 minAbsPar,hessian,grad,value,errorflag)
  ! Check if error appears during the approximation
  if(any(errorflag/=0))   then
    message(2)=1
    return
  end if
  ! Now inverse the hessian matrix to estimate the
  ! standard errors for calculated parameters 
  call inverse(hessian,npars,ifault,tol)
  ! Check if error presents during inversing
  message(2)=ifault
  if(message(2)==1)  return
  !
  ! Extract diagnal elements of the 
  ! inversed hessian matrix
  call diag(hessian,npars,spars)
  !
  ! Check if any diagnal elment of inversed hessian is below zero
  if(any(spars<0.0D+00))  message(3)=1
  if(message(3)==1)  return
  ! 
  spars=dsqrt(spars)
  !
  ! Check if any parameter is near the boundary
  if (any(dabs(pars-l)<1.0D-04) .or. &
      any(dabs(pars-u)<1.0D-04)  )  message(4)=1
  !
  ! Check if gama>mu in MAM4 with four parameters
  if( npars==4 .and. &
      pars(2)>pars(3) )  message(5)=1
  ! 
  return 
  !
end subroutine MamED
