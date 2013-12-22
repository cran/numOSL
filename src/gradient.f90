subroutine gradient(pars,ED,Error,npars,iter,&
                    nED,grad,value,errorflag)
!--------------------------------------------------------------------------------------
! Subroutine gradient() is used to calculate the gradients of
! specified parameters that belong to a function of interest
! while this subroutine constraints the function to be fun34
! =====================================================================================
!
! npars,         input:: integer, the length of the parameters to be analyzed/
!
! nED,           input:: integer, the length of equivalent doses.
!
! iter,          input:: integer,  the number of Richardson improvement iterations.
!
! pars(npars),   input:: real values, the parameters to be analyzed.
!
! ED(nED),       input:: real values, the equivalent dose.
!
! Error(nED),    input:: real values, the standard errors of equivalent dose.
!
! grad(npars),  output:: real values, the calculated gradient.
!
! value,        output:: real value, the caluculated value for the specified function.
!
! errorflag,    output:: integer, error message:
!                        1.1) if both value and gradient can be calculated,errorflag=0;
!                        1.2) if value or gradient can not be calculated, errorflag=1.
! ====================================================================================
! Author:: Peng Jun, 2013.04.01, revised in 2013.04.21.
!
! Dependence:: inner function fun34; subroutine pnorm; function alnorm.
!
! Reference::  Paul Gilbert and Ravi Varadhan (2012). numDeriv: Accurate Numerical
!              Derivatives. R package version 2012.9-1.
!--------------------------------------------------------------------------------------              
  implicit none
  integer(kind=4),intent(in)::npars                 
  integer(kind=4),intent(in)::nED 
  integer(kind=4),intent(in)::iter                   
  integer(kind=4),intent(out)::errorflag            
  real   (kind=8),intent(in)::pars(npars) 
  real   (kind=8),intent(in)::ED(nED)                
  real   (kind=8),intent(in)::Error(nED)                        
  real   (kind=8),intent(out)::grad(npars)       
  real   (kind=8),intent(out)::value    
  !
  ! Local variables
  real(kind=8),parameter::eps=1.0D-04
  real(kind=8),parameter:: d=0.0001D+00
  real(kind=8),parameter:: ZeroTol=1.781029D-05
  integer(kind=4),parameter:: v=2
  real(kind=8),dimension(npars)::vtol,h
  real(kind=8),dimension(iter,npars)::A
  real(kind=8),dimension(npars,npars)::Diag
  integer(kind=4)::i,j
  !
  ! Initializing some values
  A=0.0D+00
  vtol=0.0D+00
  Diag=0.0D+00
  errorflag=0
  grad=0.0D+00
  ! Calculate function's value
  ! and check for Infs and NaNs
  value=fun34(pars)
  if( value.ne.value .or. &
      value+1.0D+00==value )     then
    errorflag=1
    return
  end if
  ! That a diagonal matrix
  do i=1,npars
    Diag(i,i)=1.0D+00
  end do
  ! Special tackling with 
  ! small values in pars
  do i=1,npars
    if(dabs(pars(i))<zerotol)  vtol(i)=1.0D+00
  end do
  ! Set initial h value
  h=dabs(d*pars)+eps*vtol
  !
  ! Create matrix A with iter rows, npars columns.
  ! in each iteration, generating a new row by reduce
  ! h by 1/v.
  !.................................................
  do j=1,npars
    a(1,j)=(fun34(pars+h*Diag(j,:))-&
            fun34(pars-h*Diag(j,:)))/(2.0D+00*h(j))
  end do
  ! Check if any Inf or NaN presents in a(i,:)
  if( any( a(1,:) .ne. a(1,:) ) .or. &
      any( a(1,:)+1.0D+00==a(1,:) )  )  then
    errorflag=1
    return
  end if
  ! Successively reduce h by 1/v.
  h=h/real(v,kind=8)
  !.................................................
  do i=2,iter
    do j=1,npars
      if( i/=1 .and. dabs(a(i-1,j))<1.0D-20 )  then
        a(i,j)=0.0D+00
      else 
        a(i,j)=(fun34(pars+h*Diag(j,:))-&
                fun34(pars-h*Diag(j,:)))/(2.0D+00*h(j))
      end if
    end do
    ! Check if any Inf or NaN presents in a(i,:)
    if( any( a(i,:) .ne. a(i,:) ) .or. &
        any( a(i,:)+1.0D+00==a(i,:) )  )  then
      errorflag=1
      return
    end if
    ! Successively reduce h by 1/v.
    h=h/real(v,kind=8)
  end do
  !.................................................
  ! Apply Richardson Extrapolation 
  ! to improve the accuracy of gradient
  do i=1,iter-1
    a(1:(iter-i),:)=(a(2:(iter+1-i),:)*(4.0D+00**i)-&
                     a(1:(iter-i),:))/&
                    (4.0D+00**i-1.0D+00)
  end do
  !
  grad=a(1,:)
  !  
  return
  !
  contains
!----------------------
  function fun34(x)
  !------------------------------------------------------------------------------------------------------------
  ! fun34() is a inner function contained in subroutine hessian its used for calculating minus 
  ! logged maximum likelihood value of minimum age model of three or four parameters 
  !
  ! x(npars), input     :: real values, the parameters used for calculating
  ! fun34,   output     :: real value, the value for the specified function
  !
  !  Author :: Peng Jun, 2013.03.16
  !
  ! Reference:: Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of
  !             single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design
  !             and statistical models. Archaeometry, 41, pp. 339-364.
  !
  ! Dependence:: subroutine pnorm, function alnorm, also the ED data (log-scale) provided in subroutine hessian
  !
  !--------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8)::x(npars)
    real(kind=8)::fun34
    real(kind=8)::pnorm1(nED)
    real(kind=8)::alnorm
    real(kind=8),parameter::pi=3.141592653589793238462643383279502884197D+00
    logical,parameter::upper=.false.   
    !
    if (npars==3)    then
      !
      pnorm1=(x(2)-(x(2)/x(3)**2+ED/Error**2)/(1.0D+00/x(3)**2+1.0D+00/Error**2))*dsqrt(1.0D+00/Error**2+1.0D+00/x(3)**2)
      !
      call pnorm(pnorm1,nED,upper)
      !
      fun34=-sum(dlog(x(1)/dsqrt(2.0D+00*pi*Error**2)*dexp(-(ED-x(2))**2/(2.0D+00*Error**2))+&
	        (1.0D+00-x(1))/dsqrt(2.0D+00*pi*(Error**2+x(3)**2))*dexp(-(ED-x(2))**2/(2.0D+00*(Error**2+x(3)**2)))*&
	        (1.0D+00-pnorm1)/(0.5D+00)))  
    elseif(npars==4)  then
      !
      pnorm1=(x(2)-(x(3)/x(4)**2+ED/Error**2)/(1.0D+00/x(4)**2+1.0D+00/Error**2))*dsqrt(1.0D+00/Error**2+1.0D+00/x(4)**2)
      !
      call pnorm(pnorm1,nED,upper)	 
      !
      fun34=-sum(dlog(x(1)/dsqrt(2.0D+00*pi*Error**2)*dexp(-(ED-x(2))**2/(2.0D+00*Error**2))+&
	        (1.0D+00-x(1))/dsqrt(2.0D+00*pi*(Error**2+x(4)**2))*dexp(-(ED-x(3))**2/(2.0D+00*(Error**2+x(4)**2)))*&
	        (1.0D+00-pnorm1)/(1.0D+00-alnorm((x(2)-x(3))/x(4),upper))))   
      ! 
    end if
    !
    return
    !
  end function fun34 
  !
end subroutine gradient  
