subroutine lmfunc(xdat,ydat,m,n,x,fvec,fjac,ldfjac,iflag)
! ---------------------------------------------------------------------------------------------------------------------------
! For formula I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t), K=1:7,
! lmfunc is used to calculate vectors fevc(i)=I(t,i)-ydat(i), i=1:length(ydat), 
! so is the jacobian matrix fjac, they will be passed to subroutine lmfit
! to perform the Levenberg-Marquadt optimization.
! ===========================================================================================================================
!
! m,                   input:: integer, length of xdat or ydat.
!
! n,                   input:: integer, dimension of the problem (length of x).
!
! ldfjac,              input:: integer,the leading dimension of FJAC, which must be no less than m.
!
! xdat(m),             input:: real values, times values.
!
! ydat(m),             input:: real values, signal values.
!
! x(n),                input:: real values, initial guess values.
!
! fvec(m),            output:: real values, differences between predicted values and obeserved values.
!
! fjac(ldfjac,n),     output:: real values, the jacobian matrix.
!
! iflag,               input:: integer, If IFLAG = 1 on intput, FCN should calculate the functions at X and return 
!                              this vector in FVEC. If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!                              return this matrix in FJAC. To terminate the algorithm, the user may set IFLAG negative.
! ===========================================================================================================================
! Author:: Peng Jun, 2013.05.21.
!
! Dependence:: No
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::m
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::ldfjac
  integer(kind=4),intent(in)::iflag
  real   (kind=8),dimension(m),intent(in)::xdat,ydat
  real   (kind=8),dimension(n),intent(in)::x
  real   (kind=8),dimension(m),intent(out)::fvec
  real   (kind=8),dimension(ldfjac,n),intent(out)::fjac
  ! local variables
  integer(kind=4)::i
  !
  ! calculate values for targeted function and 
  ! store them in vector fvec
  if(iflag==1) then
    fvec=0.0D+00
    do i=1,n/2
      fvec=fvec+x(i)*dexp(-x(i+n/2)*xdat)
    end do
    fvec=fvec-ydat
  ! calculate matrix jacobian in a column by column order
  else if(iflag==2)  then
    do i=1,n/2
      fjac(:,i)=dexp(-x(i+n/2)*xdat)
    end do
    do i=n/2+1,n
      fjac(:,i)=-x(i-n/2)*xdat*dexp(-x(i)*xdat)
    end do
  end if
  !
  return
end subroutine lmfunc  
