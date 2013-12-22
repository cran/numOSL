subroutine lmfunc1(xdat,ydat,m,n,x, &
                  fvec,fjac,ldfjac,iflag)
! ----------------------------------------------------------------------------------------------------------------------------
! For formula I(t)=a1*(t/max(t))*exp(-b1*t^2)/2/max(t))+...+ak*(t/max(t))*exp(-bk*t^2)/2/max(t)), where k=1:7
! lmfunc1() is used to calculate vectors fevc(i)=I(t,i)-ydat(i), i=1:length(ydat), or the jacobian matrix fjac, 
! both them will be passed to subroutine lmfit to perform the Levenberg-Marquadt optimization.
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
! fjac(ldfjac,n),     output:: real values, the jacobian matrix
!
! iflag,               input:: integer, If IFLAG = 1 on intput, FCN should calculate the functions at X and return 
!                              this vector in FVEC. If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!                              return this matrix in FJAC. To terminate the algorithm, the user may set IFLAG negative.
! ===========================================================================================================================
! Author:: Peng Jun, 2013.07.24, revised in 2013.10.05.
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
  ! Local variables
  integer(kind=4)::i
  real   (kind=8)::maxx
  !
  ! Calculate values for targeted function and 
  ! store them in vector fvec
  maxx=maxval(xdat)
  if(iflag==1) then
    fvec=0.0D+00
    do i=1,n/2
      fvec=fvec+x(i)*(xdat/maxx)*&
           dexp(-x(i+n/2)*xdat**2/2.0D+00/maxx)
    end do
    fvec=fvec-ydat
  ! Calculate matrix jacobian in a column by column order
  else if(iflag==2)  then
    do i=1,n/2
      fjac(:,i)=(xdat/maxx)*dexp(-x(i+n/2)*xdat**2/2.0D+00/maxx)
    end do
    do i=n/2+1,n
      fjac(:,i)=-xdat**2/2.0D+00/maxx*&
                dexp(-x(i)*xdat**2/2.0D+00/maxx)*&
                x(i-n/2)*(xdat/maxx)
    end do
  end if
  ! now return
  return
end subroutine lmfunc1  
