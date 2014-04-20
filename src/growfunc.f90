subroutine growfunc(xdat,ydat,m,n,x,fvec,fjac,ldfjac,iflag)
!----------------------------------------------------------------------------------------
! Subroutine growfunc is used for calculating the fvec or a jacobian matrix, those values 
! will be used in subroutine lmder1 for performing a Levenberg-Marquadt optimization.
! =======================================================================================
!
! m,         input:: integer, number of paired observations.
!
! n,         input:: integer, the dimension of the problem:
!                    3 for an exponential model y=k*(1-exp(-b*x))+c;
!                    4 for a linear+Exponential model y=k*(1-exp(-b*x))+c*x+d.
!
! ldfjac,    input:: integer, the leading dimension of FJAC, which needs to be consistent with m.
!
! iflag,     input:: integer, calculate the functions (iflag=1) or the jacobian matrix (iflag=2) at X.
!
! xdat(m),   input:: real values, regenerative dose values.
!
! ydat(m),   input:: real values, standardlized signal values.
!
! x(n),      input:: real values, initials from which either fvec or jacobian will be estimated.
!
! fvec(m),  output:: real values, differences between observations and fitted values.
!
! fjac(ldfjac,n),output:: real values, the estimated jacobian matrix.
! ===========================================================================================
! Author:: Peng Jun, 2013.05.27; revised in 2014.04.02.
!
! Dependence::No
!--------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                    intent(in):: m
  integer(kind=4),                    intent(in):: n
  integer(kind=4),                    intent(in):: ldfjac
  integer(kind=4),                    intent(in):: iflag
  real   (kind=8),dimension(m),       intent(in):: xdat,ydat
  real   (kind=8),dimension(n),       intent(in):: x
  real   (kind=8),dimension(m),       intent(out):: fvec
  real   (kind=8),dimension(ldfjac,n),intent(out):: fjac
  ! Local variables.
  real   (kind=8),dimension(4):: cx
  real   (kind=8),dimension(ldfjac,4):: cfjac
  ! 
  ! Store x in cx.
  cx=0.0D+00
  cx(1:n)=x
  cfjac=0.0D+00
  !
  if(iflag==1) then
    ! Calculate fvec.
    if(n==3) then
      fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)-ydat  ! an exponential model (do not pass the origin).
    else if(n==4) then
      fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)*xdat+cx(4)-ydat ! a linear+Exponential model (do not pass the origin).
    end if
  else if(iflag==2) then
    ! Calculate the jacobian matrix.
    if(n==3) then
      cfjac(:,1)=1.0D+00-dexp(-x(2)*xdat) 
      cfjac(:,2)=x(1)*xdat*dexp(-x(2)*xdat)  ! an exponential model (do not pass the origin).
      cfjac(:,3)=1.0D+00
    else if(n==4) then
      cfjac(:,1)=1.0D+00-dexp(-x(2)*xdat)
      cfjac(:,2)=x(1)*xdat*dexp(-x(2)*xdat)   ! a linear+Exponential model (do not pass the origin).
      cfjac(:,3)=xdat
      cfjac(:,4)=1.0D+00
    end if
    !
    fjac(:,1:n)=cfjac(:,1:n)
  end if
  return
end subroutine growfunc
