subroutine growfunc(xdat,ydat,m,n,x,fvec,fjac,ldfjac,iflag)
!----------------------------------------------------------------------------------------------------------------------------
! Subroutine for calculating fvec or jacobian matrix, those values
! will be used in subroutine lmder1 for Levenberg-Marquadt optimization.
! ===========================================================================================================================
!
! m,               input:: integer, number of paired observations.
!
! n,               input:: integer, the dimension of the problem:
!                          2 for linear model; 
!                          3 for Exponential model; 
!                          4 for linear+Exponential model.
!
! ldfjac,          input:: integer, the leading dimension of FJAC, here its need to be consistent with m.
!
! iflag,           input:: integer, calculate the functions at X (iflag=1) or calculate the jacobian matrix at X (iflag=2).
!
! xdat(m),         input:: real values, independent variables in paired fitting, here is the Dose values.
!
! ydat(m),         input:: real values, dependent variables in paried fitting, here is the standardlized signal values.
!
! x(n),            input:: real values, initial guess parameters, from which either fvec or jacobian will be estimated.
!
! fvec(n),        output:: real values, differences between the dependent variables and the fitted values, in vector format.
!
! fjac(ldfjac,n), output:: real values, the estimated jacobian matrix.
! ============================================================================================================================
! Author:: Peng Jun, 2013.05.27.
!
! Dependence::No
!-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                    intent(in)::m
  integer(kind=4),                    intent(in)::n
  integer(kind=4),                    intent(in)::ldfjac
  integer(kind=4),                    intent(in)::iflag
  real   (kind=8),dimension(m),       intent(in)::xdat,ydat
  real   (kind=8),dimension(n),       intent(in)::x
  real   (kind=8),dimension(m),       intent(out)::fvec
  real   (kind=8),dimension(ldfjac,n),intent(out)::fjac
  !
  integer(kind=4)::allerror
  real   (kind=8),dimension(4)::cx
  real   (kind=8),allocatable::cfjac(:,:)
  ! 
  ! Store x in cx
  cx=0.0D+00
  cx(1:n)=x
  !
  if(iflag==1) then
    ! Calculate fvec
    if(n==2) then                
      fvec=cx(1)*xdat+cx(2)-ydat      ! Linear model
    else if(n==3) then
      fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)-ydat  ! Exponential model
    else if(n==4) then
      fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)*xdat+cx(4)-ydat ! linear+Exponential model
    end if
  else if(iflag==2) then
    ! Calculate jacobian matrix
    if(n==2) then
      allocate( cfjac(ldfjac,2),stat=allerror)
      cfjac(:,1)=xdat       ! Linear model
      cfjac(:,2)=1.0D+00
    else if(n==3) then
      allocate( cfjac(ldfjac,3),stat=allerror)
      cfjac(:,1)=1.0D+00-dexp(-x(2)*xdat) 
      cfjac(:,2)=x(1)*xdat*dexp(-x(2)*xdat)  ! Exponential model
      cfjac(:,3)=1.0D+00
    else if(n==4) then
      allocate( cfjac(ldfjac,4),stat=allerror)
      cfjac(:,1)=1.0D+00-dexp(-x(2)*xdat)
      cfjac(:,2)=x(1)*xdat*dexp(-x(2)*xdat)   ! linear+Exponential model
      cfjac(:,3)=xdat
      cfjac(:,4)=1.0D+00
    end if
    !
    fjac=cfjac
    deallocate(cfjac)
  end if
  return
end subroutine growfunc
