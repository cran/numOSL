subroutine mcCAM(nED,nsim,ED,Error,addsigma,&
                 inis,iflog,ntry,w,m,chains,iflag)
!===========================================================================================
! Construct MCMC chains for the CAM age model (either in logged- or unlogged-scale).
!
! Author:: Peng Jun, 2014.03.15.
!
! Dependence:: subroutine SliceCAM().
!===========================================================================================
  implicit none
  integer(kind=4), intent(in):: nED                  ! Number of equivalent doses.
  integer(kind=4), intent(in):: nsim                 ! Number of simulations.
  real   (kind=8), intent(in):: ED(nED)              ! Equivalent dose values (un-logged).
  real   (kind=8), intent(in):: Error(nED)           ! Standard errors (absolute).
  real   (kind=8), intent(in):: addsigma             ! Added spread.
  real   (kind=8), intent(in):: inis(2)              ! Initials of the chains.
  integer(kind=4), intent(in):: ntry                 ! Maximum times of trials per simulation.
  real   (kind=8), intent(in):: w                    ! Size of the steps for creating interval.
  real   (kind=8), intent(in):: m                    ! Limit on steps.
  real   (kind=8), intent(out):: chains(nsim,2)      ! Simulated chains.
  integer(kind=4), intent(out):: iflag               ! Error message (0=success, 1=fail).
  logical,         intent(in):: iflog                ! Change to logged-scale or not.
  !
  ! Local variables.
  real   (kind=8):: yy(nED), xx(nED)
  real   (kind=8):: iniMu, iniSigma
  real   (kind=8):: Value
  integer(kind=4):: i, j
  real   (kind=8):: upperSigma
  real   (kind=8):: lowerMu, upperMu
  !
  ! Default return chains.
  chains=-99.0
  ! 
  ! Use logged-scale or not.
  if (iflog .eqv. .true.)  then
      xx=sqrt( (Error/ED)**2 +addsigma**2 )
      yy=log(ED)
      upperSigma=5.0
  else 
      xx=sqrt( Error**2 + addsigma**2 )
      yy=ED
      upperSigma=sum( (yy-sum(yy)/nED)**2 )/(nED-1)
  end if
  !
  ! Initials of the chains.
  if (iflog .eqv. .true.)  then
      iniMu=log( inis(1) )
  else 
      iniMu=inis(1)
  end if 
  iniSigma=inis(2)
  !
  ! Set lower and upper boundaries for Mu.
  if ( all(yy>0.0) ) then
      lowerMu=minval(yy)*0.999
      upperMu=maxval(yy)*1.001
  else if ( all(yy<=0.0) ) then
      lowerMu=minval(yy)*1.001
      upperMu=maxval(yy)*0.999
  else 
      lowerMu=minval(yy)*1.001
      upperMu=maxval(yy)*1.001
  end if
  !
  call random_seed()
  do i=1, nsim
      ! Update Mu.
      AA: do j=1, ntry
          call SliceCAM(iniMu,iniSigma,nED,yy,xx,& 
                        1,Value,iflag,w,m,lowerMu,upperMu)
          if (iflag==0) exit AA
      end do AA
      ! Error checking.
      if (iflag/=0)  return
      iniMu=Value
      chains(i,1)=Value
      !
      !
      ! Update Sigma.
      BB: do j=1, ntry
          call SliceCAM(iniMu,iniSigma,nED,yy,xx,& 
                        2,Value,iflag,w,m,0.0D+00,upperSigma)
          if (iflag==0)  exit BB
      end do BB
      ! Error checking.
      if (iflag/=0)  return
      iniSigma=Value
      chains(i,2)=Value
  end do
  !
  return
end subroutine mcCAM
!
! ------------------------------------------------------------------------------------------------------
!
subroutine SliceCAM(iniMu,iniSigma,nED,ED,Error,& 
                    which,Value,iflag,w,m,lower,upper)
! =====================================================================================================
! Update parameters in CAM age model with the Slice Sampling.
!
! Author:: Peng Jun, 2014.03.15.
!
! Dependence:: Inner function funcMu(); funcSigma().
!
! Reference :: Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
!              vol. 31, no. 3, pp. 705-767.
!
! Note:: THIS CODE IS TRANSFORMED (with minor modification) FROM THE PRIMITIVE R CODE OF RM NEAL AT:
!        http://www.cs.utoronto.ca/~radford/slice.software.html
! =====================================================================================================
  implicit none
  real   (kind=8), intent(in):: iniMu             ! Initial Mu value.
  real   (kind=8), intent(in):: iniSigma          ! Initial Sigma value.
  integer(kind=4), intent(in):: nED               ! Number of equivalent doses.
  integer(kind=4), intent(in):: which             ! Which parameter to be updated (1=Mu, 2=Sigma).
  real   (kind=8), intent(in):: ED(nED)           ! Equivalent dose values (either logged or un-logged).
  real   (kind=8), intent(in):: Error(nED)        ! Standard errors (either relative or absoulte).
  real   (kind=8), intent(in):: w                 ! Size of the steps for creating interval.
  real   (kind=8), intent(in):: m                 ! Limit on steps for creating interval.
  real   (kind=8), intent(in):: lower             ! Lower boundary of the desired distribution.
  real   (kind=8), intent(in):: upper             ! Upper boundary of the desired distribution.
  real   (kind=8), intent(out):: Value            ! Updated parameter corresponding to "which".
  integer(kind=4), intent(out):: iflag            ! Error message (0=success, 1=crashed).
  ! Local variables.
  real   (kind=8):: ran
  real   (kind=8):: gx0, gx1
  real   (kind=8):: logy
  real   (kind=8):: U, L, R, J, K
  real   (kind=8):: gL, gR
  real   (kind=8):: Error2(nED)
  !
  iflag=0
  Value=-99.0
  !
  Error2=Error**2
  !
  ! Calculate the specified function value.
  if (which==1) then
      gx0=funcMu(iniMu)
  else if (which==2) then
      gx0=funcSigma(iniSigma)
  end if
  ! Error checking.
  if(iflag/=0) return
  !
  !
  ! Transform logy in log terms.
  call random_number(ran)
  logy=gx0+log(ran)
  !
  ! Find the initial interval [L,R] to sample from.
  call random_number(ran)
  U=ran * w
  !
  if (which==1)  then
      L=iniMu - U
      R=iniMu + (w-U)
  else if (which==2) then
      L=iniSigma - U
      R=iniSigma + (w-U)
  end if
  !
  ! Expand the interval until its ends are outside the
  ! slice, or until the limit on steps is reached.
  if (m<=1.0)  then
      ! No limit on number of steps for the case m<=1.0.
      !
      ! For the left side.
      do 
          if (L<=lower) exit
          ! 
          if (which==1)  then
              gL=funcMu(L)
          else if (which==2)  then
              gL=funcSigma(L)
          end if
          ! Error checking.
          if (iflag/=0) return
          !
          if (gL<=logy) exit
          L=L-w
      end do
      ! For the right side.
      do 
          if (R>=upper) exit
          !
          if (which==1)  then
              gR=funcMu(R)
          else if (which==2) then
              gR=funcSigma(R)
          end if
          ! Error checking.
          if (iflag/=0) return
          !
          if (gR<=logy) exit
          R=R+w
      end do
      !
  else if (m>1.0) then
      ! Limit on steps for the case that m>1.0.
      call random_number(ran)
      J=floor(m*ran)
      K=(m-1.0) - J
      !
      ! For the left side.
      do while (J>0.0) 
          if (L<=lower)  exit
          ! 
          if (which==1)  then
              gL=funcMu(L)
          else if (which==2) then
              gL=funcSigma(L)
          end if 
          ! Error checking.
          if (iflag/=0) return
          !
          if (gL<=logy) exit
          L=L-w
          J=J-1.0
      end do
      !
      ! For the right side.
      do while (K>0.0) 
          if (R>=upper) exit
          !
          if (which==1)  then
              gR=funcMu(R)
          else if (which==2)  then
              gR=funcSigma(R)
          end if
          ! Error checking.
          if (iflag/=0) return
          !
          if (gR<=logy) exit
          R=R+w
          K=K-1.0
      end do
      !
  end if
  !
  ! Shrink the interval.
  if (L<lower) L=lower
  if (R>upper) R=upper
  !
  ! Sample from the interval (with shrinking).
  do 
      call random_number(ran)
      Value=L+ran*(R-L)
      !
      if (which==1) then
          gx1=funcMu(Value)
      else if (which==2)  then
          gx1=funcSigma(Value)
      end if
      ! Error checking.
      if (iflag/=0)  return
      !
      if (gx1>=logy) exit
      ! 
      if (which==1)  then
          if (Value>iniMu)  then
              R=Value
          else 
              L=Value
          end if
      else if (which==2)  then
          if (Value>iniSigma)  then
              R=Value
          else 
              L=Value
          end if
      end if
      !
  end do
  !
  return
  !
  contains 
  ! Some inner functions concerning Mu and Sigma.
  ! -----------------------------------------------------------------------------------------------
  function funcMu(x)
    implicit none
    real   (kind=8):: x
    real   (kind=8):: funcMu
    ! 
    funcMu=sum(log( 1.0/sqrt(iniSigma**2+Error2)*exp(-(ED-x)**2/2.0/(iniSigma**2+Error2)) ))
    ! Checking NAN (NA).
    if (funcMu .ne. funcMu) iflag=1
    !
    return
    !
  end function funcMu
  ! -----------------------------------------------------------------------------------------------
  function funcSigma(x)
    implicit none
    real   (kind=8):: x
    real   (kind=8):: funcSigma
    ! 
    funcSigma=sum(log( 1.0/sqrt(x**2+Error2)*exp(-(ED-iniMu)**2/2.0/(x**2+Error2)) ))
    ! Checking NAN (NA).
    if (funcSigma .ne. funcSigma) iflag=1
    !
    return
    !
  end function funcSigma
  ! -----------------------------------------------------------------------------------------------
end subroutine SliceCAM
