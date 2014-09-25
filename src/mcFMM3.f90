subroutine mcFMM3(nED,nsim,ED,Error,addsigma,&
                  inis,iflog,ntry,w,m,chains,iflag)
!===========================================================================================
! Construct MCMC chains for the 3-component FMM age model 
! (either in logged- or unlogged-scale).
!
! Author:: Peng Jun, 2014.03.16.
!
! Dependence:: subroutine SliceFMM3().
!===========================================================================================
  implicit none
  integer(kind=4), intent(in):: nED               ! Number of equivalent doses.
  integer(kind=4), intent(in):: nsim              ! Number of simulations.
  real   (kind=8), intent(in):: ED(nED)           ! Equivalent dose values (un-logged).
  real   (kind=8), intent(in):: Error(nED)        ! Standard errors (absolute).
  real   (kind=8), intent(in):: addsigma          ! Added spread.
  real   (kind=8), intent(in):: inis(2,3)         ! Initials of the chains.
  integer(kind=4), intent(in):: ntry              ! Maximum times of trials per simulation.
  real   (kind=8), intent(in):: w                 ! Size of the steps for creating interval.
  real   (kind=8), intent(in):: m                 ! Limit on steps.
  real   (kind=8), intent(out):: chains(nsim,6)   ! Simulated chains.
  integer(kind=4), intent(out):: iflag            ! Error message (0=success, 1=fail).
  logical,         intent(in):: iflog             ! Change to logged-scale or not.
  ! 
  ! Local variables.
  real   (kind=8):: yy(nED), xx(nED)
  real   (kind=8):: iniP1, iniP2, iniP3
  real   (kind=8):: iniMu1, iniMu2, iniMu3
  real   (kind=8):: Value
  integer(kind=4):: i, j
  real   (kind=8):: lowerMus, upperMus
  real   (kind=8):: sumPs
  !
  ! Default return chains.
  chains=-99.0
  !
  ! Use logged-scale or not.
  if (iflog .eqv. .true.) then
      xx=sqrt( (Error/ED)**2 + addsigma**2 )
      yy=log(ED) 
  else 
      xx=sqrt( Error**2 + addsigma**2 )
      yy=ED
  end if
  !
  ! Initials of the chains.
  iniP1=inis(1,1)/sum(inis(1,:))
  iniP2=inis(1,2)/sum(inis(1,:))
  iniP3=inis(1,3)/sum(inis(1,:)) 
  !
  if (iflog .eqv. .true.)  then
      iniMu1=log(inis(2,1))
      iniMu2=log(inis(2,2))
      iniMu3=log(inis(2,3))
  else 
      iniMu1=inis(2,1)
      iniMu2=inis(2,2)
      iniMu3=inis(2,3)
  end if
  !
  ! Set lower and upper boundaries for Mu.
  if ( all(yy>0.0) ) then
      lowerMus=minval(yy)*0.999
      upperMus=maxval(yy)*1.001
  else if ( all(yy<=0.0) ) then
      lowerMus=minval(yy)*1.001
      upperMus=maxval(yy)*0.999
  else 
      lowerMus=minval(yy)*1.001
      upperMus=maxval(yy)*1.001
  end if
  !
  call random_seed()
  do i=1, nsim
      !
      ! Update P1.
      loopP1: do j=1, ntry
          call SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                         nED,yy,xx,1,Value,iflag,w,m,0.0D+00,1.0D+00)
          if (iflag==0) exit loopP1
      end do loopP1
      ! Error checking.
      if (iflag/=0)  return
      iniP1=Value
      !
      !
      ! Update P2.
      loopP2: do j=1, ntry
          call SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                         nED,yy,xx,2,Value,iflag,w,m,0.0D+00,1.0D+00)
          if (iflag==0) exit loopP2
      end do loopP2
      ! Error checking.
      if (iflag/=0)  return
      iniP2=Value
      !
      !
      ! Update P3.
      loopP3: do j=1, ntry
          call SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                         nED,yy,xx,3,Value,iflag,w,m,0.0D+00,1.0D+00)
          if (iflag==0) exit loopP3
      end do loopP3
      ! Error checking.
      if (iflag/=0)  return
      iniP3=Value
      !
      ! Normalize Ps.
      sumPs=iniP1+iniP2+iniP3
      chains(i,1)=iniP1/sumPs
      chains(i,2)=iniP2/sumPs
      chains(i,3)=iniP3/sumPs
      !
      !
      ! Update Mu1.
      loopMu1: do j=1, ntry
          call SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                         nED,yy,xx,4,Value,iflag,w,m,lowerMus,upperMus)
          if (iflag==0) exit loopMu1
      end do loopMu1
      ! Error checking.
      if (iflag/=0) return
      iniMu1=Value
      chains(i,4)=Value
      !
      ! Update Mu2.
      loopMu2: do j=1, ntry
          call SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                         nED,yy,xx,5,Value,iflag,w,m,lowerMus,upperMus)
          if (iflag==0) exit loopMu2
      end do loopMu2
      ! Error checking.
      if (iflag/=0) return
      iniMu2=Value
      chains(i,5)=Value
      !
      ! Update Mu3.
      loopMu3: do j=1, ntry
          call SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                         nED,yy,xx,6,Value,iflag,w,m,lowerMus,upperMus)
          if (iflag==0) exit loopMu3
      end do loopMu3
      ! Error checking.
      if (iflag/=0) return
      iniMu3=Value
      chains(i,6)=Value
  end do
  !
  return
end subroutine mcFMM3
!
! -------------------------------------------------------------------------------------------------------
!
subroutine SliceFMM3(iniP1,iniP2,iniP3,iniMu1,iniMu2,iniMu3,&
                     nED,ED,Error,which,Value,iflag,w,m,lower,upper)
! =====================================================================================================
! Update parameters in 3-component FMM age model with the Slice Sampling.
!
! Author:: Peng Jun, 2014.03.16.
!
! Dependence:: Inner function funcPs(1,2,3) and funcMus(1,2,3).
!
! Reference :: Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
!              vol. 31, no. 3, pp. 705-767.
!
! Note:: THIS CODE IS TRANSFORMED (with minor modification) FROM THE PRIMITIVE R CODE OF RM NEAL AT:
!        http://www.cs.utoronto.ca/~radford/slice.software.html
! =====================================================================================================
  implicit none
  real   (kind=8), intent(in):: iniP1           ! Initial P1 value.
  real   (kind=8), intent(in):: iniP2           ! Initial P2 value.
  real   (kind=8), intent(in):: iniP3           ! Initial P3 value.
  real   (kind=8), intent(in):: iniMu1          ! Initial Mu1 value.
  real   (kind=8), intent(in):: iniMu2          ! Initial Mu2 value.
  real   (kind=8), intent(in):: iniMu3          ! Initial Mu3 value.
  integer(kind=4), intent(in):: nED             ! Number of equivalent doses.
  integer(kind=4), intent(in):: which           ! Which parameter to be updated (1,2,3 for Ps; 4,5,6 for Mus).
  real   (kind=8), intent(in):: ED(nED)         ! Equivalent dose values (either logged or un-logged).
  real   (kind=8), intent(in):: Error(nED)      ! Standard errors (either relative or absolute).
  real   (kind=8), intent(in):: w               ! Size of the steps for creating interval.
  real   (kind=8), intent(in):: m               ! Limit on steps for creating interval.
  real   (kind=8), intent(in):: lower           ! Lower boundary of the desired parameter.
  real   (kind=8), intent(in):: upper           ! Upper boundary of the desired parameter.
  real   (kind=8), intent(out):: Value          ! Updated parameter corresponding to "which".
  integer(kind=4), intent(out):: iflag          ! Error message (0=success, 1=fail).
  ! Local variables.
  real   (kind=8):: ran
  real   (kind=8):: gx0, gx1
  real   (kind=8):: logy
  real   (kind=8):: U, L, R, J, K
  real   (kind=8):: gL, gR
  real   (kind=8):: locVal
  real   (kind=8):: part1(nED), part2(nED), part3(nED)
  real   (kind=8):: Error2(nED)
  !
  iflag=0
  Value=-99.0
  !
  Error2=Error**2
  !
  ! Calculate the specified function value.
  if (which==1)  then
      gx0=funcP1(iniP1)
  else if (which==2) then
      gx0=funcP2(iniP2)
  else if (which==3) then
      gx0=funcP3(iniP3)
  else if (which==4) then
      gx0=funcMu1(iniMu1)
  else if (which==5) then
      gx0=funcMu2(iniMu2)
  else if (which==6) then
      gx0=funcMu3(iniMu3)
  end if
  ! Error checking.
  if (iflag/=0) return
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
      L=iniP1 - U
      R=iniP1 + (w-U)
  else if (which==2) then
      L=iniP2 - U
      R=iniP2 + (w-U)
  else if (which==3)  then
      L=iniP3 - U
      R=iniP3 + (w-U)
  else if (which==4) then
      L=iniMu1 - U
      R=iniMu1 + (w-U)
  else if (which==5) then
      L=iniMu2 - U
      R=iniMu2 + (w-U)
  else if (which==6) then
      L=iniMu3 - U
      R=iniMu3 + (w-U)
  end if
  !
  ! Expand the interval until its ends are outside the
  ! slice, or until the limit on steps is reached.
  if (m<=1.0)  then
      ! No limit on number of steps for the case m<=1.0.
      !
      ! For the left side.
      do 
          if (L<=lower)  exit
          !
          if (which==1) then
              gL=funcP1(L)
          else if (which==2) then
              gL=funcP2(L)
          else if (which==3) then
              gL=funcP3(L)
          else if (which==4) then
              gL=funcMu1(L)
          else if (which==5) then
              gL=funcMu2(L)
          else if (which==6) then
              gL=funcMu3(L)
          end if
          ! Error checking.
          if (iflag/=0) return
          !
          if (gL<=logy) exit
          L=L-w
      end do
      ! 
      ! For the right side.
      do
          if (R>=upper)  exit
          !
          if (which==1) then
              gR=funcP1(R)
          else if (which==2) then
              gR=funcP2(R)
          else if (which==3) then
              gR=funcP3(R)
          else if (which==4) then
              gR=funcMu1(R)
          else if (which==5)  then
              gR=funcMu2(R)
          else if (which==6)  then
              gR=funcMu3(R) 
          end if
          ! Error checking.
          if (iflag/=0)  return
          if (gR<=logy)  exit
          R=R+w
      end do
      !
  else if (m>1.0)  then
      ! Limit on steps for the case that m>1.0.
       call random_number(ran)
       J=floor(m*ran)
       K=(m-1.0) - J
       !
       ! For the left side.
       do while (J>0.0) 
           if (L<=lower)  exit
           !
           if (which==1) then
               gL=funcP1(L)
           else if (which==2) then
               gL=funcP2(L)
           else if (which==3) then
               gL=funcP3(L)
           else if (which==4) then
               gL=funcMu1(L)
           else if (which==5) then
               gL=funcMu2(L)
           else if (which==6) then
               gL=funcMu3(L)
           end if
           ! Error checking.
           if (iflag/=0) return
           !
           if (gL<=logy)  exit
           L=L-w
           J=J-1.0
       end do
       !
       ! For the right side.
       do while (K>0.0) 
           if (R>=upper) exit
           !
           if (which==1) then
               gR=funcP1(R)
           else if (which==2) then
               gR=funcP2(R)
           else if (which==3) then
               gR=funcP3(R)
           else if (which==4)  then
               gR=funcMu1(R)
           else if (which==5)  then
               gR=funcMu2(R)
           else if (which==6) then
               gR=funcMu3(R)
           end if
           ! Error checking.
           if (iflag/=0)  return
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
  ! Sample from the interval with shrinking.
  do 
      call random_number(ran)
      Value=L+ran*(R-L)
      !
      if (which==1) then
          gx1=funcP1(Value) 
      else if (which==2) then
          gx1=funcP2(Value)
      else if (which==3) then
          gx1=funcP3(Value)
      else if (which==4)  then
          gx1=funcMu1(Value) 
      else if (which==5) then
          gx1=funcMu2(Value)
      else if (which==6) then
          gx1=funcMu3(Value)
      end if
      ! Error checking.
      if (iflag/=0) return
      !
      if (gx1>=logy) exit
      !
      select case (which)
          case (1)
              locVal=iniP1
          case (2)
              locVal=iniP2
          case (3)
              locVal=iniP3
          case (4)
              locVal=iniMu1
          case (5)
              locVal=iniMu2
          case (6)
              locVal=iniMu3
      end select
      !
      if (Value>locVal) then
          R=Value
      else 
          L=Value
      end if
      !
  end do
  !
  return
  !
  contains
  ! Inner functions about Ps and Mus.
  !----------------------------------------------------------------------
  function funcP1(x)
      implicit none
      real   (kind=8):: x
      real   (kind=8):: funcP1
      ! Local variables
      real   (kind=8):: sumPs
      !
      sumPs=x+iniP2+iniP3
      part1=    x/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      part3=iniP3/sumPs/Error*exp( -(ED-iniMu3)**2/2.0/Error2 )
      funcP1=sum( log(part1+part2+part3) )
      ! Error checking.
      if (funcP1 .ne. funcP1) iflag=1
      return
  end function funcP1
  !----------------------------------------------------------------------
  function funcP2(x)
      implicit none
      real   (kind=8):: x
      real   (kind=8):: funcP2
      ! Local variables
      real   (kind=8):: sumPs
      !
      sumPs=iniP1+x+iniP3
      part1=iniP1/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=    x/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      part3=iniP3/sumPs/Error*exp( -(ED-iniMu3)**2/2.0/Error2 )
      funcP2=sum( log(part1+part2+part3) )
      ! Error checking.
      if (funcP2 .ne. funcP2) iflag=1
      return
  end function funcP2
  !----------------------------------------------------------------------
  function funcP3(x)
      implicit none
      real   (kind=8):: x
      real   (kind=8):: funcP3
      ! Local variables
      real   (kind=8):: sumPs
      !
      sumPs=iniP1+iniP2+x
      part1=iniP1/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      part3=    x/sumPs/Error*exp( -(ED-iniMu3)**2/2.0/Error2 )
      funcP3=sum( log(part1+part2+part3) )
      ! Error checking.
      if (funcP3 .ne. funcP3) iflag=1
      return
  end function funcP3
  !----------------------------------------------------------------------
  function funcMu1(x)
      implicit none
      real   (kind=8):: x
      real   (kind=8):: funcMu1
      ! Local variables
      real   (kind=8):: sumPs
      !
      sumPs=iniP1+iniP2+iniP3
      part1=iniP1/sumPs/Error*exp( -(ED-     x)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      part3=iniP3/sumPs/Error*exp( -(ED-iniMu3)**2/2.0/Error2 )
      funcMu1=sum( log(part1+part2+part3) )
      ! Error checking.
      if (funcMu1 .ne. funcMu1) iflag=1
      return
  end function funcMu1
  !----------------------------------------------------------------------
  function funcMu2(x)
      implicit none
      real   (kind=8):: x
      real   (kind=8):: funcMu2
      ! Local variables
      real   (kind=8):: sumPs
      !
      sumPs=iniP1+iniP2+iniP3
      part1=iniP1/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-     x)**2/2.0/Error2 )
      part3=iniP3/sumPs/Error*exp( -(ED-iniMu3)**2/2.0/Error2 )
      funcMu2=sum( log(part1+part2+part3) )
      ! Error checking.
      if (funcMu2 .ne. funcMu2) iflag=1
      return
  end function funcMu2
  !----------------------------------------------------------------------
  function funcMu3(x)
      implicit none
      real   (kind=8):: x
      real   (kind=8):: funcMu3
      ! Local variables
      real   (kind=8):: sumPs
      !
      sumPs=iniP1+iniP2+iniP3
      part1=iniP1/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      part3=iniP3/sumPs/Error*exp( -(ED-     x)**2/2.0/Error2 )
      funcMu3=sum( log(part1+part2+part3) )
      ! Error checking.
      if (funcMu3 .ne. funcMu3) iflag=1
      return
  end function funcMu3
  !----------------------------------------------------------------------
end subroutine SliceFMM3
