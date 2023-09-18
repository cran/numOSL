subroutine mcFMM2(nED,nsim,ED,Error,addsigma,&
                  inis,iflog,ntry,w,m,chains,iflag)
!===========================================================================================
! Construct MCMC chains for the 2-component FMM age model. 
! (either in logged- or unlogged-scale).
!
! Author:: Peng Jun, 2023.08.30.
!
! Dependence:: subroutine SliceFMM2().
!===========================================================================================
  implicit none
  integer, intent(in):: nED, nsim, ntry                                
  real(kind(1.0d0)), intent(in):: ED(nED), Error(nED),&  
                                  addsigma, inis(2,2), w, m                  
  real(kind(1.0d0)), intent(out):: chains(nsim,4)      
  integer, intent(out):: iflag               
  logical, intent(in):: iflog                
  !
  ! Local variables.
  real(kind(1.0d0)):: yy(nED), xx(nED),  iniP1, iniP2,& 
                      iniMu1, iniMu2, Value, lowerMus,& 
                      upperMus, sumPs
  integer:: i, j, seed
  !
  seed=123456789
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
  !
  if (iflog .eqv. .true.)  then
      iniMu1=log(inis(2,1))
      iniMu2=log(inis(2,2))
  else 
      iniMu1=inis(2,1)
      iniMu2=inis(2,2)
  end if
  !
  ! Set lower and upper boundaries for Mus.
  if ( all(yy>0.0) )  then
      lowerMus=minval(yy)*0.999
      upperMus=maxval(yy)*1.001
  else if ( all(yy<=0.0) )  then
      lowerMus=minval(yy)*1.001
      upperMus=maxval(yy)*0.999
  else 
      lowerMus=minval(yy)*1.001
      upperMus=maxval(yy)*1.001
  end if
  !
  !!!call random_seed()
  do i=1, nsim
      !
      ! Update P1.
      loopP1: do j=1, ntry
          call SliceFMM2(iniP1,iniP2,iniMu1,iniMu2,nED,yy,xx,&
                         1,Value,iflag,w,m,0.0D+00,1.0D+00,seed)
          if (iflag==0) exit loopP1
      end do loopP1
      ! Error checking.
      if (iflag/=0)  return
      iniP1=Value
      !
      !
      ! Update P2.
      loopP2: do j=1, ntry
          call SliceFMM2(iniP1,iniP2,iniMu1,iniMu2,nED,yy,xx,&
                         2,Value,iflag,w,m,0.0D+00,1.0D+00,seed)
          if (iflag==0)  exit loopP2
      end do loopP2
      ! Error checking.
      if (iflag/=0)  return
      iniP2=Value
      !
      ! Normalize Ps.
      sumPs=iniP1+iniP2
      chains(i,1)=iniP1/sumPs
      chains(i,2)=iniP2/sumPs
      !
      !
      ! Update Mu1.
      loopMu1: do j=1, ntry
          call SliceFMM2(iniP1,iniP2,iniMu1,iniMu2,nED,yy,xx,&
                         3,Value,iflag,w,m,lowerMus,upperMus,seed)
          if (iflag==0)  exit loopMu1
      end do loopMu1
      ! Error checking.
      if (iflag/=0)  return
      iniMu1=Value
      chains(i,3)=Value
      ! 
      ! Update Mu2.
      loopMu2: do j=1, ntry
          call SliceFMM2(iniP1,iniP2,iniMu1,iniMu2,nED,yy,xx,&
                         4,Value,iflag,w,m,lowerMus,upperMus,seed)
          if (iflag==0)  exit loopMu2
      end do loopMu2
      ! Error checking.
      if (iflag/=0)  return
      iniMu2=Value
      chains(i,4)=Value
  end do
  !
  return
end subroutine mcFMM2
!
!-----------------------------------------------------------------------------
!
subroutine SliceFMM2(iniP1,iniP2,iniMu1,iniMu2,nED,ED,Error,&
                     which,Value,iflag,w,m,lower,upper,seed)
! =====================================================================================================
! Update parameters in 2-component FMM age model with the Slice Sampling.
!
! Author:: Peng Jun, 2023.09.09.
!
! Dependence:: function r8_uniform_01, inner functions funcPs(1,2) and funcMus(1,2).
!
! Reference :: Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
!              vol. 31, no. 3, pp. 705-767.
!
! Note:: THIS CODE IS TRANSFORMED (with minor modification) FROM THE PRIMITIVE R CODE OF RM NEAL AT:
!        http://www.cs.utoronto.ca/~radford/slice.software.html
! =====================================================================================================
  implicit none
  integer, intent(in):: nED, which
  real(kind(1.0d0)), intent(in):: iniP1, iniP2, iniMu1, iniMu2,& 
                                  ED(nED), Error(nED), w, m,&
                                  lower, upper          
  real(kind(1.0d0)), intent(out):: Value         
  integer, intent(out):: iflag 
  integer, intent(inout):: seed  
  !      
  ! Local variables.
  real(kind(1.0d0)):: ran, gx0, gx1, logy, U, L, R, J, K, gL, gR,& 
                      locVal, part1(nED), part2(nED), Error2(nED),&
                      r8_uniform_01
  !
  iflag=0
  Value=-99.0
  !
  Error2=Error**2
  !
  ! Calculate the specified function value.
  if (which==1)  then
      gx0=funcP1(iniP1)
  else if (which==2)  then
      gx0=funcP2(iniP2)
  else if (which==3) then
      gx0=funcMu1(iniMu1)
  else if (which==4) then
      gx0=funcMu2(iniMu2)
  end if
  !
  ! Error checking.
  if (iflag/=0) return
  !
  ! 
  ! Transform logy in log terms.
  !!!call random_number(ran)
  ran=r8_uniform_01(seed)
  logy=gx0+log(ran)
  !
  ! Find the initial interval [L,R] to sample from.
  !!!call random_number(ran)
  ran=r8_uniform_01(seed)
  U=ran * w
  !
  if (which==1)  then
      L=iniP1 - U
      R=iniP1 + (w-U)
  else if (which==2)  then
      L=iniP2 - U
      R=iniP2 + (w-U)
  else if (which==3)  then
      L=iniMu1 - U
      R=iniMu1 + (w-U)
  else if (which==4)  then
      L=iniMu2 - U
      R=iniMu2 + (w-U)
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
          if (which==1)  then
              gL=funcP1(L)
          else if (which==2)  then
              gL=funcP2(L)
          else if (which==3)  then
              gL=funcMu1(L)
          else if (which==4)  then
              gL=funcMu2(L)
          end if
          ! Error checking.
          if (iflag/=0)  return
          !
          if (gL<=logy)  exit
          L=L-w
      end do
      !
      ! For the right side.
      do 
          if(R>=upper)  exit
          !
          if (which==1)  then
              gR=funcP1(R)
          else if (which==2)  then
              gR=funcP2(R)
          else if (which==3)  then
              gR=funcMu1(R)
          else if (which==4)  then
              gR=funcMu2(R)
          end if
          ! Error checking.
          if (iflag/=0)  return
          if (gR<=logy)  exit
          R=R+w
      end do
      !
  else if (m>1.0)  then
      ! Limit on steps for the case that m>1.0.
      !!!call random_number(ran)
      ran=r8_uniform_01(seed)
      J=floor(m*ran)
      K=(m-1.0) - J
      !
      ! For the left side.
      do while (J>0.0) 
          if (L<=lower)  exit
          !
          if (which==1)  then
              gL=funcP1(L)
          else if (which==2)  then
              gL=funcP2(L)
          else if (which==3)  then
              gL=funcMu1(L)
          else if (which==4)  then
              gL=funcMu2(L)
          end if
          ! Error checking.
          if (iflag/=0)  return
          !
          if (gL<=logy)  exit
          L=L-w
          J=J-1.0
      end do
      !
      ! For the right side.
      do while (K>0.0) 
          if (R>=upper)  exit
          !
          if (which==1)  then
              gR=funcP1(R)
          else if (which==2)  then
              gR=funcP2(R)
          else if (which==3) then
              gR=funcMu1(R)
          else if (which==4)  then
              gR=funcMu2(R)
          end if
          ! Error checking.
          if (iflag/=0)  return
          !
          if (gR<=logy)  exit
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
      !!!call random_number(ran)
      ran=r8_uniform_01(seed)
      Value=L+ran*(R-L)
      !
      if (which==1)  then
          gx1=funcP1(Value)
      else if (which==2)  then
          gx1=funcP2(Value)
      else if (which==3)  then
          gx1=funcMu1(Value)
      else if (which==4)  then
          gx1=funcMu2(Value)
      end if
      ! Error checking.
      if (iflag/=0)  return
      !
      if (gx1>=logy)  exit
      !
      select case (which)
          case (1)
              locVal=iniP1
          case (2)
              locVal=iniP2
          case (3)
              locVal=iniMu1
          case (4)
              locVal=iniMu2
      end select
      !
      if (Value>locVal)  then
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
      real(kind(1.0d0)):: x, funcP1
      !
      ! Local variables
      real(kind(1.0d0)):: sumPs
      !
      sumPs=x+iniP2
      part1=    x/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      funcP1=sum( log(part1+part2) )
      ! Error checking.
      if (funcP1 .ne. funcP1) iflag=1
      return
  end function funcP1
  !----------------------------------------------------------------------
  function funcP2(x)
      implicit none
      real(kind(1.0d0)):: x, funcP2 
      !
      ! Local variables
      real(kind(1.0d0)):: sumPs
      !
      sumPs=iniP1+x
      part1=iniP1/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=    x/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      funcP2=sum( log(part1+part2) )
      ! Error checking.
      if (funcP2 .ne. funcP2) iflag=1
      return
  end function funcP2
  !----------------------------------------------------------------------
  function funcMu1(x)
      implicit none
      real(kind(1.0d0)):: x, funcMu1 
      !
      ! Local variables
      real(kind(1.0d0)):: sumPs
      !
      sumPs=iniP1+iniP2
      part1=iniP1/sumPs/Error*exp( -(ED-     x)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-iniMu2)**2/2.0/Error2 )
      funcMu1=sum( log(part1+part2) )
      ! Error checking.
      if (funcMu1 .ne. funcMu1) iflag=1
      return
  end function funcMu1
  !----------------------------------------------------------------------
  function funcMu2(x)
      implicit none
      real(kind(1.0d0)):: x, funcMu2 
      !
      ! Local variables
      real(kind(1.0d0)):: sumPs
      !
      sumPs=iniP1+iniP2
      part1=iniP1/sumPs/Error*exp( -(ED-iniMu1)**2/2.0/Error2 )
      part2=iniP2/sumPs/Error*exp( -(ED-     x)**2/2.0/Error2 )
      funcMu2=sum( log(part1+part2) )
      ! Error checking.
      if (funcMu2 .ne. funcMu2) iflag=1
      return
  end function funcMu2
  !----------------------------------------------------------------------
end subroutine SliceFMM2
