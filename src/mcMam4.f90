subroutine mcMam4(nED,nsim,ED,Error,addsigma,&
                  inis,iflog,ntry,w,m,chains,iflag)
!===========================================================================================
! Construct MCMC chains for the MAM-4 age model (either in logged- or unlogged-scale).
!
! Author:: Peng Jun, 2014.03.13.
!
! Dependence:: subroutine SliceMam4().
!===========================================================================================
  implicit none
  integer(kind=4),intent(in):: nED              ! Number of equivalent doses.
  integer(kind=4),intent(in):: nsim             ! Number of simulations.
  real   (kind=8),intent(in):: ED(nED)          ! Equivalent dose values (un-logged).
  real   (kind=8),intent(in):: Error(nED)       ! Standard errors (absolute).
  real   (kind=8),intent(in):: addsigma         ! Added spread.
  real   (kind=8),intent(in):: inis(4)          ! Initials of the chains.
  integer(kind=4),intent(in):: ntry             ! Maximum times of trials in parameter updating in on simulation.
  real   (kind=8),intent(in):: w                ! Size of the steps for creating interval.
  real   (kind=8),intent(in):: m                ! Limit on steps.
  real   (kind=8),intent(out):: chains(nsim,4)  ! Simulated chains.
  integer(kind=4),intent(out):: iflag           ! Error message (0=success, 1=fail).
  logical,         intent(in):: iflog           ! Change to logged0scale or not.
  !
  ! Local variables.
  real   (kind=8):: yy(nED), xx(nED)
  real   (kind=8):: iniP, iniGama, iniMu, iniSigma
  real   (kind=8):: Value
  integer(kind=4):: i, j
  real   (kind=8):: upperSigma
  real   (kind=8):: lowerGama, upperGama
  !
  ! Default return chains.
  chains=-99.0
  !
  ! Use logged-scale or not.
  if (iflog .eqv. .true.)  then
      xx=sqrt( (Error/ED)**2 + addsigma**2 )
      yy=log(ED)
      upperSigma=5.0
  else 
      xx=sqrt( Error**2 + addsigma**2 )
      yy=ED
      upperSigma=sum((yy- sum(yy)/nED)**2)/(nED-1)
  end if
  !
  ! Initials of the chains.
  iniP=inis(1)
  if (iflog .eqv. .true.) then
      iniGama=log( inis(2) )
      iniMu=log( inis(3) )
  else 
      iniGama=inis(2)
      iniMu=inis(3)
  end if
  iniSigma=inis(4)
  !
  ! Set lower and upper boundaries for gama and mu.
  if ( all(yy>0.0) ) then
      lowerGama=minval(yy)*0.999
      upperGama=maxval(yy)*1.001
  else if ( all(yy<=0.0) ) then
      lowerGama=minval(yy)*1.001
      upperGama=maxval(yy)*0.999
  else 
      lowerGama=minval(yy)*1.001
      upperGama=maxval(yy)*1.001
  end if
  !
  call random_seed()
  do i=1, nsim
      ! Update P.
      AA: do j=1, ntry
          call SliceMam4(iniP,iniGama,iniMu,iniSigma,nED,yy,xx,&
                         1,Value,iflag,w,m,0.0D+00,1.0D+00)
          if (iflag==0) exit AA
      end do AA
      if (iflag/=0) return
      iniP=Value
      chains(i,1)=Value
      ! 
      !
      ! Update Gama.
      BB: do j=1, ntry
          call SliceMam4(iniP,iniGama,iniMu,iniSigma,nED,yy,xx,&
                         2,Value,iflag,w,m,lowerGama,upperGama)
          if (iflag==0) exit BB
      end do BB
      ! Error checking.
      if (iflag/=0) return
      iniGama=Value
      chains(i,2)=Value
      !
      !
      ! Update Mu.
      CC: do j=1, ntry
          call SliceMam4(iniP,iniGama,iniMu,iniSigma,nED,yy,xx,&
                         3,Value,iflag,w,m,lowerGama,upperGama)
          if (iflag==0) exit CC
      end do CC
      ! Error checking.
      if (iflag/=0) return
      iniMu=Value
      chains(i,3)=Value
      !
      !
      ! Update Sigma.
      DD: do j=1, ntry
          call SliceMam4(iniP,iniGama,iniMu,iniSigma,nED,yy,xx,&
                         4,Value,iflag,w,m,0.0D+00,upperSigma)
          if (iflag==0) exit DD
      end do DD
      ! Error checking.
      if (iflag/=0) return
      iniSigma=Value
      chains(i,4)=Value
  end do
  !
  return
end subroutine mcMam4
!
!-----------------------------------------------------------------------------------------------------------
!
subroutine SliceMam4(iniP,iniGama,iniMu,iniSigma,nED,ED,Error,&
                     which,Value,iflag,w,m,lower,upper)
! =====================================================================================================
! Update parameters in a MAM-4 age model with the Slice Sampling.
!
! Author:: Peng Jun, 2014.03.13.
!
! Dependence:: Inner function funcP(); funcGama(); funcMu; funcSigma(); 
!              external function alnorm(); external subroutine pnorm().
!
! Reference :: Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
!              vol. 31, no. 3, pp. 705-767.
!
! Note:: THIS CODE IS TRANSFORMED (with minor modification) FROM THE PRIMITIVE R CODE OF RM NEAL AT:
!        http://www.cs.utoronto.ca/~radford/slice.software.html
! =====================================================================================================
  implicit none
  real   (kind=8), intent(in):: iniP         ! Initial proportion.
  real   (kind=8), intent(in):: iniGama      ! Initial Gama value.
  real   (kind=8), intent(in):: iniMu        ! Initial Mu value.
  real   (kind=8), intent(in):: iniSigma     ! Initial Sigma value.
  integer(kind=4), intent(in):: nED          ! Number of equivalent doses.
  integer(kind=4), intent(in):: which        ! Which parameter to be updated (1=p, 2=gama, 3=mu, 4=sigma).
  real   (kind=8), intent(in):: ED(nED)      ! Equivalent dose values (either logged or un-logged).
  real   (kind=8), intent(in):: Error(nED)   ! Standard errors ( either relative or absolute)
  real   (kind=8), intent(in):: w            ! Size of the steps for creating interval.
  real   (kind=8), intent(in):: m            ! Limit on steps (m<=1.0 for no limit, m>1.0 for putting some limit).
  real   (kind=8), intent(in):: lower        ! Lower boundary of the desired distribution.
  real   (kind=8), intent(in):: upper        ! Upper boundary of the desired distribution.
  real   (kind=8), intent(out):: Value       ! Updated parameter corresponding to "which".
  integer(kind=4), intent(out):: iflag       ! Error message (0=success, 1=fail).
  ! Local variables
  real   (kind=8):: ran
  real   (kind=8):: gx0, gx1
  real   (kind=8):: logy
  real   (kind=8):: U, L, R, J, K
  real   (kind=8):: gL, gR
  !
  iflag=0
  Value=-99.0
  ! 
  ! Calculate the specified function value.
  if (which==1) then 
      gx0=funcP(iniP)
  else if (which==2) then
      gx0=funcGama(iniGama)
  else if (which==3) then
      gx0=funcMu(iniMu)
  else if (which==4) then
      gx0=funcSigma(iniSigma)
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
  if (which==1) then
      L=iniP - U
      R=iniP + (w-U)
  else if (which==2) then
      L=iniGama - U
      R=iniGama + (w-U)
  else if (which==3) then
      L=iniMu - U
      R=iniMu + (w-U)
  else if (which==4) then
      L=iniSigma - U
      R=iniSigma + (w-U)
  end if
  !
  ! Expand the interval until its ends are outside the
  ! slice, or until the limit on steps is reached.
  if (m<=1.0) then
      ! No limit on number of steps for the case m<=1.0.
      !
      ! For the left side.
      do  
          if (L<=lower) exit
          !
          if (which==1) then
              gL=funcP(L)
          else if (which==2) then
              gL=funcGama(L)
          else if (which==3) then
              gL=funcMu(L)
          else if (which==4) then
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
          if(R>=upper) exit
          !
          if (which==1) then
              gR=funcP(R)
          else if (which==2) then 
              gR=funcGama(R)
          else if (which==3) then
              gR=funcMu(R)
          else if (which==4) then
              gR=funcSigma(R)
          end if
          ! Error checking.
          if (iflag/=0) return
          !
          if (gR<=logy) exit
          R=R+w
      end do
  else if (m>1.0) then
      ! Limit on steps for the case that m>1.0.
      call random_number(ran)
      J=floor(m*ran)
      K=(m-1.0)-J
      !
      ! For the left side.
      do while (J>0.0) 
          if (L<=lower) exit
          !
          if (which==1) then
              gL=funcP(L)
          else if (which==2) then
              gL=funcGama(L)
          else if (which==3) then
              gL=funcMu(L)
          else if (which==4) then
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
          if (which==1) then 
              gR=funcP(R)
          else if (which==2) then
              gR=funcGama(R)
          else if (which==3) then
              gR=funcMu(R)
          else if (which==4) then
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
      call  random_number(ran)
      Value=L+ran*(R-L)
      !
      if (which==1) then 
          gx1=funcP(Value)
      else if (which==2) then
          gx1=funcGama(Value)  
      else if (which==3) then
          gx1=funcMu(Value)
      else if (which==4) then
          gx1=funcSigma(Value)
      end if
      ! Error checking.
      if (iflag/=0) return
      !
      if (gx1>=logy) exit
      !
      if (which==1) then
          if (Value>iniP) then
              R=Value
          else 
              L=Value
          end if
      else if (which==2) then
          if (Value>iniGama) then 
              R=Value
          else 
              L=Value
          end if
      else if (which==3) then 
          if (Value>iniMu) then
              R=Value
          else 
              L=Value
          end if 
      else if (which==4) then 
          if (Value>iniSigma) then 
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
  ! Some inner functions concerning p, gama, mu, sigma.
  !------------------------------------------------------------------------------------------------
  function funcP(x)
    implicit none
    real   (kind=8):: x
    real   (kind=8):: funcP
    ! local variables
    real   (kind=8):: alnorm
    real   (kind=8):: pnormValues(nED)
    logical, parameter:: right=.false.
    !
    pnormValues= (iniGama- (iniMu/iniSigma**2+ED/Error**2)/(1.0/iniSigma**2+1.0/Error**2) ) * &
                  sqrt(1.0/Error**2+1.0/iniSigma**2)
    ! 
    call pnorm(pnormValues, nED, right)
    !
    funcP=sum(log( x/Error*exp(-(ED-iniGama)**2/2.0/Error**2) + &
                  (1.0-x)/sqrt(Error**2+iniSigma**2)*exp(-(ED-iniMu)**2/2.0/(iniSigma**2+Error**2)) * &
                  (1.0-pnormValues)/(1.0-alnorm((iniGama-iniMu)/iniSigma, right))  ))
    ! Checking NaN (NA)
    if (funcP .ne. funcP) iflag=1
    !
    return
  end function funcP
  !------------------------------------------------------------------------------------------------
  function funcGama(x)
    implicit none
    real   (kind=8):: x
    real   (kind=8):: funcGama
    ! local variables
    real   (kind=8):: alnorm
    real   (kind=8):: pnormValues(nED)
    logical, parameter:: right=.false.
    !
    pnormValues= (x- (iniMu/iniSigma**2+ED/Error**2)/(1.0/iniSigma**2+1.0/Error**2) ) * &
                  sqrt(1.0/Error**2+1.0/iniSigma**2)
    ! 
    call pnorm(pnormValues, nED, right)
    !
    funcGama=sum(log( iniP/Error*exp(-(ED-x)**2/2.0/Error**2) + &
                     (1.0-iniP)/sqrt(Error**2+iniSigma**2)*exp(-(ED-iniMu)**2/2.0/(iniSigma**2+Error**2)) * &
                     (1.0-pnormValues)/(1.0-alnorm((x-iniMu)/iniSigma, right))  ))
    ! Checking NaN (NA)
    if (funcGama .ne. funcGama) iflag=1
    !
    return
  end function funcGama
  !------------------------------------------------------------------------------------------------
  function funcMu(x)
    implicit none
    real   (kind=8):: x
    real   (kind=8):: funcMu
    ! local variables
    real   (kind=8):: alnorm
    real   (kind=8):: pnormValues(nED)
    logical, parameter:: right=.false.
    !
    pnormValues= (iniGama- (x/iniSigma**2+ED/Error**2)/(1.0/iniSigma**2+1.0/Error**2) ) * &
                  sqrt(1.0/Error**2+1.0/iniSigma**2)
    ! 
    call pnorm(pnormValues, nED, right)
    !
    funcMu=sum(log( iniP/Error*exp(-(ED-iniGama)**2/2.0/Error**2) + &
                   (1.0-iniP)/sqrt(Error**2+iniSigma**2)*exp(-(ED-x)**2/2.0/(iniSigma**2+Error**2)) * &
                   (1.0-pnormValues)/(1.0-alnorm((iniGama-x)/iniSigma, right))  ))
    ! Checking NaN (NA)
    if (funcMu .ne. funcMu) iflag=1
    !
    return
  end function funcMu
  !------------------------------------------------------------------------------------------------
  function funcSigma(x)
    implicit none
    real   (kind=8):: x
    real   (kind=8):: funcSigma
    ! local variables
    real   (kind=8):: alnorm
    real   (kind=8):: pnormValues(nED)
    logical, parameter:: right=.false.
    !
    pnormValues= (iniGama- (iniMu/x**2+ED/Error**2)/(1.0/x**2+1.0/Error**2) ) * &
                  sqrt(1.0/Error**2+1.0/x**2)
    ! 
    call pnorm(pnormValues, nED, right)
    !
    funcSigma=sum(log( iniP/Error*exp(-(ED-iniGama)**2/2.0/Error**2) + &
                      (1.0-iniP)/sqrt(Error**2+x**2)*exp(-(ED-iniMu)**2/2.0/(x**2+Error**2)) * &
                      (1.0-pnormValues)/(1.0-alnorm((iniGama-iniMu)/x, right))  ))
    ! Checking NaN (NA)
    if (funcSigma .ne. funcSigma) iflag=1
    !
    return
  end function funcSigma
  !------------------------------------------------------------------------------------------------
end subroutine SliceMam4
