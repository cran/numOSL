subroutine diffev(npars,tim,sig,ntim,lamda,ithn,value,&
                  predict,agents,np,f,cr,maxiter,tol,errorflag)
!---------------------------------------------------------------------------------------------------------------------------------------
! Fitting OSL signal decay curve using differential evolution algorithm, 
! I=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t) will be fitted.
! ======================================================================================================================================
! npars,                input:: integer, the dimension of the problem (length of parameters).
!
! tim(ntim),            input:: real values, the time values.
!
! sig(lntim),           input:: real values, the signal values.
!
! ntim,                 input:: integer, the length of the observations (tim or signal).
!
! value,               output:: real value, minimized sum of the square of rediduals.
!
! lamda(npars),        output:: real values, estimated lamda values.
!
! ithn(npars),         output:: real values, estimated ithn values.
!
! predict(ntim),       output:: real values, fitted values correspond to signal values.
!
! agents(np,npars),    output:: real values, agents values calculated by differential evolution.
!
! np,                   input:: integer, argument used of differential evolution algorithm.
!
! f,                    input:: real value, argument used for differential evolution algorithm.
!
! cr,                   input:: real value, argument used for differential evolution algorithm.
!
! maxiter,              input:: integer, the maximum iterative number allowed for differential evolution algorithm.
!
! tol,                  input:: real value, the iteration of differential evolution will be stopped if tol has been reached.
!
! errorflag(2),        output:: integer, error message generated during the calling:
!                               1.1) if tol has been reached withn the maximum iterative number, errorflag(1)=0;
!                               1.2) if tol has not been reached withn the maximum iterative number, errorflag(1)=1;
!                               2.1) if find at least one appropriate agent, errorflag(2)=0;
!                               2.2) if fail in find at least one appropriate agent, errorflag(2)=1.
! =======================================================================================================================================
! Author:: Peng Jun, 2013.05.31; revised in 2013.06.03; revised in 2013.12.14.
! 
! Dependence:: subroutine leaveone; subroutine NonReplaceSample; 
!              subroutine targfunc; subroutine sort; subroutine sd.
!
! References:: Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
!              OSL decay curves. Radiation Measurements 41, 886-891.             
!
!              Differential evolution, http://en.wikipedia.org/wiki/Differential_evolution
!----------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                    intent(in)::npars
  integer(kind=4),                    intent(in)::ntim
  integer(kind=4),                    intent(in)::np
  integer(kind=4),                    intent(in)::maxiter
  real   (kind=8),                    intent(in)::f
  real   (kind=8),                    intent(in)::cr
  real   (kind=8),                    intent(in)::tol
  real   (kind=8),dimension(ntim),    intent(in)::tim
  real   (kind=8),dimension(ntim),    intent(in)::sig
  real   (kind=8),                    intent(out)::value
  integer(kind=4),dimension(2),       intent(out)::errorflag
  real   (kind=8),dimension(npars),   intent(out)::lamda
  real   (kind=8),dimension(npars),   intent(out)::ithn
  real   (kind=8),dimension(ntim),    intent(out)::predict
  real   (kind=8),dimension(np,npars),intent(out)::agents
  !
  ! local variables
  real   (kind=8)                    ::ran1      ! a ramdom number for comparing with critical random values
  real   (kind=8),dimension(npars)   ::ran       ! random numbers used for initializing agents
  real   (kind=8),dimension(2,npars) ::space     ! vector space
  integer(kind=4),dimension(npars)   ::colsamp   ! index of generated random offspring
  real   (kind=8),dimension(npars)   ::y         ! random offspring vectors in combaining with agents3
  integer(kind=4),dimension(1)       ::indx      ! one sampled index from colsamp
  integer(kind=4),dimension(np)      ::npsamp    ! index of row of agents
  integer(kind=4),dimension(np-1)    ::samp      ! index for (np-1) rows obtained from agents' index
  integer(kind=4),dimension(3)       ::samp3     ! three random index sampled from samp
  real   (kind=8),dimension(3,npars) ::agents3   ! three random vectors obtained from agents
  real   (kind=8),dimension(npars)   ::yithn,agentsithn
  real   (kind=8)                    ::yvalue,agentsvalue
  integer(kind=4)                    ::errorflagy, errorflagagents
  real   (kind=8),dimension(npars)   ::outsd
  real   (kind=8),dimension(np)      ::valueArary
  integer(kind=4)                    ::MinIndex
  real   (kind=8),parameter          ::targtol=1.0D-09
  real   (kind=8),dimension(2,7),parameter::tspace=reshape((/1.0D-01,10.0D+00,&        
                                                             1.0D-04, 5.0D+00,&
                                                             1.0D-04, 3.0D+00,&
                                                             1.0D-05, 1.0D+00,&  
                                                             1.0D-05, 0.5D+00,&
                                                             1.0D-06, 0.1D+00,&
                                                             1.0D-06, 0.05D+00/), (/2,7/) )      ! the total vector space 
  integer(kind=4)::i,j,iter
  integer(kind=4),parameter::typ=1
  !
  ! initialzing space using values in tspace
  space=tspace(:,1:npars)
  !
  ! generate uniform random numbers 
  ! to initialize matrix agents
  call random_seed()
  do i=1,np
      call random_number(ran)
      agents(i,:)=space(1,:)+ran*( space(2,:)-space(1,:) ) 
  end do 
  !
  ! initializing colsamp and npsamp
  colsamp=(/(i,i=1,npars)/)
  npsamp=(/(i,i=1,np)/)
  !
  ! start the process of differential evolution
  iter=0
  valueArary=1.01D+30
  do i=1, maxiter
    ! update agents by row
    do j=1,np
      !
      ! initializing samp, samp must exclude j
      call leaveone(npsamp,np,j,samp)
      !
      ! Picking out three rows form agents randomly,
      ! they must be not only distinct from each other
      ! out also distinct from the agent values in the jth row
      call NonReplaceSample(samp,samp3,np-1,3)
      !
      ! storing three random agents vector in agents3
      agents3=agents(samp3,:)
      !
      ! Picking out a random index ranging from 1 to npars from colsamp
      call NonReplaceSample(colsamp,indx,npars,1)
      !
      ! Generating a uniformly distributed number from (0,1)
      call random_number(ran1)
      !
      ! Generate y based on obtained agents3
      if(ran1<cr .or. j==indx(1) )  then
        y=agents3(1,:) + f*( agents3(2,:)-agents(3,:) )
      else 
        y=agents(j,:)
      end if
      !
      ! Updating agents values in the jth row 
      ! if all  values in y are larger than 0
      if( all(y>0.0D+00) )  then
        !
        ! calculate correspond ithn values, lamda values using y values
        call targfunc(y,npars,tim,sig,targtol,typ,&
                      ntim,yvalue,yithn,errorflagy)
        !
        ! calculate correspond ithn values, lamda values, using the ith agent values
        call targfunc(agents(j,:),npars,tim,sig,targtol,typ,&
                      ntim,agentsvalue,agentsithn,errorflagagents)
        !
        ! replace the ith agent by y if:
        ! 1) function value can be successfully calculated using Linear Algebra
        !    method and all ithn values for y are larger than 0; 
        ! 2) function value of y not less than function value for the ith agent
        ! ***
        ! At this points also note that we don't care if error appears when calling
        ! targfunc using the jth agents value, if so, agentsvalue will be 1.0D+30,
        ! then y will always be accepted as long as errorflagy=0
        if( errorflagy==0 .and. yvalue<=agentsvalue )  then
           ! sort y in descending order
           call sort(y,npars)
           agents(j,:)=y
           valueArary(j)=yvalue
        end if
        !
      end if
      !
    end do
    ! after each updation of matrix agents, we check if its flat enough for stopping.
    ! standard deviation values for final agents, calculating by cloumn order,
    ! to check if matrix agents is flat enough for output
    call sd(agents,np,npars,outsd)
    !
    ! if matrix agents is flat enough, stopping the whole iteration, else continue
    if( all(outsd<=tol) )  exit 
    ! add iterative number by 1 each time
    iter=iter+1
  end do 
  !
  ! check convergence 
  errorflag(1)=0
  if(iter>=maxiter)  errorflag(1)=1
  !
  ! calculate lamda values for output, looping throuth
  ! matrix agents row by row, to pick out a row that
  ! given the miniumu function value.
  value=1.0D+30
  ! if no appropriate agent has been found, set errorflag(2)=1
  errorflag(2)=1
  do i=1, np
    if( valueArary(i)<value )  then
      MinIndex=i
      value=valueArary(i)
      errorflag(2)=0
    end if
  end do
  !
  if(errorflag(2)==0) then
    call targfunc(agents(MinIndex,:),npars,tim,sig,targtol,typ,&
                  ntim,agentsvalue,agentsithn,errorflagagents)
    value=agentsvalue
    ithn=agentsithn
    lamda=agents(MinIndex,:)
    ! set predict
    predict=0.0D+00
    do i=1,npars
      predict=predict+ithn(i)*dexp(-lamda(i)*tim)
    end do
  else 
    ! if errorflag(2)=1, then return all values 
    ! to -99.0
    value=-99.0D+00
    ithn=-99.0D+00
    lamda=-99.0D+00
    predict=-99.0D+00
  end if
  !
  return
end subroutine diffev
