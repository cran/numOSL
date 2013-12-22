subroutine calED2(Dose,ltx,ndose,inltx,outDose,pars,npars,predtval,upb,&
                  nstart,parserrors,value,mcED,method,motoiter,errorflag)
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!  Fitting a dose-response curve with Levenberg-Marquadt method and calculate equivalent dose (standard error) by interpolation (origin).
! ========================================================================================================================================================================
! ndose,             input:: integer, length of redoses (Redose1, Redose2,...).
!
! nstart,            input:: integer, number of random trials.
!
! upb,               input:: real value, upper boundary for b value, b is generated unifiormly in (0, upb).
!
! npars,             input:: integer, model used for fitting the dose-response curve:
!                            if npars=1, a linear model of the form: y=ax will be fitted, where ndat>=1;
!                            if npars=2, a Exponential model of the form: y=a(1-exp(-bx))  will be fitted, where ndat>=2;
!                            if npars=3, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx will be fitted, where ndat>=3.
!
! dose(ndose),       input:: real values, Redose data used to build the dose-response curve.
!
! ltx(ndose,2),      input:: real values, standardlized OSL signal data used to build the dose-response curve.
!
! inltx(2),          input:: real values, standardlized OSL signal values that used for calculating EDs and EDs' standard errors.
!
! outDose(2),        output:: real values, calculated EDs and standard errors that correlated to standardlized OSL signal.
!
! pars(npars),       output:: real values, calculated characteristical parameters of a dose-response curve.
!
! parserrors(npars), output:: real values, estimated parameters' standard errors.
!
! predtval(ndose),   output:: real values, fitted values that correspond to ltx.
!
! value,             output:: real value, sum of the square of residuals.
!
! mcED(motoiter),    output:: real values, equivalent doses generated with monte carlo simulation.
!
! motoiter,           input:: integer, numbers of Moto Carlo simulations that using to estimate EDs' standard errors.
!
! method,             input:: integer, method for calculating dose's standard error, 1 for simple method; 2 for Moto Carlo method.
!
! errorflag(2),      output:: integer, error message generated when calling subroutine calED2():
!                             1.1) if successed in calling lmfit2() or linearfit() to estimate parameters, errorflag(1)=123;
!                             1.2) if fail in calling lmfit2() to estimate parameters, errorflag(1) will be one in (0,4,5,6,7);
!                             1.3) if fail in calling linearfit() to estimate parameters, errorflag(1) will be 0;
!                             2.1) if successed in calculating parameters' standard errors, errorflag(2)=0;
!                             2.2) if fail in calling linearfit() or lmfit2() to estimate parameters' standard errors, errorflag(2) will be one in (1,2,3).
! =====================================================================================================================================================================
! Author:: Peng Jun, 2013.09.20, revised in 2013.09.21.
!
! Dependence:: subroutine inipars2; subroutine lmfit2; subroutine interpolate2; subroutine linearfit; subroutine r8vec_normal.
!
!              Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot
!              regenerative dose measurements. Ancient TL 25, pp. 15-24.
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                     intent(in)::ndose
  integer(kind=4),                     intent(in)::npars
  integer(kind=4),                     intent(in)::motoiter
  integer(kind=4),                     intent(in)::method
  integer(kind=4),                     intent(in)::nstart
  real   (kind=8),                     intent(in)::upb
  real   (kind=8),dimension(ndose),    intent(in)::Dose
  real   (kind=8),dimension(ndose,2),  intent(in)::ltx
  real   (kind=8),dimension(2),        intent(in)::inltx
  real   (kind=8),dimension(2),        intent(out)::outDose
  real   (kind=8),dimension(npars),    intent(out)::pars
  real   (kind=8),dimension(npars),    intent(out)::parserrors
  real   (kind=8),dimension(ndose),    intent(out)::predtval
  real   (kind=8),dimension(motoiter), intent(out)::mcED
  real   (kind=8),                     intent(out)::value
  integer(kind=4),dimension(2),        intent(out)::errorflag
  !
  ! Local variables
  real   (kind=8),parameter::lmtol=1.0D-07                 ! toleracne for stopping iterations for subroutine lmfit2
  real   (kind=8)::minvalue                                ! minimized value of subroutine interpolate2
  real   (kind=8)::maxDose                                 ! maximum Redose value
  real   (kind=8)::averageErr                              ! average value of the sum of the squared residual errors
  real   (kind=8),dimension(2)::lowup                      ! low and up boundary dose values
  real   (kind=8)::lowltx,upltx                            ! low and up boundary ltx value
  integer(kind=4)::seed                                    ! seed for normal distribution random number generation
  real   (kind=8),dimension(ndose)::ranltx                 ! random Lx/Tx
  real   (kind=8),dimension(npars)::Ranpars,Ranparserrors  ! characteristic parameters of simulated dose-response curve 
  real   (kind=8),dimension(ndose)::Ranpredtval            ! fitted values correspond to ranltx
  real   (kind=8)::ranvalue                                ! sum of the square of residual generaged with random simulation
  integer(kind=4),dimension(2)::Ranerrorflag               ! errorflag generated during random fitting
  integer(kind=4)::i,j,mccount                             ! iterative count
  real   (kind=8)::sumdose,sumdose2,mcDose                 ! cumulative sum of ED, sum of the square of ED, mc ED of
  !                                                        ! Moto Carlo ED values
  ! Variables for subroutine inipars2()
  real   (kind=8),parameter::initol=1.0D-09
  real   (kind=8),dimension(1):: rand
  integer(kind=4),dimension(2)::info
  real   (kind=8),dimension(npars-1)::outpars
  real   (kind=8),dimension(3)::cpars
  real   (kind=8),dimension(npars)::lmpars,lmparserrors
  real   (kind=8),dimension(ndose)::lmpredtval 
  integer(kind=4),dimension(2)::lmerrorflag
  real   (kind=8)::bvalue,lmvalue,itervalue
  !
  ! Default output if encountering any error during the calculation
  outDose=-99.0D+00
  mcED=-99.0D+00
  value=-99.0D+00
  pars=-99.0D+00
  parserrors=-99.0D+00
  predtval=-99.0D+00
  ! Default returned errorflag if appearing any error
  errorflag(1)=0
  errorflag(2)=1
  !
  ! Initialize parameters
  if(npars==1) then  ! For linear model
    ! Call linearfit() to obtain the characterized parameters of the dose-response
    ! curve, return lmpars, lmparserrors, lmpredtval, lmvalue and lmerrorflag
    call linearfit(Dose,ltx(:,1),ndose,lmpars,npars,&
                   .true.,lmparserrors,lmpredtval,&
                   lmvalue,lmerrorflag) 
    pars=lmpars
    parserrors=lmparserrors
    predtval=lmpredtval
    value=lmvalue
    errorflag=lmerrorflag
  else if(npars==2 .or. npars==3) then ! For expentional or linear plus expentional model
    cpars=0.0D+00
    itervalue=1.0D+30
    call random_seed()
    ! Try to fit the dose-response curve with upb values
    loop: do i=1, nstart
      ! Generate upb 
      call random_number(rand)
      bvalue=upb*rand(1)
      ! Use generated upb to initailize parameters
      call inipars2(bvalue,npars+1,ndose,dose,&
                    ltx(:,1),outpars,initol,info)
      ! Error checking of inipars2()
      if(info(1)==0 .and. info(2)==0) then
        if(npars==2)  then
          cpars(1)=outpars(1)
          cpars(2)=bvalue
        else if(npars==3) then
          cpars(1)=outpars(1)
          cpars(2)=bvalue
          cpars(3)=outpars(2)
        end if
        ! Pass cpars to lmpars
        lmpars=cpars(1:npars)
      end if
      ! Cycle the looping if fail in parameters initialization
      if(info(1)/=0 .or. info(2)/=0) cycle loop
      ! Call lmfit2() to get the characterized parameters for Dose-Response 
      ! curve, return pars, parserrors,predtval,value, errorflag
      call lmfit2(Dose,ltx(:,1),ndose,lmpars,npars,&
                  .true.,lmparserrors,lmpredtval,&
                  lmvalue,lmtol,lmerrorflag) 
      ! Check if any error appears in lmfit2()
      if(lmerrorflag(1)==123  .and. lmvalue<itervalue) then
        ! Store improved results only
        itervalue=lmvalue
        pars=lmpars
        parserrors=lmparserrors
        predtval=lmpredtval
        value=lmvalue
        errorflag=lmerrorflag
      end if
    end do loop
  end if
  ! Error checking, if error appears in linearfit() or lmfit2(), return
  if(errorflag(1)/=123) return
  !
  maxDose=maxval(Dose)
  ! Calculate equivalent dose using inltx(1), store it in outDose(1)
  call interpolate2(outDose(1),inltx(1),pars,npars,&
                    0.0D+00,maxDose*1.3D+00,minvalue)
  ! Estimate standard error of ED
  if(method==1) then
    ! 1) Method 1, simple transformation
    averageErr=value/real(ndose,kind=8)
    lowltx=inltx(1)-dsqrt((inltx(2))**2+averageErr)
    upltx =inltx(1)+dsqrt((inltx(2))**2+averageErr)
    ! Calculate low bounded ED, store it in lowup(1)
    call interpolate2(lowup(1),lowltx,pars,npars,&
                      0.0D+00,maxDose*1.3D+00,minvalue)
    ! Calculate up bounded ED, store it in lowup(2)
    call interpolate2(lowup(2),upltx,pars,npars,&
                      0.0D+00,maxDose*1.3D+00,minvalue)
    ! Calculate standard error of ED with simple transformation
    outDose(2)=(lowup(2)-lowup(1))/2.0D+00
  else if(method==2) then
    ! 2) Method 2, Monte Carlo iteration
    ! Set random seed
    seed=332571951
    ! Initializing sum of mccount(mccount), Sum of mcdose(mcdose), Sum of the square of
    ! mcdose(sumdose2) and count number (i), they will be updated during the simulation
    mccount=0
    sumdose=0.0D+00
    sumdose2=0.0D+00
    i=0
    ! Moto Carlo simulations
    do 
      ! Set Ranpars to be pars
      Ranpars=pars
      ! Generate random ltx values with mean=ltx(1,j), 
      ! sd=ltx(2,j), they will be store in array ranltx
      do j=1, ndose
        call r8vec_normal(1,ltx(j,1),ltx(j,2),seed,ranltx(j))
      end do
      ! Fitting x(Dose) .VS. random ltx
      if(npars==1) then
        ! Linear model
        call linearfit(Dose,ranltx,ndose,Ranpars,npars,.false.,&
                       Ranparserrors,Ranpredtval,Ranvalue,Ranerrorflag)
      else if(npars==2 .or. npars==3)  then
        ! Non-linear model
        call lmfit2(Dose,ranltx,ndose,Ranpars,npars,.false.,&
                    Ranparserrors,Ranpredtval,Ranvalue,&
                    lmtol,Ranerrorflag)
      end if
      ! Error check of linearfit() or lmfit2()
      if(Ranerrorflag(1)==123) then
        ! Count really performed monte carlo numbers
        i=i+1
        !! Generate random signal value and store it
        ! in mcED(i) with mean=inltx(1), sd=inltx(2)
        call r8vec_normal(1,inltx(1),inltx(2),seed,mcED(i))
        ! Interpolate using Ranpars and random values mcED(i), 
        ! calculated ED will be stored in the ith mcED
        call interpolate2(mcDose,mcED(i),Ranpars,npars,&
                          0.0D+00,maxDose*1.5D+00,minvalue)
        mcED(i)=mcDose
        ! Updating mccount, sumdoes, sumdose2
        mccount=mccount+1
        sumdose=sumdose+mcDose
        sumdose2=sumdose2+mcDose**2
      end if
      ! Terminate the simulation or not?
      if(i==motoiter) exit
    end do
    ! Calculate standard error for ED
    outDose(2)=dsqrt((real(mccount,kind=8)*sumdose2-sumdose**2)/&
                      real(mccount,kind=8)/real(mccount-1,kind=8))
  end if
  ! Terminiate and return
  return
end subroutine calED2
