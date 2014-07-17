subroutine calED0(Dose,ltx,ndose,ninltx,inltx,outDose,mcED,pars,npars,&
                  predtval,upb,nstart,parserrors,value,method,nsim,errorflag)
!-----------------------------------------------------------------------------------------
! Calculate a number of equivalent dose values (pass the origin).
! ========================================================================================
!
! ndose,               input:: integer, length of regenerative dose values.
!
! ninltx,              input:: integer, number of wanted equivalent doses.
!
! nstart,              input:: integer, number of random trials.
!
! upb,                 input:: real value, upper boundary for b value, b generates uniformly in (0, upb).
!
! npars,               input:: integer, model used for fitting the dose-response curve:
!                        if npars=1, fit a linear model of the form: y=ax (ndat>=1);
!                        if npars=2, fit an exponential model of the form: y=a(1-exp(-bx)) (ndat>=2);
!                        if npars=3, fit a linear+Exponential model of the form: y=a(1-exp(-bx))+cx (ndat>=3).
!
! dose(ndose),         input:: real values, regenerative dose values.
!
! ltx(ndose,2),        input:: real values, standardlized OSL signal values.
!
! inltx(ninltx,2),     input:: real values, inputs used for calculating EDs and standard errors.
!
! outDose(ninltx,2),  output:: real values, calculated EDs (standard errors).
!
! mcED(nsim,ninltx),  output:: real values, simulated ED values.
!
! pars(npars),        output:: real values, calculated characteristical parameters of the growth curve.
!
! parserrors(npars),  output:: real values, parameters' standard errors.
!
! predtval(ndose),    output:: real values, fitted values that correspond to ltx.
!
! value,              output:: real value, sum of the squared residuals.
!
! nsim,                input:: integer, numbers of Moto Carlo simulations used to estimate EDs' standard errors.
!
! method,              input:: integer, method for error assessing, 1 for simple method; 2 for Moto Carlo method.
!
! errorflag(2),       output:: integer, error messages:
!                              if success in fitting growth curve, errorflag(1)=123;
!                              if parameters' standard errors can be estimated and all pars>0.0, errorflag(2)=0.
! =========================================================================================================================
! Author:: Peng Jun, 2013.09.21; revised in 2014.04.05.
!
! Dependence:: subroutine inipars0; subroutine interpolate0; subroutine growFit;  
!              subroutine linearfit; subroutine r8vec_normal
!
!              Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived 
!              from single aliquot regenerative dose measurements. Ancient TL 25, pp. 15-24.
!--------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                        intent(in):: ndose
  integer(kind=4),                        intent(in):: npars
  integer(kind=4),                        intent(in):: ninltx
  integer(kind=4),                        intent(in):: nsim
  integer(kind=4),                        intent(in):: method
  integer(kind=4),                        intent(in):: nstart
  real   (kind=8),                        intent(in):: upb
  real   (kind=8),dimension(ndose),       intent(in):: Dose
  real   (kind=8),dimension(ndose,2),     intent(in):: ltx
  real   (kind=8),dimension(ninltx,2),    intent(in):: inltx
  real   (kind=8),                       intent(out):: value
  integer(kind=4),dimension(2),          intent(out):: errorflag
  real   (kind=8),dimension(ninltx,2),   intent(out):: outDose
  real   (kind=8),dimension(nsim,ninltx),intent(out):: mcED
  real   (kind=8),dimension(npars),      intent(out):: pars
  real   (kind=8),dimension(npars),      intent(out):: parserrors
  real   (kind=8),dimension(ndose),      intent(out):: predtval
  !
  ! Local variables.
  real   (kind=8),parameter          :: lmtol=1.0D-07 
  real   (kind=8),parameter          :: initol=4.053817e-10  !.Machine$double.eps^0.6 in R 
  real   (kind=8),dimension(3)       :: cpars
  real   (kind=8),dimension(npars-1) :: outpars 
  real   (kind=8),dimension(ninltx)  :: lowltx, upltx
  real   (kind=8),dimension(ninltx,2):: lowup
  real   (kind=8),dimension(ndose)   :: lmpredtval, simltx
  real   (kind=8),dimension(npars)   :: lmpars, lmparserrors
  real   (kind=8)                    :: mcDose, sumDose, sumDose2, maxDose, minValue,&
                                        avgErr, bvalue, lmvalue, loopvalue, ran, mcSig(1)
 
  integer(kind=4):: seed, mcCount, info, i, j, lmerrorflag(2)
  !
  ! Default returned values.
  value=-99.0
  ! Initialize errorflag.
  errorflag(1)=0
  errorflag(2)=1
  outDose=-99.0
  mcED=-99.0
  pars=-99.0
  parserrors=-99.0
  predtval=-99.0
  !
  ! Fitting the growth curve.
  if (npars==1) then  
      ! for a linear model.
      call linearfit(Dose,ltx(:,1),ndose,pars,npars,&
                     parserrors,predtval,value,errorflag)
      if (any(pars<=0.0D+00))  errorflag(2)=1
  else if (npars==2 .or. npars==3) then
      ! for an expentional or a linear plus expentional model.
      cpars=0.0
      loopvalue=1.0D+30
      call random_seed()
      Loop: do i=1, nstart
          call random_number(ran)
          bvalue=upb*ran
          call inipars0(bvalue,npars,ndose,Dose,&
                        ltx(:,1),outpars,initol,info)
          ! Error checking
          if (info/=0) cycle Loop
          !
          if (info==0) then
              cpars(1)=outpars(1)
              cpars(2)=bvalue             
              if (npars==3) cpars(3)=outpars(2)
              lmpars=cpars(1:npars)
          end if
          !
          call growFit(Dose,ltx(:,1),ndose,lmpars,npars,.true.,&
                       lmparserrors,lmpredtval,lmvalue,lmtol,lmerrorflag)
          !
          ! Check if the exponential curve becomes saturation.
          if (npars==2 .and. lmpars(1)*lmpars(2)*dexp(-lmpars(2)*maxval(Dose))<=2.013409D-05) cycle Loop
          !
          ! Only store improved estimates.
          if(lmerrorflag(1)==123 .and. lmerrorflag(2)==0 .and. &
             all(lmpars>0.0D+00) .and. lmvalue<loopvalue) then
              loopvalue=lmvalue
              pars=lmpars
              parserrors=lmparserrors
              predtval=lmpredtval
              value=lmvalue
              errorflag=lmerrorflag
          end if
      end do Loop
  end if
  !
  ! Error checking, return if fails.
  if(errorflag(1)/=123 .or. errorflag(2)/=0) return
  !
  ! Set upper limits on maxDose.
  if (npars==2 .or. npars==3) then
      maxDose=1.5D+00*maxval(Dose)
  end if
  !
  ! Calculate ED values.
  if (npars==1) then
      outDose(:,1)=inltx(:,1)/pars(1)
  else if (npars==2 .or. npars==3) then
      do i=1, ninltx
          call interpolate0(outDose(i,1),inltx(i,1),pars,npars,&
                            -5.0D+00,maxDose,minValue)
      end do
  end if
  !
  !
  ! Calculate standard errors of EDs.
  if(method==1) then
      ! 1) Method 1, simple transformation.
      avgErr=value/real(ndose,kind=8)
      lowltx=inltx(:,1)-sqrt((inltx(:,2))**2+avgErr)
      upltx =inltx(:,1)+sqrt((inltx(:,2))**2+avgErr)
      !
      if (npars==1)  then
          lowup(:,1)=lowltx(:)/pars(1)
          lowup(:,2)= upltx(:)/pars(1)
          !
      else if (npars==2 .or. npars==3) then
          do i=1, ninltx
              call interpolate0(lowup(i,1),lowltx(i),pars,npars,&
                                -5.0D+00,maxDose,minValue)
              call interpolate0(lowup(i,2),upltx(i),pars,npars,&
                                -5.0D+00,maxDose,minValue)
          end do
      end if
      !
      outDose(:,2)=(lowup(:,2)-lowup(:,1))/2.0D+00
      !
  ! 2) Method 2, Monte Carlo iteration
  else if (method==2) then
      ! Set random seed. 
      seed=123456789
      ! Looping the ith inltx.
      do i=1, ninltx
          mcCount=0
          sumDose=0.0
          sumDose2=0.0
          !
          ! Moto Carlo simulations.
          Inner: do 
              ! Generate random ltx values.
              do j=1, ndose
                  call r8vec_normal(1,ltx(j,1),ltx(j,2),seed,simltx(j))
              end do
              ! Fitting simulated observations.
              if (npars==1) then
                  ! for a linear model.
                  call linearfit(Dose,simltx,ndose,lmpars,npars,&
                                 lmparserrors,lmpredtval,lmvalue,lmerrorflag)
                  !
              else if (npars==2 .or. npars==3) then
                  ! for an expentional or a linear plus expentional model.
                  call random_number(ran)
                  bvalue=upb*ran
                  call inipars0(bvalue,npars,ndose,Dose,&
                                simltx,outpars,initol,info)
                  ! Error checking
                  if (info/=0) cycle Inner
                  !
                  if (info==0) then
                      cpars(1)=outpars(1)
                      cpars(2)=bvalue             
                      if (npars==3) cpars(3)=outpars(2)
                      lmpars=cpars(1:npars)
                  end if
                  !
                  call growFit(Dose,simltx,ndose,lmpars,npars,.true.,&
                               lmparserrors,lmpredtval,lmvalue,lmtol,lmerrorflag)
                  !
                  ! Check if the exponential curve becomes saturation.
                  if (npars==2 .and. lmpars(1)*lmpars(2)*dexp(-lmpars(2)*maxval(Dose))<=2.013409D-05) cycle Inner
              end if  
              !  
              ! Error check.
              if (lmerrorflag(1)==123 .and. lmerrorflag(2)==0 .and. all(lmpars>0.0D+00)) then
                  ! Simulate mcSig.
                  call r8vec_normal(1,inltx(i,1),inltx(i,2),seed,mcSig(1))
                  ! Interpolation.
                  if (npars==1) then
                      mcDose=mcSig(1)/lmpars(1)
                  else if (npars==2 .or. npars==3) then
                      call interpolate0(mcDose,mcSig(1),lmpars,npars,&
                                        -5.0D+00,maxDose,minValue)
                      if ( (minValue .ne. minValue) .or. (minValue>1.0D-03) ) cycle Inner
                  end if
                  !
                  ! Accumulate McCount, SumDose, Sum of EDs.
                  mcCount=mcCount+1
                  sumDose=sumDose+mcDose
                  sumDose2=sumDose2+mcDose**2
                  mcED(mcCount,i)=mcDose
              end if
              ! Checking termination.
              if(mcCount==nsim) exit Inner
          end do Inner
          !
          ! Calculate standard error for the ith ED
          outDose(i,2)=sqrt((real(mcCount,kind=8)*sumDose2-sumDose**2)/ &
                             real(mcCount,kind=8)/real(mcCount-1,kind=8))
      end do
  end if
  !
  return
end subroutine calED0
