subroutine sgcED2(Dose,ltx,ndose,ninltx,inltx,outDose,pars,npars,predtval,&
                  upb,nstart,parserrors,value,method,motoiter,errorflag)
!--------------------------------------------------------------------------------------------------------------------------------------------
! Calculate equivalent dose values using Standardised growth curves (SGC) method (origin).
! ===========================================================================================================================================
!
! ndose,               input:: integer, length of reDoses (Redose1, Redose2,...) or Ltx(Lx1/Tx1, Lx2/Tx2,...).
!
! ninltx,              input:: integer, number of standardlised OSL signal values from which to calculate EDs.
!
! nstart,              input:: integer, number of random trials.
!
! upb,                 input:: real value, upper boundary for b value, b is generated in (0, upb).
!
! npars,               input:: integer, model used for fitting the dose-response curve:
!                              if npars=1, a linear model of the form: y=ax will be fitted, where ndat>=1;
!                              if npars=2, a Exponential model of the form: y=a(1-exp(-bx))  will be fitted, where ndat>=2;
!                              if npars=3, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx will be fitted, where ndat>=3.
!
! dose(ndose),         input:: real values, Redose data used to build a dose-response curve.
!
! ltx(ndose,2),        input:: real values, standardlized OSL signal data used to build a dose-response curve.
!
! inltx(ninltx,2),     input:: real values, standardlized OSL signal values that used for calculating EDs and standard errors.
!
! outDose(ninltx,2),  output:: real values, calculated EDs and standard errors that correlated to standardlized OSL signal.
!
! pars(npars),        output:: real values, calculated characteristical parameters of a dose-response curve.
!
! parserrors(npars),  output:: real values, parameters' standard errors.
!
! predtval(ndose),    output:: real values, fitted values that correspond to ltx.
!
! value,              output:: real value, sum of the square of residual.
!
! motoiter,            input:: integer, numbers of Moto Carlo simulations that using to estimate EDs' standard errors.
!
! method,              input:: integer, method for calculating EDs' standard error, 1 for simple method; 2 for Moto Carlo method.
!
! errorflag(2),       output:: integer, error message generated when calling subroutine SGC:
!                              1.1) if no error appears when calling linearfit() or lmfit1() to estimate parameters, errorflag(1)=123;
!                              1.2) if error appears when calling lmfit1() to estimate parameters, errorflag(1) will be one in (0,4,5,6,7,8,9,10);
!                              1.3) if error appears when calling linearfit() to estimate parameters, errorflag(1) will be 1;
!                              2.1) if no error appears when calling linearfit() or lmfit1() to estimate parameters' std.error, errorflag(2)=0;
!                              2.2) if error appears when calling linearfit() or lmfit1() to estimate parameters' std.error, errorflag(2)=1.
! ================================================================================================================================================
! Author:: Peng Jun, 2013.09.21.
!
! Dependence:: subroutine inipars; subroutine lmfit2; subroutine interpolate2; subroutine linearfit; subroutine r8vec_normal
!
! References:: Roberts,H.M. and Duller,G.A.T., 2004. Standardised growth curves for optical dating of sediment
!              using multiple-grain aliquots. Radiation Measurements 38, pp. 241-252.
!
!              Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot
!              regenerative dose measurements. Ancient TL 25, pp. 15-24.
!----------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                    intent(in)::ndose
  integer(kind=4),                    intent(in)::npars
  integer(kind=4),                    intent(in)::ninltx
  integer(kind=4),                    intent(in)::motoiter
  integer(kind=4),                    intent(in)::method
  integer(kind=4),                    intent(in)::nstart
  real   (kind=8),                    intent(in)::upb
  real   (kind=8),dimension(ndose),   intent(in)::Dose
  real   (kind=8),dimension(ndose,2), intent(in)::ltx
  real   (kind=8),dimension(ninltx,2),intent(in)::inltx
  real   (kind=8),                    intent(out)::value
  integer(kind=4),                    intent(out)::errorflag(2)
  real   (kind=8),dimension(ninltx,2),intent(out)::outDose
  real   (kind=8),dimension(npars),   intent(out)::pars
  real   (kind=8),dimension(npars),   intent(out)::parserrors
  real   (kind=8),dimension(ndose),   intent(out)::predtval
  ! local variables
  real   (kind=8),parameter::lmtol=1.0D-07 
  real   (kind=8),parameter::initol=1.0D-09
  real   (kind=8)::maxDose,minvalue
  real   (kind=8),dimension(1)::rand
  integer(kind=4),dimension(2)::info
  real   (kind=8),dimension(3)::cpars
  real   (kind=8),dimension(npars-1)::outpars 
  real   (kind=8),dimension(ninltx)::lowltx,upltx
  real   (kind=8),dimension(ninltx,2)::lowup
  real   (kind=8),dimension(nDose)::ranltx, Ranpredtval
  real   (kind=8),dimension(npars)::Ranpars, Ranparserrors
  integer(kind=4),dimension(2)::Ranerrorflag
  real   (kind=8)::Ranvalue,mcDose,sumdose,sumdose2
  real   (kind=8)::averageerr
  real   (kind=8)::mcSig(1)
  integer(kind=4)::seed
  integer(kind=4)::mccount
  integer(kind=4)::i, j
  real   (kind=8),dimension(npars)::lmpars,lmparserrors
  real   (kind=8),dimension(ndose)::lmpredtval
  integer(kind=4),dimension(2)::lmerrorflag
  real   (kind=8)::bvalue, lmvalue, itervalue
  ! Default returned values if any error appears
  outDose=-99.0D+00
  value=-99.0D+00
  pars=-99.0D+00
  parserrors=-99.0D+00
  predtval=-99.0D+00
  ! Initialize errorflag
  errorflag(1)=0
  errorflag(2)=1
  if(npars==1) then  ! for linear model
    ! Call linearfit to get the characterized parameters for Dose-Response 
    ! curve, return lmpars, lmparserrors,lmpredtval,lmvalue,lmerrorflag
    call linearfit(Dose,ltx(:,1),ndose,lmpars,npars,&
                   .true.,lmparserrors,lmpredtval,&
                   lmvalue,lmerrorflag) 
    pars=lmpars
    parserrors=lmparserrors
    predtval=lmpredtval
    value=lmvalue
    errorflag=lmerrorflag
  else if(npars==2 .or. npars==3) then  ! for expentional or linear plus expentional model
    cpars=0.0D+00
    itervalue=1.0D+30
    call random_seed()
    loop: do i=1, nstart
      call random_number(rand)
      bvalue=upb*rand(1)
      call inipars2(bvalue,npars+1,ndose,dose,&
                    ltx(:,1),outpars,initol,info)
      ! Error checking
      if(info(1)==0 .and. info(2)==0) then
        if(npars==2) then
          cpars(1)=outpars(1)
          cpars(2)=bvalue
        else if(npars==3) then
          cpars(1)=outpars(1)
          cpars(2)=bvalue
          cpars(3)=outpars(2)
        end if
        lmpars=cpars(1:npars)
      end if
      if(info(1)/=0 .or. info(2)/=0) cycle loop
      ! Call lmfit2 to get the characterized parameters for dose-response 
      ! curve, return lmpars, lmparserrors, lmpredtval, lmvalue, lmerrorflag
      call lmfit2(Dose,ltx(:,1),ndose,lmpars,npars,&
                  .true.,lmparserrors,lmpredtval,&
                  lmvalue,lmtol,lmerrorflag) 
      ! Error checking
      ! Store improved results only
      if(lmerrorflag(1)==123 .and. lmvalue<itervalue) then
        itervalue=lmvalue
        pars=lmpars
        parserrors=lmparserrors
        predtval=lmpredtval
        value=lmvalue
        errorflag=lmerrorflag
      end if
    end do loop
  end if
  ! Error checking, if linearfit or lmfit2 fails, return
  if(errorflag(1)/=123) return
  ! Calculate ED values corresponding to inltx
  maxDose=maxval(Dose)
  do i=1,ninltx
    call interpolate2(outDose(i,1),inltx(i,1),pars,npars,&
                      0.0D+00,maxDose*1.3D+00,minvalue)
  end do
  ! Calculate standard errors of ED values
  if(method==1) then
    ! 1) Method 1, simple transformation
    averageErr=value/real(ndose,kind=8)
    lowltx=inltx(:,1)-sqrt((inltx(:,2))**2+averageErr)
    upltx =inltx(:,1)+sqrt((inltx(:,2))**2+averageErr)
    do i=1,ninltx
      ! Calculate low bounded Equivalent Dose
      call interpolate2(lowup(i,1),lowltx(i),pars,npars,&
                        0.0D+00,maxDose*1.3D+00,minvalue)
      ! Calculate up bounded Equivalent Dose
      call interpolate2(lowup(i,2),upltx(i),pars,npars,&
                        0.0D+00,maxDose*1.3D+00,minvalue)
    end do
    ! Calculate standard errors for EDs with simple tansformation
    outDose(:,2)=(lowup(:,2)-lowup(:,1))/2.0D+00
  ! 2) Method 2, Monte Carlo iteration
  else if(method==2) then
    ! Looping the ith inltx 
    do i=1, ninltx
      ! Set random seed 
      seed=332571951
      ! Initializing McCount numbers, Sum of Mc Dose, 
      ! Sum of the square of Mc Dose, i, they will be reused
      mccount=0
      sumdose=0.0D+00
      sumdose2=0.0D+00
      ! Moto Carlo simulations
      Inner: do 
        ! Set Ranpars to be pars
        Ranpars=pars
        ! Generate random ltx values with mean=ltx(1,j), 
        ! sd=ltx(2,j), store them in ranltx
        do j=1, ndose
          call r8vec_normal(1,ltx(j,1),ltx(j,2),seed,ranltx(j))
        end do
        ! Fitting x (Dose) .VS. random ltx (ranltx)
        if(npars==1) then
          call linearfit(Dose,ranltx,ndose,Ranpars,npars,.false.,&
                         Ranparserrors,Ranpredtval,Ranvalue,Ranerrorflag)
        else if(npars==2 .or. npars==3) then
          call lmfit2(Dose,ranltx,ndose,Ranpars,npars,&
                      .false.,Ranparserrors,Ranpredtval,&
                      Ranvalue,lmtol,Ranerrorflag)  
        end if    
        ! Error check of linearfit or lmfit2
        if(Ranerrorflag(1)==123) then
          ! Generate random values mcSig with mean=inltx(1), sd=inltx(2)
          call r8vec_normal(1,inltx(i,1),inltx(i,2),seed,mcSig(1))
          ! Interpolation using Ranpars and random values mcSig
          call interpolate2(mcDose,mcSig(1),Ranpars,npars,&
                            0.0D+00,maxDose*1.5D+00,minvalue)
          ! Accumulate McCount, SumDose, Sum of EDs
          mccount=mccount+1
          sumdose=sumdose+mcDose
          sumdose2=sumdose2+mcDose**2
        end if
        ! Checking terminate
        if(mccount==motoiter) exit Inner
      end do Inner
      ! Calculate standard error for the ith ED
      outDose(i,2)=sqrt((real(mccount,kind=8)*sumdose2-sumdose**2)/ &
                         real(mccount,kind=8)/real(mccount-1,kind=8))
    end do
  end if
  return
end subroutine sgcED2
