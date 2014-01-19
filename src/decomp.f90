subroutine decomp(ncomp,tim,sig,ntim,pars,Stdpars,value,transf,&
                  predtval,factor,f,cr,maxiter,tol,errorflag)
!-------------------------------------------------------------------------------------------------------------------------------
! Subroutine decomp() is used to decompose CW OSL decay curve.
! ==============================================================================================================================
!
! ncomp,             input:: integer, the number of components wanted to be decomposed.
!
! tim(ntim),         input:: real values, time values.
!
! sig(ntim),         input:: real values, signal values.
!
! ntim,              input:: integer, length of signal values.
!
! transf,            input:: logical value, whether transform estimated values or not.
!
! pars(2*ncomp),    output:: real values, calculated parameters.
!
! Stdpars(2*ncomp), output:: real values, calculated standard errors of parameters.
!
! value,            output:: real value, minimized sum of square of residuals.
!
! predtval(ntim),   output:: real values, fitted values correspond to signal values.
!
! factor,            input:: integer, arugument for scaling np in differential evolution, np=ncomp*factor.
!
! f   ,              input:: real value, arugument of differential evolution.
!
! cr,                input:: real value, arugument of differential evolution.
!
! maxiter,           input:: integer, maximum allowed iteration of differential evolution iteration.
!
! tol,               input:: real value, tolerance for stoping iteration of differential evolution.
!
!
! errorflag(3),     output:: integer values, error message generated during the running:
!                            1.1) if parameters can be initialized using differential evolution, errorflag(1)=0;
!                            1.2) if parameters can not be initialized using differential evolution,errorflag(1)=1;
!
!                            2.1) if parameters have been optimized successfully using Levenberg-Marquadt method, errorflag(2)=0;
!                            2.2) if parameters can not be (or not have been) optimized using Levenberg-Marquadt method, errorflag(2)=1;
!
!                            3.1) if success in simple trails, errorflag(3)=0;
!                            3.2) if fail in simple trails, errorflag(3)=1.
! =====================================================================================================================================
! Author:: Peng Jun, 2013.06.05, revised in 2013.10.05; revised in 2013.12.16; revised in 2014.01.02.
!
! Dependence:: subroutine diffev; subroutine lmfit; subroutine comb_next; subroutine targfunc.
!
! References:: Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
!              applications. Geochronometria 13, 135–141.
!
!              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
!              OSL decay curves. Radiation Measurements 41, 886-891.
!     
!              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
!              luminescence components in different quartz samples: implications for dose measurement. Radiation
!              Measurements, 37 (4-5), pp. 441-449.
!            
!              Jorge More, Burton Garbow, Kenneth Hillstrom,User Guide for MINPACK-1,
!              Technical Report ANL-80-74, Argonne National Laboratory, 1980.   
!---------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                   intent(in)::ncomp
  integer(kind=4),                   intent(in)::ntim 
  integer(kind=4),                   intent(in)::factor
  integer(kind=4),                   intent(in)::maxiter
  real   (kind=8),                   intent(in)::f
  real   (kind=8),                   intent(in)::cr
  real   (kind=8),                   intent(in)::tol
  logical,                           intent(in)::transf
  real   (kind=8),dimension(ntim),   intent(in)::tim
  real   (kind=8),dimension(ntim),   intent(in)::sig
  real   (kind=8),dimension(ntim),   intent(out)::predtval
  real   (kind=8),dimension(2*ncomp),intent(out)::pars
  real   (kind=8),dimension(2*ncomp),intent(out)::Stdpars
  integer(kind=4),dimension(3),      intent(out)::errorflag
  real   (kind=8),                   intent(out)::value
  !
  ! Variables for subroutine diffev
  real   (kind=8),dimension(factor*ncomp,ncomp)::agents
  real   (kind=8),dimension(ncomp)::dlamda,dithn
  real   (kind=8),dimension(ntim)::dpredict
  integer(kind=4),dimension(2)::dErr
  real   (kind=8)::dvalue
  !
  ! Variables for subroutine lmfit
  real   (kind=8),dimension(2*ncomp)::lmpars
  real   (kind=8),dimension(2*ncomp)::lmStdpars
  real   (kind=8),dimension(ntim)::lmpredtval
  real   (kind=8),parameter::lmtol=1.0D-07
  real   (kind=8)::lmvalue,lmcond
  integer(kind=4)::lmErr
  !
  ! Variables for subroutine targfunc and subroutine comb_next
  logical:: done
  real   (kind=8),parameter::targtol=1.0D-09
  integer(kind=4)::targErr
  real   (kind=8)::tvalue
  integer(kind=4),dimension(ncomp)::iarray
  real   (kind=8),dimension(ncomp)::tlamda,tithn
  integer(kind=4),dimension(7),parameter::permdex=(/7,21,35,35,21,7,1/)
  real   (kind=8),dimension(7),parameter::initry=(/32.0D+00,  &      ! UF
                                                   2.5D+00 ,  &      ! Fast 
                                                   0.62D+00,  &      ! Medium 
                                                   0.15D+00,  &      ! Slow1 
                                                   0.023D+00, &      ! Slow2 
                                                   0.0022D+00,&      ! Slow3 
                                                   0.0003D+00 /)     ! Slow4
  ! Local variables
  integer(kind=4)::i
  integer(kind=4),parameter::typ=1
  real   (kind=8),dimension(2*ncomp)::cpars,cStdpars
  real   (kind=8),dimension(ntim)::cpredtval
  real   (kind=8)::cvalue
  logical        ::saving
  real   (kind=8)::loopcond,mincond
  !
  ! A saving mark
  saving=.false.
  ! Default output values, i.e. if wantted output can not be
  ! updating by some of: 
  ! 1) differential evolution; 
  ! 2) Levenberg-Marquadt; 
  ! 3) simple trials. 
  ! Then all output will be -99.0D+00
  pars=-99.0D+00
  Stdpars=-99.0D+00
  predtval=-99.0D+00
  value=-99.0D+00
  ! 
  ! Call diffev() to initialize parameters
  call diffev(ncomp,tim,sig,ntim,dlamda,dithn,dvalue,dpredict,&
              agents,factor*ncomp,f,cr,maxiter,tol,dErr)
  !
  ! Update outputted values if possible
  if(dErr(2)==0) then
    pars=(/dithn,dlamda/)
    ! Stdpars can not be estimated by differential 
    ! evolution, so need not update here
    ! Stdpars=-99.0D+00
    predtval=dpredict
    value=dvalue 
    ! If diffev successed, set errorflag(1)=0
    errorflag(1)=0
  else 
    ! If diffev failed, set errorflag(1)=1
    errorflag(1)=1
  end if 
  !
  errorflag(2)=1
  ! Call lmfit() if diffev() successed
  if(errorflag(1)==0) then
    lmpars=(/dithn,dlamda/)
    call lmfit(tim,sig,ntim,lmpars,2*ncomp,typ,transf,&
               lmcond,lmStdpars,lmpredtval,lmvalue,lmtol,lmErr)
    ! Update output [pars, stdpars, predtval and value] if possible
    if(lmErr==1) then
      pars=lmpars
      Stdpars=lmStdpars
      predtval=lmpredtval
      value=lmvalue
      ! If lmfit() successed, set errorflag(2)=0
      errorflag(2)=0
    else if(lmErr==4 .or. lmErr==5 .or. lmErr==6) then
      ! Save those values that from which parameters' std.errors can not be estimated, too
      cpars=lmpars
      cStdpars=lmStdpars
      cpredtval=lmpredtval
      cvalue=lmvalue
      ! Set saving mark to be true
      saving=.true.
    end if
    !
  end if
  !
  errorflag(3)=1
  ! Set loopcond
  if(errorflag(2)==0) then
    loopcond=lmcond
  else 
    loopcond=1.0D+30
  end if
  mincond=1.0D+30
  !
  ! Call simple trials to try to further improve the optimization
  done=.true.
  Loop: do i=1,permdex(ncomp)
    ! Obtain random set of decay rates
    call comb_next(done,7,ncomp,iarray)
    !
    tlamda=initry(iarray)
    ! Obtain tithn
    call targfunc(tlamda,ncomp,tim,sig,targtol,typ,&
                  ntim,tvalue,tithn,targErr)
    if(targErr==0) then
      lmpars=(/tithn,tlamda/)
    else 
      lmpars(1:ncomp)=1.0D+05*tlamda
      lmpars(ncomp+1:)=tlamda
    end if
    !
    call lmfit(tim,sig,ntim,lmpars,2*ncomp,typ,transf,&
               lmcond,lmStdpars,lmpredtval,lmvalue,lmtol,lmErr)
    ! Update output [pars, stdpars, predtval and value] if possible
    if(lmErr==1 .and. lmcond<loopcond) then
      pars=lmpars
      Stdpars=lmStdpars
      predtval=lmpredtval
      value=lmvalue
      loopcond=lmcond
      ! If its a successful simple trail, set errorflag(3)=0
      errorflag(3)=0
    else if ((lmErr==4 .or. lmErr==5 .or. lmErr==6) .and. lmcond<mincond) then
      ! Save those values that from which pars' std.errors can not be estimated, too
      cpars=lmpars
      cStdpars=lmStdpars
      cpredtval=lmpredtval
      cvalue=lmvalue
      mincond=lmcond
      ! Set saving to be true
      saving=.true.
    end if
    !
  end do Loop
  !
  if(saving .eqv. .true.) then
    if(errorflag(2)/=0 .and. errorflag(3)/=0) then
      pars=cpars
      ! Stdpars are not avialble here
      Stdpars=-99.0D+00
      predtval=cpredtval
      value=cvalue
    end if
  end if
  !
  ! If transf=true, then ithn (a[1],a[2],...) need to be resetted
  ! i.e. transform a[i] to a[i]/b[i]
  if(transf .eqv. .true.) then
    ! Do not transform the solution that pars=-99.9D+00 !!!
    if(all(pars>0.0D+00)) pars(1:ncomp)=pars(1:ncomp)/pars(ncomp+1:)
  end if
  !
  return
end subroutine decomp
