subroutine fitlm(ncomp,tim,sig,ntim,pars,Stdpars,&
                 value,predtval,transf,errorflag)
!----------------------------------------------------------------------------------------------------------------------
! Subroutine fitlm() is used for fitting OSL signal of type "lm" using Levenberg-Marquadt method,  a series combination
! of initial parameters will be given to call subroutine lmfit() to perform the fitting process.
! =====================================================================================================================
!
! ncomp,             input:: integer, number of components to be decomposed.
!
! ntim,              input:: integer, length of tim or signal
!
! tim(ntim),         input:: real values, time values.
!
! sig(ntim),         input:: real values, signal values.
!
! transf,            input:: logical value, whether transform estimated parameters or not.
!
! pars(2*ncomp),    output:: real values, estimated parameters.
!
! Stdpars(2*ncomp), output:: real values, estimated standard errors of parameters.
!
! value,            output:: real value, minimized objective function value.
!
! predtval(ntim),   output:: real values, predicted signal values that corresponded to signal.
!
! errorfalg,        output:: integer, error message generated during the calculation:
!                            1.1) a successful work given errorflag=0; 
!                            1.2) a totally fail work given errorflag=1.
! =====================================================================================================================
!     Author:: Peng Jun, 2013.07.24, revised in 2013.08.03, revised in 2013.10.05; revised in 2013.12.16.
!
! Dependence:: subroutine comb_next; subroutine targfunc; subroutine lmfit.
! 
! References:: Bluszcz, A., 1996. Exponential function fitting to TL growth data 
!              and similar applications. Geochronometria 13, 135–141.
!              
!              Bulur, E., 2000. A simple transformation for converting CW-OSL 
!              curves to LM-OSL curves. Radiation Measurements 32, 141-145.
!
!              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
!              luminescence components in different quartz samples: implications for dose measurement. Radiation
!              Measurements, 37 (4-5), pp. 441-449.
!---------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                   intent(in)::ncomp
  integer(kind=4),                   intent(in)::ntim 
  real   (kind=8),dimension(ntim),   intent(in)::tim
  real   (kind=8),dimension(ntim),   intent(in)::sig
  logical,                           intent(in)::transf
  real   (kind=8),dimension(ntim),   intent(out)::predtval
  real   (kind=8),dimension(2*ncomp),intent(out)::pars
  real   (kind=8),dimension(2*ncomp),intent(out)::Stdpars
  integer(kind=4),                   intent(out)::errorflag
  real   (kind=8),                   intent(out)::value
  !
  ! Variables for subroutine lmfit()
  real   (kind=8),dimension(2*ncomp)::lmpars,cpars
  real   (kind=8),dimension(2*ncomp)::lmStdpars,cStdpars
  real   (kind=8),dimension(ntim)::lmpredtval,cpredtval
  real   (kind=8),parameter::lmtol=1.0D-07
  real   (kind=8)::lmvalue,cvalue
  integer(kind=4)::lmErr
  ! Variables for subroutine targfunc() and subroutine comb_next()
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
  integer(kind=4), parameter::typ=2
  integer(kind=4):: i
  logical:: saving
  real   (kind=8)::loopvalue,minvalue
  ! 
  saving=.false.
  loopvalue=1.0D+30
  minvalue=1.0D+30
  errorflag=1
  done=.true.
  Loop: do i=1,permdex(ncomp)
    ! Obtain random set
    call comb_next(done,7,ncomp,iarray)
    !
    tlamda=initry(iarray)
    ! Obtain tithn
    call targfunc(tlamda,ncomp,tim,sig,targtol,&
                  typ,ntim,tvalue,tithn,targErr)
    !
    if(targErr==0) then
      lmpars=(/tithn,tlamda/)
    else 
      lmpars(1:ncomp)=100000.0D+00*tlamda
      lmpars(ncomp+1:)=tlamda
    end if
    call lmfit(tim,sig,ntim,lmpars,2*ncomp,typ,transf,&
               lmStdpars,lmpredtval,lmvalue,lmtol,lmErr)
    ! Set pars, stdpars, predtval and value if possible
    if(lmErr==1 .and. lmvalue<loopvalue) then
      pars=lmpars
      Stdpars=lmStdpars
      predtval=lmpredtval
      value=lmvalue
      loopvalue=lmvalue
      ! Successful simple trails given errorflag=0 
      errorflag=0
    end if
    ! If parameters' standard errors can not be estimated, save it too
    if((lmErr==4 .or. lmErr==5 .or. lmErr==6) .and. lmvalue<minvalue) then
      cpars=lmpars
      cStdpars=lmStdpars
      cpredtval=lmpredtval
      cvalue=lmvalue
      minvalue=lmvalue
      ! Set saving to be true
      saving=.true.
    end if
    !??? Exit if at least one successful simple tries has been achieved
    !??? if(errorflag==0) exit Loop
  end do Loop
  ! If errorflag=0, then return values will be areadly 
  ! saved, or else return values need to be resetted
  if(errorflag/=0) then
    if(saving .eqv. .true.) then
      ! Return a solution that Std.pars can not obtained
      pars=cpars
      Stdpars=-99.0D+00
      predtval=cpredtval
      value=cvalue
    else if(saving .eqv. .false.) then
      ! Return a solution that is totally not possible
      ! if simple trials are failures
      pars=    -99.0D+00
      Stdpars= -99.0D+00
      predtval=-99.0D+00
      value=   -99.0D+00
    end if
  end if
  ! Whether transform a[i] to a[i]/b[i] or not ?
  if(transf .eqv. .true.) then
    ! Need not to transform the solution that pars=-99.0D+00 !!!
    if(all(pars>0.0D+00)) pars(1:ncomp)=pars(1:ncomp)/pars(ncomp+1:)
  end if
  return
end subroutine fitlm
