! ***************************************************************************************************************************
! FILE "decomp_Constant" consists a series of subroutines rewritten from "SRC" file of R package RadialPlotter(version 2.2).
! Those subroutines are used for decompose decay curve plus constant subtracted.
! 1).   decomp_C(): decompose CW-OSL decay curves.
! 2).   deffev_C(): initilize parameters in with differential evolution method.
! 3).    fitlm_C(): decompose LM-OSL decay curves.
! 4).    lmfit_C(): optimize parameters with Levenberg-Marquadt method.
! 5).  lmfunc1_C(): a subroutine for lmfit_C() (LM-OSL).
! 6).   lmfunc_C(): a subroutine for lmfit_C() (CW-OSL).
! 7).   lmhess_C(): estimate parameters' standard errors with difference-approximation method.
! 8). targfunc_C(): calculate magnitudes from decay rates with Linear Algebra method.
! Edited in 2013.12.16. Peng Jun.
! ****************************************************************************************************************************
!
subroutine decomp_C(ncomp,tim,sig,ntim,pars,Stdpars,value,transf,&
                    predtval,factor,f,cr,maxiter,tol,errorflag)
!-------------------------------------------------------------------------------------------------------------------------------
! Subroutine decomp_C() is used to decompose CW OSL decay curve (plus a constant).
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
! pars(2*ncomp+1),   output:: real values, calculated parameters.
!
! Stdpars(2*ncomp+1),output:: real values, calculated standard errors of parameters.
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
!                            3.2) if simple trails have not been performed (or fail in simple trials), errorflag(3)=1.
! =====================================================================================================================================
! Author:: Peng Jun, 2013.12.15.
!
! Dependence:: subroutine diffev_C; subroutine lmfit_C; 
!              subroutine comb_next; subroutine targfunc_C.
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
  real   (kind=8),dimension(2*ncomp+1),intent(out)::pars
  real   (kind=8),dimension(2*ncomp+1),intent(out)::Stdpars
  integer(kind=4),dimension(3),      intent(out)::errorflag
  real   (kind=8),                   intent(out)::value
  !
  ! Variables for subroutine diffev_C
  real   (kind=8),dimension(factor*ncomp,ncomp)::agents
  real   (kind=8),dimension(ncomp)::dlamda,dithn
  real   (kind=8)                 ::dconstant
  real   (kind=8),dimension(ntim)::dpredict
  integer(kind=4),dimension(2)::dErr
  real   (kind=8)::dvalue
  !
  ! Variables for subroutine lmfit_C
  real   (kind=8),dimension(2*ncomp+1)::lmpars
  real   (kind=8),dimension(2*ncomp+1)::lmStdpars
  real   (kind=8),dimension(ntim)::lmpredtval
  real   (kind=8),parameter::lmtol=1.0D-07
  real   (kind=8)::lmvalue
  integer(kind=4)::lmErr
  !
  ! Variables for subroutine targfunc_C and subroutine comb_next
  logical:: done
  real   (kind=8),parameter::targtol=1.0D-09
  integer(kind=4)::targErr
  real   (kind=8)::tvalue
  integer(kind=4),dimension(ncomp)::iarray
  real   (kind=8),dimension(ncomp)::tlamda,tithn
  real   (kind=8)                 ::tconstant
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
  real   (kind=8),dimension(2*ncomp+1)::cpars,cStdpars
  real   (kind=8),dimension(ntim)::cpredtval
  real   (kind=8)::cvalue
  logical        ::saving
  real   (kind=8)::loopvalue,minvalue
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
  ! Call diffev_C() to initialize parameters
  call diffev_C(ncomp,tim,sig,ntim,dlamda,dithn,dconstant,dvalue,&
                dpredict,agents,factor*ncomp,f,cr,maxiter,tol,dErr)
  !
  ! Update outputted values if possible
  if(dErr(2)==0) then
    pars=(/dithn,dlamda,dconstant/)
    ! Stdpars can not be estimated by differential 
    ! evolution, so need not update here
    ! Stdpars=-99.0D+00
    predtval=dpredict
    value=dvalue 
    ! If diffev_C successed, set errorflag(1)=0
    errorflag(1)=0
  else 
    ! If diffev_C failed, set errorflag(1)=1
    errorflag(1)=1
  end if 
  !
  errorflag(2)=1
  ! Call lmfit()_C if diffev()_C successed
  if(errorflag(1)==0) then
    lmpars=(/dithn,dlamda,dconstant/)
    call lmfit_C(tim,sig,ntim,lmpars,2*ncomp+1,typ,transf,&
                 lmStdpars,lmpredtval,lmvalue,lmtol,lmErr)
    ! Update output [pars, stdpars, predtval and value] if possible
    if(lmErr==1) then
      pars=lmpars
      Stdpars=lmStdpars
      predtval=lmpredtval
      value=lmvalue
      ! If lmfit_C() successed, set errorflag(2)=0
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
  ! If diffev()_C or lmfit()_C fails, then simple trials will be conducted
  if(errorflag(1)/=0 .or. errorflag(2)/=0) then
    loopvalue=1.0D+30
    minvalue=1.0D+30
    done=.true.
    Loop: do i=1,permdex(ncomp)
      ! Obtain random set of decay rates
      call comb_next(done,7,ncomp,iarray)
      !
      tlamda=initry(iarray)
      ! Obtain tithn
      call targfunc_C(tlamda,ncomp,tim,sig,targtol,typ,&
                      ntim,tvalue,tithn,tconstant,targErr)
      if(targErr==0) then
        lmpars=(/tithn,tlamda,tconstant/)
      else 
        lmpars(1:ncomp)=10000.0D+00*tlamda
        lmpars(ncomp+1:2*ncomp)=tlamda
        lmpars(2*ncomp+1)=200.0
      end if
      !
      call lmfit_C(tim,sig,ntim,lmpars,2*ncomp+1,typ,transf,&
                   lmStdpars,lmpredtval,lmvalue,lmtol,lmErr)
      ! Update output [pars, stdpars, predtval and value] if possible
      if(lmErr==1 .and. lmvalue<loopvalue) then
        pars=lmpars
        Stdpars=lmStdpars
        predtval=lmpredtval
        value=lmvalue
        loopvalue=lmvalue
        ! If its a successful simple trail, set errorflag(3)=0
        errorflag(3)=0
      else if ((lmErr==4 .or. lmErr==5 .or. lmErr==6) .and. lmvalue<minvalue) then
        ! Save those values that from which pars' std.errors can not be estimated, too
        cpars=lmpars
        cStdpars=lmStdpars
        cpredtval=lmpredtval
        cvalue=lmvalue
        minvalue=lmvalue
        ! Set saving to be true
        saving=.true.
      end if
      !??? Exit if at least one successful simple trail has been achieved
      !??? if(errorflag(3)==0) exit Loop
      !
    end do Loop
    !
  end if
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
    if(all(pars>0.0D+00)) pars(1:ncomp)=pars(1:ncomp)/pars(ncomp+1:2*ncomp)
  end if
  !
  return
end subroutine decomp_C

subroutine diffev_C(npars,tim,sig,ntim,lamda,ithn,constant,value,&
                    predict,agents,np,f,cr,maxiter,tol,errorflag)
!---------------------------------------------------------------------------------------------------------------------------------------
! Fitting OSL signal decay curve using differential evolution algorithm, I=a1*exp(-b1*t)+a2*exp(-b2*t)+...+
!                                                                          ak*exp(-bk*t)+c will be fitted.
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
! constant,            output:: real value, estimated constant.
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
! Author:: Peng Jun, 2013.12.14.
! 
! Dependence:: subroutine leaveone; subroutine NonReplaceSample; 
!              subroutine targfunc_C; subroutine sort; subroutine sd.
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
  real   (kind=8),                    intent(out)::constant
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
  real   (kind=8)                    ::yconstant,agentsconstant
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
                                                             1.0D-06, 0.05D+00/), (/2,7/) )  ! the whole vector space 
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
        call targfunc_C(y,npars,tim,sig,targtol,typ,ntim,&
                        yvalue,yithn,yconstant,errorflagy)
        !
        ! calculate correspond ithn values, lamda values, using the ith agent values
        call targfunc_C(agents(j,:),npars,tim,sig,targtol,typ,ntim,&
                        agentsvalue,agentsithn,agentsconstant,errorflagagents)
        !
        ! replace the ith agent by y if:
        ! 1) function value can be successfully calculated using Linear Algebra
        !    method and all ithn values for y are larger than 0; 
        ! 2) function value of y not less than function value for the ith agent
        ! ***
        ! At this points also note that we don't care if error appears when calling
        ! targfunc_C using the jth agents value, if so, agentsvalue will be 1.0D+30,
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
    call targfunc_C(agents(MinIndex,:),npars,tim,sig,targtol,typ,ntim,&
                    agentsvalue,agentsithn,agentsconstant,errorflagagents)
    value=agentsvalue
    ithn=agentsithn
    lamda=agents(MinIndex,:)
    constant=agentsconstant
    ! set predict
    predict=constant
    do i=1,npars
      predict=predict+ithn(i)*dexp(-lamda(i)*tim)
    end do
  else 
    value=-99.0D+00
    ithn=-99.0D+00
    lamda=-99.0D+00
    constant=-99.0D+00
    predict=-99.0D+00
  end if
  !
  return
end subroutine diffev_C

subroutine fitlm_C(ncomp,tim,sig,ntim,pars,Stdpars,&
                   value,predtval,transf,errorflag)
!----------------------------------------------------------------------------------------------------------------------
! Subroutine fitlm_C() is used for fitting OSL signal of type "lm" (plus a constant) using Levenberg-Marquadt method,  
! a series combination of initial parameters will be given to call subroutine lmfit_C() to perform the fitting process.
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
! pars(2*ncomp+1),  output:: real values, estimated parameters.
!
! Stdpars(2*ncomp+1),output:: real values, estimated standard errors of parameters.
!
! value,            output:: real value, minimized objective function value.
!
! predtval(ntim),   output:: real values, predicted signal values that corresponded to signal.
!
! errorfalg,        output:: integer, error message generated during the calculation:
!                            1.1) a successful work given errorflag=0; 
!                            1.2) a totally fail work given errorflag=1.
! =====================================================================================================================
!     Author:: Peng Jun, 2013.07.24, revised in 2013.08.03, revised in 2013.10.05.
!
! Dependence:: subroutine comb_next; subroutine targfunc_C; subroutine lmfit_C.
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
  real   (kind=8),dimension(2*ncomp+1),intent(out)::pars
  real   (kind=8),dimension(2*ncomp+1),intent(out)::Stdpars
  integer(kind=4),                   intent(out)::errorflag
  real   (kind=8),                   intent(out)::value
  !
  ! Variables for subroutine lmfit_C()
  real   (kind=8),dimension(2*ncomp+1)::lmpars,cpars
  real   (kind=8),dimension(2*ncomp+1)::lmStdpars,cStdpars
  real   (kind=8),dimension(ntim)::lmpredtval,cpredtval
  real   (kind=8),parameter::lmtol=1.0D-07
  real   (kind=8)::lmvalue,cvalue
  integer(kind=4)::lmErr
  ! Variables for subroutine targfunc_C() and subroutine comb_next()
  logical:: done
  real   (kind=8),parameter::targtol=1.0D-09
  integer(kind=4)::targErr
  real   (kind=8)::tvalue
  integer(kind=4),dimension(ncomp)::iarray
  real   (kind=8),dimension(ncomp)::tlamda,tithn
  real   (kind=8)                 ::tconstant
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
    call targfunc_C(tlamda,ncomp,tim,sig,targtol,&
                    typ,ntim,tvalue,tithn,tconstant,targErr)
    !
    if(targErr==0) then
      lmpars=(/tithn,tlamda,tconstant/)
    else 
      lmpars(1:ncomp)=100000.0D+00*tlamda
      lmpars(ncomp+1:2*ncomp)=tlamda
      lmpars(2*ncomp+1)=1500.0D+00
    end if
    call lmfit_C(tim,sig,ntim,lmpars,2*ncomp+1,typ,transf,&
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
    if(all(pars>0.0D+00)) pars(1:ncomp)=pars(1:ncomp)/pars(ncomp+1:2*ncomp)
  end if
  return
end subroutine fitlm_C

subroutine lmfit_C(xdat,ydat,ndat,pars,npars,typ,transf,&
                   stderror,predtval,value,tol,info)
!--------------------------------------------------------------------------------------------------------------------------------
! Specifying a fitting model, lmfit_C() will fit the model using Levenberg-Marquadt method.
! typ=1 has the formula I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t)+c, where k=1:7;
! typ=2 has the formula I(t)=a1*(t/max(t))*exp(-b1*t^2)/2/max(t))+...+
!                            ak*(t/max(t))*exp(-bk*t^2)/2/max(t))+c*(t/max(t)), where k=1:7.
! ===============================================================================================================================
!
! ndat,                input:: integer, length of xdat or y dat.
!
! npars,               input:: integer, dimension of parameters (or length of pars).
!
! typ,                 input:: integer, 1 for fitting 'cw', 2 for fitting 'lm'.
!
! transf,              output:: logical, whether transform se(a[i]) to se(a[i]/b[i]) or not.
!
! xdat(ndat),          input:: real values, xdat (or independent variables x).
!
! ydat(ndat),          input:: real values, ydat (or dependent variables y).
!
! pars(npars),  input/output:: real values, initial guess values for parameters to be optimized,
!                              overwritten to be final results for output.
!
! stderror(npars),    output:: real values, estimated standard errors for parameters.
!
! predtval(ndat),     output:: real values, fited values correspond to ydat.
!
! value,              output:: real value, final total residual error.
!
! tol,                 input:: real value,  termination occurs when the algorithm estimates either that the relative error in the  
!                              sum of squares is at most TOL or that the relative error between X and the solution is at most TOL.
!
! info,               output:: error message generated during the calculation:
!                              1.1) successful work, info=1;
!                              1.2) error when calling lmder1, info=2;
!                              1.3) at least one estimated parameter is below zero, info=3;
!                              1.4) error when calling lmhess to approximate model's hessian matrix, info=4;
!                              1.5) error when attempt to inverse the approximated hessian matrix, info=5;
!                              1.6) if diagnal elements of inversed hessian matrix has at least one minus value, info=6.
! ================================================================================================================================
! Author:: Peng Jun, 2013.12.15.
!
! Dependence:: subroutine lmfunc_C; subroutine lmfunc1_C; subroutine lmder1; 
!              subroutine lmhess_C; subroutine inverse; subroutine diag.
!---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ndat
  integer(kind=4),intent(in)::npars
  integer(kind=4),intent(in)::typ
  real   (kind=8),dimension(ndat),intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat),intent(out)::predtval
  real   (kind=8),intent(in)::tol
  logical,        intent(in)::transf
  real   (kind=8),intent(out)::value
  integer(kind=4),intent(out)::info
  !
  ! Variables for subroutine lmhess_C
  real   (kind=8),parameter::lmtol=2.220446D-16   ! used for singular matrix diagnose in subroutine (.Machine$double.eps in R)
                                                  ! lmhess and subroutine inverse
  real   (kind=8),parameter::minAbsPar=0.0D+00    ! for lmhess use
  real   (kind=8),dimension(npars,npars)::hessian ! hessian matrix by finite-difference approximation
  real   (kind=8),dimension(npars)::gradient      ! gradient by finite-difference approximation
  integer(kind=4),dimension(5)::hesserror         ! error message generated in lmhess
  ! Variables for subroutine inverse
  integer(kind=4)::inverror                       ! for singular matrix diagnose in subroutine inverse
  ! Local variables
  integer(kind=4)::ldfjac,iflag
  real   (kind=8),dimension(ndat)::fvec           ! fitted residual in vector form
  real   (kind=8),dimension(ndat,npars)::fjac     ! jacobian matrix
  integer(kind=4)::i
  !
  ldfjac=ndat
  !
  ! Default returned stderrors if error appears
  stderror=-99.0D+00
  ! Specify iflag to be 1 to calculate fvec
  iflag=1
  !
  ! Using initial pars to caculate fvec, which will be used in subroutine lmder1
  if(typ==1) then
    ! For 'cw'
    call lmfunc_C(xdat,ydat,ndat,npars,pars,fvec,fjac,ldfjac,iflag)
  else if(typ==2) then
    ! For 'lm'
    call lmfunc1_C(xdat,ydat,ndat,npars,pars,fvec,fjac,ldfjac,iflag)
  end if
  !
  ! Optimizing initial pars using Levenberg-Marquadt method
  ! and return pars and info for output
  if(typ==1) then
    ! For 'cw'
    call lmder1(lmfunc_C,ndat,npars,pars,fvec,&
                fjac,ldfjac,tol,info,xdat,ydat)
  else if(typ==2) then
    ! For 'lm'
    call lmder1(lmfunc1_C,ndat,npars,pars,fvec,&
                fjac,ldfjac,tol,info,xdat,ydat)
  end if
  !
  ! Calculate sum of squre of residual 
  value=sum((fvec)**2)
  ! Calculate fitted values correspond to ydat
  predtval=fvec+ydat
  !
  ! Check if any error appears when calling lmder1, and reset info 
  ! to 2 and return if so, else reset info to 1 and continue
  if(info==1 .or. info==2 .or. info==3) then
    info=1
  else 
    info=2
    return
  end if  
  !
  ! Check if any estimated parameter is below zero,
  ! if so, reset info to be 3 and return
  if(any(pars<0.0D+00)) then
    info=3
    return
  end if
  ! Estimate pars' standard errors with subroutine lmhess_C
  if(typ==1) then
    call lmhess_C(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
                  hessian,gradient,value,hesserror,1)
  else if(typ==2) then
    call lmhess_C(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
                  hessian,gradient,value,hesserror,2)
  end if
  ! Reset value after calling subroutine lmhess_C
  value=value**2
  ! Check if any error appears when calling lmhess_C, 
  ! if so, reset info to be 4 and return
  if(any(hesserror/=0))  then
    info=4
    return
  end if
  !
  ! Set hessian to be its inverse form
  call inverse(hessian,npars,inverror,lmtol)
  ! Check if any error appears when calling inverse, 
  ! if so, reset info to be 5 and return
  if(inverror==1)  then
    info=5
    return
  end if
  ! Extract diagnal elements from inversed hessian 
  ! matrix and save it in arrary stderror
  call diag(hessian,npars,stderror)
  ! check if any elements in stderror is 
  ! below zero, if so, set info to be 6 and return
  if(any(stderror<0)) then
    info=6
    return
  end if
  ! Rest stderror
  stderror=dsqrt(stderror)
  !
  ! transform se(a[i]) to se(a[i]/b[i]) or not?
  if(transf .eqv. .true.) then
    ! If transf=true, then std.err of ithn (a1,a2,...) need to be resetted
    do i=1,(npars-1)/2
      stderror(i)=pars(i)/pars(i+(npars-1)/2)*dsqrt(hessian(i,i)/pars(i)**2+&
                  hessian(i+(npars-1)/2,i+(npars-1)/2)/pars(i+(npars-1)/2)**2-&
                   2.0D+00*hessian(i,i+(npars-1)/2)/pars(i)/pars(i+(npars-1)/2))
    end do
  end if
  !
  return
end subroutine lmfit_C

subroutine lmfunc1_C(xdat,ydat,m,n,x, &
                     fvec,fjac,ldfjac,iflag)
! ----------------------------------------------------------------------------------------------------------------------------
! For formula I(t)=a1*(t/max(t))*exp(-b1*t^2)/2/max(t))+...+ak*(t/max(t))*exp(-bk*t^2)/2/max(t))+c*(t/max(t)), where k=1:7
! lmfunc1_C() is used to calculate vectors fevc(i)=I(t,i)-ydat(i), i=1:length(ydat), or the jacobian matrix fjac, 
! both them will be passed to subroutine lmfit_C to perform the Levenberg-Marquadt optimization.
! ===========================================================================================================================
!
! m,                   input:: integer, length of xdat or ydat.
!
! n,                   input:: integer, dimension of the problem (length of x).
!
! ldfjac,              input:: integer,the leading dimension of FJAC, which must be no less than m.
!
! xdat(m),             input:: real values, times values.
!
! ydat(m),             input:: real values, signal values.
!
! x(n),                input:: real values, initial guess values.
!
! fvec(m),            output:: real values, differences between predicted values and obeserved values.
!
! fjac(ldfjac,n),     output:: real values, the jacobian matrix
!
! iflag,               input:: integer, If IFLAG = 1 on intput, FCN should calculate the functions at X and return 
!                              this vector in FVEC. If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!                              return this matrix in FJAC. To terminate the algorithm, the user may set IFLAG negative.
! ===========================================================================================================================
! Author:: Peng Jun, 2013.12.15.
!
! Dependence:: No
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::m
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::ldfjac
  integer(kind=4),intent(in)::iflag
  real   (kind=8),dimension(m),intent(in)::xdat,ydat
  real   (kind=8),dimension(n),intent(in)::x
  real   (kind=8),dimension(m),intent(out)::fvec
  real   (kind=8),dimension(ldfjac,n),intent(out)::fjac
  ! Local variables
  integer(kind=4)::i
  real   (kind=8)::maxx
  !
  ! Calculate values for targeted function and 
  ! store them in vector fvec
  maxx=maxval(xdat)
  if(iflag==1) then
    fvec=x(n)*xdat/maxx
    do i=1,(n-1)/2
      fvec=fvec+x(i)*(xdat/maxx)*&
           dexp(-x(i+(n-1)/2)*xdat**2/2.0D+00/maxx)
    end do
    fvec=fvec-ydat
  ! Calculate matrix jacobian in a column by column order
  else if(iflag==2)  then
    do i=1,(n-1)/2
      fjac(:,i)=(xdat/maxx)*dexp(-x(i+(n-1)/2)*xdat**2/2.0D+00/maxx)
    end do
    do i=(n-1)/2+1,n-1
      fjac(:,i)=-xdat**2/2.0D+00/maxx*&
                dexp(-x(i)*xdat**2/2.0D+00/maxx)*&
                x(i-(n-1)/2)*(xdat/maxx)
    end do
    fjac(:,n)=xdat/maxx
  end if
  ! now return
  return
end subroutine lmfunc1_C  

subroutine lmfunc_C(xdat,ydat,m,n,x,fvec,fjac,ldfjac,iflag)
! ---------------------------------------------------------------------------------------------------------------------------
! For formula I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t)+c, K=1:7,
! lmfunc_C is used to calculate vectors fevc(i)=I(t,i)-ydat(i), i=1:length(ydat), 
! so is the jacobian matrix fjac, they will be passed to subroutine lmfit_C
! to perform the Levenberg-Marquadt optimization.
! ===========================================================================================================================
!
! m,                   input:: integer, length of xdat or ydat.
!
! n,                   input:: integer, dimension of the problem (length of x).
!
! ldfjac,              input:: integer,the leading dimension of FJAC, which must be no less than m.
!
! xdat(m),             input:: real values, times values.
!
! ydat(m),             input:: real values, signal values.
!
! x(n),                input:: real values, initial guess values.
!
! fvec(m),            output:: real values, differences between predicted values and obeserved values.
!
! fjac(ldfjac,n),     output:: real values, the jacobian matrix.
!
! iflag,               input:: integer, If IFLAG = 1 on intput, FCN should calculate the functions at X and return 
!                              this vector in FVEC. If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!                              return this matrix in FJAC. To terminate the algorithm, the user may set IFLAG negative.
! ===========================================================================================================================
! Author:: Peng Jun, 2013.12.15.
!
! Dependence:: No
!----------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::m
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::ldfjac
  integer(kind=4),intent(in)::iflag
  real   (kind=8),dimension(m),intent(in)::xdat,ydat
  real   (kind=8),dimension(n),intent(in)::x
  real   (kind=8),dimension(m),intent(out)::fvec
  real   (kind=8),dimension(ldfjac,n),intent(out)::fjac
  ! local variables
  integer(kind=4)::i
  !
  ! calculate values for targeted function and 
  ! store them in vector fvec
  if(iflag==1) then
    fvec=x(n)
    do i=1,(n-1)/2
      fvec=fvec+x(i)*dexp(-x(i+(n-1)/2)*xdat)
    end do
    fvec=fvec-ydat
  ! calculate matrix jacobian in a column by column order
  else if(iflag==2)  then
    do i=1,(n-1)/2
      fjac(:,i)=dexp(-x(i+(n-1)/2)*xdat)
    end do
    do i=(n-1)/2+1,n-1
      fjac(:,i)=-x(i-(n-1)/2)*xdat*dexp(-x(i)*xdat)
    end do
    fjac(:,n)=1.0D+00
  end if
  !
  return
end subroutine lmfunc_C  

subroutine lmhess_C(pars,xdat,ydat,npars,ndat,tol,minAbsPar,&
                    hessian,gradient,value,errorflag,model)
!------------------------------------------------------------------------------------------------------
! Subroutine lmhess_C() is used to calculate the gradient, hessian matrix of a 
! given function contained in its inner, that is fun34 in decay curve fitting, 
! using finite-difference approximation (cw or lm decay curve).
! ====================================================================================================
!
! pars(npars)          :: input, real values, the parameters of the function.
!
! xdat (ndat)          :: input, real values, the xdat values.
!
! ydat (ndat)          :: input, real values, the ydat values.
!
! npars                :: input, integer, the length of tim ( or signal).
!
! ndat                 :: input, integer, the dimension of the parameters.
!
! tol                  :: input, real value, tolerance value for diagnosing sigular matrix.
!
! minAbspar            :: input, real value, the allowed minimum absolute parameter.
!
! hessian(npars,npars) :: output, real values, the hessian matrix.
!
! gradient(npars)      :: output, real values, the gradient of the parameters.
!
! value                :: output, real value, the correspond function value for the specified parameters.
! 
! errorflag(5)         :: output, integer values, error message generated during the calling:
!                         1) if subroutine GJordan is called sucessfully, errorflag(1)=0, otherwise 1;
!                         2) if function value can be calculated,         errorflag(2)=0, otherwise 1;
!                         3) if gradient can be calculated,               errorflag(3)=0, otherwise 1; 
!                         4) if hessian can be calculated,                errorflag(4)=0, otherwise 1; 
!                         5) if no error appears in arrary allocation,    errorflag(5)=0, otherwise 1.
!
! model                :: input, integer, a model to be used for approximation:
!                         1) y=I(1)*exp(-lamda(1)*x)+
!                              I(2)*exp(-lamda(2)*x)+...+
!                              I(k)*exp(-lamda(k)*x)+c,  fitting a "cw" signal curve
!
!                         2) y=a1*(x/max(x))*exp(-b1*x^2)/2/max(x))+
!                              a2*(x/max(x))*exp(-b2*x^2)/2/max(x))+...+
!                              ak*(x/max(x))*exp(-bk*x^2)/2/max(x))+c*(x/max(x)), fitting a "lm" signal curve
! ======================================================================================================
! Dependence:: subroutine GJordan, inter function fun34.
!
! Author:: Peng Jun, 2013.12.15.
!
! Reference:  Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the
!             R Development Core Team (2013). nlme: Linear and Nonlinear Mixed
!             Effects Models. R package version 3.1-108.
!-------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::npars                 
  integer(kind=4),intent(in)::ndat 
  integer(kind=4),intent(in)::model                   
  integer(kind=4),intent(out)::errorflag(5)             
  real   (kind=8),intent(in)::pars(npars)            
  real   (kind=8),intent(in)::xdat(ndat)                
  real   (kind=8),intent(in)::ydat(ndat)            
  real   (kind=8),intent(in)::tol                   
  real   (kind=8),intent(in)::minAbsPar              
  real   (kind=8),intent(out)::gradient(npars)       
  real   (kind=8),intent(out)::hessian(npars,npars)  
  real   (kind=8),intent(out)::value                
  ! local variables
  integer(kind=4)::i,j,p
  integer(kind=4)::ncols
  integer(kind=4)::solerror
  real  (kind=8)::incr(npars)
  real  (kind=8)::diagpar(npars,npars)
  real  (kind=8),allocatable::frac(:),cfrac(:)
  real  (kind=8)::ffrac(1+2*npars)
  real  (kind=8),allocatable::cols(:,:),ccols(:,:),shifted(:,:),transcols(:,:)
  real  (kind=8)::pcols(npars,2*npars+1)
  real  (kind=8),allocatable::xcols(:,:),cxcols(:,:)
  real  (kind=8),allocatable::pxcols(:,:)
  real  (kind=8),parameter::eps=6.055454e-06       
  real  (kind=8)::maxx           
  !
  ! Decide the incr values
  ! for each initial par in pars, check if it
  ! is smaller than minabspar, the incr will
  ! be decided through this check
  do i=1,npars
    if(dabs(pars(i))<=minAbsPar)  then
      incr(i)=minAbsPar*eps
    else
      incr(i)=dabs(pars(i))*eps
    end if
  end do
  !
  ! build a diagal matirx diagpar
  diagpar=0.0D+00
  do i=1,npars
    diagpar(i,i)=1.0D+00
  end do
  ! creat ffrac with length of 2*npars+1
  ! and specify values for it
  ffrac(1)=1.0D+00
  ffrac(2:npars+1)=incr
  ffrac(npars+2:2*npars+1)=incr**2
  ! create pcols to be npars rows
  ! and 2*npars columns, then storing
  ! diagpar in it
  pcols(:,1)=0.0D+00
  pcols(:,2:npars+1)=diagpar
  pcols(:,npars+2:2*npars+1)=-diagpar
  !
  ! default return values if error appears
  value=-99.0D+00
  gradient=-99.0D+00
  hessian=-99.0D+00
  errorflag=0
  !
  ! in a total of (npars-1) times looping
  ! add  (npars-i) columns to cols in each loop number
  ! hence a toal of npars*(npars-1)/2 columns for cols
  allocate( cols(1:npars,1:npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! in a total of (npars-1) times looping
  ! add (npars-i) number to frac in each loop number
  ! hence a total of npars*(npars-1)/2 values for frac
  allocate( frac(1:npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  !
  ncols=0
  do i=1,npars-1
    ! in each loop(i), generate ccols and cfrac with 
    ! npars rows (npars-i) columns 
    allocate( ccols(1:npars,1:npars-i), stat=errorflag(5))
    if(errorflag(5)/=0) return
    allocate( cfrac(1:npars-i), stat=errorflag(5))
    if(errorflag(5)/=0) return
    ! fill ccols and cfrac in each loop
    do j=i+1,npars
      ccols(:,j-i)=diagpar(:,i)+diagpar(:,j)
      cfrac(j-i)=incr(i)*incr(j)       
    end do 
    ! now store ccols and cfrac to cols
    ! and frac respectively
    cols(:,ncols+1:ncols+npars-i)=ccols
    frac(ncols+1:ncols+npars-i)=cfrac
    ! delocate ccols and cfrac, preparing
    ! for new allocations
    deallocate(ccols)
    deallocate(cfrac)
    ! update started filling index ncols
    ncols=ncols+npars-i 
    !
  end do
  ! now add up pcols and cols togeteher to be ccols
  ! ccols has a total of ( 2*npars+1 + npars*(npars-1)/2 ) columns
  allocate( ccols(1:npars,1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  ccols(:,1:2*npars+1)=pcols
  ccols(:,2*npars+2:npars*(npars-1)/2+2*npars+1)=cols
  ! not need cols presently, so release it
  deallocate(cols)
  ! allocate cols again to store ccols, and release ccols
  allocate(cols(1:npars,1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  cols=ccols
  deallocate(ccols)
  ! now add up ffrac and frac together to be cfrac,
  ! cfrac has a total of ( 2*npars+1 + npars*(npars-1)/2 ) values
  allocate( cfrac(1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  cfrac(1:2*npars+1)=ffrac
  cfrac(2*npars+2:npars*(npars-1)/2+2*npars+1)=frac
  ! not need frac so release it
  deallocate(frac)
  ! allocate frac again, store cfrac in it then release cfrac
  allocate(frac(1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  frac=cfrac
  deallocate(cfrac) 
  ! specify p to be the column number of cols
  ! that is npars*(npars-1)/2+2*npars+1
  p=npars*(npars-1)/2+2*npars+1
  ! now allocate shifted to be the same shape 
  ! with cols and store some values in it
  allocate(shifted(1:npars,p), stat=errorflag(5))
  if(errorflag(5)/=0) return
  do i=1,p
    shifted(:,i)=pars+incr*cols(:,i)
  end do
  ! allocate transcols to store transposed cols
  ! transcols has p rows, npars columns
  allocate( transcols(1:p,1:npars), stat=errorflag(5))
  if(errorflag(5)/=0) return
  transcols=transpose(cols)
  ! now cols will not be needed, release it
  deallocate(cols)
  ! allocate pxcols to store transcols 
  ! in differ style
  allocate(pxcols(p,1+2*npars), stat=errorflag(5))
  if(errorflag(5)/=0) return
  pxcols(:,1)=1.0D+00
  pxcols(:,2:npars+1)=transcols
  pxcols(:,npars+2:2*npars+1)=(transcols)**2
  !
  ncols=0
  ! now allocate xcols with p rows and 
  ! npars*(npars-1)/2 columns
  allocate( xcols(p,1:npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! 
  do i=1,npars-1
    ! in each loop, allocate cxcols 
    allocate(cxcols(p,1:npars-i), stat=errorflag(5))
    if(errorflag(5)/=0) return
    do j=i+1,npars
      cxcols(:,j-i)=transcols(:,i)*transcols(:,j)
    end do
    ! store cxcols to xcols and release it
    xcols(:,ncols+1:ncols+npars-i)=cxcols
    deallocate(cxcols)
    ! update started filling index
    ncols=ncols+npars-i  
  end do
  ! now store pxcols and xcols together to cxcols
  allocate( cxcols(p,1:1+2*npars+npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  cxcols(:,1:1+2*npars)=pxcols
  cxcols(:,2+2*npars:1+2*npars+npars*(npars-1)/2)=xcols
  ! release xcols and pxcols
  deallocate(xcols)
  deallocate(pxcols)
  ! allocate xcols again to store cxcols
  allocate(xcols(p,1:1+2*npars+npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  xcols=cxcols
  ! release cxcols
  deallocate(cxcols)
  ! allocate pxcols again and store some
  ! values to it then release shifted
  ! now pxcols has p rows and 1 column
  allocate(pxcols(1:p,1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  do i=1,p
    pxcols(i,:)=fun34(shifted(:,i))
  end do
  deallocate(shifted)
  ! solve xcols %*% X = pxcols
  ! store solved X values in pxcols
  call GJordan(xcols,pxcols,p,1,solerror,tol)
  if(solerror==1)  errorflag(1)=1
  ! scale pxcols with frac and release frac,
  ! release xcols
  pxcols(:,1)=pxcols(:,1)/frac
  deallocate(xcols)
  deallocate(frac)
  ncols=2*npars+2 
  ! fill diagpar with some new values,  note that
  ! non-diagnal part of digpar are zeros
  ! these value are comes from pxcols, and will
  ! be used to constructe hessian matrix
  do j=1,npars
    do i=1,npars
        if (i==j) diagpar(i,j)=pxcols(1+npars+i,1)
        if (i>j)  diagpar(j+1:npars,j)=pxcols(ncols:ncols+npars-j-1,1)  
    end do
    ncols=ncols+npars-j
  end do
  ! estimate gradients
  do i=1,npars
    gradient(i)=pxcols(1+i,1)
  end do
  ! estimate fun(pars)
  value=pxcols(1,1)
  ! now pxcols will not be needed, release it
  deallocate(pxcols)
  ! estimate hessian matrix
  hessian=diagpar+transpose(diagpar)
  !
  ! check Inf and NaN for value
  if(value .ne. value .or. &
     value+1.0D+00==value)              errorflag(2)=1
  ! check Inf and NaN for gradient
  if(any(gradient .ne. gradient) .or. &
     any(gradient+1.0D+00==gradient))   errorflag(3)=1
  ! check Inf and NaN for hessian
  if(any(hessian .ne. hessian) .or. &
     any(hessian+1.0D+00==hessian))     errorflag(4)=1
  ! if no error appears, return 
  return
  ! INNER FUNCTION FOR HESSION APPROXIMATION
  contains
    function fun34(x)
      implicit none
      real(kind=8)::fun34
      real(kind=8),dimension(npars)::x
      real(kind=8),dimension(ndat)::fvec
      integer(kind=4)::k
      !
      if(model==1) then
        ! for fitting 'cw' with a constant
        fvec=x(npars)
        do k=1,(npars-1)/2
          fvec=fvec+x(k)*dexp(-x(k+(npars-1)/2)*xdat)   
        end do
      else if(model==2) then
        ! for fitting 'lm' with a constant
        maxx=maxval(xdat)
        fvec=x(npars)*xdat/maxx
        do k=1,(npars-1)/2
          fvec=fvec+x(k)*(xdat/maxx)*&
               dexp(-x(k+(npars-1)/2)*xdat**2/2.0D+00/maxx)  
        end do
      end if
      fun34=dsqrt(sum((fvec-ydat)**2))
      return
    end function fun34
end subroutine lmhess_C

subroutine targfunc_C(lamda,nlamda,tim,sig,tol,typ,&
                      ntim,value,ithn,constant,errorflag)
!-----------------------------------------------------------------------------------------------------------------------------------
! targfunc_C() is a subroutine for calculating ithn plus a constant;
! initial lamda values must be provided, then ithn values
! will be calculated using Linear Algebra method according to Bluszcz A (1996).
! For type "cw": I(t)=ithn1*exp(-lamda1*t)+ithn2*exp(-lamda2*t)+...+
!                     ithnk*exp(-lamdak*t)+c, where k=1,2,...,7;
! For type "lm": I(t)=ithn1*(t/max(t))*exp(-lamda1*t^2)/2/max(t))+...+
!                     ithnk*(t/max(t))*exp(-lamdak*t^2)/2/max(t))+c*(t/max(t)), where k=1,2,...,7.
! ===================================================================================================================================
!
! lamda(nlamda),     input:: real values, lamda values.
!
! nlamda,            input:: integer, length of lamda values.
!
! tim(ntim),         input:: real values, time values.
!
! sig(ntim),         input:: real values, signal values.
!
! tol,               input:: real value, maximum tolerance for identifying linear independence.
!
! typ,               input:: integer, fitting type, 1 for type 'cw, 2 for type 'lm'.
!
! ntim,              input:: integer, length of tim and sig.
!
! value,            output:: cal caculated function value.
!
! ithn(nlamda),     output:: estimated I0i values (i=1,2,...maxcomp).
!
! constant,         output:: estimated constant.
!
! errorflag,        output:: integer, error message generated during the calculation:
!                            1.1) if no error appears in Gjordan and no estimated ithn is below zero, errorflag=0;
!                            1.2) if either error appears in Gjordan or any estimated ithn<0, errorflag=1.
! ====================================================================================================================================
! Author:: Peng Jun, 2013.12.14.
!
! Dependence:: subroutine GJordan.
!
! References:: Bluszcz, A., 1996. Exponential function fitting to TL growth data 
!              and similar applications. Geochronometria 13, 135–141.
!
!              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution 
!              to fitting OSL decay curves. Radiation Measurements 41, 886-891.
!--------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                  intent(in)::nlamda      ! length of lamda values
  integer(kind=4),                  intent(in)::ntim        ! length of tim and sig
  real   (kind=8),                  intent(in)::tol         ! maximum tolerance for linear independence
  integer(kind=4),                  intent(in)::typ         ! fitting type ('cw' or 'lm')
  real   (kind=8),dimension(nlamda),intent(in)::lamda       ! lamda values
  real   (kind=8),dimension(ntim),  intent(in)::tim,sig     ! time and signal values
  real   (kind=8),                  intent(out)::value      ! targeted function value
  real   (kind=8),dimension(nlamda),intent(out)::ithn       ! I0i values (i=1,2,...maxcomp)
  real   (kind=8),                  intent(out)::constant   ! estimated constant
  integer(kind=4),                  intent(out)::errorflag  ! error message generaged during the calculation
  ! Local variables
  real(kind=8),dimension(ntim,nlamda+1)::coef               ! coefficent matrix 
  real(kind=8),dimension(ntim,1)::ssig                      ! for storing sig values using matrix
  real(kind=8),dimension(nlamda+1,1)::iithn                 ! for storing ithn values using matrix
  real(kind=8),dimension(ntim)::rowsumcoef                  ! sum values in matrix coef by row
  real(kind=8),dimension(nlamda+1,nlamda+1)::ccoef          ! t(coef)%*%coef 
  integer(kind=4)::maxt                                     ! total stimulation time (P)
  integer(kind=4)::i                                        ! iterative lopping index
  !
  ! Initializing coefficents 
  if(typ==1) then
    ! For type 'cw'
    do i=1,nlamda
      coef(:,i)=dexp(-lamda(i)*tim(:))
    end do
    coef(:,nlamda+1)=1.0D+00
  else if(typ==2) then
    ! For type 'lm'
    maxt=maxval(tim)
    do i=1,nlamda
      coef(:,i)=(tim(:)/maxt)*dexp(-lamda(i)*(tim(:))**2/2.0D+00/maxt)
    end do
    coef(:,nlamda+1)=tim(:)/maxt
  end if
  !
  ! Store sig[array: ntim] values in matrix ssig[matrix: ntim*1]
  ssig(:,1)=sig
  !
  ! Calculate ithn values using Linear Algebra method
  ! initialize iithn and ccoef for Linear Algebra using
  iithn=matmul(transpose(coef),ssig)
  ccoef=matmul(transpose(coef),coef)
  call GJordan(ccoef,iithn,nlamda+1,1,errorflag,tol)
  ! At this point, if no error appears when calling
  ! GJordan, errorflag will be 0
  !
  ! Pass iithn[nlamda*1] to ithn[nlamda]
  ithn=iithn(1:nlamda,1)
  ! Pass iithn[nlamda+1] to constant
  constant=iithn(nlamda+1,1)
  !
  ! Return if error appears. error situations including:
  ! 1) matrix ccoef is singular; 2) any obtained ithn value is below zero
  ! then errorflag=1, ithn=ithn(wrong values), value=1.00D+30
  if(errorflag==1 .or. any(ithn<1.0D-10) .or. constant<1.0D-10)  then
    value=1.00D+30
    errorflag=1
    return
  end if
  !
  ! Reset matirx coef with calculated ithn
  if(typ==1) then
    ! For type 'cw'
    do i=1,nlamda
      coef(:,i)=ithn(i)*dexp(-lamda(i)*tim(:))
    end do
    coef(:,nlamda+1)=constant
  else if(typ==2) then
    ! For type 'lm'
    do i=1,nlamda
      coef(:,i)=ithn(i)*(tim(:)/maxt)*&
                dexp(-lamda(i)*(tim(:))**2/2.0D+00/maxt)
    end do
    coef(:,nlamda+1)=constant*(tim(:)/maxt)
  end if
  !
  ! Calculate rowsums of matrix coef
  rowsumcoef=sum(coef,2)
  !
  ! Calculate residual value
  value=sum((sig-rowsumcoef)**2) 
  ! Now return
  return
end subroutine targfunc_C
