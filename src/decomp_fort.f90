subroutine decomp_fort(tim,sig,ntim,pars,stdp,n2,uw,addc,typ,&
                       factor,f,cr,maxiter,tol,fvec1,fmin,message)
!------------------------------------------------------------------------
! Subroutine decomp_fort is used for OSL decay curve decomposition.
!------------------------------------------------------------------------
!   tim(ntim):: input, real values, time values.
!   sig(ntim):: input, real values, decay signal values.
!        ntim:: input, integer, number of data points.
!    pars(n2):: output, real values, estimated parameters.
!    stdp(n2):: output, real values, estimated std of pars.
!          n2:: input, integer, dimension of the problem (>=2).
!          uw:: input, integer, 0=un-weighted, 1=weighted.
!        addc:: input, integer, 0=non-constant, 1=constant.
!         typ:: input, integer, type of OSL, 1=CW-OSL, 2=LM-OSL.
!      factor:: input, integer, NP=factor*(n2-addc)/2.
!           f:: input, real value, the differential weight.
!          cr:: input, real value, the crossover probability.
!     maxiter:: input, integer, the maximum number of iterations.
!         tol:: input, real value, tolerance for stopping iteration.
! fvec1(ntim):: output, real values, predicted decay signal values.
!        fmin:: output, real value, minimumized objective.
!     message:: output, integer, 0=success, 1=fail.
!------------------------------------------------------------------------
! Author:: Peng Jun, 2023.08.30. 
!------------------------------------------------------------------------
! Dependence:: subroutine diffev; 
!              subroutine lmfit; 
!              subroutine comb_next.
!------------------------------------------------------------------------
! Reference::   Bluszcz, A., Adamiec, G., 2006. Application of 
!               differential evolution to fitting OSL decay curves. 
!               Radiation Measurements 41, 886-891.
!-----------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer, intent(in):: ntim, n2, uw, addc,&
                          typ, factor, maxiter
    real(kind(1.0d0)), intent(in):: tim(ntim), sig(ntim),&
                                    f, cr, tol
    real(kind(1.0d0)), intent(out):: pars(n2), stdp(n2),&
                                     fvec1(ntim), fmin
    integer, intent(out):: message
    ! Local variables.
    integer:: np, n1, iflag, i,& 
              ivec((n2-addc)/2),& 
              lmInfo, model
    real(kind(1.0d0)):: constant, lamda((n2-addc)/2),ithn((n2-addc)/2),&
                        agents(factor*(n2-addc)/2,(n2-addc)/2), cpars(n2),&
                        cstdp(n2), cfvec1(ntim), cfmin, minValue, wght1(ntim)
    integer, parameter:: nperm(7)=(/7,21,35,35,21,7,1/)
    real(kind(1.0d0)), parameter:: initry(7)=(/32.0,2.5,0.62,0.15,&
                                             0.023,0.0022,0.0003/)
    logical:: done
    !
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    message = 1
    !
    n1 = (n2-addc)/2
    np = factor*n1  
    !
    if (uw==0) then
        wght1 = 1.0
    else if (uw==1) then
        wght1 = sqrt(sig)
    end if
    !
    if (typ==1) then
        model = 4
    else if (typ==2) then
        model = 5
    end if
    !
    call diffev(tim,sig,wght1,ntim,np,f,cr,maxiter,tol,typ,addc,&
                n1,lamda,ithn,constant,agents,fvec1,fmin,iflag)
    if (iflag==0) then
        cpars(1:n1) = ithn
        cpars(n1+1:2*n1) = lamda
        if (addc==1) cpars(n2) = constant
        call lmfit(tim,sig,wght1,ntim,cpars,cstdp,n2,&
                   model,cfvec1,cfmin,lmInfo)
    end if
    !
    minValue = 1.0D+20
    if (iflag==0 .and. lmInfo==0) then
        pars = cpars
        stdp = cstdp
        fvec1 = cfvec1
        fmin = cfmin
        minValue = cfmin
        message = 0
    end if
    !
    done = .true.
    loopA: do i=1, nperm(n1)
        if (iflag==0) then
            cpars(1:n1) = sum(ithn)/real(n1)
        else 
            cpars(1:n1) = 5000.0
        end if
        !
        call comb_next(7,n1,ivec,done)
        cpars(n1+1:2*n1) = initry(ivec)
        !
        if (addc==1) then
            if (typ==1) cpars(n2) = sum(sig(ntim-2:ntim))/3.0 
            if (typ==2) cpars(n2) = tim(ntim)*sum(sig(ntim-2:ntim)*&
                                    tim(ntim-2:ntim))/sum((tim(ntim-2:ntim))**2)
        end if
        !
        call lmfit(tim,sig,wght1,ntim,cpars,cstdp,n2,&
                   model,cfvec1,cfmin,lmInfo)
        if (lmInfo==0 .and. cfmin<minValue) then
            pars = cpars
            stdp = cstdp
            fvec1 = cfvec1
            fmin = cfmin
            minValue = cfmin
            message = 0
        end if
    end do loopA
    !
    return
end subroutine decomp_fort
