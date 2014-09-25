subroutine decomp(tim,sig,ntim,pars,stdp,n2,addc,typ,factor,&
                  f,cr,maxiter,tol,fvec1,fvalue,message)
!----------------------------------------------------------------------------
! Subroutine decomp() is used for OSL decay curve decomposition.
!----------------------------------------------------------------------------
!   tim(ntim):: input, real values, the time values.
!   sig(ntim):: input, real values, the decay signal values.
!        ntim:: input, integer, number of data points.
!    pars(n2):: output, real values, estimated parameters.
!    stdp(n2):: output, real values, estimated std of pars.
!          n2:: input, integer, dimension of the problem.
!        addc:: input, integer, add a constant (1) or not (0).
!         typ:: input, integer, type of OSL, 1=CW-OSL, 2=LM-OSL.
!      factor:: input, integer, NP=factor*(n2-addc)/2.
!           f:: input, real value, the differential weight.
!          cr:: input, real value, the crossover probability.
!     maxiter:: input, integer, the maximum number of iterations.
!         tol:: input, real value, a tolerance for stopping the iteration.
! fvec1(ntim):: output, real values, predicted decay signal values.
!      fvalue:: output, real value, the minimumized sum of sqaured residuals.
!     message:: output, integer, 0 means a successful work, else 1.
!---------------------------------------------------------------------------
! Author:: Peng Jun, 2014.09.03.
!---------------------------------------------------------------------------
! Dependence:: subroutine diffev; subroutine lmdcyve; subroutine comb_next.-
!---------------------------------------------------------------------------
! Reference::   Bluszcz, A., Adamiec, G., 2006. Application of 
!               differential evolution to fitting OSL decay curves. 
!               Radiation Measurements 41, 886-891.
!--------------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ntim, n2, addc, typ,&
                                  factor, maxiter
    real   (kind=8), intent(in):: tim(ntim), sig(ntim),&
                                  f, cr, tol
    real   (kind=8), intent(out):: pars(n2), stdp(n2),&
                                   fvec1(ntim), fvalue
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: np, n1, iflag, i, ivec((n2-addc)/2), lmInfo
    real   (kind=8):: constant, lamda((n2-addc)/2), ithn((n2-addc)/2),&
                      agents(factor*(n2-addc)/2,(n2-addc)/2), cpars(n2),&
                      cstdp(n2), cfvec1(ntim), cfvalue, cond, minCond
    integer(kind=4), parameter:: nperm(7)=(/7,21,35,35,21,7,1/)
    real   (kind=8), parameter:: initry(7)=(/32.0,2.5,0.62,0.15,&
                                             0.023,0.0022,0.0003/)
    logical:: done
    !
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fvalue = -99.0
    message = 1
    !
    n1 = (n2-addc)/2
    np = factor*n1  
    !
    call diffev(tim,sig,ntim,np,f,cr,maxiter,tol,typ,addc,&
                n1,lamda,ithn,constant,agents,fvec1,fvalue,iflag)
    !
    if (iflag==0) then
        cpars(1:n1) = ithn
        cpars(n1+1:2*n1) = lamda
        if (addc==1) cpars(n2) = constant
        call lmdcycve(tim,sig,ntim,cpars,cstdp,n2,&
                      addc,typ,cfvec1,cfvalue,cond,lmInfo)
    end if
    !
    minCond = 1.0D+20
    if (iflag==0 .and. lmInfo==0) then
        pars = cpars
        stdp = cstdp
        fvec1 = cfvec1
        fvalue = cfvalue
        minCond = cond
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
                                                  tim(ntim-2:ntim))/&
                                              sum((tim(ntim-2:ntim))**2)
        end if
        !
        call lmdcycve(tim,sig,ntim,cpars,cstdp,n2,&
                      addc,typ,cfvec1,cfvalue,cond,lmInfo)
        if (lmInfo==0 .and. cond<minCond) then
            pars = cpars
            stdp = cstdp
            fvec1 = cfvec1
            fvalue = cfvalue
            minCond = cond
            message = 0
        end if
    end do loopA
    !
    return
end subroutine decomp
