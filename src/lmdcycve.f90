subroutine lmdcycve(tim,sig,ntim,pars,stdp,n2,addc,&
                    typ,fvec1,fvalue,cond,message)
!--------------------------------------------------------------------------
! Subroutine lmdcycve() is used for estimating the parameters
! of an OSL decay curve using the levenberg-marquardt algorithm.
!--------------------------------------------------------------------------
!   tim(ntim):: input, real values, the time values.
!   sig(ntim):: input, real values, the decay signal values.
!        ntim:: input, integer, number of data points.
!    pars(n2):: input/output, real values, the parameters, overwritten.
!    stdp(n2):: output, real values, the standard errors of parameters.
!          n2:: input, integer, the dimension of the problem.
!        addc:: input, integer, add a constant (1) not (0).
!         typ:: input, integer, type of OSL, 1=CW-OSL, 2=LM-OSL.
! fvec1(ntim):: output, real values, fitted signal values.
!      fvalue:: output, real value, minimized sum of squared residuals.
!        cond:: output, real vlaue, log-condition (Inf) of hessian matrix.
!     message:: output, integer, 0 means a successful work, else 1.
!--------------------------------------------------------------------------
! Author:: Peng Jun, 2014.08.31.
!--------------------------------------------------------------------------
! Dependence:: subroutine lmdif1; subroutine numHess; subroutine inverse.--
!              Inner subroutine func; inner function fcn.------------------
!--------------------------------------------------------------------------
! Reference:: J.J. Moré, "The Levenberg-Marquardt algorithm: implementation 
!             and theory," in Lecture Notes in Mathematics 630: Numerical  
!             Analysis, G.A. Watson (Ed.), Springer-Verlag: Berlin, 
!             1978, pp.105-116.
!--------------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ntim, n2, addc, typ
    real   (kind=8), intent(in):: tim(ntim), sig(ntim)                     
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: stdp(n2), fvec1(ntim), fvalue, cond
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, IWA(n2), LWA, errorflag, singular, i
    real   (kind=8):: WA(ntim*n2+5*n2+ntim), fvec(ntim), DPMPAR, TOL, cond1, cond2
    real   (kind=8):: timmaxt(ntim), tim2maxt2(ntim), hess(n2,n2), diag(n2)
    !
    stdp = -99.0
    fvec1 = -99.0
    fvalue = -99.0
    cond = -99.0
    !
    if (typ==2) then
        timmaxt = tim/maxval(tim)
        tim2maxt2 = tim**2/maxval(tim)/2.0
    end if
    !
    TOL = dsqrt(DPMPAR(1))
    LWA = ntim*n2+5*n2+ntim
    !
    call lmdif1(func,ntim,n2,pars,fvec,&
                TOL,info,IWA,WA,LWA)
    if (info==1 .or. info==2 .or. info==3) then
        message = 0
    else 
        message = 1
        return
    end if 
    !
    if (any(pars<=0.0)) then
        message = 1
        return
    end if
    !
    fvec1 = sig + fvec
    fvalue = sum(fvec**2)   
    !
    call  numHess(fcn,pars,n2,hess,errorflag)
    if (errorflag/=0) then
        message = 1
        return
    end if
    cond1 = maxval(sum(abs(hess),dim=2))
    !
    call inverse(hess,n2,singular) 
    if (singular/=0) then
        message = 1
        return
    end if
    cond2 = maxval(sum(abs(hess),dim=2))
    !
    cond = log(cond1*cond2)
    !
    do i=1, n2
        diag(i) = hess(i,i)
    end do
    if (any(diag<0.0)) then
        message = 1
        return
    end if
    !
    stdp = sqrt(diag)
    !
    return
    !
    contains 
        ! Calculate a vector that contains residuals.
        subroutine func(m,n,x,fvec,iflag)
            implicit none 
            integer(kind=4):: m, n, iflag, i
            real   (kind=8):: x(n), fvec(m)                          
            !
            if (typ==1) then
                if (addc==0) then
                    fvec = 0.0
                    do i=1, n/2
                        fvec = fvec + x(i)*x(i+n/2)*&
                               exp(-x(i+n/2)*tim)
                    end do
                else if (addc==1) then
                    fvec = x(n)
                    do i=1, (n-1)/2
                        fvec = fvec + x(i)*x(i+(n-1)/2)*&
                               exp(-x(i+(n-1)/2)*tim)
                    end do
                end if
            else if (typ==2) then
                if (addc==0) then
                    fvec = 0.0
                    do i=1, n/2
                        fvec = fvec + x(i)*timmaxt*x(i+n/2)*&
                               exp(-x(i+n/2)*tim2maxt2)
                    end do
                else if (addc==1) then
                    fvec = x(n)*timmaxt
                    do i=1, (n-1)/2
                        fvec = fvec + x(i)*timmaxt*x(i+(n-1)/2)*&
                               exp(-x(i+(n-1)/2)*tim2maxt2) 
                    end do
                end if
            end if
            fvec = fvec - sig
            !
            return
        end subroutine func
        ! Calculate the sqrt of the sum of the squared residuals.
        function fcn(x)
            implicit none 
            integer(kind=4):: i
            real   (kind=8):: x(n2), vec(ntim), fcn                       
            !
            if (typ==1) then
                if (addc==0) then
                    vec = 0.0
                    do i=1, n2/2
                        vec = vec + x(i)*x(i+n2/2)*&
                              exp(-x(i+n2/2)*tim)
                    end do
                else if (addc==1) then
                    vec = x(n2)
                    do i=1, (n2-1)/2
                        vec = vec + x(i)*x(i+(n2-1)/2)*&
                              exp(-x(i+(n2-1)/2)*tim)
                    end do
                end if
            else if (typ==2) then
                if (addc==0) then
                    vec = 0.0
                    do i=1, n2/2
                        vec = vec + x(i)*timmaxt*x(i+n2/2)*&
                              exp(-x(i+n2/2)*tim2maxt2)
                    end do
                else if (addc==1) then
                    vec = x(n2)*timmaxt
                    do i=1, (n2-1)/2
                        vec = vec + x(i)*timmaxt*x(i+(n2-1)/2)*&
                              exp(-x(i+(n2-1)/2)*tim2maxt2) 
                    end do
                end if
            end if
            vec = vec - sig
            fcn = sqrt(sum(vec**2))
            !
            return
        end function fcn
end subroutine lmdcycve
