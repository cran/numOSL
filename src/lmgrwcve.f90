subroutine lmgrwcve(dose,ltx,sltx,ndat,pars,stdp,n2,&
                    model,origin,fvec1,fvalue,cond,message)
!-------------------------------------------------------------------
! Subroutine lmgrwcve() is used for fitting a growth curve.
!-------------------------------------------------------------------
!  dose(ndat):: input, rea values, the dose values.
!   ltx(ndat):: input, real values, the Lx/Tx values.
!  sltx(ndat):: input, real values, the std of Lx/Tx values.
!        ndat:: input, integer, number of data points.
!    pars(n2):: input/output, real values, estimated parameters.
!    stdp(n2):: output, real values, std of parameters.
!          n2:: input, integer, dimension of the problem.
!       model:: input, integer, fitting model, 1=exp, 2=lexp, 3=dexp.
!      origin:: input, integer, 0=origin, 1=non-origin.
! fvec1(ndat):: output, real values, fitted Lx/Tx values.
!        cond:: output, real vlaue, the logged condition number.
!     message:: output, integer, 0=success, 1=fail.
!-------------------------------------------------------------------
! Author:: Peng Jun, 2014.09.05.
!-------------------------------------------------------------------
! Dependence:: subroutine lmdif1; subroutine numHess;---------------
!              subroutine inverse; inner function func, fcn.--------
!-------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2, model, origin
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat), sltx(ndat)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: stdp(n2), fvec1(ndat),& 
                                   fvalue, cond
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, IWA(n2), LWA, errorflag, singular, i
    real   (kind=8):: WA(ndat*n2+5*n2+ndat), fvec(ndat),& 
                      DPMPAR, TOL, cond1, cond2
    real   (kind=8):: hess(n2,n2), diag(n2)
    !
    stdp = -99.0
    fvec1 = -99.0
    fvalue = -99.0
    cond = -99.0
    !
    TOL = dsqrt(DPMPAR(1))
    LWA = ndat*n2+5*n2+ndat
    !
    call lmdif1(func,ndat,n2,pars,fvec,&
                TOL,info,IWA,WA,LWA)
    if (info==1 .or. info==2 .or. info==3) then
        message = 0 
    else 
        message = 1
        return
    end if
    !
    if (any(pars(1:n2-origin)<=0.0)) then
        message = 1
        return
    end if
    !
    fvec1 = ltx + fvec*sltx
    fvalue = sum(fvec**2)
    !
    call numHess(fcn,pars,n2,hess,errorflag)
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
        ! Calculate a vector that contains weighted residuals.
        subroutine func(m,n,x,fvec,iflag)
            implicit none
            integer(kind=4):: m, n, iflag
            real   (kind=8):: x(n), fvec(m), locx(5)
            !
            locx = 0.0
            locx(1:n) = x(1:n)
            !
            if (model==1) then
                fvec = locx(1)*(1.0-exp(-locx(2)*dose))+&
                       locx(3)
            else if (model==2) then
                fvec = locx(1)*(1.0-exp(-locx(2)*dose))+&
                       locx(3)*dose+locx(4)
            else if (model==3) then
                fvec = locx(1)*(1.0-exp(-locx(2)*dose))+&
                       locx(3)*(1.0-exp(-locx(4)*dose))+&
                       locx(5)
            end if
            fvec = (fvec - ltx) / sltx
            !
            return
        end subroutine func
        ! Calculate the sqrt of the sum 
        ! of the squared weighted residuals.
        function fcn(x)
            implicit none
            real   (kind=8):: x(n2), vec(ndat),&
                              locx(5), fcn
            !
            locx = 0.0
            locx(1:n2) = x(1:n2)
            !
            if (model==1) then
                vec = locx(1)*(1.0-exp(-locx(2)*dose))+&
                      locx(3)
            else if (model==2) then
                vec = locx(1)*(1.0-exp(-locx(2)*dose))+&
                      locx(3)*dose+locx(4)
            else if (model==3) then
                vec = locx(1)*(1.0-exp(-locx(2)*dose))+&
                      locx(3)*(1.0-exp(-locx(4)*dose))+&
                      locx(5)
            end if
            vec = (vec - ltx) / sltx
            fcn = sqrt(sum(vec**2))
            !
            return
        end function fcn
end subroutine lmgrwcve 
