subroutine fitGrowth(dose,ltx,sltx,ndat,n2,pars,stdp,upb,&
                     model,origin,nstart,fvec1,fvalue,message)
!-----------------------------------------------------------------
! Subroutine fitGrowth() is used to fitting a growth curve.
!-----------------------------------------------------------------
! dose(ndat):: input, real values, the dose values.
! ltx(ndat) :: input, real values, the standardised OSLs.
! sltx(ndat):: input, real values, errors of OSLs.
!       ndat:: input, integer, number of data points.
!         n2:: input, integer, dimension of the problem.
!   pars(n2):: output, real values, estimated pars.
!   stdp(n2):: output, real values, errors of pars.
!        upb:: input, real value, upper limit on b.
!      model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp.
!     origin:: input, integer, 0=origin, 1=non-origin.
!     nstart:: input, integer, number of simulations.
!      fvec1:: output, real value, the minimized chi-square.
!    message:: output, integer, 0=success, 1=fail.
!-----------------------------------------------------------------
! Author:: Peng Jun, 2014.09.24.
!-----------------------------------------------------------------
! Dependence:: subroutine linearFit; subroutine inipars;----------
!              subroutine lmgrwcve.
!-----------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2, model,&
                                  origin, nstart
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat),&
                                  sltx(ndat), upb
    real   (kind=8), intent(out):: pars(n2), stdp(n2),&
                                   fvec1(ndat), fvalue
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, i
    real   (kind=8):: minCond, maxDose, ran(2), outp(3),&
                      locp(5), cpars(n2), cstdp(n2),&
                      cfvec1(ndat), cfvalue, cond, grad
    !
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fvalue = -99.0
    message = 1
    !
    if (model==0) then
        call linearFit(dose,ltx,sltx,ndat,pars,&
                       stdp,n2,fvec1,fvalue,info)
        if (info==0) message = 0
        if (message/=0) return
    else if (model==1 .or. model==2 .or. model==3) then
        minCond = 1.0D+20
        maxDose = maxval(dose)
        !
        loopA: do i=1, nstart
            call random_number(ran)
            ran = ran*upb
            call inipars(ran(1),ran(2),model,origin,&
                         dose,ltx,sltx,ndat,outp,info)
            if (info/=0) cycle loopA
            !
            if (info==0) then
                locp = 0.0
                if (model==1 .or. model==2) then
                    locp(1) = outp(1)
                    locp(2) = ran(1)
                    locp(3) = outp(2)
                    locp(4) = outp(3)
                else if (model==3) then
                    locp(1) = outp(1)
                    locp(2) = ran(1)
                    locp(3) = outp(2)
                    locp(4) = ran(2)
                    locp(5) = outp(3)
                end if
                cpars = locp(1:n2)
            end if
            !
            call lmgrwcve(dose,ltx,sltx,ndat,cpars,cstdp,n2,&
                          model,origin,cfvec1,cfvalue,cond,info)
            !
            if (info/=0) cycle loopA
            !
            if (model==1 .or. model==3) then
                locp = 0.0
                locp(1:n2) = cpars
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                       locp(3)*locp(4)*exp(-locp(4)*maxDose)
                if (grad<1.0D-13) cycle loopA
            end if
            !
            if (cond<minCond) then
                pars = cpars
                stdp = cstdp
                fvec1 = cfvec1
                fvalue = cfvalue
                minCond = cond
                message = 0
            end if
        end do loopA
        !
        if (message/=0) return
    end if
    !
    return
end subroutine fitGrowth                                         
