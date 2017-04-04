subroutine fitGrowth_fort(dose,ltx,sltx,ndat,n2,pars,stdp,&
                          model,uw,fvec1,fmin,message)
!--------------------------------------------------------------
! Subroutine fitGrowth_fort is used to fitting a growth curve.
!--------------------------------------------------------------
! dose(ndat):: input, real values, dose values.
! ltx(ndat) :: input, real values, standardised OSLs.
! sltx(ndat):: input, real values, errors of OSLs.
!       ndat:: input, integer, number of data points.
!         n2:: input, integer, number of pars (>=1).
!   pars(n2):: output, real values, estimated pars.
!   stdp(n2):: output, real values, errors of pars.
!      model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp.
!         uw:: input, integer, 0=un-weighted, 1=weighted.
!      fvec1:: output, real value, minimized objective.
!    message:: output, integer, 0=success, 1=fail.
!--------------------------------------------------------------
! Author:: Peng Jun, 2017.03.30. 
!--------------------------------------------------------------
! Dependence:: subroutine linefit; 
!              subroutine inipars;
!              subroutine lmfit.
!--------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2,&
                                  model, uw
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat),&
                                  sltx(ndat)
    real   (kind=8), intent(out):: pars(n2), stdp(n2),&
                                   fvec1(ndat), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, i, j
    real   (kind=8):: maxDose, ran(2), outp(3),&
                      locp(5), cpars(n2), cstdp(n2),&
                      cfvec1(ndat), cfmin, wght1(ndat),&
                      minValue, inib(24)
    !
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fmin = -99.0
    message = 1
    !
    if (uw==0) then
        ! un-Weighted.
        wght1 = 1.0
    else if (uw==1) then
        ! Weighted.
        wght1 = sltx
    end if
    !
    !
    if (model==0) then
        !
        ! Linear model.
        call linefit(dose,ltx,wght1,ndat,pars,&
                     stdp,n2,fvec1,fmin,info)
        if (info==0) message = 0
        if (message/=0) return
        !
    else
        !
        ! Non-linear model.
        minValue = 1.0D+20
        maxDose = maxval(dose)
        !
        do i=1, 12
            inib(2*i-1:2*i) = (10.0)**(i-11)*(/1.0, 5.0/)
        end do
        !
        ! Initialization.
        if (model==1 .or. model==2) then
            !
            loopA: do i=1, 24
                !
                ran(1) = inib(i)
                ran(2) = 0.0
                !
                call inipars(ran(1),ran(2),model,n2,&
                             dose,ltx,wght1,ndat,outp,info)
                if (info/=0) cycle loopA
                !
                locp = 0.0
                locp(1) = outp(1)
                locp(2) = ran(1)
                locp(3) = outp(2)
                locp(4) = outp(3)
                !
                cpars = locp(1:n2)
                !
                call lmfit(dose,ltx,wght1,ndat,cpars,cstdp,&
                           n2,model,cfvec1,cfmin,info)
                if (info/=0) cycle loopA
                !
                if (cfmin<minValue) then
                    pars = cpars
                    stdp = cstdp
                    fvec1 = cfvec1
                    fmin = cfmin
                    minValue = cfmin
                    message = 0
                end if
            end do loopA
            !
        else if (model==3) then
            loopB: do i=1, 24
                loopC: do j=i, 24
                    !
                    ran(1) = inib(i)
                    ran(2) = inib(j)
                    !
                    call inipars(ran(1),ran(2),model,n2,&
                                 dose,ltx,wght1,ndat,outp,info)
                    if (info/=0) cycle loopC
                    !
                    locp = 0.0
                    locp(1) = outp(1)
                    locp(2) = ran(1)
                    locp(3) = outp(2)
                    locp(4) = ran(2)
                    locp(5) = outp(3)
                    !
                    cpars = locp(1:n2)
                    !
                    call lmfit(dose,ltx,wght1,ndat,cpars,cstdp,&
                               n2,model,cfvec1,cfmin,info)
                    if (info/=0) cycle loopC
                    !
                    if (cfmin<minValue) then
                        pars = cpars
                        stdp = cstdp
                        fvec1 = cfvec1
                        fmin = cfmin
                        minValue = cfmin
                        message = 0
                    end if
                end do loopC
            end do loopB
            !
        end if
        !
    end if
    !
    return
end subroutine fitGrowth_fort                                         
