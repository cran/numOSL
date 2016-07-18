subroutine fitGOK(dose,ltx,sltx,ndat,n2,pars,stdp,&
                  uw,fvec1,fmin,message)
!--------------------------------------------------------------
! Subroutine fitGOK is used for fitting a GOK model.
!--------------------------------------------------------------
! dose(ndat):: input, real values, dose values.
! ltx(ndat) :: input, real values, standardised OSLs.
! sltx(ndat):: input, real values, errors of OSLs.
!       ndat:: input, integer, number of data points.
!         n2:: input, integer, number of pars (>=1).
!   pars(n2):: output, real values, estimated pars.
!   stdp(n2):: output, real values, errors of pars.
!         uw:: input, integer, 0=un-weighted, 1=weighted.
!      fvec1:: output, real value, minimized objective.
!    message:: output, integer, 0=success, 1=fail.
!--------------------------------------------------------------
! Author:: Peng Jun, 2016.07.08.
!--------------------------------------------------------------
! Dependence:: subroutine inipars;
!              subroutine lmfit1.
!--------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, n2,&
                                  uw
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat),&
                                  sltx(ndat)
    real   (kind=8), intent(out):: pars(n2), stdp(n2),&
                                   fvec1(ndat), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: info, i, j
    integer(kind=4), parameter:: model=7
    real   (kind=8):: maxDose, ran(2), outp(3),&
                      locp(4), cpars(n2), cstdp(n2),&
                      cfvec1(ndat), cfmin, grad,& 
                      wght1(ndat), minValue, inib(24), rcv(11)
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
    minValue = 1.0D+20
    maxDose = maxval(dose)
    !
    do i=1, 12
        inib(2*i-1:2*i) = (10.0)**(i-11)*(/1.0, 5.0/)
    end do
    !
    rcv = (/((10.0)**(i-6), i=1, 11)/)
    !
    loopA: do i=1, 24
        loopB: do j=1, 11
            !
            ran(1) = inib(i)
            ran(2) = rcv(j)
            !
            call inipars(ran(1),ran(2),model,n2,&
                         dose,ltx,wght1,ndat,outp,info)
            if (info/=0) cycle loopB
            !
            locp = 0.0
            locp(1) = outp(1)
            locp(2) = ran(1)
            locp(3) = ran(2)
            locp(4) = outp(2)
            !
            cpars = locp(1:n2)
            !
            call lmfit1(dose,ltx,wght1,ndat,cpars,cstdp,&
                        n2,cfvec1,cfmin,info)
            if (info/=0) cycle loopB
            !
            ! Check the saturating level.
            locp = 0.0
            locp(1:n2) = cpars
            grad = locp(1)*locp(2)*(1.0+locp(2)*locp(3)*maxDose)**&
                   (-1.0/locp(3)-1.0)
            if (grad<1.0D-13) cycle loopB
            !
            if (cfmin<minValue) then
                pars = cpars
                stdp = cstdp
                fvec1 = cfvec1
                fmin = cfmin
                minValue = cfmin
                message = 0
            end if
        end do loopB
    end do loopA
    !
    return
    !
end subroutine fitGOK                                       
