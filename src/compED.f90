subroutine compED(ed1,sed1,ndat,ncomp,addsigma,&
                  pars,stdp,maxlik,bic,info)
!--------------------------------------------------------
! Subroutine compED() is used for estimating a
! finite mixture with a given number of components.
!--------------------------------------------------------
!     ed1(ndat):: input, real values, unlogged ED values.
!    sed1(ndat):: input, real values, absolute errors.
!          ndat:: input, integer, number of data points.
!      addsigma:: input, real value, additional error.
! pars(2,ncomp):: output, real values, estimated pars.
! stdp(2,ncomp):: output, real values, estimated stdpars.
!        maxlik:: output, real value, logged likelihood.
!           bic:: output, real value, BIC value.
!          info:: output, integer, 0=success, 1=fail.
!--------------------------------------------------------
! Author:: Peng Jun, 2014.10.01.
!--------------------------------------------------------
! Dependence:: subroutine fmmED; Subroutine camED.-------
!--------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, ncomp
    real   (kind=8), intent(in):: ed1(ndat), sed1(ndat), addsigma
    real   (kind=8), intent(out):: pars(2,ncomp), stdp(2,ncomp),&
                                   maxlik, bic
    integer(kind=4), intent(out):: info
    ! Local variables.
    integer(kind=4):: i, message
    real   (kind=8):: ed(ndat), sed(ndat), maxMaxlik,&
                      cpars(2,ncomp), cstdp(2,ncomp),&
                      cmaxlik, cbic, inimu(ncomp+4)
    !
    ed = log(ed1)
    sed = sed1/ed1
    ! 
    if (ncomp==1) then
        ! CAM.
        info = 0
        call camED(ed,sed,ndat,addsigma,&
                   pars,stdp,maxlik,bic)
        return
    else 
        ! FMM.
        info = 1
        maxMaxlik = -1.0D+20
        !
        do i=1, ncomp+4
            inimu(i) = minval(ed) +&
                      (maxval(ed)-minval(ed))/(ncomp+3)*&
                       real(i-1)
        end do
        !           
        do i=1, 5
            pars(1,:) = 1.0/real(ncomp)
            pars(2,:) = inimu(i:i+ncomp-1)
            !
            call fmmED(ed,sed,ndat,ncomp,addsigma,&
                       pars,stdp,maxlik,bic,message)
            if (message==0 .and. maxlik>maxMaxlik) then
                cpars = pars
                cstdp = stdp
                cmaxlik = maxlik 
                cbic = bic
                maxMaxlik = maxlik
                info = 0
            end if
        end do
        !
        if (info/=0) return
        !
        pars = cpars
        stdp = cstdp
        maxlik = cmaxlik
        bic = cbic
    end if
    !
    return
end subroutine compED                      
