subroutine goodComp(ed1,sed1,ndat,maxcomp,&
                    gcomp,addsigma,message)
!----------------------------------------------------
! Subroutine goodComp() is used for searching
! the optimal number of components with BIC.
!----------------------------------------------------
!  ed1(ndat):: input, real values, unlogged EDs.
! sed1(ndat):: input, real vlaues, errors of EDs.
!       ndat:: input, integer, number of data points.
!    maxcomp:: input, integer, maximum component.
!      gcomp:: output, integer, the optimal component.
!   addsigma:: input, real value, additional error.
!    message:: output, integer, 0=success, 1=fail.
!----------------------------------------------------
! Author: Peng Jun, 2014.09.16.
!----------------------------------------------------
! Dependence:: subroutine compED.--------------------
!----------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, maxcomp
    real   (kind=8), intent(in):: ed1(ndat), sed1(ndat), addsigma
    integer(kind=4), intent(out):: gcomp, message
    ! Local variables.
    real   (kind=8):: minBic, maxlik, bic
    real   (kind=8), allocatable:: pars(:,:), stdp(:,:)
    integer(kind=4):: i, info
    !
    minBic = 1.0D+20 
    gcomp = -99
    message = 1
    !
    do i=2, maxcomp
        allocate(pars(2,i),stdp(2,i),stat=info)
        if (info/=0) then
            message = 1
            return
        end if
        !
        call compED(ed1,sed1,ndat,i,addsigma,&
                    pars,stdp,maxlik,bic,info)
        if (info==0 .and. bic<minBic) then
            minBIc = bic
            gcomp = i
            message = 0
        end if
        !
        deallocate(pars,stdp,stat=info)
        if (info/=0) then
            message = 1
            return
        end if
        !
    end do
    !
    return
end subroutine goodcomp                              
