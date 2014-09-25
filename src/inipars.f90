subroutine inipars(b1,b2,model,origin,dose,&
                   ltx,sltx,ndat,outp,info)
!------------------------------------------------------------
! Subroutine inipars() is used to estimate other
! parameters of a growth curve using fixed b values.
!------------------------------------------------------------
!         b1:: input, real value, the first b value.
!         b2:: input, real value, the second b value.
!      model:: input, integer, 1=exp, 2=lexp, 3=dexp.
!     origin:: input, integer, 0=origin, 1=non-origin.
! dose(ndat):: input, real values, dose values.
!  ltx(ndat):: input, real values, Lx/Tx values.
! sltx(ndat):: input, real values, std of Lx/Tx values.
!       ndat:: input, integer, number of data points.
!    outp(3):: output, real values, the estimated parameters.
!       info:: output, integer, 0=success, 1=fail.
!------------------------------------------------------------
! Author:: Peng Jun, 2014.09.05.
!------------------------------------------------------------
! Dependence:: subroutine gjordan.---------------------------
!------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: model, origin, ndat
    real   (kind=8), intent(in):: b1, b2, dose(ndat),&
                                  ltx(ndat), sltx(ndat)
    real   (kind=8), intent(out):: outp(3)
    integer(kind=4), intent(out):: info
    ! Local variables.
    real   (kind=8):: coefMat(ndat,3), ltxMat(ndat,1)
    real   (kind=8), allocatable:: aMat(:,:), bMat(:,:)
    integer(kind=4):: ndim, alc, singular
    ! 
    outp = 0.0
    info = 0
    !
    if (model==1) then
        if (origin==0) then
            ! y=a*(1-exp(-b*x)).
            coefMat(:,1) = (1.0-exp(-b1*dose)) / sltx
            ndim = 1
        else if (origin==1) then
            ! y=a*(1-exp(-b*x))+c.
            coefMat(:,1) = (1.0-exp(-b1*dose)) / sltx
            coefMat(:,2) = 1.0 / sltx
            ndim = 2
        end if
    else if (model==2) then
        if (origin==0) then
            ! y=a*(1-exp(-b*x))+c*x.
            coefMat(:,1) = (1.0-exp(-b1*dose)) / sltx
            coefMat(:,2) = dose / sltx
            ndim = 2
        else if (origin==1) then
            ! y=a*(1-exp(-b*x))+c*x+d.
            coefMat(:,1) = (1.0-exp(-b1*dose)) / sltx
            coefMat(:,2) = dose / sltx
            coefMat(:,3) = 1.0 / sltx
            ndim = 3
        end if
    else if (model==3) then
        if (origin==0) then
            ! y=a*(1-exp(-b*x))+c*(1-exp(-d*x)).
            coefMat(:,1) = (1.0-exp(-b1*dose)) / sltx
            coefMat(:,2) = (1.0-exp(-b2*dose)) / sltx
            ndim = 2
        else if (origin==1) then
            ! y=a*(1-exp(-b*x))+c*(1-exp(-d*x))+e.
            coefMat(:,1) = (1.0-exp(-b1*dose)) / sltx
            coefMat(:,2) = (1.0-exp(-b2*dose)) / sltx
            coefMat(:,3) = 1.0 / sltx
            ndim = 3
        end if
    end if
    !
    ltxMat(:,1) = ltx / sltx
    allocate(aMat(ndim,ndim),bMat(ndim,1),stat=alc) 
    if (alc/=0) then
        info=1
        return
    end if
    !
    aMat = matmul(transpose(coefMat(:,1:ndim)),coefMat(:,1:ndim))
    bMat = matmul(transpose(coefMat(:,1:ndim)),ltxMat)
    !
    call gjordan(aMat,bMat,ndim,1,singular)
    outp(1:ndim) = bMat(1:ndim,1)
    deallocate(aMat,bMat,stat=alc)
    if (any(outp(1:ndim-origin)<=0.0)) then
        info = 1
        return
    end if
    !
    if (singular/=0 .or. alc/=0) then
        info=1
        return
    end if
    !
    return
end subroutine inipars
