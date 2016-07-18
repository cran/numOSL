subroutine inipars(b1,b2,model,n2,dose,&
                   ltx,sltx,ndat,outp,info)
!--------------------------------------------------------
! Subroutine inipars is used to estimate other
! parameters of a growth curve using fixed b values.
!--------------------------------------------------------
!         b1:: input, real value, first b value.
!         b2:: input, real value, second b value.
!      model:: input, integer, model for estimation,
!              1=exp, 2=lexp, 3=dexp, 7=gok.
!         n2:: input, integer, number of pars (>=2).
! dose(ndat):: input, real values, dose values.
!  ltx(ndat):: input, real values, Lx/Tx values.
! sltx(ndat):: input, real values, std of Lx/Tx values.
!       ndat:: input, integer, number of data points.
!    outp(3):: output, real values, estimated parameters.
!       info:: output, integer, 0=success, 1=fail.
!--------------------------------------------------------
! Author:: Peng Jun, 2016.06.29.
!--------------------------------------------------------
! Dependence:: subroutine gjordan.
!--------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: model, n2, ndat
    real   (kind=8), intent(in):: b1, b2, dose(ndat),&
                                  ltx(ndat), sltx(ndat)
    real   (kind=8), intent(out):: outp(3)
    integer(kind=4), intent(out):: info
    ! Local variables.
    real   (kind=8):: coefMat(ndat,3), ltxMat(ndat,1)
    real   (kind=8), allocatable:: aMat(:,:), bMat(:,:)
    integer(kind=4):: ndim, alc, singular, origin
    ! 
    outp = 0.0
    info = 0
    !
    if (model==1) then
        if (n2==2) then
            ! y=a*(1-exp(-b*x)).
            coefMat(:,1) = (1.0-exp(-b1*dose))/sltx
            ndim = 1
            origin = 0
        else if (n2==3) then
            ! y=a*(1-exp(-b*x))+c.
            coefMat(:,1) = (1.0-exp(-b1*dose))/sltx
            coefMat(:,2) = 1.0/sltx
            ndim = 2
            origin = 1
        end if
    else if (model==2) then
        if (n2==3) then
            ! y=a*(1-exp(-b*x))+c*x.
            coefMat(:,1) = (1.0-exp(-b1*dose))/sltx 
            coefMat(:,2) = dose/sltx 
            ndim = 2
            origin = 0
        else if (n2==4) then
            ! y=a*(1-exp(-b*x))+c*x+d.
            coefMat(:,1) = (1.0-exp(-b1*dose))/sltx 
            coefMat(:,2) = dose/sltx 
            coefMat(:,3) = 1.0/sltx 
            ndim = 3
            origin = 1
        end if
    else if (model==3) then
        if (n2==4) then
            ! y=a*(1-exp(-b*x))+c*(1-exp(-d*x)).
            coefMat(:,1) = (1.0-exp(-b1*dose))/sltx 
            coefMat(:,2) = (1.0-exp(-b2*dose))/sltx 
            ndim = 2
            origin = 0
        else if (n2==5) then
            ! y=a*(1-exp(-b*x))+c*(1-exp(-d*x))+e.
            coefMat(:,1) = (1.0-exp(-b1*dose))/sltx 
            coefMat(:,2) = (1.0-exp(-b2*dose))/sltx 
            coefMat(:,3) = 1.0/sltx 
            ndim = 3
            origin = 1
        end if
    else if (model==7) then
        if (n2==3) then
            ! y=a*(1-(1+b*c*x)^(-1/c)).
            coefMat(:,1) = (1.0-(1.0+b1*b2*dose)**(-1.0/b2))/sltx 
            ndim = 1
            origin = 0
        else if (n2==4) then
            ! y=a*(1-(1+b*c*x)^(-1/c))+d.
            coefMat(:,1) = (1.0-(1.0+b1*b2*dose)**(-1.0/b2))/sltx
            coefMat(:,2) = 1.0/sltx 
            ndim = 2
            origin = 1
        end if
    end if
    !
    ltxMat(:,1) = ltx/sltx 
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
