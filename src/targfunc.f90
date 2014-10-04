subroutine targfunc(lamda,n1,tim,sig,wght,ntim,typ,&
                    addc,ithn,constant,fmin,iflag)
!----------------------------------------------------------------
! Subroutine targfunc() is used for calculating the
! number of trapped electrons with given decay rates
! using a Linear Algebra method (Bluszcz A, 1996).
!----------------------------------------------------------------
!  lamda(n1):: input, real values, decay rates.
!         n1:: input, integer, number of decay rates (>=1).
!  tim(ntim):: input, real vlaues, time values.
!  sig(ntim):: input, real values, decay signal values.
! wght(ntim):: input, real values, weights of signal values.
!       ntim:: input, integer, number of data points.
!        typ:: input, integer, type of OSL, 1=CW-OSL, 2=LM-OSL.
!       addc:: input, integer, 0=non-constant, 1=constant.
!   ithn(n1):: output, real vlaues, number of trapped electrons.
!   constant:: output, real vlaue, constant, constant=0 if addc=0.
!       fmin:: output, real value, the objective value.
!      iflag:: output, integer, 0=success, 1=fail.
!----------------------------------------------------------------
! Author:: Peng Jun, 2014.09.28.
!----------------------------------------------------------------
! Dependence:: subroutine gjordan().-----------------------------
!----------------------------------------------------------------
! Reference:: Bluszcz, A., 1996. Exponential function fitting to 
!             TL growth data and similar applications.
!             Geochronometria 13, 135–141.
!----------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n1, ntim, typ, addc
    real   (kind=8), intent(in):: lamda(n1), tim(ntim),& 
                                  sig(ntim), wght(ntim)
    real   (kind=8), intent(out):: ithn(n1), constant, fmin
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    real   (kind=8):: coefMat1(ntim,n1+1), sigMat(ntim,1),& 
                      timmaxt(ntim),tim2maxt2(ntim),&
                      bMat(n1,1), bMat1(n1+1,1),& 
                      aMat(n1,n1), aMat1(n1+1,n1+1)
    integer(kind=4):: i, singular
    !
    iflag = 0
    fmin = 1.0D+20
    !
    if (typ==1) then
        do i=1, n1
            coefMat1(:,i) = lamda(i)*&
                            exp(-lamda(i)*tim)/wght
        end do
        if (addc==1) coefMat1(:,n1+1) = 1.0/wght
    else if (typ==2) then
        timmaxt = tim/tim(ntim)
        tim2maxt2 = tim**2/tim(ntim)/2.0
        do i=1, n1
            coefMat1(:,i) = timmaxt*lamda(i)*&
                            exp(-lamda(i)*tim2maxt2)/wght
        end do
        if (addc==1) coefMat1(:,n1+1) = timmaxt/wght
    end if
    !
    sigMat(:,1) = sig/wght
    !
    if (addc==1) then
         bMat1 = matmul(transpose(coefMat1), sigMat)
         aMat1 = matmul(transpose(coefMat1), coefMat1)
         call gjordan(amat1,bmat1,n1+1,1,singular)
         ithn = bMat1(1:n1,1)
         constant = bMat1(n1+1,1)
    else if (addc==0) then
         bMat = matmul(transpose(coefMat1(:,1:n1)), sigMat)
         aMat = matmul(transpose(coefMat1(:,1:n1)), coefMat1(:,1:n1))
         call gjordan(amat,bmat,n1,1,singular)
         ithn = bMat(1:n1,1)
         constant = 0.0
    end if
    !
    if (singular/=0 .or.&
        any(ithn<=0.0) .or.&
        constant<0.0) then
        iflag = 1
        return
    end if
    !
    if (typ==1) then
        do i=1, n1
            coefMat1(:,i) = ithn(i)*lamda(i)*&
                            exp(-lamda(i)*tim)
        end do
        coefMat1(:,n1+1) = constant
    else if (typ==2) then
        do i=1, n1
            coefMat1(:,i) = ithn(i)*timmaxt*lamda(i)*&
                            exp(-lamda(i)*tim2maxt2)
        end do
        coefMat1(:,n1+1) = constant*timmaxt
    end if
    !
    fmin = sum(((sig-sum(coefMat1,dim=2))/wght)**2)
    !
    return
end subroutine targfunc   
