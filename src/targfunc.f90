subroutine targfunc(lamda,n1,tim,sig,ntim,typ,&
                    addc,ithn,constant,fvalue,iflag)
!---------------------------------------------------------------------------------
! Subroutine targfunc() is used for calculating the number of trapped electrons
! with given decay rates using a Linear Algebra method (Bluszcz A, 1996).
!---------------------------------------------------------------------------------
! lamda(n1):: input, real values, the given decay rates.
!        n1:: input, integer, dimension of the problem.
! tim(ntim):: input, real vlaues, time values.
! sig(ntim):: input, real values, stimulated decay signal values.
!      ntim:: input, integer, number of data points.
!       typ:: input, integer, type of decay OSL, 1 means CW-OSL, 2 means LM-OSL.
!      addc:: input, integer, do a background subtraction (1) or not (0).
!  ithn(n1):: output, real vlaues, calculated number of trapped electrons.
!  constant:: output, real vlaue, estimated constant, constant=0 if addc=0.
!    fvalue:: output, real value, calculated sum of squared residuals.
!     iflag:: output, integer, 0 means a successful work, 1 means
!                     any(ithn<0) or constant<0 or singular model.
!--------------------------------------------------------------------------------
! Author:: Peng Jun, 2014.08.29.
!--------------------------------------------------------------------------------
! Dependence:: subroutine gjordan().---------------------------------------------
!--------------------------------------------------------------------------------
! Reference:: Bluszcz, A., 1996. Exponential function fitting to TL growth data 
!             and similar applications. Geochronometria 13, 135–141.
!--------------------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n1, ntim, typ, addc
    real   (kind=8), intent(in):: lamda(n1), tim(ntim), sig(ntim)
    real   (kind=8), intent(out):: ithn(n1), constant, fvalue
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    real   (kind=8):: coefMat1(ntim,n1+1), sigMat(ntim,1), timmaxt(ntim),&
                      bMat(n1,1), bMat1(n1+1,1), tim2maxt2(ntim),&
                      aMat(n1,n1), aMat1(n1+1,n1+1), maxt
    integer(kind=4):: i, singular
    !
    iflag = 0
    fvalue = 1.0D+20
    !
    if (typ==1) then
        do i=1, n1
            coefMat1(:,i) = lamda(i)*&
                            exp(-lamda(i)*tim)
        end do
        coefMat1(:,n1+1) = 1.0
    else if (typ==2) then
        maxt = maxval(tim)
        timmaxt = tim/maxt
        tim2maxt2 = tim**2/maxt/2.0
        do i=1, n1
            coefMat1(:,i) = timmaxt*lamda(i)*&
                            exp(-lamda(i)*tim2maxt2)
        end do
        coefMat1(:,n1+1) = timmaxt
    end if
    !
    sigMat(:,1) = sig
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
    fvalue = sum((sig-sum(coefMat1,dim=2))**2)
    !
    return
end subroutine targfunc   
