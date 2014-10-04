subroutine diffev(tim,sig,wght,ntim,np,f,cr,maxiter,tol,typ,addc,&
                  n1,lamda,ithn,constant,agents,fvec1,fmin,iflag)
!------------------------------------------------------------------------------
! Subroutine diffev() is used for estimating parameters of an 
! OSL decay curve by the differential evolution method.
!------------------------------------------------------------------------------
!     tim(ntim):: input, real values, time values.
!     sig(ntim):: input, real values, decay signal values.
!    wght(ntim):: input, real values, weigths od signal values.
!          ntim:: input, integer, number of data points.
!            np:: input, integer, number of differential populations (np>=4).
!             f:: input, real value, differential weight.
!            cr:: input, real value, crossover probability.
!       maxiter:: input, integer, maximum number of iterations.
!           tol:: input, real value, tolerance for stopping the iteration.
!           typ:: input, integer, type of decay curve, 1=CW-OSL, 2=LM-OSL.
!          addc:: input, integer, 0=non-constant, 1=constant.
!            n1:: input, integer, dimension of the problem (>=1).
!     lamda(n1):: output, real values, estimated decay rates.
!      ithn(n1):: output, real values, estimated number of trapped electrons.
!      constant:: output, real value, estimated constant background.
! agents(np,n1):: output, real values, differential populations.
!   fvec1(ntim):: output, real values, predicted decay signal values.
!          fmin:: output, real value, best (minimum) objective.
!         iflag:: output, integer, 0=success, 1=fail.
!------------------------------------------------------------------------------
! Author:: Peng Jun, 2014.09.28.
!------------------------------------------------------------------------------
! Dependence:: subroutine leaveOne; subroutine sample;-------------------------
!              subroutine targfunc; subroutine hpsort; subroutine relSD.-------
!------------------------------------------------------------------------------
! Reference:: Bluszcz, A., Adamiec, G., 2006. Application of 
!             differential evolution  to fitting OSL decay 
!             curves. Radiation Measurements 41, 886-891.
! Differential evolution, http://en.wikipedia.org/wiki/Differential_evolution
!------------------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ntim, np, maxiter, typ, addc, n1
    real   (kind=8), intent(in):: tim(ntim), sig(ntim), wght(ntim),&
                                  f, cr, tol
    real   (kind=8), intent(out):: lamda(n1), ithn(n1), constant,&
                                   agents(np,n1), fvec1(ntim), fmin                              
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    real   (kind=8), dimension(2,7), parameter::&
              tSpace=reshape((/1.0D+00,10.0D+00,&        
                               1.0D-01, 5.0D+00,&
                               1.0D-02, 3.0D+00,&
                               1.0D-03, 1.0D+00,&  
                               1.0D-04, 0.5D+00,&
                               1.0D-05, 0.1D+00,&
                               1.0D-06, 0.05D+00/), (/2,7/) )    
    real   (kind=8):: ran(n1), space(2,n1), ran1(1), yy(n1),&
                      fagents(np), agents3(3,n1), constants(np),&
                      rsd(n1), ithn1(n1), ithn2(n1), ithns(np,n1)
    integer(kind=4):: npsamp(np), n1samp(n1), npsamp1(np-1),&
                      samp3(3), order(n1), indx(1)
    integer(kind=4):: i, j, k, iflag1, iflag2
    real   (kind=8):: fmin1, fmin2, minValue,&
                      constant1, constant2, maxt
    !
    space = log(tSpace(:,1:n1))
    !
    do i=1, np
        call random_number(ran)
        agents(i,:) = exp(space(1,:) + ran*(space(2,:)-space(1,:)))
    end do
    !
    npsamp = (/(i,i=1,np)/)
    n1samp = (/(i,i=1,n1)/)
    !
    fagents = 1.1D+20
    minValue = 1.0D+20
    !
    loopA: do i=1, maxiter
        do j=1, np
            call leaveOne(npsamp,np,j,npsamp1)
            !
            call sample(npsamp1,np-1,3,samp3)
            !
            agents3 = agents(samp3,:)
            !
            call sample(n1samp,n1,1,indx)
            ! 
            do k=1, n1
                call random_number(ran1)
                ! 
                if (ran1(1)<cr .or. indx(1)==k) then
                    yy(k) = agents3(1,k) + f*(agents3(2,k)-agents(3,k)) 
                else 
                    yy(k) = agents(j,k)
                end if
            end do
            !
            if (all(yy>0.0)) then
                call targfunc(yy, n1, tim, sig, wght, ntim, typ, addc,&
                              ithn1, constant1, fmin1, iflag1)
                !
                call targfunc(agents(j,:), n1, tim, sig, wght, ntim, typ,&
                              addc, ithn2, constant2, fmin2, iflag2)
                if (iflag1==0 .and. fmin1<fmin2)  then
                    agents(j,:) = yy
                    ithns(j,:) =  ithn1
                    fagents(j) = fmin1
                    constants(j) = constant1
                end if
            end if
            call hpSort(agents(j,:),n1,order) 
            !
        end do
        !
        call relSD(agents,np,n1,rsd)
        if (maxval(rsd)<=tol) exit loopA
        !
    end do loopA
    !
    iflag = 1
    fmin = -99.0
    constant = -99.0 
    fvec1 = -99.0
    lamda = -99.0
    ithn = -99.0
    do i=1, np
        if (fagents(i)<minValue) then
            minValue = fagents(i)
            fmin = fagents(i)
            lamda = agents(i,:)
            ithn = ithns(i,:)
            constant = constants(i)
            iflag = 0 
        end if
    end do
    !
    if (iflag==1) return
    !   
    if (typ==1) then
        fvec1 = constant
        do i=1, n1
            fvec1 = fvec1 + ithn(i)*lamda(i)*&
                    exp(-lamda(i)*tim)
        end do
    else if (typ==2) then
        maxt = tim(ntim)
        fvec1 = constant*tim/maxt
        do i=1, n1
            fvec1 = fvec1 + ithn(i)*tim/maxt*lamda(i)*&
                    exp(-lamda(i)*tim**2/maxt/2.0)
        end do
    end if
    !
    return
end subroutine diffev
!       
subroutine sample(vec,n,k,vec1)
!--------------------------------------------------------------
! Subroutine sample() is used to draw k integers from
! integers ranging from 1 to nmax without replacement.
!--------------------------------------------------------------
!  vec(n):: input, integers, the total samples.
!       n:: input, integer, the upper integer.
!       k:: input, integer, number of integers to be sampled.
! vec1(k):: output, integers, the sampled integers.
!--------------------------------------------------------------
! Author:: Peng Jun, 2014.08.28.
!--------------------------------------------------------------
! Reference:: Bert F and Green JR, 1977. Fortran subroutines
!             for random sampling without replacement. Behavior
!             Research Methodt and Instrumentation, pp 559.
!--------------------------------------------------------------
! Dependence:: NO.---------------------------------------------
!--------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n, k, vec(n)
    integer(kind=4), intent(out):: vec1(k)
    ! Local variables.
    integer(kind=4):: m, j, ll
    real   (kind=8):: ran(1)
    !
    m = 0
    loopA: do j=1, n
        call random_number(ran)
        ll = int((float(n-j+1))*ran(1)) + 1
        if (ll>k-m) cycle loopA
        m = m + 1
        vec1(m) = vec(j)
        if (m>=k) exit loopA
    end do loopA
    !
    return
end subroutine sample
!
subroutine relSD(mat,m,n,sd)
!-----------------------------------------------------------
! Subroutine relSD() is used for calculating the
! relative standard errors of a matrix column by column. 
!-----------------------------------------------------------
! mat(m,n):: input, real values, a matrix.
!        m:: input, integer, row number of the matrix.
!        n:: input, integer, column number of the matrix.
!    sd(n):: output, real value, the resulting relative sds.
!-----------------------------------------------------------
! Author:: Peng Jun, 2014.08.29.
!-----------------------------------------------------------
! Dependence:: NO.------------------------------------------
!-----------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: m, n
    real   (kind=8), intent(in):: mat(m,n)
    real   (kind=8), intent(out):: sd(n)
    ! Local variables.
    integer(kind=4):: i
    real   (kind=8):: mean(n)
    !
    mean = sum(mat,dim=1) / real(m)
    do i=1, n
        sd(i) = sqrt(sum((mat(:,i)-mean(i))**2)/real(m-1)) / mean(i)
    end do
    !
    return
end subroutine relSD
!
subroutine leaveOne(vec,n,which,vec1)
!--------------------------------------------------------
! Subroutine leaveOne is used for removing an integer
! from current vector that contains a number of integers.
!--------------------------------------------------------
!    vec(n):: input, integers, a vector.
!         n:: input, integer, dimension of vec.
!     which:: input, integer, which one to be removed.
! vec1(n-1):: output, integers, resulting vector.
!--------------------------------------------------------
! Author:: Peng Jun, 2013.05.20.
!--------------------------------------------------------
! Dependence:: NO.---------------------------------------
!--------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n, which
    integer(kind=4), intent(in):: vec(n)
    integer(kind=4), intent(out):: vec1(n-1)
    !
    if (which==1) then
        vec1 = vec(2:n)
    else if (which==n) then
        vec1 = vec(1:n-1)
    else 
        vec1 = (/vec(1:which-1),vec(which+1:n)/)
    end if
    !
    return
end subroutine leaveOne
