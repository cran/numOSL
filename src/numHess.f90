subroutine numHess(xd,yd,syd,nd,model,&
                   pars,np,hess,iflag)
!--------------------------------------------------------
! Subroutine numHess is used for approximating the 
! Hessian matrix of a function at given parameters using 
! a numerical differention method called "Richardson".
!--------------------------------------------------------
!      xd(nd):: input, real values, observations x.
!      yd(nd):: input, real values, observations y.
!     syd(nd):: input, real values, observations sy.
!          nd:: input, integer, number of data points.
!       model:: input, integer, 0=linear model;
!                               6=MAM model.                             
!    pars(np):: input, real values, parameters.
!          np:: input, integer, number of pars.
! hess(np,np):: output, real values, the Hessian matrix.
!       iflag:: output, integer, 0=success; 1=fail.
!--------------------------------------------------------
! Author:: Peng Jun, 2016.07.06.
!--------------------------------------------------------
! Dependence:: inner function func; 
!              subroutine pnorm; 
!              function alnorm.
!--------------------------------------------------------
! Reference:: Paul Gilbert and Ravi Varadhan (2012). 
!             numDeriv: Accurate Numerical
!             Derivatives. R package version 2012.9-1. 
!             Availiable at:
!             http://CRAN.R-project.org/package=numDeriv
!--------------------------------------------------------
! NOTE: THIS SUBROUTINE IS BASED ON FUNCTION 
!       "genD" IN R PACKAGE numDeriv.
!--------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: nd, model, np
    real   (kind=8), intent(in):: xd(nd), yd(nd),& 
                                  syd(nd), pars(np)
    real   (kind=8), intent(out):: hess(np,np)
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    real   (kind=8), parameter:: eps=1.0D-04,&
                                 d=0.1D+00,&
                                 zTol=1.781029D-05
    integer(kind=4), parameter:: r=4
    real   (kind=8):: h0(np), h(np), Hdiag(np),& 
                      incr1(np), incr2(np),&
                      Dd(np*(np+3)/2), Daprox(r),& 
                      Haprox(r)
    real   (kind=8):: f0, f1, f2, m4
    integer(kind=4):: i, j, k, m, u
    real   (kind=8), parameter:: pi=&
    3.141592653589793238462643383279502884197D+00
    real   (kind=8):: xx(15)
    !
    iflag = 0
    hess = 0.0D+00
    f0 = func(pars)
    if(f0 .ne. f0 .or. f0+1.0D+00==f0) then
        iflag = 1
        return
    end if
    !
    do i=1, np
        h0(i) = dabs(d*pars(i))
        if (dabs(pars(i))<zTol) h0(i) = h0(i) + eps
    end do
    !
    Dd = 0.0D+00
    Hdiag = 0.0D+00
    Daprox = 0.0D+00
    Haprox = 0.0D+00
    !
    do i=1, np
        h = h0
        incr1 = 0.0D+00 
        do k=1, r
            incr1(i) = h(i)
            f1 = func(pars+incr1) 
            f2 = func(pars-incr1)
            Daprox(k) = (f1-f2) / (2.0D+00*h(i))
            Haprox(k) = (f1-2.0D+00*f0+f2) / (h(i))**2
            h = h / 2.0D+00
        end do
        !
        do m=1, r-1
            m4 = (4.0D+00)**m
            do k=1, r-m
                Daprox(k) = (Daprox(k+1)*m4-Daprox(k))/&
                            (m4-1.0D+00)
                Haprox(k) = (Haprox(k+1)*m4-Haprox(k))/&
                            (m4-1.0D+00)
            end do
        end do
        Dd(i) = Daprox(1)
        Hdiag(i) = Haprox(1)
    end do
    !
    u = np
    do i=1, np
        do j=1, i
            u = u + 1
            if (i==j) then 
                Dd(u) = Hdiag(i)
            else 
                h = h0
                incr1 = 0.0D+00 
                incr2 = 0.0D+00
                do k=1, r
                    incr1(i) = h(i)
                    incr2(j) = h(j)
                    f1 = func(pars+incr1+incr2)
                    f2 = func(pars-incr1-incr2)
                    Daprox(k) = (f1-2.0D+00*f0+f2-&
                                 Hdiag(i)*(h(i))**2-&
                                 Hdiag(j)*(h(j))**2)/&
                                (2.0D+00*h(i)*h(j))
                    h = h / 2.0D+00
                end do
                !
                do m=1, r-1
                    m4 = (4.0D+00)**m
                    do k=1, r-m
                        Daprox(k) = (Daprox(k+1)*m4-Daprox(k))/& 
                                    (m4-1.0D+00)
                    end do 
                end do
                Dd(u) = Daprox(1)
            end if
        end do
    end do
    ! 
    u = np
    do i=1, np
        do j=1, i
            u = u + 1
            hess(i,j) = Dd(u)
        end do
    end do
    hess = hess + transpose(hess)
    do i=1, np
        hess(i,i) = hess(i,i) / 2.0D+00
    end do
    if(any(hess .ne. hess) .or.&
       any(hess .gt. huge(0.0D+00)) .or.&
       any(hess .lt. -huge(0.0D+00)))  iflag=1
    !
    return
    !
    ! Inner function func.
    !---------------------
    contains
        function func(x)
            implicit none
            real   (kind=8):: x(np), func
            real   (kind=8):: vec(nd), alnorm
            integer(kind=4):: i
            !
            xx = 0.0D+00
            xx(1:np) = x
            !
            ! Linear model.
            if (model==0) then
                vec = (xx(1)*xd+xx(2)-yd)/syd
                func = sum(vec**2)
                return
            end if
            !
            ! Minimum age model.
            if (model==6) then
                if (np==3) then
                    vec = (xx(2)-(xx(2)/(xx(3))**2+yd/xd**2)/&
                          (1.0D+00/(xx(3))**2+1.0D+00/xd**2))*&
                           dsqrt(1.0D+00/xd**2+1.0D+00/(xx(3))**2)
                    call pnorm(vec,nd,.false.)
                    func = -sum(dlog(xx(1)/dsqrt(2.0D+00*pi*xd**2)*&
                            dexp(-(yd-xx(2))**2/(2.0D+00*xd**2))+&
                           (1.0D+00-xx(1))/dsqrt(2.0D+00*pi*(xd**2+(xx(3))**2))*&
                            dexp(-(yd-xx(2))**2/(2.0D+00*(xd**2+(xx(3))**2)))*&
                            2.0D+00*(1.0D+00-vec)))  
                else if (np==4) then
                    vec = (xx(2)-(xx(3)/(xx(4))**2+yd/xd**2)/&
                          (1.0D+00/(xx(4))**2+1.0D+00/xd**2))*&
                           dsqrt(1.0D+00/xd**2+1.0D+00/(xx(4))**2)
                    call pnorm(vec,nd,.false.)
                    func = -sum(dlog(xx(1)/dsqrt(2.0D+00*pi*xd**2)*&
                            dexp(-(yd-xx(2))**2/(2.0D+00*xd**2))+&
                           (1.0D+00-xx(1))/dsqrt(2.0D+00*pi*(xd**2+(xx(4))**2))*&
                            dexp(-(yd-xx(3))**2/(2.0D+00*(xd**2+(xx(4))**2)))*&
                           (1.0D+00-vec)/(1.0D+00-alnorm((xx(2)-xx(3))/xx(4),.false.))))   
                end if  
                return
            end if
            !
        end function func 
        !----------------      
end subroutine numHess
