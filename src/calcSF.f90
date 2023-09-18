subroutine calcSF(ax,bx,dose,ltx,pars,&
                  nd,n2,model,SF,fmin)
!-----------------------------------------------------
! Subroutine calcSF() is used for 
! calculating the scaling factor.
!-----------------------------------------------------
!       ax:: input, real value, lower limit.
!       bx:: input, real value, upper limit.
! dose(nd):: input, real values, regenerative doses.
!  ltx(nd):: input, real value, OSL signal values.
! pars(n2):: input, real vlaues, parameters.
!       nd:: input, integer, number of data points.
!       n2:: input, integer, number ([2,4]) of pars.
!    model:: input, integer, 0=line, 1=exp, 2=line+exp,
!                            3=dexp, 7=gok
!       SF:: output, real value, calculated SF.
!     fmin:: output, real value, minimized objective.
!-----------------------------------------------------
! Author:: Peng Jun, 2023.08.30.
!-----------------------------------------------------
! Dependence:: inner function fcn.
!-----------------------------------------------------
! NOTE: THIS ROUTINE IS BASED ON THE FREE FORTRAN 77 
!       CODE AT: http://www.netlib.org/fmm/fmin.f .
! Part of the R package, http://www.R-project.org .
!-----------------------------------------------------
    implicit none
    ! Arguments.
    integer, intent(in):: nd, n2, model
    real(kind(1.0d0)), intent(in):: ax, bx, dose(nd),&
                                    ltx(nd), pars(n2)
    real(kind(1.0d0)), intent(out):: SF, fmin
    ! Local variables.
    real(kind(1.0d0)):: a, b, c, d, e, p, q, r, u, v, w, x
    real(kind(1.0d0)):: t2, fu, fv, fw, fx, xm, eps, tol1, tol3
    
    real(kind(1.0d0)):: xx(5)
    real(kind(1.0d0)), parameter:: tol=1.490116D-08
    !
    c = (3.0D+00 - dsqrt(5.0D+00)) * 0.5D+00
    eps = EPSILON(0.0D+00) 
    tol1 = eps + 1.0D+00
    eps = dsqrt(eps)

    a = ax
    b= bx
    v = a + c * (b - a)
    w = v
    x = v

    d = 0.0D+00
    e = 0.0D+00
    fx = fcn(x)
    fv = fx
    fw = fx
    tol3 = tol/3.0D+00

    do 
        xm = (a + b) * 0.5D+00
        tol1 = eps * dabs(x) + tol3
        t2 = tol1 * 2.0D+00

        if ( dabs(x - xm) <= t2 - (b - a) * 0.5D+00 ) exit

        p = 0.0D+00
        q = 0.0D+00
        r = 0.0D+00
         
        if (dabs(e) > tol1) then
            r = (x - w) * (fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0D+00

            if (q > 0.0D+00) then
                p = -p
            else 
                q = -q
            end if
            r = e
            e = d
        end if
        

        if ((dabs(p) >= dabs(q * 0.5D+00 * r)) .or. &
            (p <= q * (a - x)) .or. (p >= q * (b - x))) then
            if (x < xm) then
                e = b - x
            else 
                e = a - x
            end if
            d = c * e
        else 
            d = p / q
            u = x + d
            
            if ((u - a < t2) .or. (b - u < t2)) then
                d = tol1
                if (x >= xm) d = -d
            end if
        end if

        if (dabs(d) >= tol1) then
            u = x + d
        else if (d > 0.0D+00) then
            u = x + tol1
        else 
            u = x - tol1
        end if

        fu = fcn(u)

        if (fu <= fx) then
            if (u < x) then
                b = x
            else 
                a = x
            end if

            v = w; w = x; x = u
            fv = fw; fw = fx; fx = fu

        else 
            if (u < x) then
                a = u
            else 
                b = u
            end if
                 
            if ((fu <= fw) .or. (w == x)) then
                v = w; fv = fw
                w = u; fw = fu
            else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
                v = u; fv = fu
            end if
        end if

    end do
    
    SF = x
    fmin = fcn(x)
    
    return 

    contains
    ! Inner function fcn.
    function fcn(x)
        implicit none
        real(kind(1.0d0)):: fcn, x
        !
        xx = 0.0D+00 
        xx(1:n2) = pars
        !
        if (model==0) then 
            fcn = sum((xx(1)*dose+xx(2)-x*ltx)**2)
        else if (model==1) then
            fcn = sum((xx(1)*(1.0-exp(-xx(2)*dose))+&
                       xx(3)-x*ltx)**2)
        else if (model==2) then
            fcn = sum((xx(1)*(1.0-exp(-xx(2)*dose))+&
                       xx(3)*dose+xx(4)-x*ltx)**2)
        else if (model==3) then
            fcn = sum((xx(1)*(1.0-exp(-xx(2)*dose))+&
                       xx(3)*(1.0-exp(-xx(4)*dose))+&
                       xx(5)-x*ltx)**2)
        else if (model==7) then                     
            fcn = sum((xx(1)*(1.0-(1.0+xx(2)*xx(3)*dose)**&
                      (-1.0/xx(3)))+xx(4)-x*ltx)**2)
        end if
        !
        ! Test for Inf or NaN.
        if ((fcn .ne. fcn) .or. (fcn .gt. huge(0.0D+00))) fcn = huge(0.0D+00)
        !
        return
    end function fcn
end subroutine calcSF
