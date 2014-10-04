subroutine interpolate(ax,bx,ltx,pars,&
                       n2,model,dose,fmin)
!-----------------------------------------------------
! Subroutine interpolate() is used for 
! interpolating an equivalent dose.
!-----------------------------------------------------
!       ax:: input, real value, lower limit.
!       bx:: input, real value, upper limit.
!      ltx:: input, real value, OSL value.
! pars(n2):: input, real vlaues, parameters.
!       n2:: input, integer, number ([2,5]) of pars.
!    model:: input, integer, 1="exp",2="lexp",3="dexp".
!     dose:: output, real value, resulting value.
!     fmin:: output, real value, minimized objective.
!-----------------------------------------------------
! Author:: Peng Jun, 2014.09.26.
!-----------------------------------------------------
! Dependence:: inner function fcn.
!-----------------------------------------------------
! NOTE: THIS ROUTINE IS BASED ON THE FREE FORTRAN 77 
!       CODE AT: http://www.netlib.org/fmm/fmin.f .
!-----------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n2, model
    real   (kind=8), intent(in):: ax, bx, ltx,&
                                  pars(n2)
    real   (kind=8), intent(out):: dose, fmin
    ! Local variables.
    real   (kind=8):: a, b, c, d, e, eps,&
                      xm, p, q, r, tol1,&
                      tol2, u, v, w
    real   (kind=8):: fu, fv, fw, fx, x
    real   (kind=8):: dabs, dsqrt, dsign
    real   (kind=8):: xx(5)
    real   (kind=8), parameter:: tol=1.490116e-08
    !
        c = 0.5D+00*(3.0D+00-dsqrt(5.0D+00))
    !
        eps = 1.0D+00
    100 eps = eps/2.0D+00
        tol1 = 1.0D+00 + eps
        if (tol1 .gt. 1.0D+00) goto 100
        eps = dsqrt(eps)
    !
        a = ax
        b = bx
        v = a + c*(b-a)
        w = v
        x = v
        e = 0.0D+00
        fx = fcn(x)
        fv = fx
        fw = fx
    !
    200 xm = 0.5D+00*(a+b)
        tol1 = eps*dabs(x) + tol/3.0D+00
        tol2 = 2.0D+00*tol1
    !
        if (dabs(x-xm) .le. (tol2-0.5D+00*(b-a))) goto 900
    !
        if (dabs(e) .le. tol1) goto 400
    !
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        p = (x-v)*q - (x-w)*r
        q = 2.0D+00*(q-r)
        if (q .gt. 0.0D+00) p=-p
        q = dabs(q)
        r = e
        e = d
    !
        if (dabs(p) .ge. dabs(0.5D+00*q*r)) goto 400
        if (p .le. q*(a-x)) goto 400
        if (p .ge. q*(b-x)) goto 400
    !
        d = p/q
        u = x+d
    !
        if ((u-a) .lt. tol2) d = dsign(tol1,xm-x)
        if ((b-u) .lt. tol2) d = dsign(tol1,xm-x)
        goto 500
    !
    400 if (x .ge. xm) e = a-x
        if (x .lt. xm) e = b-x
        d = c*e
    !
    500 if (dabs(d) .ge. tol1) u = x+d
        if (dabs(d) .lt. tol1) u = x+dsign(tol1,d)
        fu = fcn(u)
    !
        if (fu .gt. fx) goto 600
        if (u .ge. x) a = x
        if (u .lt. x) b = x
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
        goto 200
    !
    600 if (u .lt. x) a = u
        if (u .ge. x) b = u
        if (fu .le. fw) goto 700
        if (w .eq. x) goto 700
        if (fu .le. fv) goto 800
        if (v .eq. x) goto 800
        if (v .eq. w) goto 800
        goto 200
    !
    700 v = w
        fv = fw
        w = u
        fw = fu
        goto 200
    !
    800 v = u
        fv = fu
        goto 200
    !
    900 dose = x
        fmin = fcn(x)
    !
        return
    !
    contains
    ! Inner function fcn.
    function fcn(x)
        implicit none
        real   (kind=8):: fcn, x
        !
        xx = 0.0D+00 
        xx(1:n2) = pars
        !
        if (model==1) then
            ! Exp model (n2 = 2 or 3)
            fcn = (xx(1)*(1.0-dexp(-xx(2)*x))+&
                   xx(3)-ltx)**2
        else if (model==2) then
            ! Linear plus exp model (n2 = 3 or 4).
            fcn = (xx(1)*(1.0-dexp(-xx(2)*x))+&
                   xx(3)*x+xx(4)-ltx)**2
        else if (model==3) then
            ! Double exp model (n2 = 4 or 5).
            fcn = (xx(1)*(1.0-dexp(-xx(2)*x))+&
                   xx(3)*(1.0-dexp(-xx(4)*x))+&
                   xx(5)-ltx)**2
        end if
        !
        return
    end function fcn
end subroutine interpolate
