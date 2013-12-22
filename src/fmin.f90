function fmin(ax,bx,tol,ltx,cpars,npars)
!------------------------------------------------------------------------------------------------
! Approximation x to the point where a given function attains a minimum value on the interval 
! (ax,bx)  is determined. The method used is a combination of  golden  section  search  and
! successive parabolic interpolation. 
! ===============================================================================================
!
! ax,       input:: real value, left endpoint of initial interval.
!
! bx,       input:: real value, right endpoint of initial interval.
!
! tol,      input:: real value, desired length of the interval of uncertainty of the final result.
!
! ltx,      input:: real value, standardlised OSL signal from which equivalent dose to be calculated.
!
! cpars(4), input:: real values, characteristic parameters of the dose-response curve.
!
! npars,    input:: integer, dimension of the fitting model.
!
! fmin,    output:: real value, calculated equivalent dose correspond to ltx.
! ================================================================================================
! Author:: Peng Jun, 2013.06.22.
!
! Dependence:: inner function line; inner function exper; inner function linexp.
!
! Reference:: http://www.netlib.org/fmm/fmin.f
!-------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4)::npars
  real   (kind=8)::ax,bx,tol
  real   (kind=8),dimension(4)::cpars
  real   (kind=8)::ltx
  real   (kind=8)::fmin
  ! Local variables
  real   (kind=8):: a,b,c,d,e,eps,xm,p,q,&
                    r,tol1,tol2,u,v,w
  real   (kind=8):: fu,fv,fw,fx,x
  !
  ! Squared inverse of the golden ratio
  c=0.5D+00*(3.0D+00-dsqrt(5.0d+00))
  !
  ! eps is approximately the square root of 
  ! the relative machine precision
  eps=1.0D+00
  10 eps=eps/2.0D+00
     tol1=1.0D+00+eps
     if(tol1 .gt. 1.0D+00) goto 10
     eps=dsqrt(eps)
  !
  ! Initialization
  a=ax 
  b=bx
  v=a+c*(b-a)
  w=v
  x=v
  e=0.0D+00
  !
  ! Calculate a model dependent fx
  if(npars==2)  then
    fx=line(x)
  else if(npars==3) then
    fx=exper(x)
  else if(npars==4) then
    fx=linexp(x)
  end if
  fv=fx
  fw=fx
  !
  ! Main iterations
  20 xm=0.5D+00*(a+b)
     tol1=eps*dabs(x)+tol/3.0D+00
     tol2=2.0D+00*tol1
  !
  ! Check converge (90)
  if(dabs(x-xm) .le. (tol2-0.5D+00*(b-a)) ) goto 90
  !
  ! Check if golden-section is necessary (40)
  if(dabs(e).le. tol1) goto 40
  !
  ! Fit parabola
  r=(x-w)*(fx-fv)
  q=(x-v)*(fx-fw)
  p=(x-v)*q-(x-w)*r
  q=2.0D+00*(q-r)
  if(q .gt. 0.0D+00) p=-p
  q=dabs(q)
  r=e
  e=d
  !
  ! Check if parabola can be accepted
  30 if(dabs(p) .ge. dabs(0.5D+00*q*r) ) goto 40
     if(p .le. q*(a-x)) goto 40
     if(p .ge. q*(b-x)) goto 40
  !
  ! Parabolic interpolation
  d=p/q
  u=x+d
  !
  ! Check if targeted function value f(x) 
  ! is too close to ax or bx
  if( (u-a) .lt. tol2) d=dsign(tol1,xm-x)
  if( (b-u) .lt. tol2) d=dsign(tol1,xm-x)
  goto 50
  !
  ! Golden-section step
  40 if(x .ge. xm) e=a-x
     if(x .lt. xm) e=b-x
     d=c*e
  !
  ! Check if f(x) is too close to x
  50 if(dabs(d) .ge. tol1) u=x+d
     if(dabs(d) .lt. tol1) u=x+dsign(tol1,d)
  ! Calculate a model dependent fu
  if(npars==2) then
    fu=line(u)
  else if(npars==3) then
    fu=exper(u)
  else if(npars==4) then
    fu=linexp(u)
  end if
  !
  ! Update a, v, v, w, x
  if(fu .gt. fx) goto 60
  if(u .ge. x) a=x
  if(u .lt. x) b=x
  v=w
  fv=fw
  w=x
  fw=fx
  x=u
  fx=fu
  goto 20
  !
  60 if(u .lt. x) a=u
     if(u .ge. x) b=u
     if(fu .le. fw) goto 70
     if(w .eq. x) goto 70
     if(fu .le. fv) goto 80
     if(v .eq. x) goto 80
     if(v .eq. w) goto 80
     goto 20
  !
  70 v=w
     fv=fw
     w=u
     fw=fu
     goto 20
  80 v=u
     fv=fu
     goto 20
  !
  ! END
  90 fmin=x
  return
  !
  contains
  ! ***************
  ! 
  ! 1) Linear function y=a*x+b
  function line(x)
    implicit none
    real(kind=8)::line,x
    line=(cpars(1)*x+cpars(2)-ltx)**2
    return
  end function line
  ! 2) Exponential function y=a*(1-exp(-b*x))+c
  function exper(x)
    implicit none
    real(kind=8)::exper,x
    exper=(cpars(1)*(1.0D+00-dexp(-cpars(2)*x))+cpars(3)-ltx)**2
    return
  end function exper
  ! 3) Exponential plus linear function y=a*(1-exp(-b*x))+c*x+d
  function linexp(x)
    implicit none
    real(kind=8):: linexp, x
    linexp=(cpars(1)*(1.0D+00-dexp(-cpars(2)*x))+cpars(3)*x+cpars(4)-ltx)**2
    return
  end function linexp
  !
  ! ****************
end function fmin
