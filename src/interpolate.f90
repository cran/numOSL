subroutine interpolate(ltx,pars,n2,model,&
                       ax,bx,dose,fvalue)
!------------------------------------------------------------
! Subroutine interpolate() is used for estimating a dose
! value from a given growth curve via interpolation.
!------------------------------------------------------------
!      ltx:: input, real value, the standardised Lx/Tx value.
! pars(n2):: input, real values, the parameters of the curve.
!       n2:: input, integer, the number of parameters.
!    model:: input, integer, 1=exp, 2=lexp, 3=dexp.
!       ax:: input, real value, the lower bound of dose.
!       bx:: input, real value, the upper bound of dose.
!     dose:: output, real value, the estimated dose value.
!   fvalue:: output, real value, the minimized objection.
!------------------------------------------------------------
! Author:: Peng Jun, 2014.09.10.
!------------------------------------------------------------
! Dependence:: external function fmin; inner function f.-----
!------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n2, model
    real   (kind=8), intent(in):: pars(n2), ltx, ax, bx
    real   (kind=8), intent(out):: dose, fvalue
    ! Local variables.
    real   (kind=8), parameter:: tol=1.490116e-08
    real   (kind=8):: locp(5), fmin
    !
    dose = fmin(ax,bx,f,tol)
    !
    locp = 0.0
    locp(1:n2) = pars
    !
     if (model==1) then
         fvalue = (locp(1)*(1.0-exp(-locp(2)*dose))+&
                   locp(3)-ltx)**2
     else if (model==2) then
         fvalue = (locp(1)*(1.0-exp(-locp(2)*dose))+&
                   locp(3)*dose+locp(4)-ltx)**2
     else if (model==3) then
         fvalue = (locp(1)*(1.0-exp(-locp(2)*dose))+&
                   locp(3)*(1.0-exp(-locp(4)*dose))+&
                   locp(5)-ltx)**2
    end if
    !
    return  
    contains
        ! Calculate the objection to be minimized. 
        function f(x) 
            implicit none
            real   (kind=8):: x, f, locp1(5)
            !
            locp1 = 0.0
            locp1(1:n2) = pars
            !
            if (model==1) then
                ! Exp model.
                f = (locp1(1)*(1.0-exp(-locp1(2)*x))+&
                     locp1(3)-ltx)**2
            else if (model==2) then
                ! Linear plus exp model.
                f = (locp1(1)*(1.0-exp(-locp1(2)*x))+&
                     locp1(3)*x+locp1(4)-ltx)**2
            else if (model==3) then
                ! Double exp model.
                f = (locp1(1)*(1.0-exp(-locp1(2)*x))+&
                     locp1(3)*(1.0-exp(-locp1(4)*x))+&
                     locp1(5)-ltx)**2
            end if
            !
            return
        end function f
end subroutine interpolate
