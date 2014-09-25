subroutine apmamstd(ed,sed,ndat,pars,&
                    stdp,np,iflag)
!---------------------------------------------------
! Subroutine apmamstd() is used for estimating
! the standard errors of a minimum age model 
! using finite difference approximation.
!---------------------------------------------------
! ed(ndat):: input, real values, logged EDs.
! sed(ndat:: input, real values, relative stdEDs.
!     ndat:: input, integer, number of data points.
! pars(np):: input, real vlaues, parameters.
! stdp(np):: output, estimated standard errors.
!       np:: input, integer, dimension of the model.
!    iflag:: output, integer, 0=success, 1=fail.
!---------------------------------------------------
! Author:: Peng Jun ,2014.09.18.
!---------------------------------------------------
! Dependence:: subroutine numHess;------------------
!              subroutine inverse;------------------ 
!              inner function fun34.---------------- 
!              function alnorm; subroutine pnorm.---
!---------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, np
    real   (kind=8), intent(in):: ed(ndat), sed(ndat), pars(np)
    real   (kind=8), intent(out):: stdp(np)
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    integer(kind=4):: i, errorflag, singular
    real   (kind=8):: hess(np,np), diag(np)
    real   (kind=8), parameter:: pi=&
    3.141592653589793238462643383279502884197D+00
    !
    iflag = 0
    stdp = -99.0
    !
    call numHess(fun34,pars,np,hess,errorflag)
    if (errorflag/=0) then
        iflag = 1
        return
    end if
    !
    call inverse(hess,np,singular)
    if (singular/=0) then
        iflag = 1
        return
    end if
    !
    do i=1, np
        diag(i) = hess(i,i)
    end do
    if (any(diag<0.0)) then
        iflag = 1
        return
    end if
    !
    stdp = sqrt(diag)
    !
    return
    ! 
    contains
        !
        function fun34(x)
            implicit none
            real   (kind=8):: x(np), x1(4), fun34,&
                              alnorm, pnorm1(ndat)
            logical,parameter:: upper=.false.   
            !
            x1 = 0.0
            x1(1:np) = x
            !
            if (np==3)    then
                !
                pnorm1 = (x1(2)-(x1(2)/x1(3)**2+ed/sed**2)/&
                         (1.0/x1(3)**2+1.0/sed**2))*&
                          sqrt(1.0/sed**2+1.0/x1(3)**2)
                !
                call pnorm(pnorm1,ndat,upper)
                !
                fun34 = -sum(log(x1(1)/sqrt(2.0*pi*sed**2)*&
                             exp(-(ed-x1(2))**2/(2.0*sed**2))+&
	                    (1.0-x1(1))/sqrt(2.0*pi*(sed**2+x1(3)**2))*&
                             exp(-(ed-x1(2))**2/(2.0*(sed**2+x1(3)**2)))*&
	                     2.0*(1.0-pnorm1)))  
            else if (np==4)  then
                !
                pnorm1 = (x1(2)-(x1(3)/x1(4)**2+ed/sed**2)/&
                         (1.0/x1(4)**2+1.0/sed**2))*&
                          sqrt(1.0/sed**2+1.0/x1(4)**2)
                !
                call pnorm(pnorm1,ndat,upper)	 
                !
                fun34 = -sum(log(x1(1)/sqrt(2.0*pi*sed**2)*&
                             exp(-(ed-x1(2))**2/(2.0*sed**2))+&
	                    (1.0-x1(1))/sqrt(2.0*pi*(sed**2+x1(4)**2))*&
                             exp(-(ed-x1(3))**2/(2.0*(sed**2+x1(4)**2)))*&
	                    (1.0-pnorm1)/(1.0-alnorm((x1(2)-x1(3))/x1(4),upper))))   
                ! 
            end if
            !
            return
        end function fun34   
end subroutine apmamstd      
