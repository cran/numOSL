subroutine fmmED(ed,sed,ndat,ncomp,addsigma,&
                 pars,stdp,maxlik,bic,info)
!----------------------------------------------------------------
! Subroutine fmmED() is used for estimating 
! parameters of a finite mixture model.
!----------------------------------------------------------------
!      ed(ndat):: input, real values, the logged ED values.
!     sed(ndat):: input, real values, relative errors of EDs.
!          ndat:: input, integer, number of data points.
!         ncomp:: integer, the number of components.
!      addsigma:: input, real value, additional uncertainty.
! pars(2,ncomp):: input/output, the patameters.
! stdp(2,ncomp):: output, rela values, the std of parameters.
!        maxlik:: output, real value, logged maximum likelihood.
!           bic:: output, real value, the BIC value.
!          info:: output, integer, 0=success, 1=fail.
!----------------------------------------------------------------
! Author:: Peng Jun, 2014.09.15.
!----------------------------------------------------------------
! Dependence:: subroutine apfmmstd.------------------------------
!----------------------------------------------------------------
! Reference:: Galbraith, RF, Green PF, 1990. Estimating the 
!             component ages in a finite mixture. Nuclear Tracks
!             and Radiation Measurements, 17, page 197-206.
!----------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, ncomp
    real   (kind=8), intent(in):: ed(ndat), sed(ndat), addsigma
    real   (kind=8), intent(inout):: pars(2,ncomp)
    real   (kind=8), intent(out):: stdp(2,ncomp), maxlik, bic
    integer(kind=4), intent(out):: info
    ! Local variables.
    real   (kind=8), parameter:: eps=1.0D-8,&
    PI=3.141592653589793238462643383279502884197D+00
    integer(kind=4), parameter:: maxiter=5000
    real   (kind=8):: w(ndat), pf(ndat,ncomp), rsumpf(ndat),&
                      p(ndat,ncomp), pp(ndat,ncomp), wp(ndat,ncomp),&
                      newp(2,ncomp)
    integer(kind=4):: i, j, message
    !
    stdp = -99.0
    info = 0
    !
    w = 1.0/(addsigma**2+sed**2)
    !
    loopA: do i=1, maxiter
        do j=1, ncomp
            pf(:,j) = pars(1,j)*sqrt(w)*&
                      exp(-0.5*w*(ed-pars(2,j))**2)
        end do
        !
        rsumpf = sum(pf, dim=2)
        !
        do j=1, ndat
            pp(j,:) = pf(j,:) / rsumpf(j)
            wp(j,:) = w(j)*pp(j,:)
            p(j,:) = pp(j,:)*w(j)*ed(j)
        end do
        !
        newp(1,:) = sum(pp,dim=1)/real(ndat)
        newp(2,:) = sum(p,dim=1)/sum(wp,dim=1)
        !
        if (sum(abs(pars(1,:)-newp(1,:)))+&
            sum(abs(pars(2,:)-newp(2,:)))<=eps) then
            exit loopA
        else 
            pars(1,:) = newp(1,:)
            pars(2,:) = newp(2,:)
        end if
    end do loopA
    !
    maxlik = sum(log(1.0/sqrt(2.0*PI)*sum(pf,dim=2)))
    bic = -2.0*maxlik + (2.0*real(ncomp)-1.0)*log(real(ndat))
    !
    if ( (maxlik .ne. maxlik) .or.&
         (maxlik+1.0==maxlik) ) then
        info = 1
        return
    end if
    call apfmmstd(pars,ncomp,ed,sed,ndat,&
                  addsigma,stdp,message)
    if (message==0) then
        pars(2,:) = exp(pars(2,:))
        stdp(2,:) = stdp(2,:)*pars(2,:)
    else 
        info = 1
    end if
    !
    return
end subroutine fmmED    
