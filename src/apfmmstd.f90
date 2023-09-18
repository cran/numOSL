subroutine apfmmstd(pars,ncomp,ed,sed,ndat,&
                    addsigma,stdp,info)
!----------------------------------------------------------------
! Subroutine apfmmstd is used for estimating the
! standard errors of parameters in a finite mixture
! using covariance matrix approximation.
!----------------------------------------------------------------
! pars(2,ncomp):: input, real values, the parameters.
!         ncomp:: input, integer, the number of components.
!      ed(ndat):: input, real values, the logged ED values.
!     sed(ndat):: input, real values, the relative errors of EDs.
!          ndat:: input, integer, number of data points.
!      addsigma:: input, real value, the addtional uncertainty.
! stdp(2,ncomp):: real vlaues, the estimated standard errors.
!          info:: output, integer, 0=success, 1=fail.
!----------------------------------------------------------------
! Author:: Peng Jun, 2023.08.30.
!----------------------------------------------------------------
! Dependence:: subroutine inverse_ger.
!----------------------------------------------------------------
! Reference:: Galbraith RF, 1988. Graphical Display of  
!             Estimates Having Differing Standard Errors. 
!             Techno-metrics, 30, page 271-281.
!----------------------------------------------------------------
    implicit none
    ! Arguments.
    integer, intent(in):: ncomp, ndat
    real(kind(1.0d0)), intent(in):: pars(2,ncomp), ed(ndat),&  
                                    sed(ndat), addsigma
    real(kind(1.0d0)), intent(out):: stdp(2,ncomp)
    integer, intent(out):: info
    ! Local variables.
    real(kind(1.0d0)):: z(ndat), s(ndat), w(ndat),&
                        aaMat(ndat,ncomp), bbMat(ndat,ncomp),& 
                        aMat(ncomp-1,ncomp-1), bMat(ncomp-1,ncomp),&
                        cMat(ncomp,ncomp), pdf(ndat,ncomp), rsumpdf(ndat),&
                        prob(ndat,ncomp), iMat(ncomp,ncomp),scnp1,&
                        covar(2*ncomp-1,2*ncomp-1), diag(2*ncomp-1)
    integer:: i, j, singular
    !
    z = ed
    s = sed
    w = 1.0/(addsigma**2+s**2)
    !
    aaMat = 0.0
    bbMat = 0.0
    aMat = 0.0
    bMat = 0.0
    cMat = 0.0
    iMat = 0.0
    info = 0
    stdp = -99.0
    !
    do i=1, ncomp
        pdf(:,i) = pars(1,i)*sqrt(w)*&
                   exp(-0.5*w*(z-pars(2,i))**2)
    end do
    !
    rsumpdf = sum(pdf, dim=2)
    !
    do i=1, ndat
        prob(i,:) = pdf(i,:)/rsumpdf(i)
    end do
    !
    do i=1, ncomp
        aaMat(:,i) = w*(z-pars(2,i))
        bbMat(:,i) = -w + (w*(z-pars(2,i)))**2
        iMat(i,i) = 1.0
    end do
    !
    do i=1, ncomp-1
        do j=1, ncomp-1
            aMat(i,j) = sum((prob(:,i)/pars(1,i)-&
                             prob(:,ncomp)/pars(1,ncomp))*&
                            (prob(:,j)/pars(1,j)-&
                             prob(:,ncomp)/pars(1,ncomp)))
        end do
    end do
    !
    do i=1, ncomp-1
        do j=1, ncomp
            bMat(i,j) = sum(prob(:,j)*aaMat(:,j)*&
                            (prob(:,i)/pars(1,i)-&
                             prob(:,ncomp)/pars(1,ncomp)-&
                             iMat(i,j)/pars(1,i)+&
                             iMat(ncomp,j)/pars(1,ncomp)))
        end do
    end do
    !
    do i=1, ncomp
        do j=1, ncomp
            cMat(i,j) = sum(prob(:,i)*prob(:,j)*&
                            aaMat(:,i)*aaMat(:,j)-&
                            iMat(i,j)*bbMat(:,i)*prob(:,i))
        end do
    end do
    !
    covar(1:ncomp-1,1:ncomp-1) = aMat
    covar(1:ncomp-1,ncomp:) = bMat
    covar(ncomp:,1:ncomp-1) = transpose(bMat)
    covar(ncomp:,ncomp:) = cMat
    !
    call inverse_ger(covar,2*ncomp-1,singular)
    if (singular/=0) then
        info = 1
        return
    end if
    !
    do i=1, 2*ncomp-1
        diag(i) = covar(i,i)
    end do
    if (any(diag<0.0)) then
         info = 1
         return
    end if                      
    !
    stdp(1,1:ncomp-1) = sqrt(diag(1:ncomp-1))
    !
    scnp1 = sum(covar(1:ncomp-1,1:ncomp-1))
    if (scnp1<0.0) then
        info = 1
        return
    end if
    stdp(1,ncomp) = sqrt(scnp1)
    !
    stdp(2,:) = sqrt(diag(ncomp:))
    !
    return
end subroutine apfmmstd
