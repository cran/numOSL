subroutine FineED(ED,Error,n,ncomp,spreadsigma,&
                  pars,maxlik,BIC,maxiter,eps)
!----------------------------------------------------------------------------------------------------------
! FineED() is used to perform the finite mixture age model fitting.
! It first change the data to log-scale and initialize the parameters automatically with various
! ranges, then call subroutine FMMED() to work out the parameter that gives the minimum BIC value.
! =========================================================================================================
!
! ED(n),input                 :: real values, the Equivalent Dose, must be unlogged.
!
! Error(n), input             :: real values, the assocaited absolute error of Equivalent Dose.
!                        
! n, input                    :: integer, the size of ED.
!
! ncomp, input                :: integer, the number of component want to be decomposed.
!
! spreadsigma, input          :: real value, the spread dispersion to the relative error of Equivalent Dose.
!
! pars(2,ncomp), output       :: real values, the estimated parameters.
!
! maxlik, output              :: real value, the maximum logged likelihood value.
!
! BIC, output                 :: real value, the BIC value.
!
! maxiter, input              :: integer, the allowed maximum iterative number.
!
! eps, input                  :: real value, the maximum tolerance for stopping the iterative process.
! ===========================================================================================================
! Author:: Peng Jun, 2013.03.03; revised in 2014.03.30.
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!               Standard Errors.Techno-metrics, 30, page 271-281.
!
!               Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!               finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!               
! Dependence  :: subroutine FMMED; subroutine kmeans.
!------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                   intent(in)::n
  integer(kind=4),                   intent(in)::ncomp
  integer(kind=4),                   intent(in)::maxiter
  real   (kind=8),dimension(n),      intent(in)::ED
  real   (kind=8),dimension(n),      intent(in)::Error
  real   (kind=8),                   intent(in)::eps
  real   (kind=8),                   intent(in)::spreadsigma 
  real   (kind=8),                   intent(out)::maxlik
  real   (kind=8),                   intent(out)::BIC
  real   (kind=8),dimension(2,ncomp),intent(out)::pars
  !
  ! Arguments for subroutine kmeans()
  integer(kind=4),parameter::iter=10
  integer(kind=4),parameter::nstart=1000
  integer(kind=4)::belo(n),clusp(ncomp)
  real   (kind=8)::energ(ncomp),clust(ncomp)
  integer(kind=4)::ifault
  ! Local variables
  integer(kind=4)::i,j
  real(kind=8)::minBIC
  real(kind=8)::sMaxlik
  real(kind=8),dimension(2,ncomp)::cpars,kpars
  real(kind=8),dimension(n)::sED,sError
  !
  ! Transform absolute error of equivalent dose to relative errors
  sError=Error/ED
  ! Transform equivalent dose to logged-scale
  sED=dlog(ED)
  ! minBIC is very large
  minBIC=1.0D+30
  !
  do i=1,3
    ! Initialize proportions
    pars(1,:)=1.0D+00/real(ncomp,kind=8)
    ! Initialize characteristic equivalent doses 
    do j=1,ncomp
      pars(2,j)=minval(sED)+(maxval(sED)-minval(sED))*real(j+i-2,kind=8)/real(ncomp+1,kind=8)
    end do
    !
    call FMMED(sED,sError,n,ncomp,spreadsigma,&
               pars,maxlik,BIC,maxiter,eps)
    ! Store improved results only
    if(BIC<minBIC)  then
      cpars=pars
      minBIC=BIC
      sMaxlik=maxlik
    end if
    !
  end do
  !
  ! Check if k-means method initialized parameters
  ! can further improved the results (decrease BIC)
  call kmeans(1,n,ncomp,iter,nstart,sED,&
              belo,clust,clusp,energ,ifault)
  !
  !kpars(1,:)=1.0D+00/real(ncomp,kind=8)
  kpars(1,:)=real(clusp,kind=8)/sum(clusp)
  kpars(2,:)=clust
  !
  call FMMED(sED,sError,n,ncomp,spreadsigma,&
             kpars,maxlik,BIC,maxiter,eps)
  !
  if(BIC<minBIC)  then
    cpars=kpars
    minBIC=BIC
    sMaxlik=maxlik
  end if
  ! Parameters for output
  pars=cpars
  BIC=minBIC
  maxlik=sMaxlik
  !
  return
end subroutine FineED
