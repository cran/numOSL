subroutine FineED1(ED,Error,n,ncomp,spreadsigma,pars,&
                   spars,maxlik,BIC,maxiter,eps,tol,errorflag)
!----------------------------------------------------------------------------------------------------------------------------------
! FineED1() is used to perform the fitting of a finite mixture age model, it first change 
! the data to log-scale and initialize the parameters automatically with various ranges, 
! then call subroutine FMMED() to work out the parameter that gives the minimum BIC value,
! it different from subroutine FineED() that it also calculates the parameters' standard error, 
! it also allows the perform of Central Age Model analysis by specifying ncomp=1.
! =================================================================================================================================
!
! ED(n),input                 :: real values, the Equivalent Dose, must be unlogged.
!
! Error(n), input             :: real values, the assocaited absolute error of Equivalent Dose.
!                          
! n, input                    :: integer, the size of Equivalent Dose (or Error).
!
! ncomp, input                :: integer, the number of components want to be decomposed.
!
! spreadsigma, input          :: real value, the spread dispersion to the relative error of Equivalent Dose.
!
! pars(2,ncomp), output       :: real values, the estimated parameters.
!
! spars(2,ncomp), output      :: real values, the estimated standard errors of parameters.
!
! maxlik, output              :: real value, the maximum logged-likelihood value.
!
! BIC, output                 :: real value, the BIC value.
!
! maxiter, input              :: integer, the allowed maximum iterative number.
!
! eps, input                  :: real value, the maximum tolerance for stopping the iterative process.
!
! tol, input                  :: real value, tolerance for diagnosing matrix to be singular if max(abs(matrix)) is smaller than tol.
!
! errorflag, output           :: error message generated during the calculation:
!                                1.1) if successed in calculating parameters' std.errors, errorflag=0;
!                                1.2) if failed in calculating parameters' std.errors, errorflag=1.
! ==================================================================================================================================
! Author:: Peng Jun, 2013.03.05; revised in 2014.03.30.
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!               Standard Errors.Techno-metrics, 30, page 271-281.
!
!               Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!               finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!               
! Dependence  :: subroutine FMMED; subroutine kmeans; subroutine CAM; subroutine AppCovar.
!-----------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                   intent(in)::n
  integer(kind=4),                   intent(in)::ncomp
  integer(kind=4),                   intent(in)::maxiter
  real   (kind=8),dimension(n),      intent(in)::ED
  real   (kind=8),dimension(n),      intent(in)::Error
  real   (kind=8),                   intent(in)::spreadsigma 
  real   (kind=8),                   intent(in)::eps 
  real   (kind=8),                   intent(in)::tol
  real   (kind=8),dimension(2,ncomp),intent(out)::pars
  real   (kind=8),dimension(2,ncomp),intent(out)::spars
  real   (kind=8),                   intent(out)::maxlik
  real   (kind=8),                   intent(out)::BIC
  integer(kind=4),                   intent(out)::errorflag
  !
  ! Variables for subroutine kmeans
  integer(kind=4),parameter::iter=10
  integer(kind=4),parameter::nstart=1000
  integer(kind=4)::belo(n),clusp(ncomp)
  real   (kind=8)::energ(ncomp),clust(ncomp)
  integer(kind=4)::ifault
  !
  ! Local variables
  integer(kind=4)::i,j
  real(kind=8)::minBIC
  real(kind=8)::sMaxlik
  real(kind=8),dimension(2,ncomp)::cpars,kpars
  real(kind=8),dimension(n)::sED,sError
  !
  ! Transform absolute std.errors of equivalent doses to relative errors
  sError=Error/ED
  ! Transform equivalent doses to logged-scale
  sED=dlog(ED)
  !
  ! For a central age model
  if (ncomp==1)  then 
    errorflag=0
    call CAM(ED,Error,n,spreadsigma,maxiter,&
             eps,pars,spars,maxlik,BIC)
    return
  else 
    ! For finite mixtrue age model of various number of components
    !
    ! Initialize cpars
    cpars(1,:)=pars(1,:)
    do j=1,ncomp
      cpars(2,j)=minval(sED)+(maxval(sED)-minval(sED))*real(j,kind=8)/real(ncomp+1,kind=8)
    end do
    ! Set minBIC to be a large number
    minBIC=1.0D+30
    ! 
    do i=1,3
      ! Initialize pars
      ! 1) proportions (1:ncomp)
      pars(1,:)=1.0D+00/real(ncomp,kind=8)
      ! 2) Characteristic equivalent doses (1:ncomp)
      do j=1,ncomp
        pars(2,j)=minval(sED)+(maxval(sED)-minval(sED))*real(j+i-2,kind=8)/real(ncomp+1,kind=8)
      end do
      ! Call FMMED() to perfect the initialized parameters (pars)
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
    ! Approximate standard errors of parameters
    call AppCovar(pars,ncomp,ED,Error,n, &
                  spreadsigma,spars,errorflag,tol)
    !
    ! Transform equivalent dose to un-logged scale
    pars(2,:)=dexp(pars(2,:))
    spars(2,:)=spars(2,:)*pars(2,:)
    !
    return
  end if
  !
end subroutine FineED1
