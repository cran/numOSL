subroutine MAM(ED,Error,nED,pars,spars,value,BIC,&
               npars,sigma,maxiter,tol,errorflag)
!--------------------------------------------------------------------------------------------------------------------
! MAM() attempts to estimate the values of the parameters in Minimum Age Models
! by using the L-BFGS-B method, insteading for a globle minimum, the L-BFGS-B routine
! results in a local minimum, it calls subroutine MamED() with various intial guess values
! and loopping in turn to find a most appropriate one.
! Note that this subroutine do logged equivalent dose analysis, which means that all the
! ED data must be larger than zero!
! ==================================================================================================================
! 
! The results rely strongly on the initial guess values, but in this subroutine, by
! adjusting the sigma, converge can be achieved in most situations.
!
! ED(nED), input            :: real values, the equivalent dose, must be of type un-logged.
!
! Error(nED), input         :: real values, the associated absolute errors for the equivalent dose.
!
! nED, input                :: integer, the size of the ED data.
!
! npars, input              :: integer,the dimension of the problems, 3 for MAM3, 4 for MAM4.
!
! pars(npars), output       :: real values, the estimated parameters.
!
! spars(npars), output      :: real values, the estimated standard errors for parameters.
!
! value , output            :: real value, the estimated logged maximum logged likelihood value.
!
! BIC, output               :: real value, the estimated minimum BIC value.
!
! sigma, input              :: real value, the added spread to the ED data.
!
! maxiter, input            :: integer, the allowed maximum iterative numbers.
!
! tol, input                :: real value, the allowed maximum tolerance for the hessian to be non-sigular.
!
! errorflag                 :: integer values, the error message generated in the analysis:
!                              1.1) a successful work gives errorflag=0;
!                              1.2) a fail work gives errorflag=1.
!====================================================================================================================
! Author :: Peng Jun, 2013.03.16, revised in 2013.04.04.
! 
! Reference :: Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of
!              single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design
!              and statistical models. Archaeometry, 41, pp. 339-364.
!
!              Galbraith, R.F., Roberts, R.G., 2012. Statistical aspects of equivalent dose and error calculation 
!              and display in OSL dating: An overview and some recommendations. Quaternary Geochronology,11, pp. 1-27.
!
! Dependence :: subroutine kmeans; subroutine MamED.
!-------------------------------------------------------------------------------------------------------------------
 implicit none
  integer(kind=4),intent(in)::nED
  integer(kind=4),intent(in)::npars
  integer(kind=4),intent(in)::maxiter
  real   (kind=8),intent(in)::tol
  real   (kind=8),intent(in)::sigma
  real   (kind=8),intent(out)::value
  real   (kind=8),intent(out)::BIC
  integer(kind=4),intent(out)::errorflag
  real   (kind=8),dimension(nED),intent(in)::ED
  real   (kind=8),dimension(nED),intent(in)::Error
  real   (kind=8),dimension(npars),intent(out)::pars
  real   (kind=8),dimension(npars),intent(out)::spars
  ! 
  ! Variables for subroutine kmeans()
  integer(kind=4),parameter::nclu=3,iter=20,nstart=100
  integer(kind=4)::belo(nED),clusp(nclu)
  real   (kind=8)::energ(nclu),clust(nclu)
  integer(kind=4)::ifault
  !
  ! Variables for subroutine MamED()
  integer(kind=4),dimension(5)::message
  real   (kind=8),dimension(2,npars)::bound
  !
  ! Local variables
  real   (kind=8),dimension(npars)::cpars,dpars
  real   (kind=8),dimension(npars)::cspars,dspars
  real   (kind=8)::cvalue,dvalue
  real   (kind=8)::sED(nED),sError(nED)
  real   (kind=8)::minED,maxED,meanED,SdED
  real   (kind=8):: minvalue,minvalue1,change
  integer(kind=4)::i,j,k
  real   (kind=8)::mu(3),gama
  !
  ! Here only need to change ED to log-scale
  sED=dlog(ED)
  ! Some character values for logged ED data
  minED=minval(sED)
  maxED=maxval(sED)
  meanED=sum(sED)/real(nED,kind=8)
  SdED=dsqrt(sum((sED-meanED)**2)/real(nED-1,kind=8)) 
  !
  ! Initializing boundary for L-BFGS-B method
  if (npars==3)     then  
    ! low bunodary       ! up boundary
    bound(1,1)=1.0D-04 ;  bound(2,1)=0.99D+00 
    bound(1,2)=minED   ;  bound(2,2)=maxED 
    bound(1,3)=1.0D-03 ;  bound(2,3)=5.0D+00  
    !
  elseif(npars==4)  then
    ! low bunodary       ! up boundary
    bound(1,1)=1.0D-04 ;  bound(2,1)=0.99D+00 
    bound(1,2)=minED   ;  bound(2,2)=maxED 
    bound(1,3)=minED   ;  bound(2,3)=maxED 
    bound(1,4)=1.0D-03 ;  bound(2,4)=5.0D+00  
    ! 
  end if
  !
  ! Using k-means method to initialze a series logged ED values,
  ! which are used for subsequent L-BFGS-B analysis
  call  kmeans(1,nED,nclu,iter,nstart,&
               sED,belo,clust,clusp,energ,ifault)
  ! Sort clusts in ascending order
  do i=1,nclu-1
    do j=i+1,nclu
      if (clust(i)>clust(j))  then
        change=clust(i)
        clust(i)=clust(j)
        clust(j)=change
      end if
    end do
  end do
  !
  ! Set some values for L-BFGS-B loopping
  minvalue=1.0D+30; minvalue1=1.0D+30; errorflag=1
  ! Set single gama value 
  ! for MAM3
  gama=clust(1)
  ! Set three alternative mu values 
  ! for MAM4
  mu(1)=clust(2)
  mu(2)=sum(clust)/real(nclu,kind=8)
  mu(3)=meanED
  ! Start do loopping for different cases to
  ! find out the most appropriate estimation
  pars=0.0D+00; spars=0.0D+00 
  dpars=0.0D+00; dspars=0.0D+00; dvalue=0.0D+00
  !
  ! For MAM3 with three parameters
  if (npars==3)  then
    ! In this case a total of 15 looping iterations 
    ! are tried to find out a stable solution
    do i=1,3
      do j=1,5  
          ! Update initial parameters
          cpars(1)=0.01D+00*(5.0D+00)**(i-1)   ! the proportion 
          cpars(2)=gama                        ! so-called gama value 
          cpars(3)=SdED*0.4D+00*real(j,kind=8) ! so-called sigma
          !
          ! Call MamED to optimize the parameters and standard errors
          call MamED(ED,Error,nED,cpars,cspars,cvalue,npars, &
                     sigma,maxiter,tol,bound,message)
          !
          ! Pick out the appropriate parameters that give
          ! the minimum minus log-like value and satisfy the
          ! message checking
          if( cvalue<minvalue .and. all(message(1:4)==0) )  then
            minvalue=cvalue
            value=cvalue
            pars=cpars
            spars=cspars
            errorflag=0
          end if
          !
          ! If any estimated paramter is near to the boundary, picking out 
          ! one that gives a minimum minus log-like value and store it too
          if( errorflag==1 .and. &
              cvalue<minvalue1 .and. &
              all(message(1:3)==0) .and. &
              message(4)==1 )  then
              !
            minvalue1=cvalue
            dpars=cpars
            dspars=cspars
            dvalue=cvalue
            !
          end if
          !
      end do
    end do
    !
  ! For MAM4 with four parameters
  elseif(npars==4)  then
    ! In this case a total of 27 looping iterations 
    ! are tried to find out a stable solution
    do i=1,3
      do j=1,3
        do k=1,3 
            ! Update initial parameters  
            cpars(1)=0.01D+00*(5.0D+00)**(i-1)    ! the proportion  
            cpars(2)=gama                         ! so-called gama value    
            cpars(3)=mu(k)                        ! so-called gama value
            cpars(4)=SdED*0.5D+00*real(j,kind=8)  ! so-called mu value
            !
            ! Perform MamED to pick out the parameters and standard errors
            call MamED(ED,Error,nED,cpars,cspars,cvalue,npars, &
                       sigma,maxiter,tol,bound,message)
            !
            ! Pick out the appropriate parameters that give
            ! the minimum minus log-like value and satisfy the
            ! message checking
            if( cvalue<minvalue .and. all(message==0) )  then
              minvalue=cvalue
              value=cvalue
              pars=cpars
              spars=cspars
              errorflag=0
            end if
            !
            ! If any estimated paramter is near to the boundary, picking out 
            ! one that gives a minimum minus log-like value and store it too
            if( errorflag==1 .and. &
                cvalue<minvalue1 .and. &
                all(message(1:3)==0) .and. &
                message(4)==1 .and. &
                message(5)==0 )  then
              !
              minvalue1=cvalue
              dpars=cpars
              dspars=cspars
              dvalue=cvalue
            !
            end if
            !
        end do
      end do
    end do
    !
  end if
  !
  ! If the results can not satisfy the standard of appropriate parameters,
  ! use parameters which are near to the boundary as a substitution
  if(errorflag==1)  then
    pars=dpars
    spars=dspars
    value=dvalue
  end if
  !
  ! Transform the minus log-like value to a correct one
  ! and calculate the BIC value
  value=-value
  BIC=-2.0D+00*value+real(npars,kind=8)*dlog(real(nED,kind=8))
  !
  ! Change equivalent dose and errors back to unlogged scale
  pars(2:(npars-1))=dexp(pars(2:(npars-1)))
  spars(2:(npars-1))=spars(2:(npars-1))*pars(2:(npars-1))
  !
  return
end subroutine MAM
