subroutine lmhess(pars,xdat,ydat,npars,ndat,tol,minAbsPar,&
                  hessian,gradient,value,errorflag,model)
!------------------------------------------------------------------------------------------------------
! Subroutine fdhessian() is used to calculate the gradient, hessian matrix of a 
! given function contained in its inner, that is fun34 in CW-OSL fitting models, 
! using finite-difference approximation (non-origin and origin, cw, lm)
! ====================================================================================================
!
! pars(npars)          :: input, real values, the parameters of the function.
!
! xdat (ndat)          :: input, real values, the xdat values.
!
! ydat (ndat)          :: input, real values, the ydat values.
!
! npars                :: input, integer, the length of tim ( or signal).
!
! ndat                 :: input, integer, the dimension of the parameters.
!
! tol                  :: input, real value, tolerance value for diagnosing sigular matrix.
!
! minAbspar            :: input, real value, the allowed minimum absolute parameter.
!
! hessian(npars,npars) :: output, real values, the hessian matrix.
!
! gradient(npars)      :: output, real values, the gradient of the parameters.
!
! value                :: output, real value, the correspond function value for the specified parameters.
! 
! errorflag(5)         :: output, integer values, error message generated during the calling:
!                         1) if subroutine GJordan is called sucessfully, errorflag(1)=0, otherwise 1;
!                         2) if function value can be calculated,         errorflag(2)=0, otherwise 1;
!                         3) if gradient can be calculated,               errorflag(3)=0, otherwise 1; 
!                         4) if hessian can be calculated,                errorflag(4)=0, otherwise 1; 
!                         5) if no error appears in arrary allocation,    errorflag(5)=0, otherwise 1.
!
! model                :: input, integer, a model to be used for approximation:
!                         1) y=I(1)*exp(-lamda(1)*x)+
!                              I(2)*exp(-lamda(2)*x)+...+
!                              I(k)*exp(-lamda(k)*x),  fitting a 'cw' signal curve
!
!                         2) y=a*x+b, fitting linear grow curve
!
!                         3) y=a*(1-exp(-b*x))+c, fitting exponential grow curve
!
!                         4) y=a*(1-exp(-b*x)+c*x+d, fitting exponential plus linear grow curve
!
!                         5) y=a1*(x/max(x))*exp(-b1*x^2)/2/max(x))+
!                              a2*(x/max(x))*exp(-b2*x^2)/2/max(x))+...+
!                              ak*(x/max(x))*exp(-bk*x^2)/2/max(x)), fitting a 'lm' signal curve
!
!                         6) y=a*x, fitting linear grow curve (origin)
! 
!                         7) y=a*(1-exp(-b*x)), fitting exponential grow curve (origin)
!
!                         8) y=a*(1-exp(-b*x)+c*x, fitting exponential plus linear grow curve (origin)
! ======================================================================================================
! Dependence:: subroutine GJordan, inter function fun34.
!
! Author:: Peng Jun, 2013.05.21, revised in 2013.05.23, revised in 2013.07.24.
!
! Reference:  Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the
!             R Development Core Team (2013). nlme: Linear and Nonlinear Mixed
!             Effects Models. R package version 3.1-108.
!-------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::npars                 
  integer(kind=4),intent(in)::ndat 
  integer(kind=4),intent(in)::model                   
  integer(kind=4),intent(out)::errorflag(5)             
  real   (kind=8),intent(in)::pars(npars)            
  real   (kind=8),intent(in)::xdat(ndat)                
  real   (kind=8),intent(in)::ydat(ndat)            
  real   (kind=8),intent(in)::tol                   
  real   (kind=8),intent(in)::minAbsPar              
  real   (kind=8),intent(out)::gradient(npars)       
  real   (kind=8),intent(out)::hessian(npars,npars)  
  real   (kind=8),intent(out)::value                
  ! local variables
  integer(kind=4)::i,j,p
  integer(kind=4)::ncols
  integer(kind=4)::solerror
  real  (kind=8)::incr(npars)
  real  (kind=8)::diagpar(npars,npars)
  real  (kind=8),allocatable::frac(:),cfrac(:)
  real  (kind=8)::ffrac(1+2*npars)
  real  (kind=8),allocatable::cols(:,:),ccols(:,:),shifted(:,:),transcols(:,:)
  real  (kind=8)::pcols(npars,2*npars+1)
  real  (kind=8),allocatable::xcols(:,:),cxcols(:,:)
  real  (kind=8),allocatable::pxcols(:,:)
  real  (kind=8),parameter::eps=6.055454e-06       
  real  (kind=8)::maxx           
  !
  ! Decide the incr values
  ! for each initial par in pars, check if it
  ! is smaller than minabspar, the incr will
  ! be decided through this check
  do i=1,npars
    if(dabs(pars(i))<=minAbsPar)  then
      incr(i)=minAbsPar*eps
    else
      incr(i)=dabs(pars(i))*eps
    end if
  end do
  !
  ! build a diagal matirx diagpar
  diagpar=0.0D+00
  do i=1,npars
    diagpar(i,i)=1.0D+00
  end do
  ! creat ffrac with length of 2*npars+1
  ! and specify values for it
  ffrac(1)=1.0D+00
  ffrac(2:npars+1)=incr
  ffrac(npars+2:2*npars+1)=incr**2
  ! create pcols to be npars rows
  ! and 2*npars columns, then storing
  ! diagpar in it
  pcols(:,1)=0.0D+00
  pcols(:,2:npars+1)=diagpar
  pcols(:,npars+2:2*npars+1)=-diagpar
  !
  ! default return values if error appears
  value=-99.0D+00
  gradient=-99.0D+00
  hessian=-99.0D+00
  errorflag=0
  !
  ! in a total of (npars-1) times looping
  ! add  (npars-i) columns to cols in each loop number
  ! hence a toal of npars*(npars-1)/2 columns for cols
  allocate( cols(1:npars,1:npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! in a total of (npars-1) times looping
  ! add (npars-i) number to frac in each loop number
  ! hence a total of npars*(npars-1)/2 values for frac
  allocate( frac(1:npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  !
  ncols=0
  do i=1,npars-1
    ! in each loop(i), generate ccols and cfrac with 
    ! npars rows (npars-i) columns 
    allocate( ccols(1:npars,1:npars-i), stat=errorflag(5))
    if(errorflag(5)/=0) return
    allocate( cfrac(1:npars-i), stat=errorflag(5))
    if(errorflag(5)/=0) return
    ! fill ccols and cfrac in each loop
    do j=i+1,npars
      ccols(:,j-i)=diagpar(:,i)+diagpar(:,j)
      cfrac(j-i)=incr(i)*incr(j)       
    end do 
    ! now store ccols and cfrac to cols
    ! and frac respectively
    cols(:,ncols+1:ncols+npars-i)=ccols
    frac(ncols+1:ncols+npars-i)=cfrac
    ! delocate ccols and cfrac, preparing
    ! for new allocations
    deallocate(ccols, stat=errorflag(5))
    if(errorflag(5)/=0) return
    deallocate(cfrac, stat=errorflag(5))
    if(errorflag(5)/=0) return
    ! update started filling index ncols
    ncols=ncols+npars-i 
    !
  end do
  ! now add up pcols and cols togeteher to be ccols
  ! ccols has a total of ( 2*npars+1 + npars*(npars-1)/2 ) columns
  allocate( ccols(1:npars,1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  ccols(:,1:2*npars+1)=pcols
  ccols(:,2*npars+2:npars*(npars-1)/2+2*npars+1)=cols
  ! not need cols presently, so release it
  deallocate(cols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! allocate cols again to store ccols, and release ccols
  allocate(cols(1:npars,1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  cols=ccols
  deallocate(ccols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! now add up ffrac and frac together to be cfrac,
  ! cfrac has a total of ( 2*npars+1 + npars*(npars-1)/2 ) values
  allocate( cfrac(1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  cfrac(1:2*npars+1)=ffrac
  cfrac(2*npars+2:npars*(npars-1)/2+2*npars+1)=frac
  ! not need frac so release it
  deallocate(frac, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! allocate frac again, store cfrac in it then release cfrac
  allocate(frac(1:npars*(npars-1)/2+2*npars+1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  frac=cfrac
  deallocate(cfrac, stat=errorflag(5)) 
  if(errorflag(5)/=0) return
  ! specify p to be the column number of cols
  ! that is npars*(npars-1)/2+2*npars+1
  p=npars*(npars-1)/2+2*npars+1
  ! now allocate shifted to be the same shape 
  ! with cols and store some values in it
  allocate(shifted(1:npars,p), stat=errorflag(5))
  if(errorflag(5)/=0) return
  do i=1,p
    shifted(:,i)=pars+incr*cols(:,i)
  end do
  ! allocate transcols to store transposed cols
  ! transcols has p rows, npars columns
  allocate( transcols(1:p,1:npars), stat=errorflag(5))
  if(errorflag(5)/=0) return
  transcols=transpose(cols)
  ! now cols will not be needed, release it
  deallocate(cols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! allocate pxcols to store transcols 
  ! in differ style
  allocate(pxcols(p,1+2*npars), stat=errorflag(5))
  if(errorflag(5)/=0) return
  pxcols(:,1)=1.0D+00
  pxcols(:,2:npars+1)=transcols
  pxcols(:,npars+2:2*npars+1)=(transcols)**2
  !
  ncols=0
  ! now allocate xcols with p rows and 
  ! npars*(npars-1)/2 columns
  allocate( xcols(p,1:npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! 
  do i=1,npars-1
    ! in each loop, allocate cxcols 
    allocate(cxcols(p,1:npars-i), stat=errorflag(5))
    if(errorflag(5)/=0) return
    do j=i+1,npars
      cxcols(:,j-i)=transcols(:,i)*transcols(:,j)
    end do
    ! store cxcols to xcols and release it
    xcols(:,ncols+1:ncols+npars-i)=cxcols
    deallocate(cxcols, stat=errorflag(5))
    if(errorflag(5)/=0) return
    ! update started filling index
    ncols=ncols+npars-i  
  end do
  ! now store pxcols and xcols together to cxcols
  allocate( cxcols(p,1:1+2*npars+npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  cxcols(:,1:1+2*npars)=pxcols
  cxcols(:,2+2*npars:1+2*npars+npars*(npars-1)/2)=xcols
  ! release xcols and pxcols
  deallocate(xcols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  deallocate(pxcols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! allocate xcols again to store cxcols
  allocate(xcols(p,1:1+2*npars+npars*(npars-1)/2), stat=errorflag(5))
  if(errorflag(5)/=0) return
  xcols=cxcols
  ! release cxcols
  deallocate(cxcols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! allocate pxcols again and store some
  ! values to it then release shifted
  ! now pxcols has p rows and 1 column
  allocate(pxcols(1:p,1), stat=errorflag(5))
  if(errorflag(5)/=0) return
  do i=1,p
    pxcols(i,:)=fun34(shifted(:,i))
  end do
  deallocate(shifted, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! solve xcols %*% X = pxcols
  ! store solved X values in pxcols
  call GJordan(xcols,pxcols,p,1,solerror,tol)
  if(solerror==1)  errorflag(1)=1
  ! scale pxcols with frac and release frac,
  ! release xcols
  pxcols(:,1)=pxcols(:,1)/frac
  deallocate(xcols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  deallocate(frac, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ncols=2*npars+2 
  ! fill diagpar with some new values,  note that
  ! non-diagnal part of digpar are zeros
  ! these value are comes from pxcols, and will
  ! be used to constructe hessian matrix
  do j=1,npars
    do i=1,npars
        if (i==j) diagpar(i,j)=pxcols(1+npars+i,1)
        if (i>j)  diagpar(j+1:npars,j)=pxcols(ncols:ncols+npars-j-1,1)  
    end do
    ncols=ncols+npars-j
  end do
  ! estimate gradients
  do i=1,npars
    gradient(i)=pxcols(1+i,1)
  end do
  ! estimate fun(pars)
  value=pxcols(1,1)
  ! now pxcols will not be needed, release it
  deallocate(pxcols, stat=errorflag(5))
  if(errorflag(5)/=0) return
  ! estimate hessian matrix
  hessian=diagpar+transpose(diagpar)
  !
  ! check Inf and NaN for value
  if(value .ne. value .or. &
     value+1.0D+00==value)              errorflag(2)=1
  ! check Inf and NaN for gradient
  if(any(gradient .ne. gradient) .or. &
     any(gradient+1.0D+00==gradient))   errorflag(3)=1
  ! check Inf and NaN for hessian
  if(any(hessian .ne. hessian) .or. &
     any(hessian+1.0D+00==hessian))     errorflag(4)=1
  ! if no error appears, return 
  return
  ! INNER FUNCTION FOR HESSION APPROXIMATION
  contains
    function fun34(x)
      implicit none
      real(kind=8)::fun34
      real(kind=8),dimension(npars)::x
      real(kind=8),dimension(ndat)::fvec
      real(kind=8),dimension(4)::cx
      integer(kind=4)::k
      !
      if(model==1) then
        ! for fitting 'cw'
        fvec=0.0D+00
        do k=1,npars/2
          fvec=fvec+x(k)*dexp(-x(k+npars/2)*xdat)   
        end do
      else if(model>=2 .and. model<=4) then
        ! for fitting grow curve
        cx=0.0D+00
        cx(1:npars)=x
        if(model==2) then
          ! linear
          fvec=cx(1)*xdat+cx(2)
        else if(model==3) then
          ! exponential 
          fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)
        else if(model==4) then
          ! exponential plus linear
          fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)*xdat+cx(4)
        end if
      else if(model==5) then
        ! for fitting 'lm'
        maxx=maxval(xdat)
        fvec=0.0D+00
        do k=1,npars/2
          fvec=fvec+x(k)*(xdat/maxx)*&
               dexp(-x(k+npars/2)*xdat**2/2.0D+00/maxx)  
        end do
      else if(model>=6 .and. model<=8) then
        ! for fitting grow curve
        cx=0.0D+00
        cx(1:npars)=x
        if(model==6) then
          ! linear
          fvec=cx(1)*xdat
        else if(model==7) then
          ! exponential 
          fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))
        else if(model==8) then
          ! exponential plus linear
          fvec=cx(1)*(1.0D+00-dexp(-cx(2)*xdat))+cx(3)*xdat
        end if
      end if
      fun34=dsqrt(sum((fvec-ydat)**2))
      return
    end function fun34
end subroutine lmhess
