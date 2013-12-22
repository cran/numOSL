subroutine GJordan(A,B,N,M,error,tol)
!----------------------------------------------------------------------------------------------------
! GJordan is a subroutine used to do Guass-Jordan elimation by full pivoting elimination.
! ===================================================================================================
!
! A(N,N)::   input/output, real values, matrix that contains the coefficient of the linear equations.
!
! B(N,M)::   input/output, real values, the values of equations in the right side.
!
! N     ::   input, integer value, the dimension of matrix A that wanted to be solved.
!
! M     ::   input, integer value, the number of colums of matrix B.
!
! tol   ::   input, real value, diagnose A to be a singular matrix if max(A) is smaller than tol.
!
! error ::   output, integer value, if A is a singular matrix, return 1, otherwise 0.
! ==================================================================================================
! Dependence:: No
!
! Author:: Peng Jun, 2012.01.09, revised in 2013.01.28.
!
! Reference:: 徐士良， fortran 常用算法集
!---------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::N
  integer(kind=4),intent(in)::M
  integer(kind=4),intent(out)::error
  integer(kind=4)::JS(N),IS
  real   (kind=8),intent(in)::tol
  real   (kind=8),intent(inout)::A(N,N)
  real   (kind=8),intent(inout)::B(N,M)
  ! Local variables 
  real   (kind=8)::Q
  real   (kind=8),allocatable::D(:)
  integer(kind=4)::k,i,j
  !
  error=0
  !
  ! Start the major loop
  do k=1,N
  !
    Q=0.0D+00
    ! Find out the largest absolute value from row and colum index 
    ! that range from k to N, and save row, colum index that it belongs to
      do i=k,N
        do j=k,N
	  if(abs(A(i,j))>=Q) then
	    Q=abs(A(i,j))
	    JS(k)=j
	    IS=i
	  end if
	end do
      end do
    ! If the largest absolute value is zero, set error
    ! to 1 and return
    if (abs(Q)<tol)  then
      error=1
      return
      end if
      ! Swap the row that contains the largest absolute value with row k
      ! it's a row swap
      allocate(D(1:N-k+1))
      D=A(k,k:N)
      A(k,k:N)=A(IS,k:N)
      A(IS,k:N)=D
      deallocate(D)
      ! Swap the correspond values in matirx B
      allocate(D(1:M))
      D=B(k,1:M)
      B(k,1:M)=B(IS,1:M)
      B(IS,1:M)=D
      deallocate(D)
      ! Swap the colum that contains the largest absolute value with
      ! colum k
      allocate(D(1:N))
      D=A(1:N,k)
      A(1:N,k)=A(1:N,JS(k))
      A(1:N,JS(k))=D
      deallocate(D)	 
      ! Devide each values of colums that range from k+1 to N
      ! by A(k,k) 
      A(k,k+1:n)=A(k,k+1:n)/A(k,k)
      !
      ! Devide each values of  colums that range from 1 to M
      ! by A(k,k)
      B(k,1:m)=B(k,1:m)/A(k,k)
      !
      ! Substract values in all row but k by A(i,k)*A(k,j)
      ! in both matrix A and B
      do i=1,N
        if(i/=k)  then  
	  A(i,k+1:N)=A(i,k+1:N)-A(i,k)*A(k,k+1:N)
	  B(i,1:M)=B(i,1:M)-A(i,k)*B(k,1:M)
        end if 
      end do
     ! End the major loop
  end do
  ! Swap the row that contains the largest absolute value
  ! to the first row
  allocate(D(1:M))
  do k=N,1,-1
    D=B(k,1:M)
    B(k,1:M)=B(JS(k),1:M)
    B(JS(k),1:M)=D	  
  end do
  deallocate(D)	
  return
end subroutine GJordan
