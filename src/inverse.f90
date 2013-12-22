subroutine inverse(a,n,error,tol)
 !--------------------------------------------------------------------
 ! subroutine inverse is used to call
 ! GJordan to perform the inverse on 
 ! matrix a
 !
 ! A(n,n)::    input/output, real values, the matrix used to do inverse
 ! n     ::    input, integer values, the dimsension of matrix A
 ! error ::    output, integer value, if error=1, A is singular,
 !             otherwise error=0
 ! tol   ::    input, real values, the maximum absolute value in 
 !             matrix A that allowed for A to be a non-singular matrix
 !
 ! Dependence:: subroutine GJordan
 !
 ! Author:: Peng Jun, 2013.01.18
 !
 !---------------------------------------------------------------------
   implicit none
   integer(kind=4),intent(in)::n
   real(kind=8),intent(inout)::a(n,n)
   integer(kind=4),intent(out)::error
   real(kind=8),intent(in)::tol
   !local variables
   real(kind=8)::b(n,n)
   integer(kind=4)::i
   !
   b=0.0D+00
   do i=1,n
     b(i,i)=1.0D+00
   end do  
   call GJordan(a,b,n,n,error,tol)   
   a=b
   return
end subroutine inverse
