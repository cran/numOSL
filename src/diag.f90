subroutine diag(a,n,b)
 !------------------------------------------------------------------------
 ! subroutine diag is used to extract the diag elements of matrix A
 ! =======================================================================
 ! A(n,n)  :: input, real values, a matrix
 !   n     :: input, integer value, the dimension of A
 !  b(N)   :: output, real values, the diag elements of A
 ! =======================================================================
 ! Dependency:: no
 !
 ! Author :: Peng Jun, 2013.01.28
 !------------------------------------------------------------------------
   implicit none
   integer(kind=4),intent(in)::n
   real(kind=8),   intent(in)::a(n,n)
   real(kind=8),   intent(out)::b(n)
   ! local variables
   integer(kind=4)::i
   !
   do i=1,n
     b(i)=a(i,i)
   end do
   return
end subroutine diag
