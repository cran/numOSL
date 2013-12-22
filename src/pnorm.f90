subroutine pnorm(x,n,upper)
!-----------------------------------------------------------------
! A wrapper subroutine calls function alnorm to calculate a series
! culmulative probability densities of the normal distribution
! ===============================================================
! Author :: Peng Jun.
!
! Dependence: function alnorm.
!-----------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  real(kind=8),intent(inout)::x(n)
  logical,intent(in)::upper
  ! local variables
  real(kind=8)::alnorm
  integer(kind=4)::i
  do i=1,n
    x(i)=alnorm(x(i),upper)
  end do
  return
end subroutine pnorm

