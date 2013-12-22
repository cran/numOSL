subroutine sd(mat,nrow,ncol,outsd)
!----------------------------------------------------------------------------------
! sd() calculate standard deviations of matrix mat by column.
! =================================================================================
!
! mat(nrow,ncol),   input:: real values, a matrix to be calculated.
!
! nrow,             input:: integer, row numbers of mat.
!
! ncol,             input:: integer, col numbers of mat.
!
! outsd(ncol),     output:: real values, calculated standard deviations by column.
! =================================================================================
! Author:: Peng Jun.
!
! Dependence:: No
!----------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::nrow
  integer(kind=4),intent(in)::ncol
  real   (kind=8),dimension(nrow,ncol),intent(in)::mat
  real   (kind=8),dimension(ncol),intent(out)::outsd
  ! local variables
  integer(kind=4)::i
  real   (kind=8)::average
  !
  do i=1, ncol
    average=sum(mat(:,i))/real(nrow,kind=8)
    outsd(i)= sqrt( sum( (mat(:,i)-average )**2 )/real(nrow-1,kind=8) )
  end do
  !
  return
end subroutine sd
