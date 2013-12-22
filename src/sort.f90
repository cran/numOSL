subroutine sort(ary,n)
!--------------------------------------------------------------
! simple sort, sorting ary in descending order.
! =============================================================
! ary(n), input/output:: real values, an arrary to be sorted.
!
! n,             input:: integer:: numbers of elements in ary.
! =============================================================
! Author:: Peng Jun, 2013.05.20
!
! Dependence:: No
!--------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  real   (kind=8),dimension(n),intent(inout)::ary
  ! local variables
  integer(kind=4)::i,j
  real   (kind=8)::change
  do i=1,n-1
    do j=i+1,n
      if (ary(i)<ary(j))  then
        change=ary(i)
        ary(i)=ary(j)
        ary(j)=change
      end if
    end do
  end do
  return
end subroutine sort
