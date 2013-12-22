subroutine leaveone(aryin,n,indx,aryout)
! ----------------------------------------------------------------------------------------------------
! Specifying a index, leaveone is used to remove correspond value and return remaining values.
! ====================================================================================================
! n,            input:: integer, length of the whole array.
!
! indx,         input:: integer, a index ranging from 1 to n, corresopnding to a value to be removed.
!
! aryin(n),     input:: integer values, the whole array.
!
! aryout(n-1), output:: integer values, remaining values after the removing.
! ====================================================================================================
! Author:: Peng Jun, 2013.05.20.
! 
! Dependence:: No
!-----------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),               intent(in)::n
  integer(kind=4),               intent(in)::indx
  integer(kind=4),dimension(n),  intent(in)::aryin
  integer(kind=4),dimension(n-1),intent(out)::aryout
  !
  if(indx==1)  then
    aryout(:)=aryin(2:n)
  else if (indx==n)  then
    aryout(:)=aryin(1:n-1)
  else
    aryout(:)=(/aryin(1:indx-1),aryin(indx+1:n)/)
  end if
  return
end subroutine leaveone
