subroutine comb_next ( done, n, k, iarray )
!
!*******************************************************************************
!
!! COMB_NEXT computes combinations of K things out of N.
!
!
!  Discussion:
!
!    The combinations are computed one at a time, in lexicographical order.
!
!  Reference:
!
!    Charles Mifsud,
!    Combination in Lexicographic Order,
!    ACM algorithm 154,
!    March 1963.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input/output, logical DONE, indicator.
!    On input, if this is the first call, the user should set
!    DONE to FALSE.  On each subsequent call, the value of
!    DONE should be left at its output value from the previous
!    call.
!
!    On output, if DONE is TRUE, then a new combination was
!    computed and returned.  If DONE is FALSE, then the list
!    of combinations was finished on the previous call.
!
!    Input, integer N, the total number of things.
!
!    Input, integer K, the number of things in each combination.
!
!    Output, integer IARRAY(K), contains the list of elements in
!    the current combination.
!
  implicit none
!
  integer k
!
  logical done
  integer i
  integer iarray(k)
  integer j
  integer n
!
  if ( done ) then
    iarray=(/(i,i=1,k)/)

    if ( k > 1 ) then
      done = .false.
    else
      done = .true.
    end if

  else

    if ( iarray(k) < n ) then
      iarray(k) = iarray(k) + 1
      return
    end if

    do i = k, 2, -1

      if ( iarray(i-1) < n-k+i-1 ) then

        iarray(i-1) = iarray(i-1) + 1

        do j = i, k
          iarray(j) = iarray(i-1) + j - ( i-1 )
        end do

        return

      end if

    end do

    done = .true.

  end if

  return
end subroutine comb_next
