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
  integer,intent(in)::n
  real(kind(1.0d0)),intent(inout)::x(n)
  logical,intent(in)::upper
  !
  ! local variables
  real(kind(1.0d0))::alnorm
  integer::i
  !
  do i=1,n
    x(i)=alnorm(x(i),upper)
  end do
  !
  return
end subroutine pnorm
!
function alnorm ( x, upper )
!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  integer, parameter:: dbdbx=kind(1.0d0)
  
  real ( dbdbx ), parameter :: a1 = 5.75885480458D+00
  real ( dbdbx ), parameter :: a2 = 2.62433121679D+00
  real ( dbdbx ), parameter :: a3 = 5.92885724438D+00
  real ( dbdbx ) alnorm
  real ( dbdbx ), parameter :: b1 = -29.8213557807D+00
  real ( dbdbx ), parameter :: b2 = 48.6959930692D+00
  real ( dbdbx ), parameter :: c1 = -0.000000038052D+00
  real ( dbdbx ), parameter :: c2 = 0.000398064794D+00
  real ( dbdbx ), parameter :: c3 = -0.151679116635D+00
  real ( dbdbx ), parameter :: c4 = 4.8385912808D+00
  real ( dbdbx ), parameter :: c5 = 0.742380924027D+00
  real ( dbdbx ), parameter :: c6 = 3.99019417011D+00
  real ( dbdbx ), parameter :: con = 1.28D+00
  real ( dbdbx ), parameter :: d1 = 1.00000615302D+00
  real ( dbdbx ), parameter :: d2 = 1.98615381364D+00
  real ( dbdbx ), parameter :: d3 = 5.29330324926D+00
  real ( dbdbx ), parameter :: d4 = -15.1508972451D+00
  real ( dbdbx ), parameter :: d5 = 30.789933034D+00
  real ( dbdbx ), parameter :: ltone = 7.0D+00
  real ( dbdbx ), parameter :: p = 0.398942280444D+00
  real ( dbdbx ), parameter :: q = 0.39990348504D+00
  real ( dbdbx ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( dbdbx ), parameter :: utzero = 18.66D+00
  real ( dbdbx ) x
  real ( dbdbx ) y
  real ( dbdbx ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end function alnorm
