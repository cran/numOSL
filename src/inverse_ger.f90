subroutine inverse_ger(mat1,n,singular)
!---------------------------------------------------------
! Subroutine inverse_ger is used for inverting a square
! matrix using Gaussian-Jordan elimation by full pivoting.
!---------------------------------------------------------
! mat1(n,n):: input/output, real values, a matrix.
!         n:: input, integer, row number of the matrix.
!  singular::  output, integer, 0=normal, 1=singular.
!---------------------------------------------------------
! Author:: Peng Jun, 2023.08.30.
!---------------------------------------------------------
! Dependence:: subroutine gjordan.
!---------------------------------------------------------
    implicit none
    ! Arguments.
    integer, intent(in):: n
    real(kind(1.0d0)), intent(inout):: mat1(n,n)
    integer, intent(out):: singular
    ! Local variables.
    real(kind(1.0d0)):: mat2(n,1)
    !
    mat2 = 1.0
    !
    call gjordan(mat1,mat2,n,1,singular)
    !
    return
    !
end subroutine inverse_ger
!
subroutine gjordan(a,b,n,m,singular)
!---------------------------------------------------------
! Subroutine gjordan is used for performing
! Gaussian-Jordan elimation by full pivoting.
!---------------------------------------------------------
! a(n,n):: input/output, real values, a square matrix.
! b(n,m):: input/output, real values, a matrix.
!      n:: input, integer, row number for matrix a.
!      m:: input, integer, column number for matrix b.
!  singular:: output, integer, 0=normal, 1=singular.
!---------------------------------------------------------
! Author:: Peng Jun, 2023.08.30.
!---------------------------------------------------------
! Dependence:: NO.
!---------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.30 OF Press et al.
! --------------------------------------------------------------------
    implicit none
    integer, intent(in):: n, m
    real(kind(1.0d0)), intent(inout):: a(n,n), b(n,m)
    integer, intent(out):: singular
    !
    integer:: i, icol, irow, j, k, l, ll,& 
              indxc(n), indxr(n), ipiv(n)
    real(kind(1.0d0)):: big, dum, pivinv
    !
    singular = 0
    !
    ipiv = 0
    !
    do i=1, n
        big = 0.0
        do j=1, n
            if (ipiv(j) .ne. 1) then
                do k=1, n
                    if (ipiv(k) .eq. 0) then
                        if (abs(a(j,k)) .ge. big) then
                            big = abs(a(j,k))
                            irow = j
                            icol = k
                        end if
                    else if (ipiv(k) .gt. 1) then
                        singular = 1
                        return
                    end if
                end do
            end if
        end do
        !
        ipiv(icol) = ipiv(icol)+1
        !
        if (irow .ne. icol) then
            do l=1, n
                dum = a(irow,l)
                a(irow,l) = a(icol,l)
                a(icol,l) = dum
            end do
            !
            do l=1, m
                dum = b(irow,l)
                b(irow,l) = b(icol,l)
                b(icol,l) = dum
            end do
        end if
        !
        indxr(i) = irow
        indxc(i) = icol
        if (a(icol,icol) .eq. 0.0) then
            singular = 1
            return
        end if
        !
        pivinv = 1.0/a(icol,icol)
        a(icol,icol) = 1.0
        !
        a(icol,1:n) = a(icol,1:n)*pivinv
        b(icol,1:m) = b(icol,1:m)*pivinv
        !
        do ll=1, n
            if(ll .ne. icol) then
                dum = a(ll,icol)
                a(ll,icol) = 0.0
                a(ll,1:n) = a(ll,1:n) - a(icol,1:n)*dum
                b(ll,1:m) = b(ll,1:m) - b(icol,1:m)*dum
            end if
        end do
    end do
    !
    do l=n, 1, -1
        if (indxr(l) .ne. indxc(l)) then
            do k=1, n
                dum = a(k,indxr(l))
                a(k,indxr(l)) = a(k,indxc(l))
                a(k,indxc(l)) = dum
            end do
        end if
    end do
    !
    return
    !
end subroutine gjordan 
