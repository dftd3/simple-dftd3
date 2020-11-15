! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later AND BSD-3-Clause
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

!> Interface to BLAS library, in case DFT-D3 is not linked against BLAS this
!> module provides redistributed reference BLAS implementations instead.
module dftd3_blas
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: d3_gemv, blas_gemv


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface d3_gemv
      module procedure :: d3_sgemv
      module procedure :: d3_dgemv
      module procedure :: d3_sgemv312
      module procedure :: d3_sgemv321
      module procedure :: d3_dgemv312
      module procedure :: d3_dgemv321
   end interface d3_gemv


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface blas_gemv
#ifdef WITH_BLAS
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgemv
#else
      module procedure :: sgemv
      module procedure :: dgemv
#endif
   end interface blas_gemv


contains


subroutine d3_sgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call d3_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine d3_sgemv312


subroutine d3_sgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call d3_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine d3_sgemv321


subroutine d3_dgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call d3_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine d3_dgemv312


subroutine d3_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call d3_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine d3_dgemv321


pure subroutine d3_sgemv(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine d3_sgemv


pure subroutine d3_dgemv(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine d3_dgemv


#ifndef WITH_BLAS
!> SPDX-Identifer: BSD-3-Clause
!> Reference BLAS level2 routine (version 3.7.0)
!> -- Reference BLAS is a software package provided by Univ. of Tennessee,
!> -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!>    December 2016
pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
   ! Scalar Arguments
   real(sp), intent(in) :: alpha
   real(sp), intent(in) :: beta
   integer, intent(in) :: incx
   integer, intent(in) :: incy
   integer, intent(in) :: lda
   integer, intent(in) :: m
   integer, intent(in) :: n
   character, intent(in) :: trans
   ! Array Arguments
   real(sp), intent(in) :: a(lda, *)
   real(sp), intent(in) :: x(*)
   real(sp), intent(inout) :: y(*)

   ! Parameters
   real(sp) :: one, zero
   parameter (one=1.0e+0_sp, zero=0.0e+0_sp)

   ! Local Scalars
   real(sp) :: temp
   integer :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny

   ! Intrinsic Functions
   intrinsic :: max

   ! Test the input parameters.
   info = 0
   if (.not.lsame(trans, 'N') .and. .not.lsame(trans, 'T') .and. &
      &.not.lsame(trans, 'C')) then
      info = 1
   else if (m.lt.0) then
      info = 2
   else if (n.lt.0) then
      info = 3
   else if (lda.lt.max(1, m)) then
      info = 6
   else if (incx.eq.0) then
      info = 8
   else if (incy.eq.0) then
      info = 11
   end if
   if (info.ne.0) then
      error stop 'dgemv'
      return
   end if

   ! Quick return if possible.
   if ((m.eq.0) .or. (n.eq.0) .or. &
      &((alpha.eq.zero).and. (beta.eq.one))) return

   ! Set  LENX  and  LENY, the lengths of the vectors x and y, and set
   ! up the start points in  X  and  Y.
   if (lsame(trans, 'N')) then
      lenx = n
      leny = m
   else
      lenx = m
      leny = n
   end if
   if (incx.gt.0) then
      kx = 1
   else
      kx = 1 - (lenx-1)*incx
   end if
   if (incy.gt.0) then
      ky = 1
   else
      ky = 1 - (leny-1)*incy
   end if

   ! Start the operations. In this version the elements of A are
   ! accessed sequentially with one pass through A.

   ! First form  y := beta*y.
   if (beta.ne.one) then
      if (incy.eq.1) then
         if (beta.eq.zero) then
            do i = 1, leny
               y(i) = zero
            end do
         else
            do i = 1, leny
               y(i) = beta*y(i)
            end do
         end if
      else
         iy = ky
         if (beta.eq.zero) then
            do i = 1, leny
               y(iy) = zero
               iy = iy + incy
            end do
         else
            do i = 1, leny
               y(iy) = beta*y(iy)
               iy = iy + incy
            end do
         end if
      end if
   end if
   if (alpha.eq.zero) return
   if (lsame(trans, 'N')) then
      ! Form  y := alpha*A*x + y.
      jx = kx
      if (incy.eq.1) then
         do j = 1, n
            temp = alpha*x(jx)
            do i = 1, m
               y(i) = y(i) + temp*a(i, j)
            end do
            jx = jx + incx
         end do
      else
         do j = 1, n
            temp = alpha*x(jx)
            iy = ky
            do i = 1, m
               y(iy) = y(iy) + temp*a(i, j)
               iy = iy + incy
            end do
            jx = jx + incx
         end do
      end if
   else
      ! Form  y := alpha*A**T*x + y.
      jy = ky
      if (incx.eq.1) then
         do j = 1, n
            temp = zero
            do i = 1, m
               temp = temp + a(i, j)*x(i)
            end do
            y(jy) = y(jy) + alpha*temp
            jy = jy + incy
         end do
      else
         do j = 1, n
            temp = zero
            ix = kx
            do i = 1, m
               temp = temp + a(i, j)*x(ix)
               ix = ix + incx
            end do
            y(jy) = y(jy) + alpha*temp
            jy = jy + incy
         end do
      end if
   end if

end subroutine sgemv


!> SPDX-Identifer: BSD-3-Clause
!> Reference BLAS level2 routine (version 3.7.0)
!> -- Reference BLAS is a software package provided by Univ. of Tennessee,
!> -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!>    December 2016
pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
   ! Scalar Arguments
   real(dp), intent(in) :: alpha
   real(dp), intent(in) :: beta
   integer, intent(in) :: incx
   integer, intent(in) :: incy
   integer, intent(in) :: lda
   integer, intent(in) :: m
   integer, intent(in) :: n
   character, intent(in) :: trans
   ! Array Arguments
   real(dp), intent(in) :: a(lda, *)
   real(dp), intent(in) :: x(*)
   real(dp), intent(inout) :: y(*)

   ! Parameters
   real(dp) :: one, zero
   parameter (one=1.0e+0_dp, zero=0.0e+0_dp)

   ! Local Scalars
   real(dp) :: temp
   integer :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny

   ! Intrinsic Functions
   intrinsic :: max

   ! Test the input parameters.
   info = 0
   if (.not.lsame(trans, 'N') .and. .not.lsame(trans, 'T') .and. &
      &.not.lsame(trans, 'C')) then
      info = 1
   else if (m.lt.0) then
      info = 2
   else if (n.lt.0) then
      info = 3
   else if (lda.lt.max(1, m)) then
      info = 6
   else if (incx.eq.0) then
      info = 8
   else if (incy.eq.0) then
      info = 11
   end if
   if (info.ne.0) then
      error stop 'dgemv'
      return
   end if

   ! Quick return if possible.
   if ((m.eq.0) .or. (n.eq.0) .or. &
      &((alpha.eq.zero).and. (beta.eq.one))) return

   ! Set  LENX  and  LENY, the lengths of the vectors x and y, and set
   ! up the start points in  X  and  Y.
   if (lsame(trans, 'N')) then
      lenx = n
      leny = m
   else
      lenx = m
      leny = n
   end if
   if (incx.gt.0) then
      kx = 1
   else
      kx = 1 - (lenx-1)*incx
   end if
   if (incy.gt.0) then
      ky = 1
   else
      ky = 1 - (leny-1)*incy
   end if

   ! Start the operations. In this version the elements of A are
   ! accessed sequentially with one pass through A.

   ! First form  y := beta*y.
   if (beta.ne.one) then
      if (incy.eq.1) then
         if (beta.eq.zero) then
            do i = 1, leny
               y(i) = zero
            end do
         else
            do i = 1, leny
               y(i) = beta*y(i)
            end do
         end if
      else
         iy = ky
         if (beta.eq.zero) then
            do i = 1, leny
               y(iy) = zero
               iy = iy + incy
            end do
         else
            do i = 1, leny
               y(iy) = beta*y(iy)
               iy = iy + incy
            end do
         end if
      end if
   end if
   if (alpha.eq.zero) return
   if (lsame(trans, 'N')) then
      ! Form  y := alpha*A*x + y.
      jx = kx
      if (incy.eq.1) then
         do j = 1, n
            temp = alpha*x(jx)
            do i = 1, m
               y(i) = y(i) + temp*a(i, j)
            end do
            jx = jx + incx
         end do
      else
         do j = 1, n
            temp = alpha*x(jx)
            iy = ky
            do i = 1, m
               y(iy) = y(iy) + temp*a(i, j)
               iy = iy + incy
            end do
            jx = jx + incx
         end do
      end if
   else
      ! Form  y := alpha*A**T*x + y.
      jy = ky
      if (incx.eq.1) then
         do j = 1, n
            temp = zero
            do i = 1, m
               temp = temp + a(i, j)*x(i)
            end do
            y(jy) = y(jy) + alpha*temp
            jy = jy + incy
         end do
      else
         do j = 1, n
            temp = zero
            ix = kx
            do i = 1, m
               temp = temp + a(i, j)*x(ix)
               ix = ix + incx
            end do
            y(jy) = y(jy) + alpha*temp
            jy = jy + incy
         end do
      end if
   end if

end subroutine dgemv


!> SPDX-Identifer: BSD-3-Clause
!> Reference BLAS level1 routine (version 3.1)
!> -- Reference BLAS is a software package provided by Univ. of Tennessee,
!> -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!>    December 2016
pure function lsame(ca, cb)
   character, intent(in) :: ca
   character, intent(in) :: cb
   logical :: lsame

   !> Intrinsic Functions
   intrinsic :: ichar
   !> Local Scalars
   integer :: inta, intb, zcode

   ! Test if the characters are equal
   lsame = ca .eq. cb
   if (lsame) return

   ! Now test for equivalence if both characters are alphabetic.
   zcode = ichar('Z')

   ! Use 'Z' rather than 'A' so that ASCII can be detected on Prime
   ! machines, on which ICHAR returns a value with bit 8 set.
   ! ICHAR('A') on Prime machines returns 193 which is the same as
   ! ICHAR('A') on an EBCDIC machine.
   inta = ichar(ca)
   intb = ichar(cb)

   if (zcode.eq.90 .or. zcode.eq.122) then
      ! ASCII is assumed - ZCODE is the ASCII code of either lower or
      ! upper case 'Z'.
      if (inta.ge.97 .and. inta.le.122) inta = inta - 32
      if (intb.ge.97 .and. intb.le.122) intb = intb - 32

   else if (zcode.eq.233 .or. zcode.eq.169) then
      ! EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
      ! upper case 'Z'.
      if (inta.ge.129 .and. inta.le.137 .or. &
         &inta.ge.145 .and. inta.le.153 .or. &
         &inta.ge.162 .and. inta.le.169) inta = inta + 64
      if (intb.ge.129 .and. intb.le.137 .or. &
         &intb.ge.145 .and. intb.le.153 .or. &
         &intb.ge.162 .and. intb.le.169) intb = intb + 64

   else if (zcode.eq.218 .or. zcode.eq.250) then
      ! ASCII is assumed, on Prime machines - ZCODE is the ASCII code
      ! plus 128 of either lower or upper case 'Z'.
      if (inta.ge.225 .and. inta.le.250) inta = inta - 32
      if (intb.ge.225 .and. intb.le.250) intb = intb - 32
   end if
   lsame = inta .eq. intb

end function lsame
#endif


end module dftd3_blas
