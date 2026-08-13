! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
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

!> Eigensolver for small dense symmetric matrices.
!>
!> The low-rank decomposition of the reference C6 coefficients requires a single
!> eigendecomposition of a matrix with at most a few hundred rows. A self-contained
!> cyclic Jacobi implementation avoids introducing a LAPACK dependency for it.
module dftd3_fourier_jacobi
   use mctc_env, only : wp
   implicit none
   private

   public :: symmetric_eigendecomposition

   !> Maximum number of Jacobi sweeps, convergence is usually reached in less than ten
   integer, parameter :: max_sweeps = 100


contains


!> Eigendecomposition of a dense symmetric matrix by cyclic Jacobi rotations.
!>
!> Eigenpairs are returned in order of descending magnitude of the eigenvalue,
!> which is the order required for a truncated low-rank expansion.
subroutine symmetric_eigendecomposition(amat, eval, evec)

   !> Symmetric matrix to decompose
   real(wp), intent(in) :: amat(:, :)

   !> Eigenvalues, sorted by descending magnitude
   real(wp), intent(out) :: eval(:)

   !> Eigenvectors, stored column-wise in the same order as the eigenvalues
   real(wp), intent(out) :: evec(:, :)

   integer :: ndim, isweep, ip, iq, jdim, kdim, imax
   real(wp) :: anorm, off, theta, tval, cval, sval, skip
   real(wp) :: ajp, ajq, apj, aqj, vjp, vjq, tmp
   real(wp), allocatable :: work(:, :)

   ndim = size(amat, 1)

   evec(:, :) = 0.0_wp
   do ip = 1, ndim
      evec(ip, ip) = 1.0_wp
   end do
   eval(:) = 0.0_wp

   anorm = sqrt(sum(amat**2))
   if (anorm <= 0.0_wp) return

   skip = 0.01_wp * epsilon(1.0_wp) * anorm
   work = amat

   do isweep = 1, max_sweeps
      off = 0.0_wp
      do ip = 1, ndim - 1
         off = off + sum(work(ip, ip+1:)**2)
      end do
      if (sqrt(2.0_wp*off) <= epsilon(1.0_wp) * anorm) exit

      do ip = 1, ndim - 1
         do iq = ip + 1, ndim
            if (abs(work(ip, iq)) <= skip) cycle

            theta = 0.5_wp * (work(iq, iq) - work(ip, ip)) / work(ip, iq)
            tval = sign(1.0_wp, theta) / (abs(theta) + sqrt(1.0_wp + theta*theta))
            cval = 1.0_wp / sqrt(1.0_wp + tval*tval)
            sval = tval * cval

            do jdim = 1, ndim
               ajp = work(jdim, ip)
               ajq = work(jdim, iq)
               work(jdim, ip) = cval*ajp - sval*ajq
               work(jdim, iq) = sval*ajp + cval*ajq
            end do
            do jdim = 1, ndim
               apj = work(ip, jdim)
               aqj = work(iq, jdim)
               work(ip, jdim) = cval*apj - sval*aqj
               work(iq, jdim) = sval*apj + cval*aqj
            end do
            work(ip, iq) = 0.0_wp
            work(iq, ip) = 0.0_wp

            do jdim = 1, ndim
               vjp = evec(jdim, ip)
               vjq = evec(jdim, iq)
               evec(jdim, ip) = cval*vjp - sval*vjq
               evec(jdim, iq) = sval*vjp + cval*vjq
            end do
         end do
      end do
   end do

   do ip = 1, ndim
      eval(ip) = work(ip, ip)
   end do

   do ip = 1, ndim - 1
      imax = ip
      do iq = ip + 1, ndim
         if (abs(eval(iq)) > abs(eval(imax))) imax = iq
      end do
      if (imax == ip) cycle
      tmp = eval(ip)
      eval(ip) = eval(imax)
      eval(imax) = tmp
      do kdim = 1, ndim
         tmp = evec(kdim, ip)
         evec(kdim, ip) = evec(kdim, imax)
         evec(kdim, imax) = tmp
      end do
   end do

end subroutine symmetric_eigendecomposition


end module dftd3_fourier_jacobi
