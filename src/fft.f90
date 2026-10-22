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

!> Complex fast Fourier transform on a three-dimensional mesh.
!>
!> Self-contained radix-2/3 Stockham transform. Allowing products of powers of
!> two and three avoids padding the particle mesh to the next power of two.
module dftd3_fourier_fft
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: fft_mesh, new_fft_mesh, fft_3d, valid_mesh_size


   !> Twiddle factors for the transforms along the three mesh directions
   type :: fft_mesh

      !> Number of points along each direction
      integer :: mesh(3) = 0

      !> Roots of unity for each direction, whose lengths need not divide each other
      complex(wp), allocatable :: tw(:, :)

   end type fft_mesh


contains


!> Smallest supported mesh size not below the requested number of points
pure function valid_mesh_size(npoint) result(nmesh)

   !> Requested number of mesh points
   integer, intent(in) :: npoint

   !> Supported number of mesh points
   integer :: nmesh
   integer :: pow3, candidate

   nmesh = 1
   do while (nmesh < npoint)
      nmesh = 2 * nmesh
   end do

   pow3 = 3
   do while (pow3 < nmesh)
      candidate = pow3
      do while (candidate < npoint)
         candidate = 2 * candidate
      end do
      nmesh = min(nmesh, candidate)
      pow3 = 3 * pow3
   end do

end function valid_mesh_size


!> Tabulate the roots of unity for transforms on the given mesh
subroutine new_fft_mesh(self, mesh)

   !> Instance of the transform
   type(fft_mesh), intent(out) :: self

   !> Number of points along each direction
   integer, intent(in) :: mesh(3)

   integer :: ik, idir

   self%mesh(:) = mesh
   allocate(self%tw(0:maxval(mesh) - 1, 3))
   do idir = 1, 3
      do ik = 0, mesh(idir) - 1
         self%tw(ik, idir) = cmplx(cos(2*pi*ik/mesh(idir)), -sin(2*pi*ik/mesh(idir)), wp)
      end do
   end do

end subroutine new_fft_mesh


!> Transform along all three directions, isign selects forward (-1) or backward (+1)
subroutine fft_3d(self, grid, isign)

   !> Instance of the transform
   type(fft_mesh), intent(in) :: self

   !> Mesh to transform in place
   complex(wp), intent(inout) :: grid(:, :, :)

   !> Sign of the exponent
   integer, intent(in) :: isign

   integer :: i1, i2, i3
   complex(wp), allocatable :: work(:, :)

   !$omp parallel default(none) shared(self, grid, isign) private(i1, i2, i3, work)
   allocate(work(size(self%tw, 1), 2))

   !$omp do schedule(runtime) collapse(2)
   do i3 = 1, self%mesh(3)
      do i2 = 1, self%mesh(2)
         call fft_1d(grid(:, i2, i3), self%tw(:, 1), work, isign)
      end do
   end do
   !$omp end do

   !$omp do schedule(runtime) collapse(2)
   do i3 = 1, self%mesh(3)
      do i1 = 1, self%mesh(1)
         call fft_1d(grid(i1, :, i3), self%tw(:, 2), work, isign)
      end do
   end do
   !$omp end do

   !$omp do schedule(runtime) collapse(2)
   do i2 = 1, self%mesh(2)
      do i1 = 1, self%mesh(1)
         call fft_1d(grid(i1, i2, :), self%tw(:, 3), work, isign)
      end do
   end do
   !$omp end do

   deallocate(work)
   !$omp end parallel

end subroutine fft_3d


!> Autosorting transform using two contiguous work vectors shared across stages
pure subroutine fft_1d(zvec, tw, work, isign)

   !> Sequence to transform
   complex(wp), intent(inout) :: zvec(:)

   !> Roots of unity for this direction
   complex(wp), intent(in) :: tw(0:)

   !> Per-thread scratch space, reused for every mesh line
   complex(wp), intent(inout) :: work(:, :)

   !> Sign of the exponent
   integer, intent(in) :: isign

   integer :: nn, stride, radix, blocks, jj, kk, ii, io, src, dst
   complex(wp) :: wval, wval2, aa, bb, cc, dd, omega

   nn = size(zvec)
   if (nn <= 1) return

   work(:nn, 1) = zvec
   omega = cmplx(0.0_wp, isign * sqrt(3.0_wp) / 2, wp)
   src = 1
   stride = 1
   do while (stride < nn)
      radix = 3
      if (mod(nn / stride, 2) == 0) radix = 2
      blocks = nn / (radix * stride)
      dst = 3 - src
      do jj = 0, stride - 1
         wval = tw(jj * blocks)
         if (isign > 0) wval = conjg(wval)
         if (radix == 3) wval2 = wval * wval
         do kk = 0, blocks - 1
            ii = jj + kk * stride + 1
            io = jj + kk * radix * stride + 1
            aa = work(ii, src)
            bb = wval * work(ii + nn/radix, src)
            if (radix == 2) then
               work(io, dst) = aa + bb
               work(io + stride, dst) = aa - bb
            else
               cc = wval2 * work(ii + 2*nn/radix, src)
               dd = omega * (bb - cc)
               work(io, dst) = aa + bb + cc
               aa = aa - 0.5_wp * (bb + cc)
               work(io + stride, dst) = aa + dd
               work(io + 2*stride, dst) = aa - dd
            end if
         end do
      end do
      src = dst
      stride = radix * stride
   end do
   zvec = work(:nn, src)

end subroutine fft_1d


end module dftd3_fourier_fft
