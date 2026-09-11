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
!> Self-contained radix-2 implementation, sufficient for the particle mesh
!> evaluation of the dispersion energy which is free to choose its mesh size.
!> Supports powers of two only, a mixed-radix backend would lift the mesh
!> quantisation but only matters once the transform dominates the runtime.
module dftd3_fourier_fft
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: fft_mesh, new_fft_mesh, fft_3d, valid_mesh_size


   !> Twiddle factors shared by the transforms along all three mesh directions
   type :: fft_mesh

      !> Number of points along each direction
      integer :: mesh(3) = 0

      !> Length of the twiddle table, the largest mesh dimension
      integer :: nmax = 0

      !> Roots of unity exp(-2*pi*i*k/nmax) for k = 0 .. nmax/2 - 1
      complex(wp), allocatable :: tw(:)

   end type fft_mesh


contains


!> Smallest supported mesh size not below the requested number of points
pure function valid_mesh_size(npoint) result(nmesh)

   !> Requested number of mesh points
   integer, intent(in) :: npoint

   !> Supported number of mesh points
   integer :: nmesh

   nmesh = 1
   do while (nmesh < npoint)
      nmesh = 2 * nmesh
   end do

end function valid_mesh_size


!> Tabulate the roots of unity for transforms on the given mesh
subroutine new_fft_mesh(self, mesh)

   !> Instance of the transform
   type(fft_mesh), intent(out) :: self

   !> Number of points along each direction
   integer, intent(in) :: mesh(3)

   integer :: ik

   self%mesh(:) = mesh
   self%nmax = maxval(mesh)

   allocate(self%tw(0:self%nmax/2 - 1))
   do ik = 0, self%nmax/2 - 1
      self%tw(ik) = cmplx(cos(2*pi*ik/self%nmax), -sin(2*pi*ik/self%nmax), wp)
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

   !$omp parallel default(none) shared(self, grid, isign) private(i1, i2, i3)

   !$omp do schedule(runtime) collapse(2)
   do i3 = 1, self%mesh(3)
      do i2 = 1, self%mesh(2)
         call fft_1d(grid(:, i2, i3), self%tw, self%nmax, isign)
      end do
   end do
   !$omp end do

   !$omp do schedule(runtime) collapse(2)
   do i3 = 1, self%mesh(3)
      do i1 = 1, self%mesh(1)
         call fft_1d(grid(i1, :, i3), self%tw, self%nmax, isign)
      end do
   end do
   !$omp end do

   !$omp do schedule(runtime) collapse(2)
   do i2 = 1, self%mesh(2)
      do i1 = 1, self%mesh(1)
         call fft_1d(grid(i1, i2, :), self%tw, self%nmax, isign)
      end do
   end do
   !$omp end do

   !$omp end parallel

end subroutine fft_3d


!> In-place radix-2 decimation in time transform of a single sequence
pure subroutine fft_1d(zvec, tw, nmax, isign)

   !> Sequence to transform
   complex(wp), intent(inout) :: zvec(:)

   !> Roots of unity for the largest mesh dimension
   complex(wp), intent(in) :: tw(0:)

   !> Length of the twiddle table
   integer, intent(in) :: nmax

   !> Sign of the exponent
   integer, intent(in) :: isign

   integer :: nn, ii, jj, kk, mm, istep, itw
   complex(wp) :: wval, temp

   nn = size(zvec)
   if (nn <= 1) return

   jj = 1
   do ii = 1, nn
      if (jj > ii) then
         temp = zvec(jj)
         zvec(jj) = zvec(ii)
         zvec(ii) = temp
      end if
      kk = nn / 2
      do while (kk >= 2 .and. jj > kk)
         jj = jj - kk
         kk = kk / 2
      end do
      jj = jj + kk
   end do

   mm = 1
   do while (mm < nn)
      istep = 2 * mm
      do kk = 1, mm
         ! every supported length divides nmax, so the table is shared
         itw = (kk - 1) * (nmax / istep)
         wval = tw(itw)
         if (isign > 0) wval = conjg(wval)
         do ii = kk, nn, istep
            temp = wval * zvec(ii + mm)
            zvec(ii + mm) = zvec(ii) - temp
            zvec(ii) = zvec(ii) + temp
         end do
      end do
      mm = istep
   end do

end subroutine fft_1d


end module dftd3_fourier_fft
