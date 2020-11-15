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

module dftd3_cutoff
   use mctc_env, only : wp
   implicit none

   public :: realspace_cutoff, get_lattice_points


   !> Coordination number cutoff
   real(wp), parameter :: cn_default = 40.0_wp

   !> Two-body interaction cutoff
   real(wp), parameter :: disp2_default = 60.0_wp

   !> Three-body interaction cutoff
   real(wp), parameter :: disp3_default = 40.0_wp


   !> Collection of real space cutoffs
   type :: realspace_cutoff
      sequence

      !> Coordination number cutoff
      real(wp) :: cn = cn_default

      !> Two-body interaction cutoff
      real(wp) :: disp2 = disp2_default

      !> Three-body interaction cutoff
      real(wp) :: disp3 = disp3_default

   end type realspace_cutoff


contains


subroutine get_lattice_points(periodic, lat, rthr, trans)
   logical, intent(in) :: periodic(:)
   real(wp), intent(in) :: rthr
   real(wp), intent(in) :: lat(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer :: rep(3)
   integer :: itr, ix, iy, iz

   if (.not.any(periodic)) then
      allocate(trans(3, 1))
      trans(:, :) = 0.0_wp
   else
      call get_translations(lat, rthr, rep)
      allocate(trans(3, product(2*rep+1)))
      itr = 0
      do ix = -rep(1), rep(1)
         do iy = -rep(2), rep(2)
            do iz = -rep(3), rep(3)
               itr = itr + 1
               trans(:, itr) = lat(:, 1)*ix + lat(:, 2)*iy + lat(:, 3)*iz
            end do
         end do
      end do
   end if

end subroutine get_lattice_points


!> generate a supercell based on a realspace cutoff, this subroutine
!  doesn't know anything about the convergence behaviour of the
!  associated property.
pure subroutine get_translations(lat, rthr, rep)
   real(wp), intent(in)  :: rthr
   real(wp), intent(in)  :: lat(3, 3)
   integer, intent(out) :: rep(3)
   real(wp) :: normx(3), normy(3), normz(3)
   real(wp) :: cos10, cos21, cos32

   ! find normal to the plane...
   call crossproduct(lat(:, 2), lat(:, 3), normx)
   call crossproduct(lat(:, 3), lat(:, 1), normy)
   call crossproduct(lat(:, 1), lat(:, 2), normz)
   ! ...normalize it...
   normx = normx/norm2(normx)
   normy = normy/norm2(normy)
   normz = normz/norm2(normz)
   ! cos angles between normals and lattice vectors
   cos10 = sum(normx*lat(:, 1))
   cos21 = sum(normy*lat(:, 2))
   cos32 = sum(normz*lat(:, 3))
   rep(1) = ceiling(abs(rthr/cos10))
   rep(2) = ceiling(abs(rthr/cos21))
   rep(3) = ceiling(abs(rthr/cos32))

contains

   pure subroutine crossproduct(a, b, c)
      real(wp), intent(in)  :: a(3), b(3)
      real(wp), intent(out) :: c(3)
      real(wp) :: x, y, z

      x=a(2)*b(3)-b(2)*a(3)
      y=a(3)*b(1)-b(3)*a(1)
      z=a(1)*b(2)-b(1)*a(2)
      c=(/x, y, z/)
   end subroutine crossproduct

end subroutine get_translations


end module dftd3_cutoff
