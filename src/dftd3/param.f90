! This file is part of s-dftd3.
! SPDX-Identifier: LGLP-3.0-or-later
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

module dftd3_param
   use mctc_env, only : wp
   implicit none

   public :: d3_param


   type :: d3_param
      real(wp) :: s6 = 1.0_wp
      real(wp) :: s8 = 1.0_wp
      real(wp) :: s9 = 0.0_wp
      real(wp) :: rs6 = 1.0_wp
      real(wp) :: rs8 = 1.0_wp
      real(wp) :: a1 = 0.4_wp
      real(wp) :: a2 = 5.0_wp
      real(wp) :: alp = 14.0_wp
   end type d3_param


end module dftd3_param
