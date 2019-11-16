! This file is part of s-dftd3.
!
! Copyright (C) 2019 Sebastian Ehlert
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module d3def_damping_parameters
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: d3_damping_parameters
   private


   type :: d3_damping_parameters
      real(wp) :: s6 = -1.0_wp
      real(wp) :: s8 = -1.0_wp
      real(wp) :: s10 = 0.0_wp
      real(wp) :: s9 = 1.0_wp
      real(wp) :: a1 = -1.0_wp
      real(wp) :: a2 = -1.0_wp
      integer :: alp = 16
   end type d3_damping_parameters


end module d3def_damping_parameters
