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

module d3def_options
   use iso_fortran_env, only: wp => real64
   use d3def_damping_parameters
   implicit none

   type :: d3_options
      character(len=:), allocatable :: func
      type(d3_damping_parameters), allocatable :: par
      real(wp) :: weighting_factor = 4.0_wp
      real(wp) :: cutoff_disp = 64.0_wp
   end type d3_options

contains

end module d3def_options
