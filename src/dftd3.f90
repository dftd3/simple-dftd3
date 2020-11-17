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

module dftd3
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_disp, only : get_dispersion
   use dftd3_ncoord, only : get_coordination_number
   use dftd3_damping, only : damping_param
   use dftd3_damping_rational, only : rational_damping_param, new_rational_damping
   use dftd3_damping_zero, only : zero_damping_param, new_zero_damping
   use dftd3_model, only : d3_model, new_d3_model
   use dftd3_param, only : d3_param, get_rational_damping, get_zero_damping
   use dftd3_version, only : get_dftd3_version
   implicit none
   private

   public :: get_dispersion, get_coordination_number, realspace_cutoff
   public :: damping_param, d3_param
   public :: get_rational_damping, get_zero_damping
   public :: rational_damping_param, new_rational_damping
   public :: zero_damping_param, new_zero_damping
   public :: d3_model, new_d3_model
   public :: get_dftd3_version


end module dftd3
