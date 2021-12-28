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

module dftd3_app_help
   use dftd3, only : get_dftd3_version
   implicit none
   private

   public :: prog_name, header, help, version


   character(len=*), parameter :: prog_name = "s-dftd3"


contains


subroutine header(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd3_version(string=version_string)
   write(unit, '(a)') &
      "-----------------------------------", &
      " s i m p l e   D F T - D 3  v"// version_string, &
      "-----------------------------------", ""

end subroutine header


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Takes an geometry input to calculate the D3 dispersion correction.", &
      "Periodic calculations are performed automatically for periodic input formats.", &
      "Specify the functional to select the correct parameters.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-i, --input <format>", "Hint for the format of the input file", &
      "--bj <method>", "Use rational (Becke-Johnson) damping function", &
      "--bj-param <list>", "Specify parameters for rational damping,", &
      "", "expected order is s6, s8, a1, a2 (requires four arguments)", &
      "--zero <method>", "Use zero (Chai-Head-Gordon) damping function", &
      "--zero-param <list>", "Specify parameters for zero damping,", &
      "", "expected order is s6, s8, rs6 (requires three arguments)", &
      "--bjm <method>", "Use modified rational damping function", &
      "--bjm-param <list>", "Specify parameters for rational damping,", &
      "", "expected order is s6, s8, a1, a2 (requires four arguments)", &
      "--zerom <method>", "Use modified zero damping function", &
      "--zerom-param <list>", "Specify parameters for modified zero damping,", &
      "", "expected order is s6, s8, rs6, bet (requires four arguments)", &
      "--op <method>", "Use optimized power damping function", &
      "--op-param <list>", "Specify parameters for optimized power,", &
      "", "expected order is s6, s8, a1, a2, bet (requires five arguments)", &
      "--atm", "Use ATM three-body dispersion", &
      "--atm-scale <s9>", "Use scaled ATM three-body dispersion", &
      "--noedisp", "Disable writing of dispersion energy to .EDISP file", &
      "--json [file]", "Dump results to JSON output (default: dftd3.json)", &
      "--grad [file]", "Request gradient evaluation,", &
      "", "write results to file (default: dftd3.txt),", &
      "", "attempts to add to Turbomole gradient and gradlatt files", &
      "--property", "Evaluate dispersion related properties", &
      "--pair-resolved", "Calculate pairwise representation of dispersion energy", &
      "-v, --verbose", "Show more, can be used multiple times", &
      "-s, --silent", "Show less, use twice to supress all output", &
      "--version", "Print program version and exit", &
      "--help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd3_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


end module dftd3_app_help
