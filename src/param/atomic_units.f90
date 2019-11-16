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

module d3par_atomic_units
   use iso_fortran_env, only: wp => real64
   use d3par_units
   implicit none

   ! conversion factor from bohr to Ångström
   real(wp),parameter :: autoaa = bohr * 1e10_wp
   real(wp),parameter :: aatoau = 1.0_wp/autoaa

   ! conversion factor from hartree to electron volts
   real(wp),parameter :: autoeV = hartree/e
   real(wp),parameter :: eVtoau = 1.0_wp/autoeV

   ! coversion factor between calorine and joule
   real(wp),parameter :: caltoj = 4.184_wp
   real(wp),parameter :: jtocal = 1.0_wp/caltoj

   ! conversion from hartree to kJ/mol
   real(wp),parameter :: autokJ = hartree*NA*1e-3_wp
   real(wp),parameter :: kJtoau = 1.0_wp/autokJ

   ! conversion from hartree to kcal/mol
   real(wp),parameter :: autokcal = autokJ*Jtocal
   real(wp),parameter :: kcaltoau = 1.0_wp/autokcal

   ! conversion from hartree to reciprocal centimeters
   real(wp),parameter :: autorcm = hartree/(h*c)*1e-2_wp
   real(wp),parameter :: rcmtoau = 1.0_wp/autorcm

   ! conversion from hartree to nanometers (wavelength)
   real(wp),parameter :: autonm = h*c/hartree * 1e+9_wp
   real(wp),parameter :: nmtoau = 1.0_wp/autonm

   ! conversion from electron mass (a.u.) to kg
   real(wp),parameter :: autokg = me
   real(wp),parameter :: kgtoau = 1.0_wp/autokg

   ! molecular mass per mole (g/mol) to electron mass (a.u.)
   real(wp),parameter :: autogmol = me*NA*1e+3_wp
   real(wp),parameter :: gmoltoau = 1.0_wp/autogmol

   ! molecular mass per mole (g/mol) to kg
   real(wp),parameter :: gmoltokg = gmoltoau*autokg
   real(wp),parameter :: kgtogmol = 1.0_wp/gmoltokg

   ! Coulomb to atomic charge units
   real(wp),parameter :: autoc = e
   real(wp),parameter :: ctoau = 1.0_wp/autoc

end module d3par_atomic_units
