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

module d3par_units
   use iso_fortran_env, only: wp => real64
   use d3par_constants, only: pi
   implicit none

   ! ========================================================================
   ! Definition of SI base units => by defining natural contants
   ! ========================================================================
   !> Planck's constant
   real(wp),parameter :: h = 6.6260715e-34_wp ! J·s = kg·m²·s⁻¹
   real(wp),parameter :: hbar = h/(2.0_wp*pi) ! J·s = kg·m²·s⁻¹

   !> speed of light in vacuum
   real(wp),parameter :: c = 299792458.0_wp ! m·s⁻¹

   !> Boltzmann's constant
   real(wp),parameter :: kb = 1.380649e-23_wp ! J·K⁻¹ = kg·m²·s⁻²·K⁻¹

   !> Avogadro's number
   real(wp),parameter :: NA = 6.02214076e23_wp ! mol⁻¹

   !> elementary charge
   real(wp),parameter :: e = 1.602176634e-19_wp ! C

   ! ========================================================================
   !  Natural constants
   ! ========================================================================
   !> fine structure constant (CODATA2018)
   real(wp),parameter :: alpha = 1.0_wp/137.035999046_wp ! dimensionless
   !> electron rest mass
   real(wp),parameter :: me = 9.10938356e-31_wp ! kg

   ! ========================================================================
   ! derived natural constants
   ! ========================================================================
   !> vacuum permeability
   real(wp),parameter :: mu0 = 2.0_wp*alpha/e**2 * h/c ! H·m⁻¹ = kg·m·s⁻²·A⁻²
   !> vacuum permittivity
   real(wp),parameter :: epsilon0 = e**2/(2.0_wp*alpha*h*c) ! F·m⁻¹ = s⁴·A²·m⁻³·kg⁻¹

   !> Coulomb's constant
   real(wp),parameter :: ke = alpha*hbar*c/e**2 ! N·m²·C⁻²

   ! ========================================================================
   !  Hartree units
   !  Defined by setting
   !  -> electron mass (me) to 1
   !  -> elementary charge (e) to 1
   !  -> reduced Planck's constant (hbar) to 1
   !  -> Coulomb's constant (ke) to 1
   ! ========================================================================
   !> Bohr radius
   real(wp),parameter :: bohr = hbar/(me*c*alpha) ! m
   !> Hartree energy
   real(wp),parameter :: hartree = me * c**2 * alpha**2 ! J = kg·m²·s⁻²

end module d3par_units
