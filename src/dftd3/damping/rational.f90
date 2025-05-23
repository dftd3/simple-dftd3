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

module dftd3_damping_rational
   use dftd3_damping, only : damping_param
   use dftd3_damping_atm, only : get_atm_dispersion, get_atm_pairwise_dispersion
   use dftd3_param, only : d3_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: rational_damping_param, new_rational_damping


   !> Rational (Becke-Johnson) damping model
   type, extends(damping_param) :: rational_damping_param
      real(wp) :: s6
      real(wp) :: s8
      real(wp) :: s9
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: alp
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

      !> Evaluate pairwise representation of additive dispersion energy
      procedure :: get_pairwise_dispersion2

      !> Evaluate pairwise representation of non-additive dispersion energy
      procedure :: get_pairwise_dispersion3

   end type rational_damping_param


   real(wp), parameter :: rs9 = 4.0_wp/3.0_wp


contains


!> Create new rational damping model
subroutine new_rational_damping(self, param)

   !> Rational damping parameters
   type(rational_damping_param), intent(out) :: self

   !> Parameters
   type(d3_param), intent(in) :: param

   self%s6 = param%s6
   self%s8 = param%s8
   self%s9 = param%s9
   self%a1 = param%a1
   self%a2 = param%a2
   self%alp = param%alp

end subroutine new_rational_damping


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)

   if (grad) then
      call get_dispersion_derivs(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
         & energy, dEdcn, gradient, sigma)
   else
      call get_dispersion_energy(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)
   end if

end subroutine get_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_dispersion_energy(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, edisp, dE

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, &
   !$omp& t8, edisp, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = self%s6*t6 + self%s8*rrij*t8

            dE = -c6ij*edisp * 0.5_wp

            energy_local(iat) = energy_local(iat) + dE
            if (iat /= jat) then
               energy_local(jat) = energy_local(jat) + dE
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_dispersion_energy_)
   energy(:) = energy(:) + energy_local(:)
   !$omp end critical (get_dispersion_energy_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_dispersion_energy


!> Evaluation of the dispersion energy expression
subroutine get_dispersion_derivs(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in) :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, d6, d8, edisp, gdisp
   real(wp) :: dE, dG(3), dS(3, 3)

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)
   real(wp), allocatable :: dEdcn_local(:)
   real(wp), allocatable :: gradient_local(:, :)
   real(wp), allocatable :: sigma_local(:, :)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, dc6dcn, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, t8, &
   !$omp& d6, d8, edisp, gdisp, dE, dG, dS) &
   !$omp shared(energy, gradient, sigma, dEdcn) &
   !$omp private(energy_local, gradient_local, sigma_local, dEdcn_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdcn_local(size(dEdcn, 1)), source=0.0_wp)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            d6 = -6*r2**2*t6**2
            d8 = -8*r2**3*t8**2

            edisp = self%s6*t6 + self%s8*rrij*t8
            gdisp = self%s6*d6 + self%s8*rrij*d8

            dE = -c6ij*edisp * 0.5_wp
            dG(:) = -c6ij*gdisp*vec
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

            energy_local(iat) = energy_local(iat) + dE
            dEdcn_local(iat) = dEdcn_local(iat) - dc6dcn(iat, jat) * edisp
            sigma_local(:, :) = sigma_local + dS
            if (iat /= jat) then
               energy_local(jat) = energy_local(jat) + dE
               dEdcn_local(jat) = dEdcn_local(jat) - dc6dcn(jat, iat) * edisp
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local(:, :) = sigma_local + dS
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_dispersion_derivs_)
   energy(:) = energy(:) + energy_local(:)
   dEdcn(:) = dEdcn(:) + dEdcn_local(:)
   gradient(:, :) = gradient(:, :) + gradient_local(:, :)
   sigma(:, :) = sigma(:, :) + sigma_local(:, :)
   !$omp end critical (get_dispersion_derivs_)
   deallocate(energy_local)
   deallocate(dEdcn_local)
   deallocate(gradient_local)
   deallocate(sigma_local)
   !$omp end parallel

end subroutine get_dispersion_derivs


!> Evaluation of the dispersion energy expression
subroutine get_dispersion3(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   call get_atm_dispersion(mol, trans, cutoff, self%s9, rs9, self%alp+2, &
      & rvdw, c6, dc6dcn, energy, dEdcn, gradient, sigma)

end subroutine get_dispersion3


!> Evaluation of the dispersion energy expression projected on atomic pairs
subroutine get_pairwise_dispersion2(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, edisp, dE

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:, :)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, &
   !$omp& t8, edisp, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1), size(energy, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = self%s6*t6 + self%s8*rrij*t8

            dE = -c6ij*edisp * 0.5_wp

            energy_local(jat, iat) = energy_local(jat, iat) + dE
            if (iat /= jat) then
               energy_local(iat, jat) = energy_local(iat, jat) + dE
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_pairwise_dispersion2_)
   energy(:, :) = energy(:, :) + energy_local(:, :)
   !$omp end critical (get_pairwise_dispersion2_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_pairwise_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_pairwise_dispersion3(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   call get_atm_pairwise_dispersion(mol, trans, cutoff, self%s9, rs9, self%alp+2, &
      & rvdw, c6, energy)

end subroutine get_pairwise_dispersion3


end module dftd3_damping_rational
