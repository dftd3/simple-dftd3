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

module dftd3_damping_mzero
   use dftd3_cutoff, only : smooth_cutoff
   use dftd3_damping, only : damping_param
   use dftd3_damping_atm, only : get_atm_dispersion, get_atm_pairwise_dispersion, &
      & get_atm_dispersion_hessian
   use dftd3_param, only : d3_param
   use dftd3_partition, only : work_partition, owns_pair
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: mzero_damping_param, new_mzero_damping


   !> Modified zero damping model
   type, extends(damping_param) :: mzero_damping_param
      real(wp) :: s6
      real(wp) :: s8
      real(wp) :: rs6
      real(wp) :: rs8
      real(wp) :: alp
      real(wp) :: bet
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2_impl => get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3_impl => get_dispersion3

      !> Evaluate pairwise representation of additive dispersion energy
      procedure :: get_pairwise_dispersion2_impl => get_pairwise_dispersion2

      !> Evaluate pairwise representation of non-additive dispersion energy
      procedure :: get_pairwise_dispersion3_impl => get_pairwise_dispersion3

      !> Evaluate pair kernel and its derivatives
      procedure :: get_damping_kernel

      !> Evaluate ATM three-body contribution to the hessian
      procedure :: get_dispersion3_hessian

   end type mzero_damping_param


   real(wp), parameter :: rs9 = 4.0_wp/3.0_wp


contains


!> Create new modified zero damping model
subroutine new_mzero_damping(self, param)

   !> Modified zero damping parameters
   type(mzero_damping_param), intent(out) :: self

   !> Parameters
   type(d3_param), intent(in) :: param

   self%s6 = param%s6
   self%s8 = param%s8
   self%s9 = param%s9
   self%rs6 = param%rs6
   self%rs8 = param%rs8
   self%alp = param%alp
   self%bet = param%bet

end subroutine new_mzero_damping


!> Pair dispersion kernel and its derivatives w.r.t. squared distance and C6
pure subroutine get_damping_kernel(self, izp, jzp, rvdw, r4r2, r2, c6, &
      & e, er, ec, err, erc, ecc)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Species indices of the pair
   integer, intent(in) :: izp, jzp

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> Squared distance of the pair
   real(wp), intent(in) :: r2

   !> C6 coefficient of the pair
   real(wp), intent(in) :: c6

   !> Pair energy and its derivatives
   real(wp), intent(out) :: e, er, ec, err, erc, ecc

   real(wp) :: rrij, r0ij, u, r1, alp6, alp8, betr0
   real(wp) :: w6, w8, w6d, w8d, w6dd, w8dd
   real(wp) :: t6, t8, t6d, t8d, t6dd, t8dd
   real(wp) :: f6, f8, f6d, f8d, f6dd, f8dd
   real(wp) :: c6t, c6td, c6tdd, c8t, c8td, c8tdd, p, pd, pdd

   e = 0.0_wp; er = 0.0_wp; ec = 0.0_wp
   err = 0.0_wp; erc = 0.0_wp; ecc = 0.0_wp
   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return

   u = r2
   r1 = sqrt(u)
   rrij = 3*r4r2(izp)*r4r2(jzp)
   r0ij = rvdw(jzp, izp)
   alp6 = self%alp
   alp8 = self%alp + 2.0_wp
   betr0 = self%bet*r0ij

   w6 = r1/(self%rs6*r0ij) + betr0
   w8 = r1/(self%rs8*r0ij) + betr0
   w6d = 1.0_wp/(self%rs6*r0ij*2.0_wp*r1)
   w8d = 1.0_wp/(self%rs8*r0ij*2.0_wp*r1)
   w6dd = -1.0_wp/(self%rs6*r0ij*4.0_wp*u*r1)
   w8dd = -1.0_wp/(self%rs8*r0ij*4.0_wp*u*r1)

   t6 = w6**(-alp6)
   t8 = w8**(-alp8)
   t6d = -alp6*w6**(-alp6 - 1.0_wp)*w6d
   t8d = -alp8*w8**(-alp8 - 1.0_wp)*w8d
   t6dd = alp6*(alp6 + 1.0_wp)*w6**(-alp6 - 2.0_wp)*w6d*w6d &
      & - alp6*w6**(-alp6 - 1.0_wp)*w6dd
   t8dd = alp8*(alp8 + 1.0_wp)*w8**(-alp8 - 2.0_wp)*w8d*w8d &
      & - alp8*w8**(-alp8 - 1.0_wp)*w8dd

   f6 = 1.0_wp/(1.0_wp + 6.0_wp*t6)
   f8 = 1.0_wp/(1.0_wp + 6.0_wp*t8)
   f6d = -6.0_wp*t6d*f6*f6
   f8d = -6.0_wp*t8d*f8*f8
   f6dd = -6.0_wp*t6dd*f6*f6 + 72.0_wp*t6d*t6d*f6**3
   f8dd = -6.0_wp*t8dd*f8*f8 + 72.0_wp*t8d*t8d*f8**3

   c6t = f6/u**3
   c6td = f6d/u**3 - 3.0_wp*f6/u**4
   c6tdd = f6dd/u**3 - 6.0_wp*f6d/u**4 + 12.0_wp*f6/u**5
   c8t = f8/u**4
   c8td = f8d/u**4 - 4.0_wp*f8/u**5
   c8tdd = f8dd/u**4 - 8.0_wp*f8d/u**5 + 20.0_wp*f8/u**6

   p = self%s6*c6t + self%s8*rrij*c8t
   pd = self%s6*c6td + self%s8*rrij*c8td
   pdd = self%s6*c6tdd + self%s8*rrij*c8tdd

   e = -c6*p
   er = -c6*pd
   ec = -p
   err = -c6*pdd
   erc = -pd

end subroutine get_damping_kernel


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, width, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma, partition)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

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

   !> Work partition of the atom pairs, defaults to the complete work
   type(work_partition), intent(in), optional :: partition

   logical :: grad

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)

   if (grad) then
      call get_dispersion_derivs(self, mol, trans, cutoff, width, rvdw, r4r2, c6, dc6dcn, &
         & energy, dEdcn, gradient, sigma, partition)
   else
      call get_dispersion_energy(self, mol, trans, cutoff, width, rvdw, r4r2, c6, energy, &
         & partition)
   end if

end subroutine get_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_dispersion_energy(self, mol, trans, cutoff, width, rvdw, r4r2, c6, energy, &
      & partition)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Work partition of the atom pairs, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, r1, r6, r8, t6, t8, f6, f8, alp6, alp8
   real(wp) :: edisp, cutoff2, r0ij, rrij, c6ij, dE
   real(wp) :: irs6r0, irs8r0, betr0, sw, dswdr

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff
   alp6 = self%alp
   alp8 = self%alp + 2.0_wp

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, trans, cutoff2, cutoff, width, alp6, alp8, rvdw, r4r2, &
   !$omp& partition) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r1, r6, r8, t6, t8, f6, &
   !$omp& f8, edisp, r0ij, rrij, c6ij, dE, irs6r0, irs8r0, betr0, sw, dswdr) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = rvdw(jzp, izp)
         c6ij = c6(jat, iat)
         irs6r0 = 1.0_wp/(self%rs6*r0ij)
         irs8r0 = 1.0_wp/(self%rs8*r0ij)
         betr0 = self%bet*r0ij
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            r1 = sqrt(r2)
            call smooth_cutoff(r1, cutoff, width, sw, dswdr)

            r6 = r2*r2*r2
            r8 = r6*r2

            t6 = (r1*irs6r0+betr0)**(-alp6)
            t8 = (r1*irs8r0+betr0)**(-alp8)

            f6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
            f8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

            edisp = sw * (self%s6 * f6 / r6 + self%s8 * rrij * f8 / r8)

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
subroutine get_dispersion_derivs(self, mol, trans, cutoff, width, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma, partition)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

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

   !> Work partition of the atom pairs, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, jtr, ic, jc
   real(wp) :: vec(3), r2, r1, r6, r8, t6, t8, d6, d8, f6, f8, alp6, alp8
   real(wp) :: edisp0, gdisp0, edisp, gdisp, cutoff2, r0ij, rrij, c6ij, dE, dG(3), dS(3, 3)
   real(wp) :: irs6r0, irs8r0, betr0, betr02rs6, betr02rs8, sw, dswdr

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)
   real(wp), allocatable :: dEdcn_local(:)
   real(wp), allocatable :: gradient_local(:, :)
   real(wp), allocatable :: sigma_local(:, :)

   cutoff2 = cutoff*cutoff
   alp6 = self%alp
   alp8 = self%alp + 2.0_wp

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, dc6dcn, trans, cutoff2, cutoff, width, alp6, alp8, rvdw, &
   !$omp& r4r2, partition) &
   !$omp private(iat, jat, izp, jzp, jtr, ic, jc, vec, r2, r1, r6, r8, t6, t8, d6, &
   !$omp& d8, f6, f8, edisp0, gdisp0, edisp, gdisp, r0ij, rrij, c6ij, &
   !$omp& dE, dG, dS, irs6r0, irs8r0, betr0, betr02rs6, betr02rs8, sw, dswdr) &
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
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = rvdw(jzp, izp)
         c6ij = c6(jat, iat)
         irs6r0 = 1.0_wp/(self%rs6*r0ij)
         irs8r0 = 1.0_wp/(self%rs8*r0ij)
         betr0 = self%bet*r0ij
         betr02rs6 = self%rs6*self%bet*r0ij*r0ij
         betr02rs8 = self%rs8*self%bet*r0ij*r0ij
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            r1 = sqrt(r2)
            call smooth_cutoff(r1, cutoff, width, sw, dswdr)

            r6 = r2*r2*r2
            r8 = r6*r2

            t6 = (r1*irs6r0+betr0)**(-alp6)
            t8 = (r1*irs8r0+betr0)**(-alp8)

            f6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
            f8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

            d6 = -6.0_wp * f6 / r2 &
               & + 6.0_wp*alp6*t6*f6**2 / (r2+betr02rs6*r1)
            d8 = -8.0_wp * f8 / r2 &
               & + 6.0_wp*alp8*t8*f8**2 / (r2+betr02rs8*r1)

            edisp0 = self%s6 * f6 / r6 + self%s8 * rrij * f8 / r8
            gdisp0 = self%s6 * d6 / r6 + self%s8 * rrij * d8 / r8
            edisp = sw * edisp0
            gdisp = sw * gdisp0 + dswdr * edisp0 / r1

            dE = -c6ij*edisp * 0.5_wp
            dG(:) = -c6ij*gdisp*vec
            do ic = 1, 3
               do jc = 1, 3
                  dS(ic, jc) = dG(ic) * vec(jc) * 0.5_wp
               end do
            end do

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
subroutine get_dispersion3(self, mol, trans, cutoff, width, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma, partition)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

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

   !> Work partition of the atom pairs, defaults to the complete work
   type(work_partition), intent(in), optional :: partition

   call get_atm_dispersion(mol, trans, cutoff, width, self%s9, rs9, self%alp+2, &
      & rvdw, c6, dc6dcn, energy, dEdcn, gradient, sigma, partition)

end subroutine get_dispersion3


!> Evaluation of the three-body contribution to the hessian
subroutine get_dispersion3_hessian(self, mol, trans, cutoff, width, rvdw, &
      & r4r2, c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, dEdcndr, dEdcndcn, partition)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivatives of the C6 w.r.t. the coordination number
   real(wp), intent(in) :: dc6dcn(:, :), d2c6dcn2(:, :), d2c6dcnij(:, :)

   !> Second derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: hessian(:, :)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)

   !> Mixed derivative w.r.t. coordination number and Cartesian coordinates
   real(wp), intent(inout) :: dEdcndr(:, :)

   !> Second derivative w.r.t. the coordination numbers
   real(wp), intent(inout) :: dEdcndcn(:, :)

   !> Work partition of the atom pairs, defaults to the complete work
   type(work_partition), intent(in), optional :: partition

   call get_atm_dispersion_hessian(mol, trans, cutoff, width, self%s9, rs9, self%alp+2, &
      & rvdw, c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, dEdcndr, dEdcndcn, &
      & partition)

end subroutine get_dispersion3_hessian


!> Evaluation of the dispersion energy expression projected on atomic pairs
subroutine get_pairwise_dispersion2(self, mol, trans, cutoff, width, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, r1, r6, r8, t6, t8, f6, f8, alp6, alp8
   real(wp) :: edisp, cutoff2, r0ij, rrij, c6ij, dE, sw, dswdr

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:, :)

   cutoff2 = cutoff*cutoff
   alp6 = self%alp
   alp8 = self%alp + 2.0_wp

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, trans, cutoff2, cutoff, width, alp6, alp8, rvdw, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r1, r6, r8, t6, t8, f6, &
   !$omp& f8, edisp, r0ij, rrij, c6ij, dE, sw, dswdr) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1), size(energy, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = rvdw(jzp, izp)
         c6ij = c6(jat, iat)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            r1 = sqrt(r2)
            call smooth_cutoff(r1, cutoff, width, sw, dswdr)

            r6 = r2*r2*r2
            r8 = r6*r2

            t6 = (r1/(self%rs6*r0ij)+self%bet*r0ij)**(-alp6)
            t8 = (r1/(self%rs8*r0ij)+self%bet*r0ij)**(-alp8)

            f6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
            f8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

            edisp = sw * (self%s6 * f6 / r6 + self%s8 * rrij * f8 / r8)

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
subroutine get_pairwise_dispersion3(self, mol, trans, cutoff, width, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(mzero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   call get_atm_pairwise_dispersion(mol, trans, cutoff, width, self%s9, rs9, self%alp+2, &
      & rvdw, c6, energy)

end subroutine get_pairwise_dispersion3


end module dftd3_damping_mzero
