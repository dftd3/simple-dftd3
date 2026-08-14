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

module dftd3_damping_cso
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

   public :: cso_damping_param, new_cso_damping


   !> CSO (C6-Scaled Only) damping model
   type, extends(damping_param) :: cso_damping_param
      real(wp) :: s6
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: a3
      real(wp) :: a4
      real(wp) :: s9
      real(wp) :: alp
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

   end type cso_damping_param


   real(wp), parameter :: rs9 = 4.0_wp/3.0_wp


contains


!> Create new CSO damping model
subroutine new_cso_damping(self, param)

   !> CSO damping parameters
   type(cso_damping_param), intent(out) :: self

   !> Parameters
   type(d3_param), intent(in) :: param

   self%s6 = param%s6
   self%a1 = param%a1
   self%a2 = param%a2
   self%a3 = param%rs6
   self%a4 = param%rs8
   self%s9 = param%s9
   self%alp = param%alp

end subroutine new_cso_damping


!> Pair dispersion kernel and its derivatives w.r.t. squared distance and C6
pure subroutine get_damping_kernel(self, izp, jzp, rvdw, r4r2, r2, c6, &
      & e, er, ec, err, erc, ecc)

   !> Damping parameters
   class(cso_damping_param), intent(in) :: self

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

   real(wp) :: rrij, r0ij, u, r1, dd6, dd6_6, a2r0ij
   real(wp) :: ef, sf, sfd, sfdd, sig, sigd, sigdd, sigu, siguu
   real(wp) :: t6, t6d, t6dd, p, pd, pdd

   e = 0.0_wp; er = 0.0_wp; ec = 0.0_wp
   err = 0.0_wp; erc = 0.0_wp; ecc = 0.0_wp

   u = r2
   r1 = sqrt(u)
   rrij = 3*r4r2(izp)*r4r2(jzp)
   r0ij = sqrt(rrij)
   dd6 = self%a3*r0ij + self%a4
   dd6_6 = dd6**6
   a2r0ij = self%a2*r0ij

   ef = exp(r1 - a2r0ij)
   sf = 1.0_wp/(1.0_wp + ef)
   ! derivatives of the sigmoid w.r.t. the distance
   sfd = -sf*(1.0_wp - sf)
   sfdd = sf*(1.0_wp - sf)*(1.0_wp - 2.0_wp*sf)
   sig = self%s6 + self%a1*sf
   sigd = self%a1*sfd
   sigdd = self%a1*sfdd
   ! convert to derivatives w.r.t. the squared distance
   sigu = 0.5_wp*sigd/r1
   siguu = 0.25_wp*sigdd/u - 0.25_wp*sigd/(u*r1)

   t6 = 1.0_wp/(u*u*u + dd6_6)
   t6d = -3.0_wp*u*u*t6*t6
   t6dd = -6.0_wp*u*t6*t6 + 18.0_wp*u**4*t6**3

   p = sig*t6
   pd = sigu*t6 + sig*t6d
   pdd = siguu*t6 + 2.0_wp*sigu*t6d + sig*t6dd

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
   class(cso_damping_param), intent(in) :: self

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

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%a1) < epsilon(1.0_wp)) return
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
   class(cso_damping_param), intent(in) :: self

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
   real(wp) :: vec(3), r2, cutoff2, rrij, r0ij, rij, d6, ef, sf, sig, t6
   real(wp) :: c6ij, edisp, dE, sw, dswdr
   real(wp) :: d6_6, a2r0ij, r4

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, trans, cutoff2, cutoff, width, r4r2, partition) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, rrij, r0ij, rij, d6, ef, &
   !$omp& sf, sig, t6, c6ij, edisp, dE, d6_6, a2r0ij, r4, sw, dswdr) &
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
         r0ij = sqrt(rrij)
         c6ij = c6(jat, iat)
         d6 = self%a3 * r0ij + self%a4
         d6_6 = d6**6
         a2r0ij = self%a2 * r0ij
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            rij = sqrt(r2)
            call smooth_cutoff(rij, cutoff, width, sw, dswdr)
            ef = exp(rij - a2r0ij)
            sf = 1.0_wp / (1.0_wp + ef)
            sig = self%s6 + self%a1 * sf
            r4 = r2 * r2
            t6 = 1.0_wp / (r4 * r2 + d6_6)

            edisp = sw * sig * t6

            dE = -c6ij*edisp * 0.5_wp

            energy_local(iat) = energy_local(iat) + dE
            if (iat /= jat) then
               energy_local(jat) = energy_local(jat) + dE
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_dispersion_energy_cso_)
   energy(:) = energy(:) + energy_local(:)
   !$omp end critical (get_dispersion_energy_cso_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_dispersion_energy


!> Evaluation of the dispersion energy expression
subroutine get_dispersion_derivs(self, mol, trans, cutoff, width, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma, partition)

   !> Damping parameters
   class(cso_damping_param), intent(in) :: self

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
   real(wp) :: vec(3), r2, cutoff2, rrij, r0ij, rij, d6, ef, sf, sig, t6
   real(wp) :: c6ij, dsig, dt6, edisp0, gdisp0, edisp, gdisp, sw, dswdr
   real(wp) :: dE, dG(3), dS(3, 3)
   real(wp) :: d6_6, a2r0ij, r4

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)
   real(wp), allocatable :: dEdcn_local(:)
   real(wp), allocatable :: gradient_local(:, :)
   real(wp), allocatable :: sigma_local(:, :)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, dc6dcn, trans, cutoff2, cutoff, width, r4r2, partition) &
   !$omp private(iat, jat, izp, jzp, jtr, ic, jc, vec, r2, rrij, r0ij, rij, d6, ef, &
   !$omp& sf, sig, t6, c6ij, dsig, dt6, edisp0, gdisp0, edisp, gdisp, dE, &
   !$omp& dG, dS, d6_6, a2r0ij, r4, sw, dswdr) &
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
         r0ij = sqrt(rrij)
         c6ij = c6(jat, iat)
         d6 = self%a3 * r0ij + self%a4
         d6_6 = d6**6
         a2r0ij = self%a2 * r0ij
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            rij = sqrt(r2)
            call smooth_cutoff(rij, cutoff, width, sw, dswdr)
            ef = exp(rij - a2r0ij)
            sf = 1.0_wp / (1.0_wp + ef)
            sig = self%s6 + self%a1 * sf
            r4 = r2 * r2
            t6 = 1.0_wp / (r4 * r2 + d6_6)

            dsig = -self%a1 * sf * (1.0_wp - sf) / rij
            dt6 = -6*r4*t6*t6

            edisp0 = sig * t6
            gdisp0 = dsig * t6 + sig * dt6
            edisp = sw * edisp0
            gdisp = sw * gdisp0 + dswdr * edisp0 / rij

            dE = -c6ij*edisp * 0.5_wp
            dG(:) = -c6ij*gdisp*vec
            do ic = 1, 3
               do jc = 1, 3
                  dS(ic, jc) = dG(ic) * vec(jc) * 0.5_wp
               end do
            end do

            energy_local(iat) = energy_local(iat) + dE
            dEdcn_local(iat) = dEdcn_local(iat) - dc6dcn(iat, jat) * edisp
            sigma_local(:, :) = sigma_local(:, :) + dS
            if (iat /= jat) then
               energy_local(jat) = energy_local(jat) + dE
               dEdcn_local(jat) = dEdcn_local(jat) - dc6dcn(jat, iat) * edisp
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local(:, :) = sigma_local(:, :) + dS
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_dispersion_derivs_cso_)
   energy(:) = energy(:) + energy_local(:)
   dEdcn(:) = dEdcn(:) + dEdcn_local(:)
   gradient(:, :) = gradient(:, :) + gradient_local(:, :)
   sigma(:, :) = sigma(:, :) + sigma_local(:, :)
   !$omp end critical (get_dispersion_derivs_cso_)
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
   class(cso_damping_param), intent(in) :: self

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
   class(cso_damping_param), intent(in) :: self

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
   class(cso_damping_param), intent(in) :: self

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
   real(wp) :: vec(3), r2, cutoff2, rrij, r0ij, rij, d6, ef, sf, sig, t6
   real(wp) :: c6ij, edisp, dE, sw, dswdr
   real(wp) :: d6_6, a2r0ij, r4

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:, :)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, c6, trans, cutoff2, cutoff, width, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, rrij, r0ij, rij, d6, ef, &
   !$omp& sf, sig, t6, c6ij, edisp, dE, d6_6, a2r0ij, r4, sw, dswdr) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1), size(energy, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = sqrt(rrij)
         c6ij = c6(jat, iat)
         d6 = self%a3 * r0ij + self%a4
         d6_6 = d6**6
         a2r0ij = self%a2 * r0ij
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            rij = sqrt(r2)
            call smooth_cutoff(rij, cutoff, width, sw, dswdr)
            ef = exp(rij - a2r0ij)
            sf = 1.0_wp / (1.0_wp + ef)
            sig = self%s6 + self%a1 * sf
            r4 = r2 * r2
            t6 = 1.0_wp / (r4 * r2 + d6_6)

            edisp = sw * sig * t6

            dE = -c6ij*edisp * 0.5_wp

            energy_local(jat, iat) = energy_local(jat, iat) + dE
            if (iat /= jat) then
               energy_local(iat, jat) = energy_local(iat, jat) + dE
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_pairwise_dispersion2_cso_)
   energy(:, :) = energy(:, :) + energy_local(:, :)
   !$omp end critical (get_pairwise_dispersion2_cso_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_pairwise_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_pairwise_dispersion3(self, mol, trans, cutoff, width, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(cso_damping_param), intent(in) :: self

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


end module dftd3_damping_cso
