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

module dftd3_damping
   use dftd3_cutoff, only : smooth_cutoff_r2
   use dftd3_fourier_ewald, only : get_dispersion_ewald
   use dftd3_fourier_kernel, only : fourier_term, max_fourier_terms, &
      & is_supported_term, get_reciprocal_cutoff
   use dftd3_model, only : d3_model
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   implicit none
   private

   public :: damping_param, dispersion_interface
   public :: get_dispersion2_hessian


   !> Target accuracy for the automatic reciprocal space cutoff
   real(wp), parameter :: kcut_tolerance = 1.0e-8_wp


   type, abstract :: damping_param
   contains
      generic :: get_dispersion2 => get_dispersion2_impl, get_dispersion2_compat
      procedure(dispersion_interface), deferred :: get_dispersion2_impl
      procedure :: get_dispersion2_compat => get_dispersion2_compat
      generic :: get_dispersion3 => get_dispersion3_impl, get_dispersion3_compat
      procedure(dispersion_interface), deferred :: get_dispersion3_impl
      procedure :: get_dispersion3_compat => get_dispersion3_compat
      generic :: get_pairwise_dispersion2 => get_pairwise_dispersion2_impl, get_pairwise_dispersion2_compat
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion2_impl
      procedure :: get_pairwise_dispersion2_compat => get_pairwise_dispersion2_compat
      generic :: get_pairwise_dispersion3 => get_pairwise_dispersion3_impl, get_pairwise_dispersion3_compat
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion3_impl
      procedure :: get_pairwise_dispersion3_compat => get_pairwise_dispersion3_compat
      !> Pair kernel and its derivatives w.r.t. squared distance and C6 coefficient
      procedure(damping_kernel_interface), deferred :: get_damping_kernel
      !> Evaluate three-body contribution to the hessian
      procedure(dispersion3_hessian_interface), deferred :: get_dispersion3_hessian
      !> Whether the damping function can be evaluated by Ewald summation
      procedure :: supports_ewald
      !> Reciprocal space representation of the damped pair potential
      procedure :: get_fourier_terms
      !> Evaluate the two-body dispersion energy by Ewald summation
      procedure :: get_dispersion2_ewald
   end type damping_param


   abstract interface
      !> Pair dispersion kernel and its partial derivatives.
      !>
      !> Returns the pair energy contribution as a function of the squared distance
      !> and the C6 coefficient together with all first and second partial
      !> derivatives. The smooth cutoff and multiplicity are applied by the caller.
      pure subroutine damping_kernel_interface(self, izp, jzp, rvdw, r4r2, r2, c6, &
            & e, er, ec, err, erc, ecc)
         import :: damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

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
      end subroutine damping_kernel_interface


      !> Evaluation of the three-body contribution to the hessian
      subroutine dispersion3_hessian_interface(self, mol, trans, cutoff, width, rvdw, &
            & r4r2, c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, dEdcndr, dEdcndcn)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

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
      end subroutine dispersion3_hessian_interface


      !> Evaluation of the dispersion energy expression
      subroutine dispersion_interface(self, mol, trans, cutoff, width, rvdw, r4r2, c6, dc6dcn, &
            & energy, dEdcn, gradient, sigma)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

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
      end subroutine dispersion_interface


      !> Evaluation of the pairwise representation of the dispersion energy
      subroutine pairwise_dispersion_interface(self, mol, trans, cutoff, width, rvdw, r4r2, c6, &
            & energy)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

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

         !> Pairwise representation of the dispersion energy
         real(wp), intent(inout) :: energy(:, :)
      end subroutine pairwise_dispersion_interface
   end interface

contains

   !> Whether the damping function has a known reciprocal space representation.
   !>
   !> Damping functions which are not rational in the interatomic distance, or
   !> whose damping radius depends on the C6 coefficient, cannot be brought into
   !> the separable form required by the Ewald summation.
   pure function supports_ewald(self) result(supported)

      !> Damping parameters
      class(damping_param), intent(in) :: self

      !> Whether the damping function can be evaluated in reciprocal space
      logical :: supported

      supported = .false.

   end function supports_ewald


   !> Reciprocal space representation of the damped pair potential.
   !>
   !> Only meaningful for damping functions reporting support for the Ewald
   !> summation, all others contribute no terms.
   pure subroutine get_fourier_terms(self, izp, jzp, rvdw, r4r2, terms, nterm)

      !> Damping parameters
      class(damping_param), intent(in) :: self

      !> Species indices of the pair
      integer, intent(in) :: izp, jzp

      !> Van-der-Waals radii for damping function
      real(wp), intent(in) :: rvdw(:, :)

      !> Expectation values for C8 extrapolation
      real(wp), intent(in) :: r4r2(:)

      !> Terms of the damped pair potential
      type(fourier_term), intent(out) :: terms(:)

      !> Number of terms of the damped pair potential
      integer, intent(out) :: nterm

      nterm = 0

   end subroutine get_fourier_terms


   !> Evaluate the two-body dispersion energy by summation over the reciprocal
   !> lattice.
   !>
   !> Requires a separable representation of the C6 coefficients in the dispersion
   !> model and a damping function supporting the Ewald summation.
   subroutine get_dispersion2_ewald(self, mol, disp, gwvec, gwdcn, energies, &
         & dEdcn, gradient, sigma, error)

      !> Damping parameters
      class(damping_param), intent(in) :: self

      !> Molecular structure data
      class(structure_type), intent(in) :: mol

      !> Dispersion model
      class(d3_model), intent(in) :: disp

      !> Weighting function for the atomic reference systems
      real(wp), intent(in) :: gwvec(:, :)

      !> Derivative of the weighting function w.r.t. the coordination number
      real(wp), intent(in), optional :: gwdcn(:, :)

      !> Dispersion energy
      real(wp), intent(inout) :: energies(:)

      !> Derivative of the energy w.r.t. the coordination number
      real(wp), intent(inout), optional :: dEdcn(:)

      !> Dispersion gradient
      real(wp), intent(inout), optional :: gradient(:, :)

      !> Dispersion virial
      real(wp), intent(inout), optional :: sigma(:, :)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: isp, jsp, it
      real(wp) :: kcut
      integer, allocatable :: nterm(:, :)
      type(fourier_term), allocatable :: terms(:, :, :)

      allocate(terms(max_fourier_terms, mol%nid, mol%nid))
      allocate(nterm(mol%nid, mol%nid))
      kcut = 0.0_wp
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            call self%get_fourier_terms(isp, jsp, disp%rvdw, disp%r4r2, &
               & terms(:, jsp, isp), nterm(jsp, isp))
            do it = 1, nterm(jsp, isp)
               if (.not.is_supported_term(terms(it, jsp, isp))) then
                  call fatal_error(error, "Damped pair potential has no known "//&
                     & "Fourier transform")
                  return
               end if
               kcut = max(kcut, get_reciprocal_cutoff(terms(it, jsp, isp), &
                  & kcut_tolerance))
            end do
         end do
      end do

      ! all scaling factors vanish, there is nothing to add
      if (kcut <= 0.0_wp) return

      if (disp%lowrank%kcut > 0.0_wp) kcut = disp%lowrank%kcut

      call get_dispersion_ewald(mol, disp%lowrank, disp%ghost, terms, nterm, kcut, &
         & gwvec, gwdcn, energies, dEdcn, gradient, sigma)

   end subroutine get_dispersion2_ewald


   !> Two-body contribution to the analytical Hessian.
   !>
   !> Accumulates the Cartesian second derivatives at fixed coordination number
   !> together with the derivatives w.r.t. the coordination number, which are
   !> contracted with the coordination number derivatives by the caller.
   subroutine get_dispersion2_hessian(self, mol, trans, cutoff, width, rvdw, r4r2, &
         & c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, dEdcndr, dEdcndcn)

      !> Damping parameters
      class(damping_param), intent(in) :: self

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

      integer :: iat, jat, izp, jzp, jtr, ic, jc, ii, jj, ia, ib, nc
      integer :: cnat(2)
      real(wp) :: vec(3), r2, cutoff2, c6ij, fac
      real(wp) :: e0, e0r, e0c, e0rr, e0rc, e0cc
      real(wp) :: sw, swr, swrr, fr, frr, fc, frc, fcc
      real(wp) :: block(3, 3), cnd(2), cnd2(2, 2), dr2(3, 2)

      ! Kept serial on purpose: the loop is O(N^2) but bound by the scattered
      ! writes into the O(N^2) hessian, so thread-private copies only add traffic
      cutoff2 = cutoff*cutoff

      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            c6ij = c6(jat, iat)
            if (iat /= jat) then
               fac = 1.0_wp
               nc = 2
               cnat(1) = iat
               cnat(2) = jat
               cnd(1) = dc6dcn(iat, jat)
               cnd(2) = dc6dcn(jat, iat)
               cnd2(1, 1) = d2c6dcn2(iat, jat)
               cnd2(2, 2) = d2c6dcn2(jat, iat)
               cnd2(1, 2) = d2c6dcnij(iat, jat)
               cnd2(2, 1) = d2c6dcnij(iat, jat)
            else
               fac = 0.5_wp
               nc = 1
               cnat(1) = iat
               cnd(1) = 2.0_wp*dc6dcn(iat, iat)
               cnd2(1, 1) = 2.0_wp*d2c6dcn2(iat, iat) + 2.0_wp*d2c6dcnij(iat, iat)
            end if
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

               call self%get_damping_kernel(izp, jzp, rvdw, r4r2, r2, c6ij, &
                  & e0, e0r, e0c, e0rr, e0rc, e0cc)
               call smooth_cutoff_r2(r2, cutoff, width, sw, swr, swrr)

               fr = fac*(swr*e0 + sw*e0r)
               frr = fac*(swrr*e0 + 2.0_wp*swr*e0r + sw*e0rr)
               fc = fac*sw*e0c
               frc = fac*(swr*e0c + sw*e0rc)
               fcc = fac*sw*e0cc

               ! Cartesian second derivatives at fixed coordination number
               if (iat /= jat) then
                  do ic = 1, 3
                     do jc = 1, 3
                        block(ic, jc) = 4.0_wp*frr*vec(ic)*vec(jc)
                     end do
                     block(ic, ic) = block(ic, ic) + 2.0_wp*fr
                  end do
                  do ic = 1, 3
                     do jc = 1, 3
                        hessian(3*(iat-1)+ic, 3*(iat-1)+jc) = &
                           & hessian(3*(iat-1)+ic, 3*(iat-1)+jc) + block(ic, jc)
                        hessian(3*(jat-1)+ic, 3*(jat-1)+jc) = &
                           & hessian(3*(jat-1)+ic, 3*(jat-1)+jc) + block(ic, jc)
                        hessian(3*(iat-1)+ic, 3*(jat-1)+jc) = &
                           & hessian(3*(iat-1)+ic, 3*(jat-1)+jc) - block(ic, jc)
                        hessian(3*(jat-1)+ic, 3*(iat-1)+jc) = &
                           & hessian(3*(jat-1)+ic, 3*(iat-1)+jc) - block(ic, jc)
                     end do
                  end do
                  dr2(:, 1) = 2.0_wp*vec
                  dr2(:, 2) = -2.0_wp*vec
               else
                  dr2(:, :) = 0.0_wp
               end if

               ! Coordination number dependence
               do ia = 1, nc
                  dEdcn(cnat(ia)) = dEdcn(cnat(ia)) + fc*cnd(ia)
                  do ib = 1, nc
                     dEdcndcn(cnat(ia), cnat(ib)) = dEdcndcn(cnat(ia), cnat(ib)) &
                        & + fcc*cnd(ia)*cnd(ib) + fc*cnd2(ia, ib)
                  end do
               end do

               if (iat /= jat) then
                  do ia = 1, nc
                     do ic = 1, 3
                        ii = 3*(iat-1) + ic
                        jj = 3*(jat-1) + ic
                        dEdcndr(ii, cnat(ia)) = dEdcndr(ii, cnat(ia)) &
                           & + frc*dr2(ic, 1)*cnd(ia)
                        dEdcndr(jj, cnat(ia)) = dEdcndr(jj, cnat(ia)) &
                           & + frc*dr2(ic, 2)*cnd(ia)
                     end do
                  end do
               end if
            end do
         end do
      end do

   end subroutine get_dispersion2_hessian


   !> Evaluation of the dispersion energy expression
   subroutine get_dispersion2_compat(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
         & energy, dEdcn, gradient, sigma)

      !> Damping parameters
      class(damping_param), intent(in) :: self

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

      call self%get_dispersion2(mol, trans, cutoff, 0.0_wp, rvdw, r4r2, c6, &
         & dc6dcn, energy, dEdcn, gradient, sigma)
   end subroutine get_dispersion2_compat


   !> Evaluation of the dispersion energy expression
   subroutine get_dispersion3_compat(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
         & energy, dEdcn, gradient, sigma)

      !> Damping parameters
      class(damping_param), intent(in) :: self

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

      call self%get_dispersion3(mol, trans, cutoff, 0.0_wp, rvdw, r4r2, c6, &
         & dc6dcn, energy, dEdcn, gradient, sigma)
   end subroutine get_dispersion3_compat


   !> Evaluation of the pairwise representation of the dispersion energy
   subroutine get_pairwise_dispersion2_compat(self, mol, trans, cutoff, rvdw, r4r2, c6, &
         & energy)

      !> Damping parameters
      class(damping_param), intent(in) :: self

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

      !> Pairwise representation of the dispersion energy
      real(wp), intent(inout) :: energy(:, :)

      call self%get_pairwise_dispersion2(mol, trans, cutoff, 0.0_wp, rvdw, r4r2, c6, energy)
   end subroutine get_pairwise_dispersion2_compat


   !> Evaluation of the pairwise representation of the dispersion energy
   subroutine get_pairwise_dispersion3_compat(self, mol, trans, cutoff, rvdw, r4r2, c6, &
         & energy)

      !> Damping parameters
      class(damping_param), intent(in) :: self

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

      !> Pairwise representation of the dispersion energy
      real(wp), intent(inout) :: energy(:, :)

      call self%get_pairwise_dispersion3(mol, trans, cutoff, 0.0_wp, rvdw, r4r2, c6, energy)
   end subroutine get_pairwise_dispersion3_compat

end module dftd3_damping
