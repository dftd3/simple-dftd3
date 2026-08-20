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

!> Separable low-rank representation of the environment dependent C6 coefficients.
!>
!> The reference C6 coefficients are indexed by a species and a reference system
!> and form a symmetric matrix, which is eigendecomposed once for a given system as
!>
!>    C6ref(Zi, p, Zj, q) = sum_l lambda_l v_l(Zi, p) v_l(Zj, q)
!>
!> Contracting the eigenvectors with the Gaussian weights of the reference systems
!> yields atom centered factors
!>
!>    C6(l, i) = sum_p v_l(Zi, p) w(i, p)
!>
!> such that the pair coefficients become separable
!>
!>    C6(i, j) = sum_l lambda_l C6(l, i) C6(l, j)
!>
!> This separability is what allows a reciprocal space evaluation of the dispersion
!> energy, since the atomic factors can be accumulated into structure factors.
module dftd3_fourier_decomposition
   use dftd3_fourier_jacobi, only : symmetric_eigendecomposition
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: d3_lowrank_c6, new_lowrank_c6, d3_lowrank_config


   !> Setup of the separable representation of the reference C6 coefficients
   type :: d3_lowrank_config

      !> Requested rank of the expansion, zero selects the rank from the tolerance
      integer :: rank = 0

      !> Requested maximum relative error of the reference C6 coefficients
      real(wp) :: tolerance = 1.0e-4_wp

      !> Reciprocal space cutoff, zero selects the cutoff from the damping radii
      real(wp) :: kcut = 0.0_wp

      !> Number of mesh points for the particle mesh evaluation, zero derives the
      !> mesh from the reciprocal space cutoff and a negative value sums over the
      !> reciprocal lattice directly
      integer :: mesh = 0

   end type d3_lowrank_config


   !> Separable representation of the reference C6 coefficients
   type :: d3_lowrank_c6

      !> Number of terms in the separable expansion
      integer :: rank = 0

      !> Maximum relative error of the reconstructed reference C6 coefficients
      real(wp) :: error = 0.0_wp

      !> Reciprocal space cutoff, zero selects the cutoff from the damping radii
      real(wp) :: kcut = 0.0_wp

      !> Number of mesh points for the particle mesh evaluation
      integer :: mesh = 0

      !> Eigenvalues of the reference C6 matrix
      real(wp), allocatable :: lambda(:)

      !> Eigenvectors resolved by reference system and species
      real(wp), allocatable :: vec(:, :, :)

   contains

      !> Contract the eigenvectors with the reference system weights
      procedure :: get_weights

      !> Assemble the pair C6 coefficients from the separable factors
      procedure :: get_atomic_c6

   end type d3_lowrank_c6


contains


!> Decompose the reference C6 coefficients into a separable expansion.
!>
!> Either a fixed rank or a target accuracy can be requested, the latter selects
!> the smallest rank reproducing all reference coefficients within the tolerance.
subroutine new_lowrank_c6(self, ref, c6ref, config)

   !> Instance of the separable representation
   type(d3_lowrank_c6), intent(out) :: self

   !> Number of reference systems for each species
   integer, intent(in) :: ref(:)

   !> Reference C6 coefficients
   real(wp), intent(in) :: c6ref(:, :, :, :)

   !> Setup of the separable representation
   type(d3_lowrank_config), intent(in) :: config

   integer :: nid, mref, ndim, isp, iref, ii, jj, il, nrank
   integer, allocatable :: spmap(:), refmap(:)
   real(wp), allocatable :: amat(:, :), eval(:), evec(:, :), res(:, :)

   nid = size(ref)
   mref = maxval(ref)
   ndim = sum(ref)

   allocate(spmap(ndim), refmap(ndim))
   ii = 0
   do isp = 1, nid
      do iref = 1, ref(isp)
         ii = ii + 1
         spmap(ii) = isp
         refmap(ii) = iref
      end do
   end do

   allocate(amat(ndim, ndim))
   do ii = 1, ndim
      do jj = 1, ndim
         amat(jj, ii) = c6ref(refmap(jj), refmap(ii), spmap(jj), spmap(ii))
      end do
   end do

   allocate(eval(ndim), evec(ndim, ndim))
   call symmetric_eigendecomposition(amat, eval, evec)

   allocate(res(ndim, ndim), source=amat)
   nrank = ndim
   do il = 1, ndim
      do ii = 1, ndim
         res(:, ii) = res(:, ii) - eval(il) * evec(:, il) * evec(ii, il)
      end do
      self%error = max_relative_error(res, amat)
      if (config%rank > 0) then
         if (il >= min(config%rank, ndim)) then
            nrank = il
            exit
         end if
      else if (self%error <= config%tolerance) then
         nrank = il
         exit
      end if
   end do

   self%rank = nrank
   self%kcut = config%kcut
   self%mesh = config%mesh
   self%lambda = eval(:nrank)
   allocate(self%vec(mref, nid, nrank), source=0.0_wp)
   do il = 1, nrank
      do ii = 1, ndim
         self%vec(refmap(ii), spmap(ii), il) = evec(ii, il)
      end do
   end do

end subroutine new_lowrank_c6


!> Maximum relative deviation of a residual matrix from the original entries
pure function max_relative_error(res, amat) result(err)

   !> Residual of the truncated expansion
   real(wp), intent(in) :: res(:, :)

   !> Original matrix
   real(wp), intent(in) :: amat(:, :)

   !> Maximum relative error
   real(wp) :: err

   real(wp) :: floor

   ! Entries far below the largest coefficient are irrelevant for the pair C6
   floor = 1.0e-6_wp * maxval(abs(amat))

   err = maxval(abs(res) / max(abs(amat), floor))

end function max_relative_error


!> Contract the eigenvectors with the Gaussian weights of the reference systems
!> to obtain the atom centered factors of the separable expansion.
subroutine get_weights(self, mol, ghost, gwvec, c6l, gwdcn, dc6ldcn)

   !> Instance of the separable representation
   class(d3_lowrank_c6), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Atoms excluded from the dispersion calculation
   logical, intent(in) :: ghost(:)

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Atom centered factors of the separable expansion
   real(wp), intent(out) :: c6l(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> Derivative of the atomic factors w.r.t. the coordination number
   real(wp), intent(out), optional :: dc6ldcn(:, :)

   integer :: iat, izp, il
   logical :: grad

   grad = present(gwdcn) .and. present(dc6ldcn)

   c6l(:, :) = 0.0_wp
   if (grad) dc6ldcn(:, :) = 0.0_wp

   do iat = 1, mol%nat
      if (ghost(iat)) cycle
      izp = mol%id(iat)
      do il = 1, self%rank
         c6l(il, iat) = sum(self%vec(:, izp, il) * gwvec(:, iat))
         if (grad) dc6ldcn(il, iat) = sum(self%vec(:, izp, il) * gwdcn(:, iat))
      end do
   end do

end subroutine get_weights


!> Assemble the pair C6 coefficients and their derivatives w.r.t. the
!> coordination number from the separable factors.
subroutine get_atomic_c6(self, mol, ghost, gwvec, gwdcn, c6, dc6dcn)

   !> Instance of the separable representation
   class(d3_lowrank_c6), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Atoms excluded from the dispersion calculation
   logical, intent(in) :: ghost(:)

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> C6 coefficients for all atom pairs
   real(wp), intent(out) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out), optional :: dc6dcn(:, :)

   integer :: iat, jat
   logical :: grad
   real(wp), allocatable :: c6l(:, :), dc6ldcn(:, :)

   grad = present(gwdcn) .and. present(dc6dcn)

   allocate(c6l(self%rank, mol%nat))
   if (grad) allocate(dc6ldcn(self%rank, mol%nat))
   call self%get_weights(mol, ghost, gwvec, c6l, gwdcn, dc6ldcn)

   c6(:, :) = 0.0_wp
   if (grad) dc6dcn(:, :) = 0.0_wp

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(c6, dc6dcn, mol, self, c6l, dc6ldcn, grad) private(iat, jat)
   do iat = 1, mol%nat
      do jat = 1, iat
         c6(iat, jat) = sum(self%lambda * c6l(:, iat) * c6l(:, jat))
         c6(jat, iat) = c6(iat, jat)
         if (grad) then
            dc6dcn(iat, jat) = sum(self%lambda * dc6ldcn(:, iat) * c6l(:, jat))
            dc6dcn(jat, iat) = sum(self%lambda * c6l(:, iat) * dc6ldcn(:, jat))
         end if
      end do
   end do

end subroutine get_atomic_c6


end module dftd3_fourier_decomposition
