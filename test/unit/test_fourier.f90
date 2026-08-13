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

module test_fourier
   use dftd3, only : get_dispersion, get_coordination_number, &
      & get_lattice_points, realspace_cutoff, get_realspace_cutoff, d3_param, &
      & rational_damping_param, new_rational_damping, zero_damping_param, &
      & new_zero_damping, mzero_damping_param, new_mzero_damping, d3_model, &
      & new_d3_model, d3_lowrank_config, get_pairwise_dispersion
   use dftd3_citation, only : citation_type, is_citation_present, doi_fourier_d3
   use dftd3_fourier_jacobi, only : symmetric_eigendecomposition
   use dftd3_fourier_kernel, only : fourier_term, get_fourier_transform, &
      & get_potential_zero
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   implicit none
   private

   public :: collect_fourier

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30.0_wp, disp2=60.0_wp, disp3=15.0_wp)

   !> Full rank is requested by asking for more terms than the expansion can hold
   type(d3_lowrank_config), parameter :: exact = d3_lowrank_config(rank=10000)

   !> Reciprocal cutoff large enough to make the truncation error negligible
   type(d3_lowrank_config), parameter :: tight = &
      & d3_lowrank_config(rank=10000, kcut=10.0_wp)

   type(d3_param), parameter :: pbe_bj = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   type(d3_param), parameter :: pbe_zero = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & rs6 = 1.217_wp, s8 = 0.722_wp, rs8 = 1.0_wp)

   !> Modified zero damping has no closed form transform and stays in real space
   type(d3_param), parameter :: pbe_mzero = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & rs6 = 1.0_wp, s8 = 1.0_wp, rs8 = 1.0_wp, bet = 0.1_wp)

   type(realspace_cutoff), parameter :: long = &
      & realspace_cutoff(cn=30.0_wp, disp2=400.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_fourier(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("jacobi-eigendecomposition", test_jacobi), &
      & new_unittest("jacobi-diagonal", test_jacobi_diagonal), &
      & new_unittest("model-setup", test_model_setup), &
      & new_unittest("model-citation", test_model_citation), &
      & new_unittest("reference-c6-mb01", test_reference_c6_mb01), &
      & new_unittest("atomic-c6-mb01", test_atomic_c6_mb01), &
      & new_unittest("atomic-c6-acetic", test_atomic_c6_acetic), &
      & new_unittest("rank-accuracy-mb01", test_rank_accuracy_mb01), &
      & new_unittest("truncation-mb01", test_truncation_mb01), &
      & new_unittest("ghost-atoms-mb01", test_ghost_atoms_mb01), &
      & new_unittest("energy-mb01", test_energy_mb01), &
      & new_unittest("energy-ewald-acetic", test_energy_acetic), &
      & new_unittest("gradient-slab", test_gradient_acetic), &
      & new_unittest("hessian-mb01", test_hessian_mb01), &
      & new_unittest("fourier-transform", test_fourier_transform), &
      & new_unittest("ewald-molecular", test_ewald_molecular), &
      & new_unittest("ewald-unsupported-damping", test_ewald_unsupported), &
      & new_unittest("ewald-hessian-rejected", test_ewald_hessian), &
      & new_unittest("ewald-pairwise-rejected", test_ewald_pairwise), &
      & new_unittest("ewald-converged-bj", test_ewald_converged_bj), &
      & new_unittest("ewald-converged-zero", test_ewald_converged_zero), &
      & new_unittest("ewald-kcut-convergence", test_ewald_kcut), &
      & new_unittest("ewald-supercell", test_ewald_supercell), &
      & new_unittest("ewald-translation", test_ewald_translation), &
      & new_unittest("ewald-gradient", test_ewald_gradient), &
      & new_unittest("ewald-sigma", test_ewald_sigma), &
      & new_unittest("cutoff-accuracy-cubic", test_cutoff_accuracy_cubic), &
      & new_unittest("cutoff-accuracy-benzene", test_cutoff_accuracy_benzene), &
      & new_unittest("cutoff-accuracy-slab", test_cutoff_accuracy_slab), &
      & new_unittest("cutoff-accuracy-molecular", test_cutoff_molecular) &
      & ]

end subroutine collect_fourier


!> Cubic test cell with two species
subroutine get_cubic(mol, alat, rep, periodic)

   !> Molecular structure data
   type(structure_type), intent(out) :: mol

   !> Lattice parameter of the primitive cell
   real(wp), intent(in) :: alat

   !> Number of repetitions along each axis
   integer, intent(in) :: rep

   !> Directions with periodic boundary conditions
   logical, intent(in), optional :: periodic(:)

   integer :: ix, iy, iz, iat
   integer, allocatable :: num(:)
   real(wp), allocatable :: xyz(:, :)
   real(wp) :: lattice(3, 3), shift(3)
   integer, parameter :: base(2) = [6, 8]
   real(wp), parameter :: frac(3, 2) = reshape(&
      & [0.05_wp, 0.10_wp, 0.15_wp, 0.55_wp, 0.35_wp, 0.45_wp], [3, 2])

   allocate(num(2*rep**3), xyz(3, 2*rep**3))
   iat = 0
   do ix = 0, rep - 1
      do iy = 0, rep - 1
         do iz = 0, rep - 1
            shift(:) = alat * [real(ix, wp), real(iy, wp), real(iz, wp)]
            num(iat+1:iat+2) = base
            xyz(:, iat+1) = alat*frac(:, 1) + shift
            xyz(:, iat+2) = alat*frac(:, 2) + shift
            iat = iat + 2
         end do
      end do
   end do

   lattice(:, :) = 0.0_wp
   lattice(1, 1) = alat*rep
   lattice(2, 2) = alat*rep
   lattice(3, 3) = alat*rep

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)

end subroutine get_cubic


!> Cyclic Jacobi rotations reproduce the input matrix and orthogonal eigenvectors
subroutine test_jacobi(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, jj
   integer, parameter :: ndim = 12
   real(wp) :: amat(ndim, ndim), eval(ndim), evec(ndim, ndim)
   real(wp) :: unity(ndim, ndim), rebuilt(ndim, ndim), val

   do ii = 1, ndim
      do jj = 1, ii
         val = sin(1.0_wp*ii*jj) + 0.5_wp*cos(2.0_wp*(ii + jj))
         amat(jj, ii) = val
         amat(ii, jj) = val
      end do
   end do

   call symmetric_eigendecomposition(amat, eval, evec)

   unity = matmul(transpose(evec), evec)
   do ii = 1, ndim
      unity(ii, ii) = unity(ii, ii) - 1.0_wp
   end do
   if (any(abs(unity) > thr)) then
      call test_failed(error, "Eigenvectors are not orthonormal")
      return
   end if

   rebuilt = 0.0_wp
   do ii = 1, ndim
      do jj = 1, ndim
         rebuilt(jj, ii) = sum(eval * evec(jj, :) * evec(ii, :))
      end do
   end do
   if (any(abs(rebuilt - amat) > thr)) then
      call test_failed(error, "Eigendecomposition does not reproduce the matrix")
      return
   end if

   do ii = 2, ndim
      if (abs(eval(ii)) > abs(eval(ii-1))) then
         call test_failed(error, "Eigenvalues are not sorted by magnitude")
         return
      end if
   end do

end subroutine test_jacobi


!> An already diagonal matrix must be returned unchanged up to the sorting
subroutine test_jacobi_diagonal(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   integer, parameter :: ndim = 4
   real(wp) :: amat(ndim, ndim), eval(ndim), evec(ndim, ndim)
   real(wp), parameter :: ref(ndim) = [3.0_wp, -2.0_wp, 1.0_wp, 0.5_wp]

   amat(:, :) = 0.0_wp
   amat(1, 1) = 1.0_wp
   amat(2, 2) = -2.0_wp
   amat(3, 3) = 3.0_wp
   amat(4, 4) = 0.5_wp

   call symmetric_eigendecomposition(amat, eval, evec)

   do ii = 1, ndim
      call check(error, eval(ii), ref(ii), thr=thr)
      if (allocated(error)) return
   end do

end subroutine test_jacobi_diagonal


!> The expansion is only set up when the optional configuration is provided
subroutine test_model_setup(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3

   call get_structure(mol, "MB16-43", "01")

   call new_d3_model(d3, mol)
   call check(error, allocated(d3%lowrank), .false.)
   if (allocated(error)) return

   call new_d3_model(d3, mol, lowrank=exact)
   call check(error, allocated(d3%lowrank), .true.)

end subroutine test_model_setup


!> Setting up the expansion suggests the reference for the fast summation method
subroutine test_model_citation(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(citation_type) :: citation

   call get_structure(mol, "MB16-43", "01")

   call new_d3_model(d3, mol, citation=citation)
   call check(error, is_citation_present(citation), .false.)
   if (allocated(error)) return

   call new_d3_model(d3, mol, lowrank=exact, citation=citation)
   call check(error, is_citation_present(citation), .true.)
   if (allocated(error)) return

   call check(error, citation%doi, doi_fourier_d3)

end subroutine test_model_citation


!> Full rank expansion reproduces the reference C6 coefficients
subroutine test_reference_c6_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   integer :: isp, jsp, iref, jref, ndim
   real(wp) :: refc6, val

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol, lowrank=exact)

   ndim = sum(d3%ref)
   call check(error, d3%lowrank%rank, ndim)
   if (allocated(error)) return

   do isp = 1, mol%nid
      do iref = 1, d3%ref(isp)
         do jsp = 1, mol%nid
            do jref = 1, d3%ref(jsp)
               refc6 = d3%c6(iref, jref, isp, jsp)
               val = sum(d3%lowrank%lambda * d3%lowrank%vec(iref, isp, :) &
                  & * d3%lowrank%vec(jref, jsp, :))
               if (abs(val - refc6) > thr2*max(abs(refc6), 1.0_wp)) then
                  call test_failed(error, "Reference C6 coefficients are not reproduced")
                  return
               end if
            end do
         end do
      end do
   end do

end subroutine test_reference_c6_mb01


!> Pair C6 coefficients and their derivatives agree with the direct evaluation
subroutine test_atomic_c6_gen(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   type(d3_model) :: d3, d3lr
   integer :: mref
   real(wp), allocatable :: cn(:), lattr(:, :)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwd2cn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), d2c6dcn2(:, :), d2c6dcnij(:, :)
   real(wp), allocatable :: c6r(:, :), dc6dcnr(:, :), d2c6dcn2r(:, :), d2c6dcnijr(:, :)

   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=exact)

   mref = maxval(d3%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat), gwdcn(mref, mol%nat), &
      & gwd2cn(mref, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn)
   call d3%weight_references(mol, cn, gwvec, gwdcn, gwd2cn)

   allocate(c6(mol%nat, mol%nat), dc6dcn(mol%nat, mol%nat), &
      & d2c6dcn2(mol%nat, mol%nat), d2c6dcnij(mol%nat, mol%nat))
   allocate(c6r(mol%nat, mol%nat), dc6dcnr(mol%nat, mol%nat), &
      & d2c6dcn2r(mol%nat, mol%nat), d2c6dcnijr(mol%nat, mol%nat))

   call d3%get_atomic_c6(mol, gwvec, gwdcn, c6r, dc6dcnr, gwd2cn, d2c6dcn2r, d2c6dcnijr)
   call d3lr%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn, gwd2cn, d2c6dcn2, d2c6dcnij)

   if (any(abs(c6 - c6r) > thr2*max(maxval(abs(c6r)), 1.0_wp))) then
      call test_failed(error, "Pair C6 coefficients do not match")
      return
   end if

   if (any(abs(dc6dcn - dc6dcnr) > thr2*max(maxval(abs(dc6dcnr)), 1.0_wp))) then
      call test_failed(error, "Derivatives of the C6 coefficients do not match")
      return
   end if

   if (any(abs(d2c6dcn2 - d2c6dcn2r) > thr2*max(maxval(abs(d2c6dcn2r)), 1.0_wp))) then
      call test_failed(error, "Second derivatives of the C6 coefficients do not match")
      return
   end if

   if (any(abs(d2c6dcnij - d2c6dcnijr) > thr2*max(maxval(abs(d2c6dcnijr)), 1.0_wp))) then
      call test_failed(error, "Mixed derivatives of the C6 coefficients do not match")
      return
   end if

end subroutine test_atomic_c6_gen


subroutine test_atomic_c6_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_atomic_c6_gen(error, mol)

end subroutine test_atomic_c6_mb01


subroutine test_atomic_c6_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "acetic")
   call test_atomic_c6_gen(error, mol)

end subroutine test_atomic_c6_acetic


!> The reconstruction error of the reference coefficients vanishes at full rank
subroutine test_rank_accuracy_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   integer :: ndim
   real(wp) :: first, half, last

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   ndim = sum(d3%ref)

   call new_d3_model(d3, mol, lowrank=d3_lowrank_config(rank=1))
   call check(error, d3%lowrank%rank, 1)
   if (allocated(error)) return
   first = d3%lowrank%error

   call new_d3_model(d3, mol, lowrank=d3_lowrank_config(rank=ndim/2))
   call check(error, d3%lowrank%rank, ndim/2)
   if (allocated(error)) return
   half = d3%lowrank%error

   call new_d3_model(d3, mol, lowrank=exact)
   call check(error, d3%lowrank%rank, ndim)
   if (allocated(error)) return
   last = d3%lowrank%error

   if (.not.(first > half .and. half > last)) then
      call test_failed(error, "Reconstruction error does not decrease with the rank")
      return
   end if

   if (last > thr2) then
      call test_failed(error, "Full rank expansion is not exact")
      return
   end if

end subroutine test_rank_accuracy_mb01


!> A truncated expansion honours the requested accuracy of the pair coefficients
subroutine test_truncation_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   integer :: mref, ndim
   real(wp), allocatable :: cn(:), lattr(:, :)
   real(wp), allocatable :: gwvec(:, :), c6(:, :), c6ref(:, :)
   real(wp), parameter :: tolerance = 1.0e-4_wp

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=d3_lowrank_config(tolerance=tolerance))

   ndim = sum(d3%ref)
   if (d3lr%lowrank%rank >= ndim) then
      call test_failed(error, "Truncation does not reduce the rank")
      return
   end if
   if (d3lr%lowrank%error > tolerance) then
      call test_failed(error, "Requested accuracy is not reached")
      return
   end if

   mref = maxval(d3%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn)
   call d3%weight_references(mol, cn, gwvec)

   allocate(c6ref(mol%nat, mol%nat), c6(mol%nat, mol%nat))
   call d3%get_atomic_c6(mol, gwvec, c6=c6ref)
   call d3lr%get_atomic_c6(mol, gwvec, c6=c6)

   if (any(abs(c6 - c6ref) > 10.0_wp*tolerance*abs(c6ref))) then
      call test_failed(error, "Pair C6 coefficients exceed the requested accuracy")
      return
   end if

end subroutine test_truncation_mb01


!> Ghost atoms must not contribute to the pair coefficients
subroutine test_ghost_atoms_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   integer :: mref
   real(wp), allocatable :: cn(:), lattr(:, :)
   real(wp), allocatable :: gwvec(:, :), c6(:, :)
   integer, parameter :: ghost(*) = [2, 5, 7]

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol, ghost=ghost, lowrank=exact)

   mref = maxval(d3%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn)
   call d3%weight_references(mol, cn, gwvec)

   allocate(c6(mol%nat, mol%nat))
   call d3%get_atomic_c6(mol, gwvec, c6=c6)

   if (any(abs(c6(:, ghost)) > 0.0_wp) .or. any(abs(c6(ghost, :)) > 0.0_wp)) then
      call test_failed(error, "Ghost atoms contribute to the pair coefficients")
      return
   end if

end subroutine test_ghost_atoms_mb01


!> The full rank expansion reproduces the dispersion energy exactly
subroutine test_energy_gen(error, mol, config, thr_energy)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Setup of the separable representation
   type(d3_lowrank_config), intent(in) :: config

   !> Threshold for the energy comparison
   real(wp), intent(in) :: thr_energy

   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   real(wp), allocatable :: energies(:), energies_ref(:)

   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=config)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), energies_ref(mol%nat))
   call get_dispersion(mol, d3, param, cutoff, energies_ref)
   call get_dispersion(mol, d3lr, param, cutoff, energies)

   call check(error, sum(energies), sum(energies_ref), thr=thr_energy)
   if (allocated(error)) then
      print "(2es24.16)", sum(energies), sum(energies_ref)
      return
   end if

   if (any(abs(energies - energies_ref) > thr_energy)) then
      call test_failed(error, "Atom-resolved energies do not match")
      return
   end if

end subroutine test_energy_gen


subroutine test_energy_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_energy_gen(error, mol, exact, 1.0e-10_wp)
   if (allocated(error)) return

   call test_energy_gen(error, mol, d3_lowrank_config(tolerance=1.0e-6_wp), 1.0e-8_wp)

end subroutine test_energy_mb01


!> On a periodic system the low-rank model takes the Ewald path, so the accuracy
!> of the expansion is measured against a full rank reciprocal summation
subroutine test_energy_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   real(wp), allocatable :: energies(:), energies_ref(:)

   call get_structure(mol, "X23", "acetic")
   call new_d3_model(d3, mol, lowrank=tight)
   call new_d3_model(d3lr, mol, &
      & lowrank=d3_lowrank_config(tolerance=1.0e-6_wp, kcut=10.0_wp))
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), energies_ref(mol%nat))
   call get_dispersion(mol, d3, param, cutoff, energies_ref)
   call get_dispersion(mol, d3lr, param, cutoff, energies)

   call check(error, sum(energies), sum(energies_ref), thr=1.0e-8_wp)
   if (allocated(error)) then
      print "(2es24.16)", sum(energies), sum(energies_ref)
   end if

end subroutine test_energy_acetic


!> Gradient and virial are unaffected by a full rank expansion.
!>
!> A slab is not fully periodic and therefore stays on the real space path, so
!> the comparison isolates the expansion from the summation technique.
subroutine test_gradient_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   real(wp) :: sigma(3, 3), sigma_ref(3, 3)
   real(wp), allocatable :: energies(:), energies_ref(:)
   real(wp), allocatable :: gradient(:, :), gradient_ref(:, :)

   call get_cubic(mol, 8.0_wp, 2, periodic=[.true., .true., .false.])
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=exact)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), energies_ref(mol%nat))
   allocate(gradient(3, mol%nat), gradient_ref(3, mol%nat))

   call get_dispersion(mol, d3, param, cutoff, energies_ref, gradient_ref, sigma_ref)
   call get_dispersion(mol, d3lr, param, cutoff, energies, gradient, sigma)

   call check(error, sum(energies), sum(energies_ref), thr=1.0e-10_wp)
   if (allocated(error)) return

   if (any(abs(gradient - gradient_ref) > thr2)) then
      call test_failed(error, "Gradients do not match")
      return
   end if

   if (any(abs(sigma - sigma_ref) > thr2)) then
      call test_failed(error, "Virials do not match")
      return
   end if

end subroutine test_gradient_acetic


!> The second derivatives of the C6 coefficients reach the hessian unchanged
subroutine test_hessian_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   real(wp), allocatable :: energies(:), hessian(:, :), hessian_ref(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=exact)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat))
   allocate(hessian(3*mol%nat, 3*mol%nat), hessian_ref(3*mol%nat, 3*mol%nat))

   call get_dispersion(mol, d3, param, cutoff, energies, hessian=hessian_ref)
   call get_dispersion(mol, d3lr, param, cutoff, energies, hessian=hessian)

   if (any(abs(hessian - hessian_ref) > thr2)) then
      call test_failed(error, "Hessians do not match")
      return
   end if

end subroutine test_hessian_mb01


!> The closed form transform agrees with the independently derived expression for
!> the leading Becke-Johnson term and with its limit at vanishing wave number
subroutine test_fourier_transform(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ik
   real(wp) :: kval, bval, phi, dphi, ref, kb
   type(fourier_term) :: term
   real(wp), parameter :: pi = 3.1415926535897932384626433832795029_wp

   do ik = 0, 20
      kval = 0.05_wp * ik
      bval = 4.5_wp
      term = fourier_term(1.0_wp, 0, 6, bval)
      call get_fourier_transform(term, kval, phi, dphi)

      kb = kval * bval
      if (kval > 0.0_wp) then
         ref = 2.0_wp*pi**2/(3.0_wp*kval*bval**4) &
            & * (2.0_wp*sin(sqrt(3.0_wp)*kb/2 - pi/6)*exp(-kb/2) + exp(-kb))
      else
         ref = 2.0_wp*pi**2/(3.0_wp*bval**3)
      end if

      call check(error, phi, ref, thr=1.0e-10_wp*abs(ref))
      if (allocated(error)) return
   end do

   ! numerical derivative of the transform
   term = fourier_term(1.0_wp, 0, 6, 4.5_wp)
   call get_fourier_transform(term, 1.0_wp, phi, dphi)
   call get_fourier_transform(term, 1.0_wp + 1.0e-6_wp, ref, kb)
   call get_fourier_transform(term, 1.0_wp - 1.0e-6_wp, kval, kb)
   call check(error, dphi, 0.5_wp*(ref - kval)/1.0e-6_wp, thr=1.0e-6_wp)
   if (allocated(error)) return

   ! bounded potential at the origin
   call check(error, get_potential_zero(term), 1.0_wp/4.5_wp**6, thr=thr)
   if (allocated(error)) return
   call check(error, get_potential_zero(fourier_term(1.0_wp, 8, 14, 4.5_wp)), &
      & 0.0_wp, thr=thr)

end subroutine test_fourier_transform


!> Without periodicity the Ewald path must not be taken
subroutine test_ewald_molecular(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   real(wp), allocatable :: energies(:), energies_ref(:)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), energies_ref(mol%nat))
   call get_dispersion(mol, d3, param, cutoff, energies_ref)
   call get_dispersion(mol, d3lr, param, cutoff, energies)

   call check(error, sum(energies), sum(energies_ref), thr=1.0e-10_wp)

end subroutine test_ewald_molecular


!> A low-rank model on a periodic system requires a damping function supporting
!> the reciprocal space summation
subroutine test_ewald_unsupported(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(mzero_damping_param) :: param
   type(error_type), allocatable :: setup_error
   real(wp), allocatable :: energies(:), energies_ref(:)

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_mzero_damping(param, pbe_mzero)

   allocate(energies(mol%nat), energies_ref(mol%nat))
   call get_dispersion(setup_error, mol, d3lr, param, cutoff, energies)
   if (.not.allocated(setup_error)) then
      call test_failed(error, "Unsupported damping function is not reported")
      return
   end if

   ! without the low-rank expansion the same setup is evaluated in real space
   call get_dispersion(setup_error, mol, d3, param, cutoff, energies)
   if (allocated(setup_error)) then
      call test_failed(error, "Real space evaluation is rejected")
      return
   end if

   call get_dispersion(mol, d3, param, cutoff, energies_ref)
   call check(error, sum(energies), sum(energies_ref), thr=thr)

end subroutine test_ewald_unsupported


!> The second derivatives are only available from the real space summation
subroutine test_ewald_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   type(error_type), allocatable :: setup_error
   real(wp), allocatable :: energies(:), hessian(:, :)

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), hessian(3*mol%nat, 3*mol%nat))
   call get_dispersion(setup_error, mol, d3lr, param, cutoff, energies, &
      & hessian=hessian)
   if (.not.allocated(setup_error)) then
      call test_failed(error, "Hessian of a low-rank model is not reported")
      return
   end if

   ! the energy alone is available from the same model
   call get_dispersion(setup_error, mol, d3lr, param, cutoff, energies)
   if (allocated(setup_error)) then
      call test_failed(error, "Ewald energy is rejected")
      return
   end if

   ! without the low-rank expansion the hessian is evaluated in real space
   call get_dispersion(setup_error, mol, d3, param, cutoff, energies, &
      & hessian=hessian)
   if (allocated(setup_error)) then
      call test_failed(error, "Real space hessian is rejected")
      return
   end if

end subroutine test_ewald_hessian


!> The pairwise decomposition would not add up to the reciprocal space energy
subroutine test_ewald_pairwise(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   type(error_type), allocatable :: setup_error
   real(wp), allocatable :: energy2(:, :), energy3(:, :), energies(:)

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   allocate(energy2(mol%nat, mol%nat), energy3(mol%nat, mol%nat), energies(mol%nat))
   call get_pairwise_dispersion(setup_error, mol, d3lr, param, cutoff, &
      & energy2, energy3)
   if (.not.allocated(setup_error)) then
      call test_failed(error, "Pairwise analysis of a low-rank model is not reported")
      return
   end if

   ! without the low-rank expansion the decomposition is available
   call get_pairwise_dispersion(setup_error, mol, d3, param, cutoff, &
      & energy2, energy3)
   if (allocated(setup_error)) then
      call test_failed(error, "Real space pairwise analysis is rejected")
      return
   end if

   call get_dispersion(mol, d3, param, cutoff, energies)
   call check(error, sum(energies), sum(energy2) + sum(energy3), thr=thr)

end subroutine test_ewald_pairwise


!> The reciprocal sum reproduces a converged real space summation
subroutine test_ewald_converged_bj(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   real(wp) :: ewald, converged, truncated
   real(wp), allocatable :: energies(:)

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   ! the low-rank model is given the short cutoff, which the Ewald path ignores
   allocate(energies(mol%nat))
   call get_dispersion(mol, d3lr, param, cutoff, energies)
   ewald = sum(energies)
   call get_dispersion(mol, d3, param, long, energies)
   converged = sum(energies)
   call get_dispersion(mol, d3, param, cutoff, energies)
   truncated = sum(energies)

   if (abs(ewald - converged) > 5.0e-8_wp) then
      call test_failed(error, "Ewald summation does not match converged real space")
      print "(3es24.16)", ewald, converged, truncated
      return
   end if

   if (abs(truncated - converged) < 1.0e-6_wp) then
      call test_failed(error, "Real space truncation error is unexpectedly small")
      print "(3es24.16)", ewald, converged, truncated
      return
   end if

end subroutine test_ewald_converged_bj


subroutine test_ewald_converged_zero(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3, d3lr
   type(zero_damping_param) :: param
   real(wp) :: ewald, converged
   real(wp), allocatable :: energies(:)
   type(d3_lowrank_config), parameter :: config = &
      & d3_lowrank_config(rank=10000, kcut=20.0_wp)

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=config)
   call new_zero_damping(param, pbe_zero)

   allocate(energies(mol%nat))
   call get_dispersion(mol, d3lr, param, long, energies)
   ewald = sum(energies)
   call get_dispersion(mol, d3, param, long, energies)
   converged = sum(energies)

   if (abs(ewald - converged) > 1.0e-8_wp) then
      call test_failed(error, "Ewald summation does not match converged real space")
      print "(2es24.16)", ewald, converged
      return
   end if

end subroutine test_ewald_converged_zero


!> Increasing the reciprocal cutoff converges the energy
subroutine test_ewald_kcut(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3lr
   type(rational_damping_param) :: param
   integer :: istep
   real(wp) :: ecur, eref, diff, last
   real(wp), allocatable :: energies(:)

   call get_cubic(mol, 8.0_wp, 1)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat))
   call new_d3_model(d3lr, mol, lowrank=d3_lowrank_config(rank=10000, kcut=14.0_wp))
   call get_dispersion(mol, d3lr, param, cutoff, energies)
   eref = sum(energies)

   last = huge(1.0_wp)
   do istep = 1, 4
      call new_d3_model(d3lr, mol, &
         & lowrank=d3_lowrank_config(rank=10000, kcut=2.0_wp*istep))
      call get_dispersion(mol, d3lr, param, cutoff, energies)
      ecur = sum(energies)
      diff = abs(ecur - eref)
      if (diff > last) then
         call test_failed(error, "Energy does not converge with the reciprocal cutoff")
         return
      end if
      last = diff
   end do

   if (last > 1.0e-9_wp) then
      call test_failed(error, "Reciprocal summation is not converged")
      return
   end if

end subroutine test_ewald_kcut


!> The energy per atom is independent of the chosen unit cell
subroutine test_ewald_supercell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3lr
   type(rational_damping_param) :: param
   real(wp) :: eprim, esuper
   real(wp), allocatable :: energies(:)

   call new_rational_damping(param, pbe_bj)

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3lr, mol, lowrank=tight)
   allocate(energies(mol%nat))
   call get_dispersion(mol, d3lr, param, cutoff, energies)
   eprim = sum(energies) / mol%nat

   call get_cubic(mol, 8.0_wp, 2)
   call new_d3_model(d3lr, mol, lowrank=tight)
   deallocate(energies)
   allocate(energies(mol%nat))
   call get_dispersion(mol, d3lr, param, cutoff, energies)
   esuper = sum(energies) / mol%nat

   call check(error, esuper, eprim, thr=1.0e-10_wp)
   if (allocated(error)) then
      print "(2es24.16)", eprim, esuper
   end if

end subroutine test_ewald_supercell


!> A rigid shift of all atoms leaves the energy unchanged
subroutine test_ewald_translation(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3lr
   type(rational_damping_param) :: param
   integer :: iat
   real(wp) :: eref, eshift
   real(wp), allocatable :: energies(:)
   real(wp), parameter :: shift(3) = [1.3_wp, -2.7_wp, 4.1_wp]

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat))
   call get_dispersion(mol, d3lr, param, cutoff, energies)
   eref = sum(energies)

   do iat = 1, mol%nat
      mol%xyz(:, iat) = mol%xyz(:, iat) + shift
   end do
   call get_dispersion(mol, d3lr, param, cutoff, energies)
   eshift = sum(energies)

   call check(error, eshift, eref, thr=1.0e-12_wp)

end subroutine test_ewald_translation


!> Analytical gradient of the reciprocal summation
subroutine test_ewald_gradient(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3lr
   type(rational_damping_param) :: param
   integer :: iat, ic
   real(wp) :: sigma(3, 3), er, el
   real(wp), allocatable :: energies(:), gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-4_wp

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), gradient(3, mol%nat), numgrad(3, mol%nat))

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_dispersion(mol, d3lr, param, cutoff, energies)
         er = sum(energies)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_dispersion(mol, d3lr, param, cutoff, energies)
         el = sum(energies)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d3lr, param, cutoff, energies, gradient, sigma)

   if (any(abs(gradient - numgrad) > 1.0e-8_wp)) then
      call test_failed(error, "Gradient does not match numerical differentiation")
      print "(3es21.13)", gradient - numgrad
      return
   end if

   if (any(abs(sum(gradient, 2)) > 1.0e-10_wp)) then
      call test_failed(error, "Sum of forces does not vanish")
      return
   end if

end subroutine test_ewald_gradient


!> Analytical virial of the reciprocal summation
subroutine test_ewald_sigma(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3lr
   type(rational_damping_param) :: param
   integer :: ic, jc
   real(wp) :: sigma(3, 3), numsigma(3, 3), er, el
   real(wp) :: eps(3, 3), lattice(3, 3)
   real(wp), allocatable :: energies(:), gradient(:, :), xyz(:, :)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   call get_cubic(mol, 8.0_wp, 1)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   allocate(energies(mol%nat), gradient(3, mol%nat))
   xyz = mol%xyz
   lattice(:, :) = mol%lattice

   do ic = 1, 3
      do jc = 1, 3
         eps(:, :) = unity
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         call get_dispersion(mol, d3lr, param, cutoff, energies)
         er = sum(energies)

         eps(:, :) = unity
         eps(jc, ic) = eps(jc, ic) - step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         call get_dispersion(mol, d3lr, param, cutoff, energies)
         el = sum(energies)

         numsigma(ic, jc) = 0.5_wp*(er - el)/step
      end do
   end do

   mol%xyz(:, :) = xyz
   mol%lattice(:, :) = lattice
   call get_dispersion(mol, d3lr, param, cutoff, energies, gradient, sigma)

   if (any(abs(sigma - numsigma) > 1.0e-7_wp)) then
      call test_failed(error, "Virial does not match numerical differentiation")
      print "(3es21.13)", sigma - numsigma
      return
   end if

end subroutine test_ewald_sigma


!> The estimated cutoff reproduces the reciprocal space energy to the requested
!> accuracy, and a tighter request produces a longer cutoff
subroutine test_cutoff_accuracy_gen(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: istep
   type(d3_model) :: d3, d3lr
   type(rational_damping_param) :: param
   type(realspace_cutoff) :: rcut, ewald_cut
   real(wp) :: eref, energy, last, achieved
   real(wp), allocatable :: energies(:)
   real(wp), parameter :: accuracy(3) = [1.0e-6_wp, 1.0e-7_wp, 1.0e-8_wp]

   call new_d3_model(d3, mol)
   call new_d3_model(d3lr, mol, lowrank=tight)
   call new_rational_damping(param, pbe_bj)

   ! the estimate keeps the default coordination number cutoff, which has to be
   ! matched by the reference to isolate the pair truncation error
   ewald_cut = realspace_cutoff()

   allocate(energies(mol%nat))
   call get_dispersion(mol, d3lr, param, ewald_cut, energies)
   eref = sum(energies)

   last = 0.0_wp
   do istep = 1, size(accuracy)
      rcut = get_realspace_cutoff(mol, d3, accuracy(istep))
      if (rcut%disp2 <= last) then
         call test_failed(error, "Cutoff does not grow with the requested accuracy")
         return
      end if
      last = rcut%disp2

      call get_dispersion(mol, d3, param, rcut, energies)
      energy = sum(energies)
      achieved = abs(energy - eref) / mol%nat
      if (achieved > accuracy(istep)) then
         call test_failed(error, "Requested accuracy is not reached")
         print "(a, es11.3, a, f9.2, a, es11.3)", "target ", accuracy(istep), &
            & "  cutoff ", rcut%disp2, "  achieved ", achieved
         return
      end if
   end do

end subroutine test_cutoff_accuracy_gen


subroutine test_cutoff_accuracy_cubic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_cubic(mol, 8.0_wp, 1)
   call test_cutoff_accuracy_gen(error, mol)

end subroutine test_cutoff_accuracy_cubic


subroutine test_cutoff_accuracy_benzene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "benzene")
   call test_cutoff_accuracy_gen(error, mol)

end subroutine test_cutoff_accuracy_benzene


!> A slab converges faster than a bulk system of the same density
subroutine test_cutoff_accuracy_slab(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   type(realspace_cutoff) :: rcut, tight_cut
   real(wp) :: eref, energy
   real(wp), allocatable :: energies(:)
   real(wp), parameter :: accuracy = 1.0e-7_wp

   call get_cubic(mol, 8.0_wp, 2, periodic=[.true., .true., .false.])
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe_bj)

   rcut = get_realspace_cutoff(mol, d3, accuracy)
   tight_cut = realspace_cutoff()
   tight_cut%disp2 = 400.0_wp

   allocate(energies(mol%nat))
   call get_dispersion(mol, d3, param, tight_cut, energies)
   eref = sum(energies)
   call get_dispersion(mol, d3, param, rcut, energies)
   energy = sum(energies)

   if (abs(energy - eref)/mol%nat > accuracy) then
      call test_failed(error, "Requested accuracy is not reached for a slab")
      print "(a, f9.2, a, es11.3)", "cutoff ", rcut%disp2, "  achieved ", &
         & abs(energy - eref)/mol%nat
      return
   end if

end subroutine test_cutoff_accuracy_slab


!> A finite system is summed exactly once all pairs are covered
subroutine test_cutoff_molecular(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   type(realspace_cutoff) :: rcut, huge_cut
   real(wp) :: energy, eref
   real(wp), allocatable :: energies(:)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe_bj)

   rcut = get_realspace_cutoff(mol, d3, 1.0e-10_wp)
   huge_cut = realspace_cutoff()
   huge_cut%disp2 = 1000.0_wp

   allocate(energies(mol%nat))
   call get_dispersion(mol, d3, param, rcut, energies)
   energy = sum(energies)
   call get_dispersion(mol, d3, param, huge_cut, energies)
   eref = sum(energies)

   call check(error, energy, eref, thr=thr)

end subroutine test_cutoff_molecular


end module test_fourier
