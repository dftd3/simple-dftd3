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

module test_gcp_hessian
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_gcp, only : gcp_param, get_gcp_param, get_geometric_counterpoise, &
      & get_geometric_counterpoise_hessian
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   implicit none
   private

   public :: collect_gcp_hessian

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp)) * 10


contains


!> Collect all exported unit tests
subroutine collect_gcp_hessian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("hf/dz", test_hf_dz_hess), &
      & new_unittest("hf/minix", test_hf_minix_hess), &
      & new_unittest("dft/lanl2", test_dft_lanl2_hess), &
      & new_unittest("hf3c", test_hf3c_hess), &
      & new_unittest("pbeh3c", test_pbeh3c_hess), &
      & new_unittest("hse3c", test_hse3c_hess), &
      & new_unittest("b973c", test_b973c_hess), &
      & new_unittest("r2scan3c", test_r2scan3c_hess), &
      & new_unittest("pbc:hf3c", test_hf3c_pbc_hess) &
      & ]

end subroutine collect_gcp_hessian


subroutine test_hf_dz_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numhess(error, mol, "hf", "dz")

end subroutine test_hf_dz_hess


subroutine test_hf_minix_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_numhess(error, mol, "hf", "minix")

end subroutine test_hf_minix_hess


subroutine test_dft_lanl2_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_numhess(error, mol, "dft", "lanl")

end subroutine test_dft_lanl2_hess


subroutine test_hf3c_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_numhess(error, mol, "hf3c")

end subroutine test_hf3c_hess


subroutine test_pbeh3c_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_numhess(error, mol, "pbeh3c")

end subroutine test_pbeh3c_hess


subroutine test_hse3c_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_numhess(error, mol, "hse3c")

end subroutine test_hse3c_hess


subroutine test_b973c_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numhess(error, mol, "b973c")

end subroutine test_b973c_hess


subroutine test_r2scan3c_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_numhess(error, mol, "r2scan3c")

end subroutine test_r2scan3c_hess


subroutine test_hf3c_pbc_hess(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_numhess(error, mol, "hf3c")

end subroutine test_hf3c_pbc_hess


!> Compare the analytical hessian against finite differences of the gradient
subroutine test_numhess(error, mol, method, basis)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Method name
   character(len=*), intent(in) :: method

   !> Basis set name
   character(len=*), intent(in), optional :: basis

   type(gcp_param) :: param
   integer :: iat, jat, ic, jc, ii, jj
   real(wp) :: energy, sigma(3, 3)
   real(wp), allocatable :: gradp(:, :), gradm(:, :), xyz(:, :)
   real(wp), allocatable :: hessian(:, :), numhessian(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   call get_gcp_param(param, mol, method=method, basis=basis)

   allocate(gradp(3, mol%nat), gradm(3, mol%nat), xyz(3, mol%nat), &
      & hessian(3*mol%nat, 3*mol%nat), numhessian(3*mol%nat, 3*mol%nat))
   xyz(:, :) = mol%xyz(:, :)
   numhessian(:, :) = 0.0_wp

   do iat = 1, mol%nat
      do ic = 1, 3
         ii = 3*(iat - 1) + ic
         mol%xyz(:, :) = xyz(:, :)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         energy = 0.0_wp
         gradp(:, :) = 0.0_wp
         call get_geometric_counterpoise(mol, param, realspace_cutoff(), energy, &
            & gradp, sigma)
         mol%xyz(:, :) = xyz(:, :)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         energy = 0.0_wp
         gradm(:, :) = 0.0_wp
         call get_geometric_counterpoise(mol, param, realspace_cutoff(), energy, &
            & gradm, sigma)
         do jat = 1, mol%nat
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               numhessian(ii, jj) = 0.5_wp*(gradp(jc, jat) - gradm(jc, jat))/step
            end do
         end do
      end do
   end do
   numhessian(:, :) = 0.5_wp*(numhessian + transpose(numhessian))
   mol%xyz(:, :) = xyz(:, :)

   call get_geometric_counterpoise_hessian(mol, param, realspace_cutoff(), hessian)

   ! the hessian is not symmetrized, it has to come out symmetric by construction
   if (any(abs(hessian - transpose(hessian)) > thr)) then
      call test_failed(error, "Hessian of counter-poise energy is not symmetric")
      print"(3es21.14)", hessian - transpose(hessian)
      return
   end if

   if (any(abs(hessian - numhessian) > thr2)) then
      call test_failed(error, &
         & "Hessian of counter-poise energy does not match finite differences")
      print"(3es21.14)", hessian - numhessian
   end if

end subroutine test_numhess


end module test_gcp_hessian
