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

module test_hessian
   use dftd3, only : get_dispersion, realspace_cutoff, damping_param, d3_param, &
      & cso_damping_param, new_cso_damping, mzero_damping_param, new_mzero_damping, &
      & optimizedpower_damping_param, new_optimizedpower_damping, &
      & rational_damping_param, new_rational_damping, zero_damping_param, &
      & new_zero_damping, z_damping_param, new_z_damping, d3_model, new_d3_model
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   implicit none
   private

   public :: collect_hessian

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_hessian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D3(BJ) hessian", test_pbed3bj_hessian), &
      & new_unittest("RPBE-D3(0) hessian", test_rpbed3zero_hessian), &
      & new_unittest("PBE-D3(0M) hessian", test_pbed3zerom_hessian), &
      & new_unittest("B97h-D3(op) hessian", test_b97hd3op_hessian), &
      & new_unittest("PBE-D3(Z) hessian", test_pbed3z_hessian), &
      & new_unittest("DSD-BLYP-D3(BJ)-ATM hessian", test_dsdblypd3bjatm_hessian), &
      & new_unittest("B3LYP-D3(CSO) hessian", test_b3lypd3cso_hessian), &
      & new_unittest("PBE-D3(CSO)-ATM hessian", test_pbed3csoatm_hessian), &
      & new_unittest("BLYP-D3(0)-ATM hessian", test_blypd3zeroatm_hessian), &
      & new_unittest("PBE-D3(BJ)-ATM smooth cutoff hessian", &
      &   test_pbed3bjatm_smooth_cutoff_hessian) &
      & ]

end subroutine collect_hessian


!> Compare the analytical hessian against finite differences of the gradient
subroutine test_numhessian(error, mol, param, cutoff)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoff
   type(realspace_cutoff), intent(in), optional :: cutoff

   integer :: iat, jat, ic, jc, ii, jj
   type(d3_model) :: d3
   type(realspace_cutoff) :: rcut
   real(wp) :: energy, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), gradp(:, :), gradm(:, :)
   real(wp), allocatable :: hessian(:, :), numhessian(:, :), xyz(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), gradp(3, mol%nat), gradm(3, mol%nat), &
      & hessian(3*mol%nat, 3*mol%nat), numhessian(3*mol%nat, 3*mol%nat), xyz(3, mol%nat))
   call new_d3_model(d3, mol)
   rcut = realspace_cutoff()
   if (present(cutoff)) rcut = cutoff
   xyz(:, :) = mol%xyz(:, :)
   hessian(:, :) = 0.0_wp
   numhessian(:, :) = 0.0_wp

   do iat = 1, mol%nat
      do ic = 1, 3
         ii = 3*(iat - 1) + ic
         mol%xyz(:, :) = xyz(:, :)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_dispersion(mol, d3, param, rcut, energy, gradp, sigma)
         mol%xyz(:, :) = xyz(:, :)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call get_dispersion(mol, d3, param, rcut, energy, gradm, sigma)
         mol%xyz(:, :) = xyz(:, :)
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

   call get_dispersion(mol, d3, param, rcut, energy, gradient, sigma, hessian=hessian)

   ! the hessian is not symmetrized, it has to come out symmetric by construction
   if (any(abs(hessian - transpose(hessian)) > thr)) then
      call test_failed(error, "Hessian of dispersion energy is not symmetric")
      print"(3es21.14)", hessian - transpose(hessian)
      return
   end if

   if (any(abs(hessian - numhessian) > thr2)) then
      call test_failed(error, "Hessian of dispersion energy does not match finite differences")
      print"(3es21.14)", hessian-numhessian
   end if

end subroutine test_numhessian


subroutine test_pbed3bj_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_rational_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_pbed3bj_hessian


subroutine test_rpbed3zero_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 0.872_wp, s8 = 0.514_wp)

   call get_structure(mol, "MB16-43", "09")
   call new_zero_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_rpbed3zero_hessian


subroutine test_pbed3zerom_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 2.340218_wp, s8 = 0.0_wp, bet = 0.129434_wp)

   call get_structure(mol, "MB16-43", "33")
   call new_mzero_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_pbed3zerom_hessian


subroutine test_b97hd3op_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 0.97388_wp, s9 = 1.0_wp, alp = 14.0_wp, bet = 6.0_wp, &
      & a1 = 0.150_wp, s8 = 0.0_wp, a2 = 4.25_wp)

   call get_structure(mol, "MB16-43", "37")
   call new_optimizedpower_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_b97hd3op_hessian


subroutine test_pbed3z_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(z_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 200770.0_wp, a2 = 0.0_wp, rs6 = 0.0_wp, rs8 = 6.25_wp)

   call get_structure(mol, "MB16-43", "39")
   call new_z_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_pbed3z_hessian


subroutine test_dsdblypd3bjatm_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 0.5_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.0_wp, s8 = 0.2130_wp, a2 = 6.0519_wp)

   call get_structure(mol, "MB16-43", "17")
   call new_rational_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_dsdblypd3bjatm_hessian


subroutine test_b3lypd3cso_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cso_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.86_wp, a2 = 2.5_wp, rs6 = 0.0_wp, rs8 = 6.25_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_cso_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_b3lypd3cso_hessian


subroutine test_pbed3csoatm_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cso_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.24_wp, a2 = 2.5_wp, rs6 = 0.0_wp, rs8 = 6.25_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_cso_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_pbed3csoatm_hessian


subroutine test_blypd3zeroatm_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.094_wp, s8 = 1.682_wp)

   call get_structure(mol, "MB16-43", "25")
   call new_zero_damping(param, inp)
   call test_numhessian(error, mol, param)

end subroutine test_blypd3zeroatm_hessian


subroutine test_pbed3bjatm_smooth_cutoff_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param), parameter :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)
   type(realspace_cutoff), parameter :: cutoff = realspace_cutoff(&
      & cn=40.0_wp, disp2=8.0_wp, disp3=8.0_wp, width2=4.0_wp, width3=4.0_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_rational_damping(param, inp)
   call test_numhessian(error, mol, param, cutoff)

end subroutine test_pbed3bjatm_smooth_cutoff_hessian


end module test_hessian
