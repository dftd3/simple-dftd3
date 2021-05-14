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

module test_periodic
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_periodic

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 100*sqrt(epsilon(1.0_wp))
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_periodic(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D3(BJ)", test_pbed3bj_acetic), &
      & new_unittest("PBEsol-D3(BJ)", test_pbesold3bj_adaman), &
      & new_unittest("TPSS-D3(BJ)", test_tpssd3bj_ammonia), &
      & new_unittest("HSE06-D3(BJ)", test_hse06d3bj_anthracene), &
      & new_unittest("BLYP-D3(0)", test_blypd3zero_benzene), &
      & new_unittest("M06L-D3(0)", test_m06ld3zero_cyanamide), &
      & new_unittest("rPW86PBE-D3(0)", test_rpw86pbed3zero_co2), &
      & new_unittest("revSSB-D3(0)", test_revssbd3zero_cytosine) &
      ! new_unittest("HSEsol-D3(BJ)-ATM", test_hsesold3bjatm_oxacb), &
      ! new_unittest("PWGGA-D3(BJ)-ATM", test_pwggad3bjatm_pyrazine), &
      ! new_unittest("B3PW91-D3(0)-ATM", test_b3pw91d3zeroatm_urea), &
      ! new_unittest("RPBE-D3(0)-ATM", test_rpbed3zeroatm_hexamine) &
      & ]

end subroutine collect_periodic


subroutine test_dftd3_gen(error, mol, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   type(d3_model) :: d3
   real(wp) :: energy

   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, cutoff, energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_dftd3_gen


subroutine test_numgrad(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   integer :: iat, ic
   type(d3_model) :: d3
   real(wp) :: energy, er, el, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   call new_d3_model(d3, mol)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_dispersion(mol, d3, param, cutoff, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_dispersion(mol, d3, param, cutoff, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d3, param, cutoff, energy, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   integer :: ic, jc
   type(d3_model) :: d3
   real(wp) :: energy, er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3), lattice(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-7_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   call new_d3_model(d3, mol)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         call get_dispersion(mol, d3, param, cutoff, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         call get_dispersion(mol, d3, param, cutoff, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d3, param, cutoff, energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr3)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_pbed3bj_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   call get_structure(mol, "X23", "acetic")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -6.673283761139015E-2_wp)

end subroutine test_pbed3bj_acetic


subroutine test_pbesold3bj_adaman(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4466_wp, s8 = 2.9491_wp, a2 = 6.1742_wp)

   call get_structure(mol, "X23", "adaman")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -7.931351552127060E-2_wp)

end subroutine test_pbesold3bj_adaman


subroutine test_tpssd3bj_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4535_wp, s8 = 1.9435_wp, a2 = 4.4752_wp)

   call get_structure(mol, "X23", "ammonia")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_tpssd3bj_ammonia


subroutine test_hse06d3bj_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.383_wp, s8 = 2.310_wp, a2 = 5.685_wp)

   call get_structure(mol, "X23", "anthracene")
   call new_rational_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_hse06d3bj_anthracene


subroutine test_blypd3zero_benzene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.094_wp, s8 = 1.682_wp)

   call get_structure(mol, "X23", "benzene")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -0.16647369651118463_wp)

end subroutine test_blypd3zero_benzene


subroutine test_m06ld3zero_cyanamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.581_wp, s8 = 0.000_wp)

   call get_structure(mol, "X23", "cyanamide")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -2.3222572072479526E-2_wp)

end subroutine test_m06ld3zero_cyanamide


subroutine test_rpw86pbed3zero_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.224_wp, s8 = 0.901_wp)

   call get_structure(mol, "X23", "CO2")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_rpw86pbed3zero_co2


subroutine test_revssbd3zero_cytosine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.221_wp, s8 = 0.560_wp)

   call get_structure(mol, "X23", "cytosine")
   call new_zero_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_revssbd3zero_cytosine


subroutine test_hsesold3bjatm_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.4650_wp, s8 = 2.9215_wp, a2 = 6.2003_wp)

   call get_structure(mol, "X23", "oxacb")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_hsesold3bjatm_oxacb


subroutine test_pwggad3bjatm_pyrazine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.2211_wp, s8 = 2.6910_wp, a2 = 6.7278_wp)

   call get_structure(mol, "X23", "pyrazine")
   call new_rational_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_pwggad3bjatm_pyrazine


subroutine test_b3pw91d3zeroatm_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.176_wp, s8 = 1.775_wp)

   call get_structure(mol, "X23", "urea")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_b3pw91d3zeroatm_urea


subroutine test_rpbed3zeroatm_hexamine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 0.872_wp, s8 = 0.514_wp)

   call get_structure(mol, "X23", "hexamine")
   call new_zero_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_rpbed3zeroatm_hexamine


end module test_periodic
