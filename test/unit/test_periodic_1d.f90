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

module test_periodic_1d
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_periodic_1d

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 1000*sqrt(epsilon(1.0_wp))
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_periodic_1d(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("revPBE-D3(BJ)", test_revpbed3bj_1d) &
      & ]

end subroutine collect_periodic_1d


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


subroutine test_revpbed3bj_1d(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 32
   integer, parameter :: num(nat) = 6
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 1.36794785746435_wp, 13.45808943446053_wp,  8.83754983226359_wp, &
      & 3.69183290816438_wp, 13.13552229161569_wp, 10.16652201690950_wp, &
      & 1.36792668081267_wp, 10.38660504434782_wp, 13.04411926632965_wp, &
      & 3.69180534781206_wp, 11.55414582295511_wp, 12.33193380846742_wp, &
      & 1.36791549262702_wp,  3.53066844289674_wp, 10.38660588677206_wp, &
      & 1.36792046664920_wp,  7.73723910626293_wp, 13.45809224934817_wp, &
      & 3.69181279359489_wp,  6.40826717723392_wp, 13.13552570942280_wp, &
      & 1.36792009865062_wp,  3.11669712338516_wp,  7.73723850632628_wp, &
      & 3.69181515738094_wp,  3.43926499914873_wp,  6.40826580885474_wp, &
      & 3.69178443989294_wp,  4.24285720771059_wp, 11.55415026712869_wp, &
      & 1.36790824853106_wp,  6.18818490375705_wp,  3.53066863732142_wp, &
      & 3.69178194163078_wp,  5.02063901427657_wp,  4.24285736953327_wp, &
      & 1.36794124909207_wp, 13.04411858182861_wp,  6.18818324080182_wp, &
      & 1.36792249732236_wp,  8.83755133592807_wp,  3.11669686076913_wp, &
      & 3.69182456413952_wp, 10.16652118921143_wp,  3.43926084011816_wp, &
      & 3.69181444966104_wp, 12.33193631088573_wp,  5.02063847821044_wp, &
      & 6.01572566324028_wp, 13.45790756713123_wp,  8.83752222635545_wp, &
      & 8.33965926123256_wp, 13.13576644753615_wp, 10.16660228658307_wp, &
      & 6.01574747573805_wp, 10.38654070512969_wp, 13.04391961251944_wp, &
      & 8.33964066450677_wp, 11.55427002850905_wp, 12.33211653730939_wp, &
      & 6.01574728097580_wp,  3.53087013230607_wp, 10.38654217813321_wp, &
      & 6.01568913853645_wp,  7.73726406411719_wp, 13.45790864082374_wp, &
      & 8.33963586549168_wp,  6.40818371470975_wp, 13.13576911116618_wp, &
      & 6.01568179676984_wp,  3.11688332536281_wp,  7.73726611148835_wp, &
      & 8.33963704688671_wp,  3.43902559351770_wp,  6.40818390180453_wp, &
      & 8.33962496288127_wp,  4.24267007149867_wp, 11.55427031066552_wp, &
      & 6.01573464280675_wp,  6.18824653544318_wp,  3.53086861480278_wp, &
      & 8.33961857277245_wp,  5.02052001792996_wp,  4.24267413625204_wp, &
      & 6.01575677304189_wp, 13.04392044501564_wp,  6.18824448603611_wp, &
      & 6.01568344836224_wp,  8.83752193432504_wp,  3.11688171781516_wp, &
      & 8.33964228963694_wp, 10.16660428027860_wp,  3.43902155668011_wp, &
      & 8.33965118613331_wp, 12.33211762632282_wp,  5.02051902430387_wp],&
      & shape(xyz))
   logical, parameter :: periodic(3) = [.true., .false., .false.]
   real(wp), parameter :: lattice(3, 3) = reshape([&
      & 9.2955628527586_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & shape(lattice))
   real(wp), parameter :: lattice_vac(3, 3) = reshape([&
      & 9.2955628527586_wp, 0.0_wp, 0.0_wp, 0.0_wp, 100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 100.0_wp], &
      & shape(lattice))
   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.5238_wp, s8 = 2.3550_wp, a2 = 3.5016_wp)

   call new(mol, num, xyz, periodic=periodic, lattice=lattice)
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -0.25436101196186600_wp)

   call new(mol, num, xyz, periodic=[.true.], lattice=lattice_vac)
   call test_dftd3_gen(error, mol, param, -0.25436101196186600_wp)

end subroutine test_revpbed3bj_1d


end module test_periodic_1d
